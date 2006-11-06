/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_insert.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix_gravity.h"
#include "fix_shear_history.h"
#include "neighbor.h"
#include "domain.h"
#include "region.h"
#include "region_block.h"
#include "region_cylinder.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

FixInsert::FixInsert(int narg, char **arg) : Fix(narg, arg)
{
  if (narg < 6) error->all("Illegal fix insert command");

  if (atom->check_style("granular") == 0)
    error->all("Must use fix insert with atom style granular");

  // required args

  ninsert = atoi(arg[3]);
  ntype = atoi(arg[4]);
  seed = atoi(arg[5]);

  PI = 4.0*atan(1.0);

  // option defaults

  int iregion = -1;
  radius_lo = radius_hi = 0.5;
  density_lo = density_hi = 1.0;
  volfrac = 0.25;
  maxattempt = 50;
  rate = 0.0;
  vxlo = vxhi = vylo = vyhi = vy = vz = 0.0;

  // optional args

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      for (iregion = 0; iregion < domain->nregion; iregion++)
	if (strcmp(arg[iarg+1],domain->regions[iregion]->id) == 0) break;
      if (iregion == domain->nregion) 
	error->all("Fix insert region ID does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"diam") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix insert command");
      radius_lo = 0.5 * atof(arg[iarg+1]);
      radius_hi = 0.5 * atof(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"dens") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix insert command");
      density_lo = atof(arg[iarg+1]);
      density_hi = atof(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"vol") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix insert command");
      volfrac = atof(arg[iarg+1]);
      maxattempt = atoi(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"rate") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      rate = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (force->dimension == 3) {
	if (iarg+6 > narg) error->all("Illegal fix insert command");
	vxlo = atof(arg[iarg+1]);
	vxhi = atof(arg[iarg+2]);
	vylo = atof(arg[iarg+3]);
	vyhi = atof(arg[iarg+4]);
	vz = atof(arg[iarg+5]);
	iarg += 6;
      } else {
	if (iarg+4 > narg) error->all("Illegal fix insert command");
	vxlo = atof(arg[iarg+1]);
	vxhi = atof(arg[iarg+2]);
	vy = atof(arg[iarg+3]);
	vz = 0.0;
	iarg += 4;
      }
    } else error->all("Illegal fix insert command");
  }

  // error check that a valid region was specified

  if (iregion == -1) error->all("Must specify a region in fix insert");

  // error checks on region

  if (domain->regions[iregion]->interior == 0)
    error->all("Must use region with side = in with fix insert");

  if (strcmp(domain->regions[iregion]->style,"block") == 0) {
    region_style = 1;
    xlo = ((RegBlock *) domain->regions[iregion])->xlo;
    xhi = ((RegBlock *) domain->regions[iregion])->xhi;
    ylo = ((RegBlock *) domain->regions[iregion])->ylo;
    yhi = ((RegBlock *) domain->regions[iregion])->yhi;
    zlo = ((RegBlock *) domain->regions[iregion])->zlo;
    zhi = ((RegBlock *) domain->regions[iregion])->zhi;
    if (xlo < domain->boxxlo || xhi > domain->boxxhi || 
	ylo < domain->boxylo || yhi > domain->boxyhi || 
	zlo < domain->boxzlo || zhi > domain->boxzhi)
      error->all("Insertion region extends outside simulation box");
  } else if (strcmp(domain->regions[iregion]->style,"cylinder") == 0) {
    region_style = 2;
    char axis = ((RegCylinder *) domain->regions[iregion])->axis;
    xc = ((RegCylinder *) domain->regions[iregion])->c1;
    yc = ((RegCylinder *) domain->regions[iregion])->c2;
    rc = ((RegCylinder *) domain->regions[iregion])->radius;
    zlo = ((RegCylinder *) domain->regions[iregion])->lo;
    zhi = ((RegCylinder *) domain->regions[iregion])->hi;
    if (axis != 'z')
      error->all("Must use a z-axis cylinder with fix insert");
    if (xc-rc < domain->boxxlo || xc+rc > domain->boxxhi || 
	yc-rc < domain->boxylo || yc+rc > domain->boxyhi || 
	zlo < domain->boxzlo || zhi > domain->boxzhi)
      error->all("Insertion region extends outside simulation box");
  } else error->all("Must use a block or cylinder region with fix insert");

  if (region_style == 2 && force->dimension == 2)
    error->all("Must use a block region with fix insert for 2d simulations");

  // random number generator, same for all procs

  random = new RanPark(seed);

  // allgather arrays

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  recvcounts = new int[nprocs];
  displs = new int[nprocs];

  // nfreq = timesteps between insertions
  // should be time for a particle to fall from top of insertion region
  //   to bottom, taking into account that the region may be moving
  // 1st insertion on next timestep

  double v_relative,delta;
  double g = 1.0;
  if (force->dimension == 3) {
    v_relative = vz - rate;
    delta = v_relative + sqrt(v_relative*v_relative + 2.0*g*(zhi-zlo)) / g;
  } else {
    v_relative = vy - rate;
    delta = v_relative + sqrt(v_relative*v_relative + 2.0*g*(yhi-ylo)) / g;
  }
  nfreq = static_cast<int> (delta/update->dt + 0.5);

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  nfirst = next_reneighbor;
  ninserted = 0;

  // nper = # to insert each time
  // depends on specified volume fraction
  // volume = volume of insertion region
  // volume_one = volume of inserted particle (with max possible radius)
  // in 3d, insure dy >= 1, for quasi-2d simulations

  double volume,volume_one;
  if (force->dimension == 3) {
    if (region_style == 1) {
      double dy = yhi - ylo;
      if (dy < 1.0) dy = 1.0;
      volume = (xhi-xlo) * dy * (zhi-zlo);
    } else volume = PI*rc*rc * (zhi-zlo);
    volume_one = 4.0/3.0 * PI * radius_hi*radius_hi*radius_hi;
  } else {
    volume = (xhi-xlo) * (yhi-ylo);
    volume_one = PI * radius_hi*radius_hi;
  }

  nper = static_cast<int> (volfrac*volume/volume_one);
  int nfinal = update->ntimestep + 1 + (ninsert-1)/nper * nfreq;

  // print stats

  if (me == 0) {
    if (screen)
      fprintf(screen,
	      "Particle insertion: %d every %d steps, %d by step %d\n",
	      nper,nfreq,ninsert,nfinal);
    if (logfile)
      fprintf(logfile,
	      "Particle insertion: %d every %d steps, %d by step %d\n",
	      nper,nfreq,ninsert,nfinal);
  }
}

/* ---------------------------------------------------------------------- */

FixInsert::~FixInsert()
{
  delete random;
  delete [] recvcounts;
  delete [] displs;
}

/* ---------------------------------------------------------------------- */

int FixInsert::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInsert::init()
{
  // insure gravity fix exists
  // for 3d must point in -z, for 2d must point in -y
  // else insertion cannot work

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"gravity") == 0) break;
  if (ifix == modify->nfix) 
    error->all("Must use fix gravity with fix insert");

  double phi = ((FixGravity *) modify->fix[ifix])->phi;
  double theta = ((FixGravity *) modify->fix[ifix])->theta;
  double PI = 2.0 * asin(1.0);
  double degree2rad = 2.0*PI / 360.0;
  double xgrav = sin(degree2rad * theta) * cos(degree2rad * phi);
  double ygrav = sin(degree2rad * theta) * sin(degree2rad * phi);
  double zgrav = cos(degree2rad * theta);

  if (force->dimension == 3) {
    if (fabs(xgrav) > EPSILON || fabs(ygrav) > EPSILON ||
	fabs(zgrav+1.0) > EPSILON)
      error->all("Gravity must point in -z to use with fix insert in 3d");
  } else {
    if (fabs(xgrav) > EPSILON || fabs(ygrav+1.0) > EPSILON ||
	fabs(zgrav) > EPSILON)
      error->all("Gravity must point in -y to use with fix insert in 2d");
  }

  // check if a shear history fix exists

  ifix_history = -1;
  if (force->pair_match("gran/history") || force->pair_match("gran/hertzian"))
    for (int i = 0; i < modify->nfix; i++)
      if (strcmp(modify->fix[i]->style,"SHEAR_HISTORY") == 0) ifix_history = i;
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixInsert::pre_exchange()
{
  int i;

  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  // nnew = # to insert this timestep

  int nnew = nper;
  if (ninserted + nnew > ninsert) nnew = ninsert - ninserted;

  // lo/hi current = z (or y) bounds of insertion region this timestep

  if (force->dimension == 3) {
    lo_current = zlo + (update->ntimestep - nfirst) * update->dt * rate;
    hi_current = zhi + (update->ntimestep - nfirst) * update->dt * rate;
  } else {
    lo_current = ylo + (update->ntimestep - nfirst) * update->dt * rate;
    hi_current = yhi + (update->ntimestep - nfirst) * update->dt * rate;
  }

  // ncount = # of my atoms that overlap the insertion region
  // nprevious = total of ncount across all procs
  
  int ncount = 0;
  for (i = 0; i < atom->nlocal; i++)
    if (overlap(i)) ncount++;

  int nprevious;
  MPI_Allreduce(&ncount,&nprevious,1,MPI_INT,MPI_SUM,world);

  // xmine is for my atoms
  // xnear is for atoms from all procs + atoms to be inserted

  double **xmine = 
    memory->create_2d_double_array(ncount,4,"fix_insert:xmine");
  double **xnear = 
    memory->create_2d_double_array(nprevious+nnew,4,"fix_insert:xnear");
  int nnear = nprevious;

  // setup for allgatherv

  int n = 4*ncount;
  MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];

  // load up xmine array
  
  double **x = atom->x;
  double *radius = atom->radius;

  ncount = 0;
  for (i = 0; i < atom->nlocal; i++)
    if (overlap(i)) {
      xmine[ncount][0] = x[i][0];
      xmine[ncount][1] = x[i][1];
      xmine[ncount][2] = x[i][2];
      xmine[ncount][3] = radius[i];
      ncount++;
    }

  // perform allgatherv to acquire list of nearby particles on all procs

  double *ptr = NULL;
  if (ncount) ptr = xmine[0];
  MPI_Allgatherv(ptr,4*ncount,MPI_DOUBLE,
		 xnear[0],recvcounts,displs,MPI_DOUBLE,world);

  // insert new atoms into xnear list, one by one
  // check against all nearby atoms and previously inserted ones
  // if there is an overlap then try again at same z (3d) or y (2d) coord
  // else insert by adding to xnear list
  // max = maximum # of insertion attempts for all particles
  // h = height, biased to give uniform distribution in time of insertion

  int success;
  double xtmp,ytmp,ztmp,radtmp,delx,dely,delz,rsq,radsum,rn,h;

  int attempt = 0;
  int max = nnew * maxattempt;
  int ntotal = nprevious+nnew;

  while (nnear < ntotal) {
    rn = random->uniform();
    h = hi_current - rn*rn * (hi_current-lo_current);
    radtmp = radius_lo + random->uniform() * (radius_hi-radius_lo);
    success = 0;
    while (attempt < max) {
      attempt++;
      xyz_random(h,xtmp,ytmp,ztmp);
      for (i = 0; i < nnear; i++) {
	delx = xtmp - xnear[i][0];
	dely = ytmp - xnear[i][1];
	delz = ztmp - xnear[i][2];
	rsq = delx*delx + dely*dely + delz*delz;
	radsum = radtmp + xnear[i][3];
	if (rsq <= radsum*radsum) break;
      }
      if (i == nnear) {
	success = 1;
	break;
      }
    }
    if (success) {
      xnear[nnear][0] = xtmp;
      xnear[nnear][1] = ytmp;
      xnear[nnear][2] = ztmp;
      xnear[nnear][3] = radtmp;
      nnear++;
    } else break;
  }

  // warn if not all insertions were performed

  ninserted += nnear-nprevious;
  if (nnear - nprevious < nnew && me == 0)
    error->warning("Less insertions than requested");

  // add new atoms in my sub-box to my arrays
  // initialize info about the atoms
  // type, diameter, density set from fix parameters
  // group mask set to "all" plus fix group
  // z velocity set to what velocity would be if particle
  //   had fallen from top of insertion region
  // this gives continuous stream of atoms
  // set npartner for new atoms to 0 (assume not touching any others)

  int m;
  double denstmp,vxtmp,vytmp,vztmp;
  double g = 1.0;

  for (i = nprevious; i < nnear; i++) {
    xtmp = xnear[i][0];
    ytmp = xnear[i][1];
    ztmp = xnear[i][2];
    radtmp = xnear[i][3];
    denstmp = density_lo + random->uniform() * (density_hi-density_lo);
    if (force->dimension == 3) {
      vxtmp = vxlo + random->uniform() * (vxhi-vxlo);
      vytmp = vylo + random->uniform() * (vyhi-vylo);
      vztmp = vz - sqrt(2.0*g*(hi_current-ztmp));
    } else {
      vxtmp = vxlo + random->uniform() * (vxhi-vxlo);
      vytmp = vy - sqrt(2.0*g*(hi_current-ytmp));
      vztmp = 0.0;
    }

    if (xtmp >= domain->subxlo && xtmp < domain->subxhi &&
	ytmp >= domain->subylo && ytmp < domain->subyhi &&
	ztmp >= domain->subzlo && ztmp < domain->subzhi) {
      atom->create_one(ntype,xtmp,ytmp,ztmp);
      m = atom->nlocal - 1;
      atom->type[m] = ntype;
      atom->radius[m] = radtmp;
      atom->density[m] = denstmp;
      if (force->dimension == 3) 
	atom->rmass[m] = 4.0*PI/3.0 * radtmp*radtmp*radtmp * denstmp;
      else
	atom->rmass[m] = PI * radtmp*radtmp * denstmp;
      atom->mask[m] = 1 | groupbit;
      atom->v[m][0] = vxtmp;
      atom->v[m][1] = vytmp;
      atom->v[m][2] = vztmp;
      if (ifix_history >= 0)
	((FixShearHistory *) modify->fix[ifix_history])->npartner[m] = 0;
    }
  }

  // tag # of new particles grow beyond all previous atoms
  // reset global natoms
  // if global map exists, reset it

  atom->tag_extend();
  atom->natoms += nnear - nprevious;
  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }

  // free local memory

  memory->destroy_2d_double_array(xmine);
  memory->destroy_2d_double_array(xnear);

  // next timestep to insert

  if (ninserted < ninsert) next_reneighbor += nfreq;
  else next_reneighbor = 0;
}

/* ----------------------------------------------------------------------
   check if particle i could overlap with a particle inserted into region
   return 1 if yes, 0 if no
   use maximum diameter for inserted particle
------------------------------------------------------------------------- */

int FixInsert::overlap(int i)
{
  double delta = radius_hi + atom->radius[i];
  double **x = atom->x;

  if (force->dimension == 3) {
    if (region_style == 1) {
      if (x[i][0] < xlo-delta || x[i][0] > xhi+delta ||
	  x[i][1] < ylo-delta || x[i][1] > yhi+delta ||
	  x[i][2] < lo_current-delta || x[i][2] > hi_current+delta) return 0;
    } else {
      if (x[i][2] < lo_current-delta || x[i][2] > hi_current+delta) return 0;
      double delx = x[i][0] - xc;
      double dely = x[i][1] - yc;
      double rsq = delx*delx + dely*dely;
      double r = rc + delta;
      if (rsq > r*r) return 0;
    }
  } else {
      if (x[i][0] < xlo-delta || x[i][0] > xhi+delta ||
	  x[i][1] < lo_current-delta || x[i][1] > hi_current+delta) return 0;
  }

  return 1;
}

/* ---------------------------------------------------------------------- */

void FixInsert::xyz_random(double h, double &x, double &y, double &z)
{
  if (force->dimension == 3) {
    if (region_style == 1) {
      x = xlo + random->uniform() * (xhi-xlo);
      y = ylo + random->uniform() * (yhi-ylo);
      z = h;
    } else {
      double r1,r2;
      while (1) {
	r1 = random->uniform() - 0.5;
	r2 = random->uniform() - 0.5;
	if (r1*r1 + r2*r2 < 0.25) break;
      }
      x = xc + 2.0*r1*rc;
      y = yc + 2.0*r2*rc;
      z = h;
    }
  } else {
    x = xlo + random->uniform() * (xhi-xlo);
    y = h;
    z = 0.0;
  }
}
