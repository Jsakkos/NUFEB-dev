/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include "atom.h"
#include <string.h>
#include "fix_ave_atom.h"
#include "fix_reaxc_species_kokkos.h"
#include "domain.h"
#include "update.h"
#include "pair_reax_c_kokkos.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "force.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "reaxc_list.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxCSpeciesKokkos::FixReaxCSpeciesKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixReaxCSpecies(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  
  // NOTE: Could improve performance if a Kokkos version of ComputeSpecAtom is added

  datamask_read = X_MASK | V_MASK | Q_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

FixReaxCSpeciesKokkos::~FixReaxCSpeciesKokkos()
{

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpeciesKokkos::init()
{
  Pair* pair_kk = force->pair_match("reax/c/kk",1);
  if (pair_kk == NULL) error->all(FLERR,"Cannot use fix reax/c/species/kk without "
                  "pair_style reax/c/kk");

  FixReaxCSpecies::init();

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->newton = 2;
  neighbor->requests[irequest]->ghost = 1;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpeciesKokkos::FindMolecule()
{
  int i,j,ii,jj,inum,itype,jtype,loop,looptot;
  int change,done,anychange;
  int *mask = atom->mask;
  int *ilist;
  double bo_tmp,bo_cut;
  double **spec_atom = f_SPECBOND->array_atom;

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      clusterID[i] = atom->tag[i];
      x0[i].x = spec_atom[i][1];
      x0[i].y = spec_atom[i][2];
      x0[i].z = spec_atom[i][3];
    }
    else clusterID[i] = 0.0;
  }

  loop = 0;
  while (1) {
    comm->forward_comm_fix(this);
    loop ++;

    change = 0;
    while (1) {
      done = 1;

      for (ii = 0; ii < inum; ii++) {
      	i = ilist[ii];
      	if (!(mask[i] & groupbit)) continue;

      	itype = atom->type[i];

        for (jj = 0; jj < MAXSPECBOND; jj++) {
      	  j = reaxc->tmpid[i][jj];

      	  if (j < i) continue;
      	  if (!(mask[j] & groupbit)) continue;

      	  if (clusterID[i] == clusterID[j] && PBCconnected[i] == PBCconnected[j]
	    && x0[i].x == x0[j].x && x0[i].y == x0[j].y && x0[i].z == x0[j].z) continue;

          jtype = atom->type[j];
      	  bo_cut = BOCut[itype][jtype];
      	  bo_tmp = spec_atom[i][jj+7];

      	  if (bo_tmp > bo_cut) {
            clusterID[i] = clusterID[j] = MIN(clusterID[i], clusterID[j]);
            PBCconnected[i] = PBCconnected[j] = MAX(PBCconnected[i], PBCconnected[j]);
            x0[i] = x0[j] = chAnchor(x0[i], x0[j]);
            if ((fabs(spec_atom[i][1] - spec_atom[j][1]) > reaxc->control->bond_cut)
             || (fabs(spec_atom[i][2] - spec_atom[j][2]) > reaxc->control->bond_cut)
             || (fabs(spec_atom[i][3] - spec_atom[j][3]) > reaxc->control->bond_cut))
              PBCconnected[i] = PBCconnected[j] = 1;
      	    done = 0;
      	  }
      	}
      }
      if (!done) change = 1;
      if (done) break;
    }
    MPI_Allreduce(&change,&anychange,1,MPI_INT,MPI_MAX,world);
    if (!anychange) break;

    MPI_Allreduce(&loop,&looptot,1,MPI_INT,MPI_SUM,world);
    if (looptot >= 400*nprocs) break;

  }
}