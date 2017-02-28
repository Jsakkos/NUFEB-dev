/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

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

#include "fix_kinetics_thermo.h"

#include <math.h>
#include <string.h>
#include <cstdio>
#include <string>
#include <sstream>

#include "atom.h"
#include "bio.h"
#include "domain.h"
#include "error.h"
#include "fix_kinetics.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "memory.h"
#include "modify.h"
#include "pointers.h"
#include "variable.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;


/* ---------------------------------------------------------------------- */

FixKineticsThermo::FixKineticsThermo(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Not enough arguments in fix kinetics/thermo command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix kinetics/thermo command");


  fixY = 0;
  if (strcmp(arg[4], "unfix") == 0) fixY = 1;

//  var = new char*[1];
//  ivar = new int[1];
//
//  for (int i = 0; i < 1; i++) {
//    int n = strlen(&arg[4+i][2]) + 1;
//    var[i] = new char[n];
//    strcpy(var[i],&arg[4+i][2]);
//  }
}

/* ---------------------------------------------------------------------- */

FixKineticsThermo::~FixKineticsThermo()
{
//  int i;
//  for (i = 0; i < 1; i++) {
//    delete [] var[i];
//  }
//  delete [] var;
//  delete [] ivar;

  memory->destroy(dG0);
  memory->destroy(khV);
  memory->destroy(liq2Gas);
}

/* ---------------------------------------------------------------------- */

int FixKineticsThermo::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsThermo::init()
{
  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  ntypes = atom->ntypes;

  bio = kinetics->bio;
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;
  ngrids = nx * ny * nz;
  temp = kinetics->temp;
  rth = kinetics->rth;
  DRGCat = kinetics->DRGCat;
  DRGAn = kinetics->DRGAn;

  nnus = bio->nnus;
  catCoeff = bio->catCoeff;
  anabCoeff = bio->anabCoeff;
  nuGCoeff = bio->nuGCoeff;
  typeGCoeff = bio->typeGCoeff;
  diss = bio->dissipation;

  khV = memory->create(khV,nnus+1,"kinetics/thermo:khV");
  liq2Gas = memory->create(liq2Gas,nnus+1,"kinetics/thermo:liq2Gas");
  dG0 = memory->create(dG0,ntypes+1,2,"kinetics/thermo:dG0");

  init_dG0();
  init_KhV();
}

/* ----------------------------------------------------------------------*/

void FixKineticsThermo::init_dG0()
{
  int *ngflag = bio->ngflag;
  int *tgflag = bio->tgflag;

  for (int i = 1; i <= ntypes; i++) {
    dG0[i][0] = 0;
    dG0[i][1] = 0;

    for (int j = 1; j <= nnus; j++) {
      int flag = ngflag[j];
      dG0[i][0] += catCoeff[i][j] * nuGCoeff[j][flag];
      dG0[i][1] += anabCoeff[i][j] * nuGCoeff[j][flag];
    }

    int flag = tgflag[i];
    dG0[i][1] += typeGCoeff[i][flag];
//    printf("dG0[%i][0] = %e \n", i, dG0[i][0]);
//    printf("dG0[%i][1] = %e \n", i, dG0[i][1]);
  }
}

/* ----------------------------------------------------------------------*/

void FixKineticsThermo::init_KhV()
{
  for (int i = 1; i < nnus + 1; i++) {
    khV[i] = 0;
    liq2Gas[i] = 0;
  }

  for (int i = 1; i < nnus + 1; i++) {
    char *nName;     // nutrient name
    nName = bio->nuName[i];

    if (bio->nuType[i] == 1) {
      char *lName = new char[strlen(nName)];     // corresponding liquid
      strncpy(lName, nName+1, strlen(nName));

      for (int j = 1; j < nnus + 1; j++) {
        char *nName2;     // nutrient name
        nName2 = bio->nuName[j];
        if (strcmp(lName, nName2) == 0) {
          liq2Gas[i] = j;
          if (nuGCoeff[j][0] > 10000) {
            khV[j] = exp((nuGCoeff[j][1] - nuGCoeff[i][1]) / (-rth * temp));
          } else {
            khV[j] = exp((nuGCoeff[j][0] - nuGCoeff[i][1]) / (-rth * temp));
          }
          break;
        }
      }
      delete [] lName;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixKineticsThermo::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  thermo();
 // output_data();
}

/* ----------------------------------------------------------------------
  thermodynamics
------------------------------------------------------------------------- */
void FixKineticsThermo::thermo()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int *type = atom->type;
  double **activity = kinetics->activity;
  double *kLa = bio->kLa;
  double **nuGas = kinetics->nuGas;
  double **nuR = kinetics->nuR;

  double vRgT = kinetics->gVol * 1000 / (kinetics->gasTrans * kinetics->temp);

  gYield = kinetics->gYield;
  nuS = kinetics->nuS;

  for (int i = 0; i < ngrids; i++) {
    // gas transfer
    for (int j = 1; j <= nnus; j++) {
      if (bio->nuType[j] == 1) {
        //get corresponding liquid ID
        int liqID = liq2Gas[j];
        double gasT = 0;
        if (liqID != 0) {
          if (nuGCoeff[liqID][0] > 10000) {
            gasT = kLa[liqID] * (activity[liqID][1] / khV[liqID] - activity[j][0]);
          } else {
            gasT = kLa[liqID] * (activity[liqID][0] / khV[liqID] - activity[j][0]);
          }
          nuGas[j][i] = gasT;
          nuGas[liqID][i] = -gasT * vRgT;
          // update nutrient consumption
          nuR[j][i] += nuGas[j][i];
          nuR[liqID][i] += nuGas[liqID][i];
        }
      }
    }

    for (int j = 1; j <= ntypes; j++) {
      double rthT = temp * rth;

      //Gibbs free energy of the reaction
      DRGCat[j][i] = dG0[j][0];  //catabolic energy values
      DRGAn[j][i] = dG0[j][1] + rthT;  //anabolic energy values

      for (int k = 1; k <= nnus; k++) {
        if (nuGCoeff[k][1] < 1e4) {
          double value = 0;
          int flag = bio->ngflag[k];
          if (activity[k][flag] == 0) value = 1e-20;
          else value = activity[k][flag];

          double dgr = rthT * log(value);
         // printf ("%e ",  value);
          DRGCat[j][i] += catCoeff[j][k] * dgr;
          DRGAn[j][i] += anabCoeff[j][k] * dgr;
        } else {
          error->all(FLERR,"nuGCoeff[1] is inf value");
        }
      }
     // printf("\n");
//      printf("DRGAn[%i][%i][0] = %e \n", i,j, DRGAn[j][i]);
//      printf("DRGCat[%i][%i][1] = %e \n", i,j, DRGCat[j][i]);

      //use catabolic and anabolic energy values to derive catabolic reaction equation
      if (fixY == 1){
        if (DRGCat[j][i] < 0) {
          gYield[j][i] = -(DRGAn[j][i] + diss[j]) / DRGCat[j][i] + bio->eD[j];
          gYield[j][i] = 1 / gYield[j][i];
        } else {
          gYield[j][i] = 0;
        }
      }
      //printf("gYield = %e \n", gYield[j][i]);
    }
  }
}

/* ----------------------------------------------------------------------
  output energy to data file
------------------------------------------------------------------------- */

//void FixKineticsThermo::output_data(){
//  std::ostringstream stm;
//  stm << update->ntimestep;
//  string str = "./DGRCat/DGRCat.csv."+ stm.str();
//  FILE *pFile = fopen (str.c_str(), "a");
//  fprintf(pFile, ",x,y,z,scalar,1,1,1,0.5\n");
//  double average = 0.0;
//
//  double xlo,xhi,ylo,yhi,zlo,zhi;
//
//  //Get computational domain size
//  if (domain->triclinic == 0) {
//    xlo = domain->boxlo[0];
//    xhi = domain->boxhi[0];
//    ylo = domain->boxlo[1];
//    yhi = domain->boxhi[1];
//    zlo = domain->boxlo[2];
//    zhi = domain->boxhi[2];
//  }
//
//  else {
//    xlo = domain->boxlo_bound[0];
//    xhi = domain->boxhi_bound[0];
//    ylo = domain->boxlo_bound[1];
//    yhi = domain->boxhi_bound[1];
//    zlo = domain->boxlo_bound[2];
//    zhi = domain->boxhi_bound[2];
//  }
//
//  double stepx = (xhi - xlo) / nx;
//  double stepy = (yhi - ylo) / ny;
//  double stepz = (zhi - zlo) / nz;
//
//  for(int i = 0; i < ngrids; i++){
//    int zpos = i/(nx * ny) + 1;
//    int ypos = (i - (zpos - 1) * (nx * ny)) / nx + 1;
//    int xpos = i - (zpos - 1) * (nx * ny) - (ypos - 1) * nx + 1;
//
//    double x = xpos * stepx - stepx/2;
//    double y = ypos * stepy - stepy/2;
//    double z = zpos * stepz - stepz/2;
//
//    average += DRGCat[2][i];
//
//    fprintf(pFile, "%i,\t%f,\t%f,\t%f,\t%f\n",i, x, y, z, DRGCat[2][i]);
//  }
//  fclose(pFile);
//
//  average = average / ngrids;
//  string str2 = "./DGRCat/Average.csv";
//  pFile = fopen (str2.c_str(), "a");
//  fprintf(pFile, "%e\n", average);
//
//}
