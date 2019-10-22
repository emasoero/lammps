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

#include "min_maske.h"
#include <mpi.h>
#include <cmath>
#include "universe.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "output.h"
#include "timer.h"

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

#define DELAYSTEP 5

/* ---------------------------------------------------------------------- */

MinMaske::MinMaske(LAMMPS *lmp) : Min(lmp) {}

/* ---------------------------------------------------------------------- */

void MinMaske::init()
{
  Min::init();

  dt = update->dt;
  last_negative = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void MinMaske::setup_style()
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    v[i][0] = v[i][1] = v[i][2] = 0.0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinMaske::reset_vectors()
{
  // atomic dof

  nvec = 3 * atom->nlocal;
  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
}

/* ----------------------------------------------------------------------
   minimization via Maske damped dynamics
------------------------------------------------------------------------- */

int MinMaske::iterate(int maxiter)
{
  bigint ntimestep;
  double vmax,vdotf,vdotfall,fdotf,fdotfloc,fdotfall,scale;
  double dtvone,dtv,dtf,dtfm;
  int flag,flagall;

  alpha_final = 0.0;

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // zero velocity if anti-parallel to force
    // else project velocity in direction of force

    double **v = atom->v;
    double **f = atom->f;
    int nlocal = atom->nlocal;

    vdotf = 0.0;
    for (int i = 0; i < nlocal; i++)
      vdotf += v[i][0]*f[i][0] + v[i][1]*f[i][1] + v[i][2]*f[i][2];
    MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);

    // Euler integration step

    double **x = atom->x;

    for (int i = 0; i < nlocal; i++) {
      if(f[i][0] > 0.0 && v[i][0] < hcx-dmax) {
	x[i][0] += dmax;
	v[i][0] += dmax;
      }
      if(f[i][1] > 0.0 && v[i][1] < hcy-dmax) {
	x[i][1] += dmax;
	v[i][1] += dmax;
      }
      if(f[i][2] > 0.0 && v[i][2] < hcz-dmax) {
	x[i][2] += dmax;
	v[i][2] += dmax;
      }
                  
      if(f[i][0] < 0.0 && v[i][0] > -hcx+dmax) {
	x[i][0] -= dmax;
	v[i][0] -= dmax;
      }
      if(f[i][1] < 0.0 && v[i][1] > -hcy+dmax) {
	x[i][1] -= dmax;
	v[i][1] -= dmax;
      }
      if(f[i][2] < 0.0 && v[i][2] > -hcz+dmax) {
	x[i][2] -= dmax;
	v[i][2] -= dmax;
      }
    }

    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;

    // energy tolerance criterion
    // only check after DELAYSTEP elapsed since velocties reset to 0
    // sync across replicas if running multi-replica minimization

    if (update->etol > 0.0 && ntimestep-last_negative > DELAYSTEP) {
      if (update->multireplica == 0) {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          return ETOL;
      } else {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) return ETOL;
      }
    }

    // force tolerance criterion
    // sync across replicas if running multi-replica minimization

    if (update->ftol > 0.0) {
      if (normstyle == MAX) fdotfloc = fnorm_max();		// max force norm
      else if (normstyle == INF) fdotfloc = fnorm_inf();	// inf force norm
      else if (normstyle == TWO) fdotfloc = fnorm_sqr();	// Euclidean force 2-norm
      else error->all(FLERR,"Illegal min_modify command");
      if (update->multireplica == 0) {
        if (fdotf < update->ftol*update->ftol) return FTOL;
      } else {
        if (fdotf < update->ftol*update->ftol) flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) return FTOL;
      }
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}
