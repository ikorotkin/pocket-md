#ifndef DATA_STRUCTURES_H_
#define DATA_STRUCTURES_H_

#include "version.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "params.h"

#ifdef CUDA
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#endif

typedef struct VecI {           /* Integer vector */
    int x, y, z;
} VecI;

typedef struct VecR {           /* Real vector */
    real x, y, z;
} VecR;

typedef struct VecD {           /* Double vector */
    double x, y, z;
} VecD;

typedef struct Atom {           /* Atom's properties: */
    VecR        c;              /* coordinates */
    VecR        vel;            /* velocities */
    VecR        f;              /* forces */
    real        m;              /* mass */
    real        invm;           /* 1/mass */
    int        *nrlist;         /* atom's neighbour list */
    int         nrlist_max;     /* atom's neighbour list length */
} Atom;


typedef struct PMD {
    Atom       *atom;           /* Atom's properties */
    VecD        box;            /* Box size, nm */
    VecD        box05;          /* Box size * 0.5 */
    VecD        cgrid_invh;     /* Cut-off grid inverse step (1/h) */
    VecI        cgrid_max;      /* Cut-off grid size */
    double      vol;            /* Volume of the box, nm^3 */
    double      t;              /* Current time, ps */
    double      dt;             /* Time step, ps */
    double      ttime;          /* Total time, ps */
    double      T;              /* Initial temperature, K */
    double      dens;           /* Density, amu/nm^3 */
    double      Ekin;           /* Total kinetic energy, kJ/mol */
    double      Vmax;           /* Velocity maximum module, nm/ps */
    double      T_curr;         /* Current temperature of the system, K */
    double      tau;            /* Berendsen thermostat tau-parameter, ps */
    double      Epot;           /* Total potential energy, kJ/mol */
    int         N;              /* Number of atoms */
    int         DOF;            /* Number of DOFs of the system */
    int         n;              /* Current time step */
    int         steps;          /* Number of time steps */
    int         threads;        /* Number of threads requested */
    int         xtc_step;       /* XTC-file output frequency */
    int         gro_step;       /* GRO-file output frequency */
    int         berendsen;      /* tau = dt_MD * berendsen, 0 - turn off thermostat */
    int         cgrid_cells;    /* Total number of the cells in the cut-off grid */
    int       **cgrid_ind;      /* Cut-off grid - atom's index [cell #][atom # in the list] */
    int        *cgrid_ind_max;  /* Cut-off grid - list length [cell #] */
    int        *cgrid_ind_w;    /* Cut-off grid - cell has already been processed (1) or not (0) [cell #] */
    int         nrlist;         /* Neighbour list update frequency */
} PMD;

#endif
