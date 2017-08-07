#ifndef PARAMS_H_
#define PARAMS_H_

#define real            double              /* Precision: float/double */

/* Defaults */
#define N_              64000               /* Number of atoms (by default) */
#define dt_             0.01                /* Time step (by default), ps */
#define steps_          100                 /* Number of time steps (by default) */
#define xtc_step_       0                   /* Write data to the XTC-file every 'xtc_step_' steps (by default) */
#define gro_step_       0                   /* Write data to the GRO-file every 'gro_step_' steps (by default) */
#define T_              300.0               /* Initial temperature (by default), K (use 367.42 K if thermostat is off) */
#define dens_           600.0               /* Density (by default), amu/nm^3 */
#define berendsen_      0                   /* Tau = dt_MD * berendsen_ (by default) */
#define nrlist_         1                   /* Update neighbour list every 'nrlist_' time steps (by default), 0 -- do not use neighbour list */

#define r_cutoff        1.0                 /* VdW cut-off radius, nm */

#define mem_max         8192                /* Maximum amount of memory for Neighbour List Method, Mb */
#define cgrid_n_max     10                  /* Assume that local density in any cgrid cell cannot be 'cgrid_n_max' times higher than mean density in the box */
#define cgrid_atom_max  3                   /* Minimum number of atoms per cgrid cell */

#define R               8.31451e-3          /* Gas constant [GROMACS units: nm - ps - a.m.u.] */
#define kB              R                   /* Boltzmann's constant [GROMACS units] */
#define Na              6.0221367e23        /* Avogadro's number, 1/mol */
#define PI              3.1415926535897932  /* PI */

/* LJ pair potential (GROMACS 4.6.7 OPLS-aa 097 AR potential) */
#define MASS_AR         39.948              /* Mass of Argon atom, amu */
#define SIGMA_AR        3.40100e-01         /* L-J Sigma parameter, nm */
#define EPSILON_AR      9.78638e-01         /* L-J Epsilon parameter, kJ/mol */


/* Derived variables */

#define S3              (SIGMA_AR*SIGMA_AR*SIGMA_AR)
#define S6              (S3*S3)
#define S12             (S6*S6)
#define ES6             (24.0*EPSILON_AR*S6)
#define ES12            (48.0*EPSILON_AR*S12)
#define ES6P            (4.0*EPSILON_AR*S6)
#define ES12P           (4.0*EPSILON_AR*S12)

#define r_cutoff2       (r_cutoff*r_cutoff)
#define r_cutoff3       (r_cutoff*r_cutoff*r_cutoff)

#endif
