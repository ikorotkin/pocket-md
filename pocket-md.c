#include "data_structures.h"
#include "pocket-md.h"
#include "estimator.h"

#ifdef XTC
#include "xtc.h"
#endif


int main(int argc, char **argv) {

    PMD pmd;

    printf("\n Pocket Molecular Dynamics for ARGON -- Serial/OpenMP VERSION %4.2f\n", VER);
    printf(  "===================================================================\n\n");

    cmd_parser(argc, argv, &pmd);               /* Read the arguments from command line */

#ifdef CUDA
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    printf("GPUs found:     %d\n", nDevices);
#endif

    pmd.n          = 0;                         /* Current time step */
    pmd.t          = 0;                         /* Current time */
    pmd.DOF        = 3*pmd.N;                   /* Number of DOFs */
    pmd.ttime      = pmd.dt*pmd.steps;          /* Total time */
    pmd.threads    = omp_get_max_threads();     /* Get number of threads */
    pmd.tau        = pmd.dt*pmd.berendsen;      /* Tau parameter for Berendsen thermostat */

    printf("OpenMP threads: %d\n", pmd.threads);

    allocator(&pmd);                            /* Allocate memory */
    initatoms(&pmd);                            /* Set atom masses */
    initbox(&pmd);                              /* Set box size */
    initcoords(&pmd);                           /* Set initial coordinates */
    initvels(&pmd);                             /* Set initial velocities */
    remove_COM_motion(&pmd);                    /* Remove initial COM motion */
    initcgrid(&pmd);                            /* Initialize cut-off grid */
    allocator_cgrid(&pmd);                      /* Allocate memory for the cut-off grid method */
    allocator_nrlist(&pmd);                     /* Allocate memory for neighbour list */

    write_current_gro(&pmd);                    /* Write initial time step to the GRO-file */

    if(pmd.nrlist <= 0) pmd.nrlist = 1;

    printf("Update nrlist:  every %d time step(s)\n", pmd.nrlist);
    printf("Thermostat:     tau = %d * dt_MD = %g ps\n", pmd.berendsen, pmd.tau);

#ifdef XTC
    if(pmd.xtc_step) xtc_open("out.xtc", "w", &pmd);      /* Open XTC-file */
#endif

    printf("\n%d atoms, %g ps -- %d steps\n\n", pmd.N, pmd.ttime, pmd.steps);
    printf("%10s %12s %10s %10s %10s %10s %10s %10s\n",
           "Step", "Time, ps", "T, K", "P, MPa", "Ekin", "Epot", "Etot", "Vmax");
    printf("-----------------------------------------------------------------------------------------\n");

    do_position(0.0, &pmd);                     /* Preparing step (with dt=0) for the MD-forces computation */
#ifndef CUDA
    do_force(&pmd);                             /* Calculate MD-forces */
#else
    do_force_cuda(&pmd);                        /* Calculate MD-forces -- CUDA version */
#endif
    do_velocity(pmd.dt*0.5, &pmd);              /* Calculate velocity for the time step dt/2 */

    double tstart = omp_get_wtime();

    /*
     **************************************** MAIN MD LOOP ****************************************
     */
    for(pmd.n = 0; pmd.n <= pmd.steps; pmd.n++) {

        calc_kin_energy(&pmd);                  /* Estimate current kinetic energy and Vmax */
        pmd.T_curr = temperature(&pmd);         /* Estimate current temperature (used in thermostat) */

        /********** Output to console **********/
        if((!(pmd.n % 10) && (pmd.n < 100)) || (!(pmd.n % 100) && (pmd.n < 1000)) || !(pmd.n % 1000)) {
            printf("%10d %12.1f %10.2f %10.2f %10.3e %10.3e %10.3e %10.3f\n",
                   pmd.n, pmd.t, temperature(&pmd), pressure(&pmd)*1.e-6, pmd.Ekin, pmd.Epot, total_energy(&pmd), pmd.Vmax);
            fflush(stdout);
        }

        /********** Output to the XTC-file **********/
#ifdef XTC
        if(pmd.xtc_step)
            if(!(pmd.n % pmd.xtc_step) && pmd.n)
                if(!xtc_write(&pmd))
                    printf("ERROR writing to the XTC-file\n");
#endif

        /********** Output to the GRO-file **********/
        if(pmd.gro_step)
            if(!(pmd.n % pmd.gro_step) && pmd.n)
                write_current_gro(&pmd);

        do_position(pmd.dt, &pmd);              /* Calculate MD-positions */

#ifndef CUDA
        if(!(pmd.n % pmd.nrlist))
            do_force(&pmd);                     /* Calculate MD-forces */
        else
            do_force_list(&pmd);                /* Calculate MD-forces using neighbour list */
#else
        do_force_cuda(&pmd);                    /* Calculate MD-forces -- CUDA version */
#endif

        do_velocity(pmd.dt, &pmd);              /* Calculate MD-velocities */

        if(pmd.berendsen)
            berendsen(&pmd);                    /* Berendsen thermostat */

        pmd.t += pmd.dt;                        /* Time step lapse */

    } /**************************************** MAIN MD LOOP ****************************************/

    double tend = omp_get_wtime();
    printf("-----------------------------------------------------------------------------------------\n\n");
    printf("Computation time:  %.4g s  --  %.4g hours  (%.3f ns/day)\n", (tend-tstart), (tend-tstart)/3600.0, pmd.t*1.e-3/(tend-tstart)*24*3600);

    write_current_gro(&pmd);                    /* Write the last time step to the GRO-file */

#ifdef XTC
    if(pmd.xtc_step) xtc_close();               /* Close XTC-file */
#endif

    printf("\n");
    return 0;

}
