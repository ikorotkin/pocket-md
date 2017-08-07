#include "data_structures.h"


void cmd_parser(int argc, char **argv, PMD *pmd) {

    int pars = 1;

    pmd->N              = N_;
    pmd->dt             = dt_;
    pmd->steps          = steps_;
    pmd->xtc_step       = xtc_step_;
    pmd->gro_step       = gro_step_;
    pmd->T              = T_;
    pmd->dens           = dens_;
    pmd->berendsen      = berendsen_;
    pmd->nrlist         = nrlist_;

    for(int i = 1; i < argc; i++) {

        if(!strcmp(argv[i], "atoms"))       {pmd->N             = atoi(argv[i+1]); pars += 2;}
        if(!strcmp(argv[i], "dt"))          {pmd->dt            = atof(argv[i+1]); pars += 2;}
        if(!strcmp(argv[i], "steps"))       {pmd->steps         = atoi(argv[i+1]); pars += 2;}
        if(!strcmp(argv[i], "xtc"))         {pmd->xtc_step      = atoi(argv[i+1]); pars += 2;}
        if(!strcmp(argv[i], "gro"))         {pmd->gro_step      = atoi(argv[i+1]); pars += 2;}
        if(!strcmp(argv[i], "T"))           {pmd->T             = atof(argv[i+1]); pars += 2;}
        if(!strcmp(argv[i], "dens"))        {pmd->dens          = atof(argv[i+1]); pars += 2;}
        if(!strcmp(argv[i], "berendsen"))   {pmd->berendsen     = atoi(argv[i+1]); pars += 2;}
        if(!strcmp(argv[i], "nrlist"))      {pmd->nrlist        = atoi(argv[i+1]); pars += 2;}

        if(!strcmp(argv[i], "help")) {

            printf("List of parameters (option - description [default value]):\n\n");

            printf("atoms     - Number of atoms [%d]\n", N_);
            printf("dt        - Time step, ps [%g]\n", dt_);
            printf("steps     - Number of time steps [%d]\n", steps_);
            printf("xtc       - Write trajectory to XTC-file every 'xtc' time steps [%d], 0 - do not write XTC\n", xtc_step_);
            printf("gro       - Write trajectory to GRO-file every 'gro' time steps [%d], 0 - do not write GRO\n", gro_step_);
            printf("T         - Initial temperature, K [%g]\n", T_);
            printf("dens      - Density, amu/nm^3 [%g]\n", dens_);
            printf("berendsen - Berendsen thermostat (tau / dt_MD) ratio [%d], 0 - turn off thermostat\n", berendsen_);
            printf("nrlist    - Update neighbour list every 'nrlist' time steps [%d], 0 - do not use neighbour list\n", nrlist_);

            printf("\nExample:\n%s atoms 64000 steps 1000\n\n", argv[0]);

            exit(1);

        }

    }

    printf("Type '%s help' for HELP\n\n", argv[0]);

    if(pars != argc) {
        printf("ERROR: Some of parameters in the command line could not be identified\n\n");
        exit(2);
    }

}
