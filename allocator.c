#include "data_structures.h"

#define Mb (1./(1024.*1024.))


void allocator(PMD *pmd) {

    pmd->atom = calloc(pmd->N, sizeof(Atom));

    if(pmd->atom == NULL) {
        printf("\nERROR: Out of memory (atom allocator)\n\n");
        exit(3);
    }

    printf("Memory (atoms): %.3f Mb\n", Mb*(double)(pmd->N*sizeof(Atom)));

}


void allocator_cgrid(PMD *pmd) {

    pmd->cgrid_ind_max = calloc(pmd->cgrid_cells, sizeof(int));
    pmd->cgrid_ind_w   = calloc(pmd->cgrid_cells, sizeof(int));

    pmd->cgrid_ind     = malloc(pmd->cgrid_cells*sizeof(int*));

    if((pmd->cgrid_ind_max == NULL) || (pmd->cgrid_ind_w == NULL) || (pmd->cgrid_ind == NULL)) {
        printf("\nERROR: Out of memory (cut-off grid allocator)\n\n");
        exit(3);
    }

    for(int i = 0; i < pmd->cgrid_cells; i++) {
        pmd->cgrid_ind[i] = malloc(pmd->N/pmd->cgrid_cells*cgrid_n_max*sizeof(int));
        if(pmd->cgrid_ind[i] == NULL) {
            printf("\nERROR: Out of memory (cut-off grid index allocator)\n\n");
            exit(3);
        }
    }

    int mem = 2*pmd->cgrid_cells*sizeof(int) + pmd->cgrid_cells*sizeof(int*) + pmd->N/pmd->cgrid_cells*cgrid_n_max*sizeof(int)*pmd->cgrid_cells;
    printf("Memory (cgrid): %.3f Mb\n", Mb*(double)(mem));

}


void allocator_nrlist(PMD *pmd) {

    /* Average number of atoms inside r_cutoff sphere */
    int neighbours = (4./3.*PI*r_cutoff3/pmd->vol)*pmd->N;

    if(neighbours < 1) neighbours = 1;

    double mem = Mb*(double)(neighbours*cgrid_n_max*sizeof(int)*pmd->N);

    if((mem > mem_max) && (pmd->nrlist > 1)) {

        printf("Memory (nlist): %.3f Mb", mem);
        printf(" <-- this is more than defined in params.h file (max_mem = %d Mb)\n", mem_max);
        printf("                Neighbour List Method will be switched OFF\n");
        printf("                Increase mem_max or decrease cgrid_n_max parameters in params.h\n");

        pmd->nrlist = 1;

    } else if(pmd->nrlist > 1) {

        for(int i = 0; i < pmd->N; i++) {
            pmd->atom[i].nrlist = malloc(neighbours*cgrid_n_max*sizeof(int));
            if(pmd->atom[i].nrlist == NULL) {
                printf("\nERROR: Out of memory (cut-off neighbour index allocator)\n\n");
                exit(3);
            }
        }

        printf("Memory (nlist): %.3f Mb\n", mem);

    } else {

        /* Neighbour List Method is OFF */

    }

}
