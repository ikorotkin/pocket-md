#include "data_structures.h"
#include "estimator.h"


void initatoms(PMD *pmd) {

    for(int i = 0; i < pmd->N; i++) {
        pmd->atom[i].m      = MASS_AR;
        pmd->atom[i].invm   = 1.0/MASS_AR;
    }

}


void initbox(PMD *pmd) {

    double mass = 0;

    for(int i = 0; i < pmd->N; i++)
        mass += pmd->atom[i].m;

    pmd->vol   = mass/pmd->dens;
    pmd->box.x = pow(pmd->vol, 1./3.);
    pmd->box.y = pmd->box.x;
    pmd->box.z = pmd->box.x;

    pmd->box05.x = pmd->box.x*0.5;
    pmd->box05.y = pmd->box.y*0.5;
    pmd->box05.z = pmd->box.z*0.5;

    printf("Box size:       %g x %g x %g nm\n", pmd->box.x, pmd->box.y, pmd->box.z);
    printf("Cut-off radius: %g nm\n", r_cutoff);

    if((pmd->box05.x < r_cutoff) || (pmd->box05.y < r_cutoff) || (pmd->box05.z < r_cutoff))
        printf("\nWARNING: Half box size is less than r_cutoff parameter\n\n");

}


void initcoords(PMD *pmd) {

    VecD h;
    int  at = 0;
    int  i_max = (int)(ceil((pow((double)(pmd->N), 1./3.))));
    int  j_max = (int)(ceil(sqrt((double)(pmd->N)/(double)(i_max))));
    int  k_max = (int)(ceil((double)(pmd->N)/(double)(i_max*j_max)));

    printf("Atom's grid:    %d x %d x %d\n", i_max, j_max, k_max);

    h.x = pmd->box.x/(double)(i_max);
    h.y = pmd->box.y/(double)(j_max);
    h.z = pmd->box.z/(double)(k_max);

    for(int k = 0; k < k_max; k++) {
        for(int j = 0; j < j_max; j++) {
            for(int i = 0; i < i_max; i++) {

                if(at == pmd->N) break;

                pmd->atom[at].c.x = h.x*(double)(i);
                pmd->atom[at].c.y = h.y*(double)(j);
                pmd->atom[at].c.z = h.z*(double)(k);

                at++;

            }
        }
    }

}


void initvels(PMD *pmd) {

    VecD   a;
    double Vmean = sqrt(3.0*R*pmd->T/MASS_AR);
    double V;

    for(int i = 0; i < pmd->N; i++) {

        a.x = (double)(rand())/(double)(RAND_MAX)*2.0 - 1.0;
        a.y = (double)(rand())/(double)(RAND_MAX)*2.0 - 1.0;
        a.z = (double)(rand())/(double)(RAND_MAX)*2.0 - 1.0;

        if((a.x*a.x + a.y*a.y + a.z*a.z) == 0.0) continue;

        V = Vmean/sqrt(a.x*a.x + a.y*a.y + a.z*a.z);

        pmd->atom[i].vel.x = a.x*V;
        pmd->atom[i].vel.y = a.y*V;
        pmd->atom[i].vel.z = a.z*V;

    }

    calc_kin_energy(pmd);
    printf("Temperature:    %g K\n", temperature(pmd));

}


void initcgrid(PMD *pmd) {

    pmd->cgrid_max.x = (int)(pmd->box.x/r_cutoff);
    pmd->cgrid_max.y = (int)(pmd->box.y/r_cutoff);
    pmd->cgrid_max.z = (int)(pmd->box.z/r_cutoff);

    pmd->cgrid_cells = pmd->cgrid_max.x*pmd->cgrid_max.y*pmd->cgrid_max.z;

    /* Restriction for gases -- number of atoms per cell should not be too small */
    if(pmd->N/pmd->cgrid_cells < cgrid_atom_max) {

        double coef = pow((double)(pmd->N)/(double)(pmd->cgrid_cells)/(double)(cgrid_atom_max), 1./3.);

        pmd->cgrid_max.x = (int)((double)(pmd->cgrid_max.x)*coef);
        pmd->cgrid_max.y = (int)((double)(pmd->cgrid_max.y)*coef);
        pmd->cgrid_max.z = (int)((double)(pmd->cgrid_max.z)*coef);

        if(pmd->cgrid_max.x < 3) pmd->cgrid_max.x = 3;
        if(pmd->cgrid_max.y < 3) pmd->cgrid_max.y = 3;
        if(pmd->cgrid_max.z < 3) pmd->cgrid_max.z = 3;

        pmd->cgrid_cells = pmd->cgrid_max.x*pmd->cgrid_max.y*pmd->cgrid_max.z;

    }

    printf("Cut-off grid:   %d x %d x %d\n", pmd->cgrid_max.x, pmd->cgrid_max.y, pmd->cgrid_max.z);

    if((pmd->cgrid_max.x < 3) || (pmd->cgrid_max.y < 3) || (pmd->cgrid_max.z < 3)) {
        printf("\nERROR: Too large cut-off radius for this box size (increase number of atoms or reduce cut-off)\n\n");
        exit(5);
    }

    for(int i = 0; i < pmd->cgrid_cells; i++) {
        pmd->cgrid_invh.x = (double)(pmd->cgrid_max.x)/pmd->box.x;
        pmd->cgrid_invh.y = (double)(pmd->cgrid_max.y)/pmd->box.y;
        pmd->cgrid_invh.z = (double)(pmd->cgrid_max.z)/pmd->box.z;
    }

}
