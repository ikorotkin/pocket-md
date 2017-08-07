#include "data_structures.h"


void write_current_gro(PMD *pmd) {

    char  fname[64];
    FILE *fout;

    sprintf(fname, "conf_%d.gro", pmd->n);
    fout = fopen(fname, "w");

    fprintf(fout, "write_current_gro(%d) :: t= %f\n%5d\n", pmd->n, pmd->t, pmd->N);

    for(int i = 0; i < pmd->N; i++) {

        fprintf(fout, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i+1, "Ar", "Ar", i+1,
                pmd->atom[i].c.x, pmd->atom[i].c.y, pmd->atom[i].c.z,
                pmd->atom[i].vel.x, pmd->atom[i].vel.y, pmd->atom[i].vel.z);

    }

    fprintf(fout, "%10.5f%10.5f%10.5f\n", pmd->box.x, pmd->box.y, pmd->box.z);
    fclose(fout);

}
