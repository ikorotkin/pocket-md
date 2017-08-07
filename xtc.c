#include "version.h"

#ifdef XTC

#include "xtcio.h"                  /* External GROMACS library */
#include "data_structures.h"

t_fileio *xtcfile;
matrix    gmxbox;
rvec     *gmxx;
real      prec;


/*
 *** Init XTC file ***
 */
void xtc_open(char *filename, char *mode, PMD *pmd) {

    xtcfile = open_xtc(filename, mode);

    gmxbox[0][1] = 0; gmxbox[0][2] = 0;
    gmxbox[1][0] = 0; gmxbox[1][2] = 0;
    gmxbox[2][0] = 0; gmxbox[2][1] = 0;

    gmxbox[0][0] = (real)(pmd->box.x);
    gmxbox[1][1] = (real)(pmd->box.y);
    gmxbox[2][2] = (real)(pmd->box.z);

    gmxx = (rvec*)calloc(pmd->N, sizeof(rvec));

    prec = 1000.0;      /* XTC file precision */

}


/*
 *** Write a frame to the XTC file ***
 */
int xtc_write(PMD *pmd) {

    for(int i = 0; i < pmd->N; i++) {
        gmxx[i][0] = (real)(pmd->atom[i].c.x);
        gmxx[i][1] = (real)(pmd->atom[i].c.y);
        gmxx[i][2] = (real)(pmd->atom[i].c.z);
    }

    return write_xtc(xtcfile, pmd->N, (pmd->n + 1), pmd->t, gmxbox, gmxx, prec);

}


/*
 *** Close the XTC file ***
 */
void xtc_close() {

    close_xtc(xtcfile);

    if(gmxx != NULL) free(gmxx);

}

#endif
