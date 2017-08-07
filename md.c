#include "data_structures.h"


void do_velocity(double t_step, PMD *pmd) {

    double dtm;

    for(int i = 0; i < pmd->N; i++) {

        dtm = t_step*pmd->atom[i].invm;

        pmd->atom[i].vel.x += dtm*pmd->atom[i].f.x;
        pmd->atom[i].vel.y += dtm*pmd->atom[i].f.y;
        pmd->atom[i].vel.z += dtm*pmd->atom[i].f.z;

    }

}


void do_position(double t_step, PMD *pmd) {

    int  pt;
    VecI ind;

    const int ind_max = pmd->N/pmd->cgrid_cells*cgrid_n_max;

    for(int i = 0; i < pmd->cgrid_cells; i++) {
        pmd->cgrid_ind_max[i] = 0;
        pmd->cgrid_ind_w[i]   = 0;
    }

    for(int i = 0; i < pmd->N; i++) {

        pmd->atom[i].c.x += t_step*pmd->atom[i].vel.x;
        pmd->atom[i].c.y += t_step*pmd->atom[i].vel.y;
        pmd->atom[i].c.z += t_step*pmd->atom[i].vel.z;

        /* PBC */
        while(pmd->atom[i].c.x < 0)           pmd->atom[i].c.x += pmd->box.x;
        while(pmd->atom[i].c.y < 0)           pmd->atom[i].c.y += pmd->box.y;
        while(pmd->atom[i].c.z < 0)           pmd->atom[i].c.z += pmd->box.z;
        while(pmd->atom[i].c.x >= pmd->box.x) pmd->atom[i].c.x -= pmd->box.x;
        while(pmd->atom[i].c.y >= pmd->box.y) pmd->atom[i].c.y -= pmd->box.y;
        while(pmd->atom[i].c.z >= pmd->box.z) pmd->atom[i].c.z -= pmd->box.z;

        /* Find atom's position in the cut-off grid method */
        ind.x = (int)(pmd->atom[i].c.x*pmd->cgrid_invh.x);
        ind.y = (int)(pmd->atom[i].c.y*pmd->cgrid_invh.y);
        ind.z = (int)(pmd->atom[i].c.z*pmd->cgrid_invh.z);

        /* Check indexes -- this error should never happen because of PBC above */
        if((ind.x < 0) || (ind.x >= pmd->cgrid_max.x) || (ind.y < 0) || (ind.y >= pmd->cgrid_max.y) || (ind.z < 0) || (ind.z >= pmd->cgrid_max.z)) {
            printf("\nERROR: Wrong indexes for the atom #%d: (%d, %d, %d) of (%d, %d, %d)\n\n",
                   i, ind.x, ind.y, ind.z, pmd->cgrid_max.x, pmd->cgrid_max.y, pmd->cgrid_max.z);
            exit(11);
        }

        /* Create neighbour list */
        pt = ind.x + ind.y*pmd->cgrid_max.x + ind.z*pmd->cgrid_max.x*pmd->cgrid_max.y;
        pmd->cgrid_ind[pt][pmd->cgrid_ind_max[pt]] = i;
        pmd->cgrid_ind_max[pt] += 1;

        /* cgrid_ind_max cannot be greater than ind_max -- if this happen, there are some regions in the box with extremely high density */
        if(pmd->cgrid_ind_max[pt] == ind_max) {
            printf("\nERROR: Too many atoms inside the cell (%d, %d, %d) in the cut-off grid method\n\n", ind.x, ind.y, ind.z);
            exit(12);
        }

        /* Reset MD-forces */
        pmd->atom[i].f.x = 0;
        pmd->atom[i].f.y = 0;
        pmd->atom[i].f.z = 0;

        /* Reset neighbour list if necessary */
        if(!(pmd->n % pmd->nrlist))
            pmd->atom[i].nrlist_max = 0;

    }

}


void do_force(PMD *pmd) {

    int  nth = pmd->threads;
    int  i, j, ji, indi, indj, pti, ptj;
    VecI ci, cj, cjn;
    VecR d, f;
    real r2, invr2, invr4, invr8, invr14, Fr;
    real Epot = 0;

#pragma omp parallel for num_threads(nth) schedule(static) private(i,j,ji,indi,indj,pti,ptj,ci,cj,cjn,d,f,r2,invr2,invr4,invr8,invr14,Fr) reduction(+:Epot)
    for(int th = 0; th < nth; th++) {

        const int cgrid_max_xy = pmd->cgrid_max.x*pmd->cgrid_max.y;

        int start_th, end_th;

        start_th = (pmd->cgrid_max.x*th)/nth;
        end_th   = (pmd->cgrid_max.x*(th + 1))/nth;

        for(ci.x = start_th; ci.x < end_th; ci.x++) {
            for(ci.y = 0; ci.y < pmd->cgrid_max.y; ci.y++) {
                for(ci.z = 0; ci.z < pmd->cgrid_max.z; ci.z++) {

                    pti = ci.x + ci.y*pmd->cgrid_max.x + ci.z*cgrid_max_xy;

                    if(nth == 1) pmd->cgrid_ind_w[pti] = 1;

                    for(cj.x = ci.x-1; cj.x <= ci.x+1; cj.x++) {

                        cjn.x = cj.x;
                        if(cjn.x < 0) cjn.x = pmd->cgrid_max.x - 1;
                        if(cjn.x > pmd->cgrid_max.x - 1) cjn.x = 0;

                        for(cj.y = ci.y-1; cj.y <= ci.y+1; cj.y++) {

                            cjn.y = cj.y;
                            if(cjn.y < 0) cjn.y = pmd->cgrid_max.y - 1;
                            if(cjn.y > pmd->cgrid_max.y - 1) cjn.y = 0;

                            for(cj.z = ci.z-1; cj.z <= ci.z+1; cj.z++) {

                                cjn.z = cj.z;
                                if(cjn.z < 0) cjn.z = pmd->cgrid_max.z - 1;
                                if(cjn.z > pmd->cgrid_max.z - 1) cjn.z = 0;

                                ptj = cjn.x + cjn.y*pmd->cgrid_max.x + cjn.z*cgrid_max_xy;

                                if((pmd->cgrid_ind_w[ptj] == 1) && (pti != ptj)) continue;

                                for(i = 0; i < pmd->cgrid_ind_max[pti]; i++) {

                                    indi = pmd->cgrid_ind[pti][i];

                                    if((pti == ptj) && (nth == 1))
                                        ji = i+1;
                                    else
                                        ji = 0;

                                    for(j = ji; j < pmd->cgrid_ind_max[ptj]; j++) {

                                        indj = pmd->cgrid_ind[ptj][j];

                                        if((nth != 1) && (indi == indj)) continue;

                                        d.x = (pmd->atom[indi].c.x - pmd->atom[indj].c.x);
                                        d.y = (pmd->atom[indi].c.y - pmd->atom[indj].c.y);
                                        d.z = (pmd->atom[indi].c.z - pmd->atom[indj].c.z);

                                        /* PBC */
                                        if(d.x >  pmd->box05.x) d.x -= pmd->box.x;
                                        if(d.y >  pmd->box05.y) d.y -= pmd->box.y;
                                        if(d.z >  pmd->box05.z) d.z -= pmd->box.z;
                                        if(d.x < -pmd->box05.x) d.x += pmd->box.x;
                                        if(d.y < -pmd->box05.y) d.y += pmd->box.y;
                                        if(d.z < -pmd->box05.z) d.z += pmd->box.z;

                                        r2 = (d.x*d.x + d.y*d.y + d.z*d.z);

                                        if(r2 > r_cutoff2) continue;

                                        invr2  = 1.0/r2;
                                        invr4  = invr2*invr2;
                                        invr8  = invr4*invr4;
                                        invr14 = invr8*invr4*invr2;

                                        Fr = (ES12*invr14 - ES6*invr8);

                                        Epot += (ES12P*invr4*invr8 - ES6P*invr2*invr4);

                                        f.x = Fr*d.x;
                                        f.y = Fr*d.y;
                                        f.z = Fr*d.z;

                                        pmd->atom[indi].f.x += f.x;
                                        pmd->atom[indi].f.y += f.y;
                                        pmd->atom[indi].f.z += f.z;

                                        if(pmd->nrlist > 1) {
                                            pmd->atom[indi].nrlist[pmd->atom[indi].nrlist_max] = indj;  // Can lead to crash if too many neighbours
                                            pmd->atom[indi].nrlist_max++;
                                        }

                                        if(nth == 1) {
                                            pmd->atom[indj].f.x -= f.x;
                                            pmd->atom[indj].f.y -= f.y;
                                            pmd->atom[indj].f.z -= f.z;
                                        }

                                    }
                                }

                            }
                        }
                    }

                }
            }
        }

    }

    if(nth == 1)
        pmd->Epot = Epot;
    else
        pmd->Epot = Epot*0.5;

}


void do_force_list(PMD *pmd) {

    int  nth = pmd->threads;
    int  jj;
    VecR d, f;
    real r2, invr2, invr4, invr8, invr14, Fr;

#pragma omp parallel for num_threads(nth) schedule(static) private(jj,d,f,r2,invr2,invr4,invr8,invr14,Fr)
    for(int th = 0; th < nth; th++) {

        int start_th, end_th;

        start_th = (pmd->N*th)/nth;
        end_th   = (pmd->N*(th + 1))/nth;

        for(int i = start_th; i < end_th; i++) {
            for(int j = 0; j < pmd->atom[i].nrlist_max; j++) {

                jj = pmd->atom[i].nrlist[j];

                d.x = (pmd->atom[i].c.x - pmd->atom[jj].c.x);
                d.y = (pmd->atom[i].c.y - pmd->atom[jj].c.y);
                d.z = (pmd->atom[i].c.z - pmd->atom[jj].c.z);

                /* PBC */
                if(d.x >  pmd->box05.x) d.x -= pmd->box.x;
                if(d.y >  pmd->box05.y) d.y -= pmd->box.y;
                if(d.z >  pmd->box05.z) d.z -= pmd->box.z;
                if(d.x < -pmd->box05.x) d.x += pmd->box.x;
                if(d.y < -pmd->box05.y) d.y += pmd->box.y;
                if(d.z < -pmd->box05.z) d.z += pmd->box.z;

                r2 = (d.x*d.x + d.y*d.y + d.z*d.z);

                invr2  = 1.0/r2;
                invr4  = invr2*invr2;
                invr8  = invr4*invr4;
                invr14 = invr8*invr4*invr2;

                Fr = (ES12*invr14 - ES6*invr8);

                f.x = Fr*d.x;
                f.y = Fr*d.y;
                f.z = Fr*d.z;

                pmd->atom[i].f.x += f.x;
                pmd->atom[i].f.y += f.y;
                pmd->atom[i].f.z += f.z;

                if(nth == 1) {
                    pmd->atom[jj].f.x -= f.x;
                    pmd->atom[jj].f.y -= f.y;
                    pmd->atom[jj].f.z -= f.z;
                }

            }
        }

    }

}

