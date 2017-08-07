#include "data_structures.h"


void calc_kin_energy(PMD *pmd) {

    double v2;

    pmd->Ekin = 0;
    pmd->Vmax = 0;

    for(int i = 0; i < pmd->N; i++) {

        v2 = (pmd->atom[i].vel.x*pmd->atom[i].vel.x +
              pmd->atom[i].vel.y*pmd->atom[i].vel.y +
              pmd->atom[i].vel.z*pmd->atom[i].vel.z);

        if(v2 > pmd->Vmax) pmd->Vmax = v2;

        pmd->Ekin += pmd->atom[i].m*v2*0.5;

    }

    pmd->Vmax = sqrt(pmd->Vmax);

}


double total_energy(PMD *pmd) {

    return (pmd->Ekin + pmd->Epot);

}


double temperature(PMD *pmd) {

    return 2.0*pmd->Ekin/((double)(pmd->DOF)*kB);

}


double pressure(PMD *pmd) {

    double P = 2.0*(pmd->Ekin - pmd->Epot)/(3.0*pmd->vol)*1.e30/Na;     /* [Pa] */

    if(P > 0)
        return P;
    else
        return 0;

}


void berendsen(PMD *pmd) {

    double T_target = pmd->T_curr + pmd->dt*(pmd->T - pmd->T_curr)/pmd->tau;
    double Vmean    = sqrt(pmd->T_curr);
    double Vtarget  = sqrt(T_target);

    for(int i = 0; i < pmd->N; i++) {

        pmd->atom[i].vel.x *= Vtarget/Vmean;
        pmd->atom[i].vel.y *= Vtarget/Vmean;
        pmd->atom[i].vel.z *= Vtarget/Vmean;

    }

}


void remove_COM_motion(PMD *pmd) {

    VecD vm = {0, 0, 0};

    for(int i = 0; i < pmd->N; i++) {
        vm.x += pmd->atom[i].vel.x;
        vm.y += pmd->atom[i].vel.y;
        vm.z += pmd->atom[i].vel.z;
    }

    vm.x /= (double)(pmd->N);
    vm.y /= (double)(pmd->N);
    vm.z /= (double)(pmd->N);

    for(int i = 0; i < pmd->N; i++) {
        pmd->atom[i].vel.x -= vm.x;
        pmd->atom[i].vel.y -= vm.y;
        pmd->atom[i].vel.z -= vm.z;
    }

}
