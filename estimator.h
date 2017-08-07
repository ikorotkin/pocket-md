#ifndef DATA_ESTIMATOR_H_
#define DATA_ESTIMATOR_H_

void   calc_kin_energy(PMD *pmd);           /* Calculate kinetic energy of the system */
double total_energy(PMD *pmd);              /* Estimate total energy of the system */
double temperature(PMD *pmd);               /* Estimate temperature of the system */
double pressure(PMD *pmd);                  /* Estimate pressure of the system */

void   berendsen(PMD *pmd);                 /* Simple Berendsen thermostat */

void   remove_COM_motion(PMD *pmd);         /* COM removal procedure */

#endif
