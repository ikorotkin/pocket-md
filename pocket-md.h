#ifndef POCKET_MD_H_
#define POCKET_MD_H_

void   cmd_parser(int argc, char **argv, PMD *pmd);     /* Read command line string */

void   allocator(PMD *pmd);                             /* Allocate memory */
void   allocator_cgrid(PMD *pmd);                       /* Allocate memory for cut-off grid method */
void   allocator_nrlist(PMD *pmd);                      /* Allocate memory for neighbour list method */

void   initatoms(PMD *pmd);                             /* Set atom masses */
void   initbox(PMD *pmd);                               /* Set box size */
void   initcoords(PMD *pmd);                            /* Set initial coordinates */
void   initvels(PMD *pmd);                              /* Set initial velocities */
void   initcgrid(PMD *pmd);                             /* Initialize cut-off grid */

void   write_current_gro(PMD *pmd);                     /* Write current time step to the GRO-file */

void   do_velocity(real t_step, PMD *pmd);              /* Compute new MD velocities */
void   do_position(real t_step, PMD *pmd);              /* Compute new MD positions */

#ifndef CUDA
void   do_force(PMD *pmd);                              /* Compute MD forces */
void   do_force_list(PMD *pmd);                         /* Compute MD forces using neighbour list */
#else
void   do_force_cuda(PMD *pmd);                         /* Compute MD forces -- CUDA version */
#endif

#endif
