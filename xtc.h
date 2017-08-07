#ifndef XTC_H_
#define XTC_H_

void xtc_open(char *filename, char *mode, PMD *pmd);
void xtc_close();
int  xtc_write(PMD *pmd);

#endif
