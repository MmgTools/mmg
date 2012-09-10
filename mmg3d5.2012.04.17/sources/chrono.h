#ifndef _CHRONO_H
#define _CHRONO_H

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifndef  ON
#define  RESET  0
#define  ON     1
#define  OFF    2
#endif

#define  BIG      1e6
#define  BIG1     1e-6
#define  TIMEMAX  12

typedef struct mytime {
	double  gini,gend,gdif,uini,uend,udif,sini,send,sdif;
  struct  timeval rutim;
	struct  rusage  ru;
  int     call;
} mytime;


/* prototypes */
void   chrono(int cmode,mytime *ptt);
void   tminit(mytime *t,int maxtim);
char  *printim(double );

#endif
