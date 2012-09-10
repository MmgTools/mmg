/* simulation of a chronograph
 * modified (08/2010) for // usage */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "chrono.h"


void  chrono(int cmode,mytime *ptt) {

  if ( cmode == RESET ) {
    ptt->call = 0;
		ptt->gini = ptt->gend = ptt->gdif = 0.0;
		ptt->sini = ptt->send = ptt->sdif = 0.0;
		ptt->uini = ptt->uend = ptt->udif = 0.0;
  }

  else {
	  gettimeofday(&(ptt->rutim), NULL);
 
    if ( cmode == ON ) {
			ptt->gini  = (double)((ptt->rutim.tv_sec ) + (ptt->rutim.tv_usec) * BIG1);

			getrusage(RUSAGE_SELF,&(ptt->ru));
			ptt->rutim = ptt->ru.ru_utime;
      ptt->uini  = (double)((ptt->rutim.tv_sec) * BIG + (ptt->rutim.tv_usec));
      ptt->rutim = ptt->ru.ru_stime;
      ptt->sini  = (double)((ptt->rutim.tv_sec )* BIG + (ptt->rutim.tv_usec));
    }
    else if ( cmode == OFF ) {
      ptt->gend  = (double)((ptt->rutim.tv_sec ) + (ptt->rutim.tv_usec) * BIG1);

      getrusage(RUSAGE_SELF,&(ptt->ru));
      ptt->rutim = ptt->ru.ru_utime;
      ptt->uend  = (double)((ptt->rutim.tv_sec ) * BIG + (ptt->rutim.tv_usec));
      ptt->rutim = ptt->ru.ru_stime;
      ptt->send  = (double)((ptt->rutim.tv_sec ) * BIG + (ptt->rutim.tv_usec));

      ptt->gdif += ptt->gend - ptt->gini;
      
      ptt->udif += (ptt->uend - ptt->uini) * BIG1;
      ptt->sdif += (ptt->send - ptt->sini) * BIG1;
    
      ptt->call++;
    }
  }
}


/* initialize time table */
void  tminit(mytime *t,int maxtim) {
	mytime  *ptt;
  int      k;

  for (k=0; k<maxtim; k++) {
		ptt = &t[k];
    ptt->call = 0;
		ptt->gini = ptt->gend = ptt->gdif = 0.0;
		ptt->sini = ptt->send = ptt->sdif = 0.0;
		ptt->uini = ptt->uend = ptt->udif = 0.0;
  }
}


/* print real time */
char *printim(double elps) {
	int    hh,mm,ss;
	char  *data;
	
	data = malloc(32*sizeof(char));
	assert(data);
  if ( elps < 60.0 )
    sprintf(data,"%5.3lfs",elps);
  else if ( elps < 3600.0 ) {
    mm = elps / 60.0;
    ss = (int)elps - mm * 60;
    sprintf(data,"%dm%ds (%7.3lfs)",mm,ss,elps);
  }
  else {
    hh = elps / 3600;
    mm = (elps - hh*3600) / 60;
    ss = elps - mm*60 - hh*3600;
    sprintf(data,"%dh%dm%ds",hh,mm,ss);
  }

	return(data);
}


