/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file common/chrono.c
 * \brief Simulation of a chronograph.
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \version 5
 * \date  08 2010
 * \copyright GNU Lesser General Public License.
 *
 * Simulation of a chronograph. Allow parallel usage.
 *
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "chrono.h"

/**
 * \fn void  chrono(int cmode,mytime *ptt)
 * \brief Function to measure time.
 * \param cmode macro that allow to reset (RESET), start (ON) or stop (OFF)
 * the chrono.
 * \param *ptt pointer toward mytime object that store the chronograph informations.
 */
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


/**
 * \fn void  tminit(mytime *t,int maxtim)
 * \brief Initialize mytime object.
 * \param *t mytime object to store the chronograph infos.
 * \param maxtim integer sepcifying the maximum number of times stored.
 */
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


/**
 * \fn void  printim(double elps,char *stim)
 * \brief Print real time.
 * \param elps elapsed time in seconds.
 * \param *stim pointer toward string containg the elapsed time at .h.m.s format.
 */
void printim(double elps,char *stim) {
  int    hh,mm,ss;

  if ( elps < 60.0 )
    sprintf(stim,"%5.3lfs",elps);
  else if ( elps < 3600.0 ) {
    mm = elps / 60.0;
    ss = (int)elps - mm * 60;
    sprintf(stim,"%dm%ds (%7.3lfs)",mm,ss,elps);
  }
  else {
    hh = elps / 3600;
    mm = (elps - hh*3600) / 60;
    ss = elps - mm*60 - hh*3600;
    sprintf(stim,"%dh%dm%ds",hh,mm,ss);
  }
}
