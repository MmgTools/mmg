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
#include "mmg2d.h"

static int MMG_inxtt[5] = {0,1,2,0,1};


/* compute iso size map */
int MMG2_doSol(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria      ptt;
  MMG5_pPoint     p1,p2;
  double     ux,uy,uz,dd;
  int        i,k,ib,ipa,ipb;
  
  sol->np = mesh->np;  
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !ptt->v[0] )  continue;

    for (i=0; i<3; i++) {
      ib  = MMG_inxtt[i+1];
      ipa = ptt->v[i];
      ipb = ptt->v[ib];
      p1  = &mesh->point[ipa];
      p2  = &mesh->point[ipb];

      ux  = p1->c[0] - p2->c[0];
      uy  = p1->c[1] - p2->c[1];
      dd  = sqrt(ux*ux + uy*uy + uz*uz);

      sol->m[ipa] += dd;
      p1->tagdel++;
      sol->m[ipb] += dd;
      p2->tagdel++;
    }
  }

  /* vertex size */
  mesh->info.hmin = 1.e20;
  mesh->info.hmax = 0.0;
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !p1->tagdel )  continue;

    sol->m[k] = sol->m[k] / (double)p1->tagdel;
 
    if ( sol->m[k] < mesh->info.hmin )
      mesh->info.hmin = sol->m[k];
    else if ( sol->m[k] > mesh->info.hmax )
      mesh->info.hmax = sol->m[k];

    p1->tagdel = 0;
  }

  if ( mesh->info.imprim < -4 )
    fprintf(stdout,"     HMIN %f   HMAX %f\n",mesh->info.hmin,mesh->info.hmax);
  return(1);
}

