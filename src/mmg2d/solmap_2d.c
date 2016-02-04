/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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
 * \file mmg2d/solmap_2d.c
 * \brief  Compute isotropic size map according to the mean of the length of the edges
 * passing through a point.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg2d.h"



/* compute iso size map */
int MMG2_doSol(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria      ptt,pt;
  MMG5_pPoint     p1,p2;
  double          ux,uy,dd;
  int             i,k,ib,ipa,ipb;
  int             MMG_inxtt[5] = {0,1,2,0,1};

  sol->np = mesh->np;
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    p1->tagdel = 0;
  }
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
      dd  = sqrt(ux*ux + uy*uy);

      sol->m[ipa] += dd;
      p1->tagdel++;
      sol->m[ipb] += dd;
      p2->tagdel++;
    }
  }

  /* vertex size */
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !p1->tagdel )  {
      sol->m[k] = mesh->info.hmax;
      continue;
    }

    sol->m[k] = MG_MIN(mesh->info.hmax,MG_MAX(mesh->info.hmin,sol->m[k] / (double)p1->tagdel));
    p1->tagdel = 0;
  }

/* compute quality */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    pt->qual = MMG2_caltri_in(mesh,sol,pt);
  }

  if ( mesh->info.imprim < -4 )
    fprintf(stdout,"     HMIN %f   HMAX %f\n",mesh->info.hmin,mesh->info.hmax);
  return(1);
}

