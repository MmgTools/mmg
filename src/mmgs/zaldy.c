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
 * \file mmgs/zaldy.c
 * \brief Memory management.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmgs.h"

/* get new point address */
int newPt(MMG5_pMesh mesh,double c[3],double n[3]) {
  MMG5_pPoint  ppt;
  int     curpt;

  if ( !mesh->npnil )  return(0);

  curpt = mesh->npnil;
  if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
  ppt   = &mesh->point[curpt];
  memcpy(ppt->c,c,3*sizeof(double));
  memcpy(ppt->n,n,3*sizeof(double));
  ppt->tag   &= ~MG_NUL;
  mesh->npnil = ppt->tmp;
  ppt->tmp    = 0;

  return(curpt);
}

void delPt(MMG5_pMesh mesh,int ip) {
  MMG5_pPoint   ppt;

  ppt = &mesh->point[ip];
  memset(ppt,0,sizeof(MMG5_Point));
  ppt->tag    = MG_NUL;
  ppt->tmp    = mesh->npnil;
  mesh->npnil = ip;
  if ( ip == mesh->np ) {
    while ( !MG_VOK((&mesh->point[mesh->np])) )  mesh->np--;
  }
}

int newElt(MMG5_pMesh mesh) {
  int     curiel;

  if ( !mesh->ntnil )  return(0);
  curiel = mesh->ntnil;

  if ( mesh->ntnil > mesh->nt )  mesh->nt = mesh->ntnil;
  mesh->ntnil = mesh->tria[curiel].v[2];
  mesh->tria[curiel].v[2] = 0;

  return(curiel);
}

void delElt(MMG5_pMesh mesh,int iel) {
  MMG5_pTria    pt;

  pt = &mesh->tria[iel];
  if ( !MG_EOK(pt) ) {
    fprintf(stdout,"  ## INVALID ELEMENT: %d.\n",iel);
    return;
  }
  memset(pt,0,sizeof(MMG5_Tria));
  pt->v[2] = mesh->ntnil;
  if ( mesh->adja )
    memset(&mesh->adja[3*(iel-1)+1],0,3*sizeof(int));
  mesh->ntnil = iel;
  if ( iel == mesh->nt ) {
    while ( !MG_EOK((&mesh->tria[mesh->nt])) )  mesh->nt--;
  }
}


int zaldy(MMG5_pMesh mesh) {
  int     million = 1048576L;
  int     k,npask,bytes;

  if ( mesh->info.mem < 0 ) {
    mesh->npmax = MG_MAX(1.5*mesh->np,NPMAX);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,NTMAX);
  }
  else {
    /* point+tria+adja */
    bytes = sizeof(MMG5_Point) + 2*sizeof(MMG5_Tria) + 3*sizeof(int);

    npask = (int)((double)mesh->info.mem / bytes * million);
    mesh->npmax = MG_MAX(1.5*mesh->np,npask);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,2*npask);
  }

  mesh->point = (MMG5_pPoint)calloc(mesh->npmax+1,sizeof(MMG5_Point));
  assert(mesh->point);
  mesh->tria  = (MMG5_pTria)calloc(mesh->ntmax+1,sizeof(MMG5_Tria));
  assert(mesh->tria);

  /* store empty links */
  mesh->npnil = mesh->np + 1;
  mesh->ntnil = mesh->nt + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  for (k=mesh->ntnil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;

  return(1);
}
