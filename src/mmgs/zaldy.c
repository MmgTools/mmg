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
int _MMG5_newPt(MMG5_pMesh mesh,double c[3],double n[3]) {
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

void _MMG5_delPt(MMG5_pMesh mesh,int ip) {
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

int _MMG5_newElt(MMG5_pMesh mesh) {
  int     curiel;

  if ( !mesh->nenil )  return(0);
  curiel = mesh->nenil;

  if ( mesh->nenil > mesh->nt )  mesh->nt = mesh->nenil;
  mesh->nenil = mesh->tria[curiel].v[2];
  mesh->tria[curiel].v[2] = 0;

  return(curiel);
}


void _MMG5_delElt(MMG5_pMesh mesh,int iel) {
  MMG5_pTria    pt;

  pt = &mesh->tria[iel];
  if ( !MG_EOK(pt) ) {
    fprintf(stdout,"  ## INVALID ELEMENT %d.\n",iel);
    exit(EXIT_FAILURE);
  }
  memset(pt,0,sizeof(MMG5_Tria));
  pt->v[2] = mesh->nenil;
  if ( mesh->adja )
    memset(&mesh->adja[3*(iel-1)+1],0,3*sizeof(int));
  mesh->nenil = iel;
  if ( iel == mesh->nt ) {
    while ( !MG_EOK((&mesh->tria[mesh->nt])) )  mesh->nt--;
  }
}

/** memory repartition for the -m option */
void _MMG5_memOption(MMG5_pMesh mesh) {
  long long  million = 1048576L;
  int        npask,bytes,memtmp;

  mesh->memMax = _MMG5_memSize();

  mesh->npmax = MG_MAX(1.5*mesh->np,_MMG5_NPMAX);
  mesh->ntmax = MG_MAX(1.5*mesh->nt,_MMG5_NTMAX);

  if ( mesh->info.mem <= 0 ) {
    if ( mesh->memMax )
      /* maximal memory = 50% of total physical memory */
      mesh->memMax = mesh->memMax*50/100;
    else {
      /* default value = 800 Mo */
      printf("  Maximum memory set to default value: %d Mo.\n",_MMG5_MEMMAX);
      mesh->memMax = _MMG5_MEMMAX*million;
    }
  }
  else {
    /* memory asked by user if possible, otherwise total physical memory */
    if ( (long long)mesh->info.mem*million > mesh->memMax && mesh->memMax ) {
      fprintf(stdout,"  ## Warning: asking for %d Mo of memory ",mesh->info.mem);
      fprintf(stdout,"when only %lld available.\n",mesh->memMax/million);
    }
    else {
      mesh->memMax= (long long)(mesh->info.mem)*million;
    }

    /* if asked memory is lower than default _MMG5_NPMAX/_MMG5_NTMAX we take lower values */

    /* point+tria+adja */
    bytes = sizeof(MMG5_Point) + sizeof(MMG5_xPoint) +
      2*sizeof(MMG5_Tria) + 3*sizeof(int) + sizeof(MMG5_Sol);

    /*init allocation need 38Mo*/
    npask = (int)((double)(mesh->info.mem-38) / bytes) * (int)million;
    mesh->npmax = MG_MIN(npask,mesh->npmax);
    mesh->ntmax = MG_MIN(2*npask,mesh->ntmax);

    /*check if the memory asked is enough to load the mesh*/
    if(mesh->np &&
       (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt )) {
      memtmp = mesh->np * bytes /(int)million + 38;
      memtmp = MG_MAX(memtmp, mesh->nt * bytes /(2 * (int)million) + 38);
      mesh->memMax = (long long) memtmp+1;
      fprintf(stdout,"  ## ERROR: asking for %d Mo of memory ",mesh->info.mem);
      fprintf(stdout,"is not enough to load mesh. You need to ask %d Mo minimum\n",
              memtmp+1);
    }
    if(mesh->info.mem < 39) {
      mesh->memMax = (long long) 39;
      fprintf(stdout,"  ## ERROR: asking for %d Mo of memory ",mesh->info.mem);
      fprintf(stdout,"is not enough to load mesh. You need to ask %d Mo minimum\n",
              39);
    }
  }

  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED (Mo)    %lld\n",
            mesh->memMax/million);

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  _MMG5_NPMAX    %d\n",mesh->npmax);
    fprintf(stdout,"  _MMG5_NTMAX    %d\n",mesh->ntmax);
  }

  return;
}

int zaldy(MMG5_pMesh mesh) {
  int     k;

  _MMG5_memOption(mesh);

  _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point);
  _MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),"initial triangles",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->tria,mesh->ntmax+1,MMG5_Tria);

  /* store empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->nt + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  for (k=mesh->nenil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;

  return(1);
}
