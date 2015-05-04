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
 * \file mmgs/hash.c
 * \brief Functions for hash tables management and triangle packing.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

/* tria packing */
static void paktri(MMG5_pMesh mesh) {
  MMG5_pTria   pt,pt1;
  int     k;

  k = 1;
  do {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) {
      pt1 = &mesh->tria[mesh->nt];
      memcpy(pt,pt1,sizeof(MMG5_Tria));
      _MMG5_delElt(mesh,mesh->nt);
    }
  }
  while ( ++k < mesh->nt );

  /* Recreate nil chain */
  mesh->nenil = mesh->nt + 1;

  for(k=mesh->nenil; k<=mesh->ntmax-1; k++){
    mesh->tria[k].v[2] = k+1;
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 if fail.
 *
 * Create adjacency table.
 *
 */
int _MMG5_hashTria(MMG5_pMesh mesh) {
  MMG5_pTria          pt;
  int                 *adja,k,i,ier;

  if ( mesh->adja )  return(1);
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING STRUCTURE\n");

  /* tassage */
  paktri(mesh);

  _MMG5_ADD_MEM(mesh,(3*mesh->ntmax+5)*sizeof(int),"adjacency table",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,int);

  ier = _MMG5_mmgHashTria(mesh, mesh->adja, 0);

//warning the following loop seems to be unused (or only in boulchknm)... to check.
  /* set seed */
  for (k=1; k<=mesh->nt; k++) {
    pt   = &mesh->tria[k];
    adja = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( !adja[i] ) {
        mesh->point[pt->v[_MMG5_inxt2[i]]].s = k;
      }
    }
  }

  return(ier);
}

/* store edges in hash table */
int assignEdge(MMG5_pMesh mesh) {
  _MMG5_Hash  hash;
  MMG5_pTria  pt;
  MMG5_pEdge  pa;
  int         k,ia;
  char        i,i1,i2;


  /* adjust hash table params */
  hash.siz  = mesh->na;
  hash.max  = 3*mesh->na+1;
  _MMG5_ADD_MEM(mesh,(hash.max+1)*sizeof(_MMG5_Hash),"hash table",return(0));
  _MMG5_SAFE_CALLOC(hash.item,hash.max+1,_MMG5_hedge);

  hash.nxt  = mesh->na;
  for (k=mesh->na; k<hash.max; k++)
    hash.item[k].nxt = k+1;

  /* hash mesh edges */
  for (k=1; k<=mesh->na; k++)
    _MMG5_hashEdge(mesh,&hash,mesh->edge[k].a,mesh->edge[k].b,k);

  /* set references to triangles */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      i1 = _MMG5_inxt2[i];
      ia = _MMG5_hashGet(&hash,pt->v[i],pt->v[i1]);
      if ( ia ) {
        i2 = _MMG5_inxt2[i1];
        pa = &mesh->edge[ia];
        pt->edg[i2] = pa->ref;
        pt->tag[i2] = pa->tag;
      }
    }
  }

  /* reset edge structure */
  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(MMG5_Edge));

  return(1);
}
