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
 * \file common/mmg_hash.c
 * \brief Functions for hash tables management and tetrahedra packing.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 * \param k index of point along the edge.
 * \return 1 if success, 0 if fail.
 *
 * Add edge \f$[a;b]\f$ to the hash table.
 *
 */
int _MMG5_hashEdge(MMG5_pMesh mesh,_MMG5_Hash *hash, int a,int b,int k) {
  _MMG5_hedge  *ph;
  int          key,ia,ib,j;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (_MMG5_KA*ia + _MMG5_KB*ib) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a == ia && ph->b == ib )
    return(1);
  else if ( ph->a ) {
    while ( ph->nxt && ph->nxt < hash->max ) {
      ph = &hash->item[ph->nxt];
      if ( ph->a == ia && ph->b == ib )  return(1);
    }
    ph->nxt   = hash->nxt;
    ph        = &hash->item[hash->nxt];
    hash->nxt = ph->nxt;

    if ( hash->nxt >= hash->max ) {
      if ( mesh->info.ddebug )
        fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",hash->max);
      _MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,0.2,_MMG5_hedge,
                         "_MMG5_edge",return(0));
      for (j=hash->nxt; j<hash->max; j++)  hash->item[j].nxt = j+1;
    }
  }
  /* insert new edge */
  ph->a = ia;
  ph->b = ib;
  ph->k = k;
  ph->nxt = 0;

  return(1);
}

/**
 * \param hash pointer toward the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 * \return the index of point stored along \f$[a;b]\f$.
 *
 * Find the index of point stored along  \f$[a;b]\f$.
 *
 */
int _MMG5_hashGet(_MMG5_Hash *hash,int a,int b) {
  _MMG5_hedge  *ph;
  int          key,ia,ib;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (_MMG5_KA*ia + _MMG5_KB*ib) % hash->siz;
  ph  = &hash->item[key];

  if ( !ph->a )  return(0);
  if ( ph->a == ia && ph->b == ib )  return(ph->k);
  while ( ph->nxt ) {
    ph = &hash->item[ph->nxt];
    if ( ph->a == ia && ph->b == ib )  return(ph->k);
  }
  return(0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward the hash table of edges.
 * \param hsiz initial size of hash table.
 * \param hmax maximal size of hash table.
 * \return 1 if success.
 *
 * Hash edges or faces.
 *
 */
int _MMG5_hashNew(MMG5_pMesh mesh,_MMG5_Hash *hash,int hsiz,int hmax) {
  int   k;

  /* adjust hash table params */
  hash->siz  = hsiz;
  hash->max  = hmax + 1;
  hash->nxt  = hsiz;

  _MMG5_ADD_MEM(mesh,(hash->max+1)*sizeof(_MMG5_hedge),"hash table",
                return(0));
  _MMG5_SAFE_CALLOC(hash->item,hmax+2,_MMG5_hedge);

  for (k=hsiz; k<hash->max; k++)
    hash->item[k].nxt = k+1;

  return(1);
}
