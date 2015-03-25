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

#define KA     7
#define KB     11

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
      delElt(mesh,mesh->nt);
    }
  }
  while ( ++k < mesh->nt );

  /* Recreate nil chain */
  mesh->ntnil = mesh->nt + 1;

  for(k=mesh->ntnil; k<=mesh->ntmax-1; k++){
    mesh->tria[k].v[2] = k+1;
  }
}

/* create adjacency */
int hashTria(MMG5_pMesh mesh) {
  MMG5_pTria          pt,pt1;
  MMG5_HGeom     hash;
  MMG5_hgeom    *ph;
  int           *adja,k,jel,hmax,dup,nmf,ia,ib;
  char           i,i1,i2,j,ok;
  unsigned int   key;

  if ( mesh->adja )  return(1);
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING STRUCTURE\n");

  /* tassage */
  paktri(mesh);

  _MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,int);

  /* adjust hash table params */
  hmax = 3.71*mesh->np;
  hash.siz  = mesh->np;
  hash.max  = hmax;
  hash.nxt  = hash.siz;
  _MMG5_SAFE_CALLOC(hash.geom,hash.max+1,MMG5_hgeom);

  for (k=hash.siz; k<hash.max-1; k++)
    hash.geom[k].nxt = k+1;

  /* hash triangles */
  mesh->base = 1;
  dup = nmf = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    pt->flag = 0;
    pt->base = mesh->base;
    adja = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      i1 = inxt[i];
      i2 = iprv[i];

      /* compute key */
      ia  = MG_MIN(pt->v[i1],pt->v[i2]);
      ib  = MG_MAX(pt->v[i1],pt->v[i2]);
      key = (KA*ia + KB*ib) % hash.siz;
      ph  = &hash.geom[key];

      /* store edge */
      if ( ph->a == 0 ) {
        ph->a = ia;
        ph->b = ib;
        ph->k = 3*k + i;
        ph->nxt = 0;
        continue;
      }
      /* update info about adjacent */
      ok = 0;
      while ( ph->a ) {
        if ( ph->a == ia && ph->b == ib ) {
          jel = ph->k / 3;
          j   = ph->k % 3;
          pt1 = &mesh->tria[jel];
          /* discard duplicate face */
          if ( pt1->v[j] == pt->v[i] ) {
            pt1->v[0] = 0;
            dup++;
          }
          /* update adjacent */
          else if ( !mesh->adja[3*(jel-1)+1+j] ) {
            adja[i] = 3*jel + j;
            mesh->adja[3*(jel-1)+1+j] = 3*k + i;
          }
          /* non-manifold case */
          else if ( adja[i] != 3*jel+j ) {
            pt->tag[i] |= MG_GEO + MG_NOM;
            pt1->tag[j]|= MG_GEO + MG_NOM;
            nmf++;
          }
          ok = 1;
          break;
        }
        else if ( !ph->nxt ) {
          ph->nxt = hash.nxt;
          ph = &hash.geom[ph->nxt];
          assert(ph);
          hash.nxt = ph->nxt;
          ph->a = ia;
          ph->b = ib;
          ph->k = 3*k + i;
          ph->nxt = 0;
          ok = 1;
          break;
        }
        else
          ph = &hash.geom[ph->nxt];
      }
      if ( !ok ) {
        ph->a = ia;
        ph->b = ib;
        ph->k = 3*k + i;
        ph->nxt = 0;
      }
    }
  }
  _MMG5_SAFE_FREE(hash.geom);

  /* set tag */
  for (k=1; k<=mesh->nt; k++) {
    pt  = &mesh->tria[k];
    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_NOM ) {
        mesh->point[pt->v[inxt[i]]].tag |= MG_NOM;
        mesh->point[pt->v[iprv[i]]].tag |= MG_NOM;
      }
    }
  }

  /* set seed */
  for (k=1; k<=mesh->nt; k++) {
    pt   = &mesh->tria[k];
    adja = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( !adja[i] )  mesh->point[pt->v[inxt[i]]].s = k;
    }
  }
  if ( nmf > 0 )  mesh->info.mani = 0;

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && dup+nmf > 0 ) {
    fprintf(stdout,"  ## ");  fflush(stdout);
    if ( nmf > 0 )  fprintf(stdout,"[non-manifold model]  ");
    if ( dup > 0 )  fprintf(stdout," %d duplicate removed\n",dup);
    fprintf(stdout,"\n");
  }
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- completed.\n");

  return(1);
}

int hashEdge(MMG5_pMesh mesh, MMG5_HGeom *hash,int a,int b,int k) {
  MMG5_hgeom  *ph;
  int         j,key,ia,ib;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (KA*ia + KB*ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( ph->a ) {
    if ( ph->a != ia || ph->b != ib ) {
      while ( ph->nxt && ph->nxt < hash->max ) {
        ph = &hash->geom[ph->nxt];
        if ( ph->a == ia && ph->b == ib )  return(1);
      }
    }
    ph->nxt = hash->nxt;
    ph      = &hash->geom[hash->nxt];
    ++hash->nxt;
    if ( hash->nxt >= hash->max ) {
      if ( mesh->info.ddebug )
        fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",hash->max);

      hash->max *= 1.2;
      _MMG5_SAFE_REALLOC(hash->geom,hash->max,MMG5_hgeom, "MMG5_hgeom");
      for (j=hash->nxt; j<hash->max; j++)
        hash->geom[j].nxt = j+1;
      return(0);
    }
  }

  /* insert new edge */
  ph->a = ia;
  ph->b = ib;
  ph->k = k;
  ph->nxt = 0;

  return(1);
}

int hashGet(MMG5_HGeom *hash,int a,int b) {
  MMG5_hgeom  *ph;
  int         key,ia,ib;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (KA*ia + KB*ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( !ph->a )  return(0);
  if ( ph->a == ia && ph->b == ib )  return(ph->k);
  while ( ph->nxt ) {
    ph = &hash->geom[ph->nxt];
    if ( ph->a == ia && ph->b == ib )  return(ph->k);
  }
  return(0);
}

/* store edges in hash table */
int assignEdge(MMG5_pMesh mesh) {
  MMG5_HGeom  hash;
  MMG5_pTria  pt;
  MMG5_pEdge  pa;
  int         k,ia;
  char        i,i1,i2;

  _MMG5_SAFE_CALLOC(hash.geom,3*mesh->na+1,MMG5_hgeom);

  /* adjust hash table params */
  hash.siz  = mesh->na;
  hash.max  = 3*mesh->na + 1;
  hash.nxt  = mesh->na;
  for (k=mesh->na; k<hash.max; k++)
    hash.geom[k].nxt = k+1;

  /* hash mesh edges */
  for (k=1; k<=mesh->na; k++)
    hashEdge(mesh,&hash,mesh->edge[k].a,mesh->edge[k].b,k);

  /* set references to triangles */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      i1 = inxt[i];
      ia = hashGet(&hash,pt->v[i],pt->v[i1]);
      if ( ia ) {
        i2 = inxt[i1];
        pa = &mesh->edge[ia];
        pt->edg[i2] = pa->ref;
        pt->tag[i2] = pa->tag;
      }
    }
  }

  /* reset edge structure */
  _MMG5_SAFE_FREE(hash.geom);
  _MMG5_SAFE_FREE(mesh->edge);

  return(1);
}

int hashNew(MMG5_HGeom *hash,int hmax) {
  int   k;

  _MMG5_SAFE_CALLOC(hash->geom,3*hmax+1,MMG5_hgeom);

  /* adjust hash table params */
  hash->siz  = hmax;
  hash->max  = 3*hmax + 1;
  hash->nxt  = hmax;
  for (k=hmax; k<hash->max; k++)
    hash->geom[k].nxt = k+1;

  return(1);
}


