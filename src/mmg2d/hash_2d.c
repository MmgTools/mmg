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
#include "mmg2d.h"

#define KTA     7
#define KTB    11

int MMG2_hashNew(HashTable *hash,int hsize,int hmax) {
  int   k;

  hash->size  = hsize;
  hash->nxtmax =hmax+1;
  hash->hnxt  = hsize;
  _MMG5_SAFE_CALLOC(hash->item,hash->nxtmax,Hedge);

  for (k=hash->size; k<hash->nxtmax; k++)
    hash->item[k].nxt = k+1;

  return(1);
}

/* Create adjacency relations between the triangles in the mesh */
int MMG2_hashTria(MMG5_pMesh mesh) {
  MMG5_pTria     pt,pt1;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1;
  int      *hcode,*link,inival,hsize,iadr;
  unsigned char   i,ii,i1,i2;
  unsigned int    key;
  
  if ( mesh->adja )  return(1);
  if ( !mesh->nt )  return(0);
  
  /* memory alloc */
  _MMG5_SAFE_CALLOC(hcode,mesh->nt+1,int);
  
  /* memory alloc */
  _MMG5_ADD_MEM(mesh,(3*mesh->ntmax+5)*sizeof(int),"adjacency table",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,int);
  
  link  = mesh->adja;
  hsize = mesh->nt;
  
  /* init */
  inival = 2147483647;
  for (k=0; k<=mesh->nt; k++)
    hcode[k] = -inival;
  
  /* build hash table */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    
    if ( !pt->v[0] )  continue;
    for (i=0; i<3; i++) {
      i1 = MMG2_idir[i+1];
      i2 = MMG2_idir[i+2];
      if ( pt->v[i1] < pt->v[i2] ) {
        mins = pt->v[i1];
        maxs = pt->v[i2];
      }
      else {
        mins = pt->v[i2];
        maxs = pt->v[i1];
      }
      
      /* compute key */
      key = KTA*mins + KTB*maxs;
      key = key % hsize + 1;
      
      /* insert */
      iadr = 3*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }
  
  /* set adjacency */
  for (l=3*mesh->nt; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 3 + 1;
    i = (l-1) % 3;
    i1 = MMG2_idir[i+1];
    i2 = MMG2_idir[i+2];
    pt = &mesh->tria[k];
    
    mins = M_MIN(pt->v[i1],pt->v[i2]);
    maxs = M_MAX(pt->v[i1],pt->v[i2]);
    
    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;

    while ( ll != inival ) {
      kk = (ll-1) / 3 + 1;
      ii = (ll-1) % 3;
      i1 = MMG2_idir[ii+1];
      i2 = MMG2_idir[ii+2];
      pt1  = &mesh->tria[kk];
      if ( pt1->v[i1] < pt1->v[i2] ) {
        mins1 = pt1->v[i1];
        maxs1 = pt1->v[i2];
      }
      else {
        mins1 = pt1->v[i2];
        maxs1 = pt1->v[i1];
      }
      
      if ( mins1 == mins  && maxs1 == maxs ) {
        /* adjacent found */
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = 3*kk + ii;
        link[ll]= 3*k + i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  _MMG5_SAFE_FREE(hcode);
  
  return(1);
}

int MMG2_hashel(MMG5_pMesh mesh) {
  MMG5_pTria     pt,pt1;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1;
  int      *hcode,*link,inival,hsize,iadr;
  unsigned char  *hvoy,i,ii,i1,i2;
  unsigned int    key;

  if ( mesh->adja )  return(1);
  if ( !mesh->nt )  return(0);

  /* memory alloc */
  _MMG5_SAFE_CALLOC(hcode,mesh->nt+1,int);

  /* memory alloc */
  _MMG5_ADD_MEM(mesh,(3*mesh->ntmax+5)*sizeof(int),"adjacency table",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,int);

  link  = mesh->adja;
  hsize = mesh->nt;
  hvoy  = (unsigned char*)hcode;

  /* init */
  inival = 2147483647;
  for (k=0; k<=mesh->nt; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !pt->v[0] )  continue;
    for (i=0; i<3; i++) {
      i1 = MMG2_idir[i+1];
      i2 = MMG2_idir[i+2];
      if ( pt->v[i1] < pt->v[i2] ) {
        mins = pt->v[i1];
        maxs = pt->v[i2];
      }
      else {
        mins = pt->v[i2];
        maxs = pt->v[i1];
      }

      /* compute key */
      key = KTA*mins + KTB*maxs;
      key = key % hsize + 1;

      /* insert */
      iadr = 3*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  for (l=3*mesh->nt; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 3 + 1;
    i = (l-1) % 3;
    i1 = MMG2_idir[i+1];
    i2 = MMG2_idir[i+2];
    pt = &mesh->tria[k];

    mins = M_MIN(pt->v[i1],pt->v[i2]);
    maxs = M_MAX(pt->v[i1],pt->v[i2]);

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    hvoy[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 3 + 1;
      ii = (ll-1) % 3;
      i1 = MMG2_idir[ii+1];
      i2 = MMG2_idir[ii+2];
      pt1  = &mesh->tria[kk];
      if ( pt1->v[i1] < pt1->v[i2] ) {
        mins1 = pt1->v[i1];
        maxs1 = pt1->v[i2];
      }
      else {
        mins1 = pt1->v[i2];
        maxs1 = pt1->v[i1];
      }

      if ( mins1 == mins  && maxs1 == maxs ) {
        /* adjacent found */
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = 3*kk + ii;
        link[ll]= 3*k + i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  _MMG5_SAFE_FREE(hcode);

  MMG2_baseBdry(mesh);
  return(1);
}

/*hash edge :
  return 1 if edge exist in the table*/
int MMG2_hashEdge(pHashTable edgeTable,int iel,int ia, int ib) {
  int       key,mins,maxs;
  Hedge     *ha;

  /* compute key */
  if ( ia < ib ) {
    mins = ia;
    maxs = ib;
  }
  else {
    mins = ib;
    maxs = ia;
  }

  key = KTA*mins + KTB*maxs;
  key = key % edgeTable->size;
  ha  = &edgeTable->item[key];
  if ( ha->min ) {
    /* edge exist*/
    if ( ha->min == mins && ha->max == maxs ) {
      return(ha->iel);
    }
    else {
      while ( ha->nxt && ha->nxt < edgeTable->nxtmax ) {
        ha = &edgeTable->item[ha->nxt];
        if ( ha->min == mins && ha->max == maxs )
          return(ha->iel);
      }
      ha->nxt = edgeTable->hnxt;
      ha      = &edgeTable->item[edgeTable->hnxt];
      ++edgeTable->hnxt;
      if ( edgeTable->hnxt == edgeTable->nxtmax ) {
        fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",edgeTable->nxtmax);
        assert(0);
        return(0);
      }
    }
  }
  /* insert */
  ha->min = mins;
  ha->max = maxs;
  ha->iel = iel;
  ha->nxt = 0;

  return(0);

}

/* Transfer some input edge data to the corresponding triangles fields */
int MMG2_assignEdge(MMG5_pMesh mesh) {
  _MMG5_Hash      hash;
  MMG5_pTria      pt;
  MMG5_pEdge      pa;
  int             k,ia;
  char            i,i1,i2;
  
  if ( !mesh->na ) return(1);
  
  /* Temporarily allocate a hash structure for storing edges */
  hash.siz = mesh->na;
  hash.max = 3*mesh->na+1;
  
  _MMG5_ADD_MEM(mesh,(hash.max+1)*sizeof(_MMG5_Hash),"hash table",return(0));
  _MMG5_SAFE_CALLOC(hash.item,hash.max+1,_MMG5_hedge);
  
  hash.nxt = mesh->na;
  
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
  
  /* Delete the hash for edges */
  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(MMG5_Edge));
  mesh->na = 0;
  
  return(1);
}

/* Create the edges in the mesh from the information stored in the triangles, or 
 by identifying the different components of the mesh 
******* Possible extension needed to take into account constrained edges *********** */
int MMG2_bdryEdge(MMG5_pMesh mesh) {
  MMG5_pTria      pt,pt1;
  MMG5_pEdge      pa;
  MMG5_pPoint     p0;
  int             k,*adja,natmp,iel;
  char            i,i1,i2;
  
  natmp = 0;
  mesh->na = 0;
  
  /* First step: Count number of boundary edges */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    
    adja = &mesh->adja[3*(k-1)+1];
    
    for (i=0; i<3; i++) {
      iel = adja[i] / 3;
      pt1 = &mesh->tria[iel];
      
      if ( iel && pt->ref <= pt1->ref ) continue;
      natmp++;
    }
  }
  
  /* Second step: Create edge mesh and store the corresponding edges */
  _MMG5_ADD_MEM(mesh,(natmp+1)*sizeof(MMG5_Edge),"edges",return(0));
  _MMG5_SAFE_CALLOC(mesh->edge,natmp+1,MMG5_Edge);
  
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    
    adja = &mesh->adja[3*(k-1)+1];
    
    for (i=0; i<3; i++) {
      iel = adja[i] / 3;
      pt1 = &mesh->tria[iel];
      
      if ( iel && pt->ref <= pt1->ref ) continue;
      
      i1 = _MMG5_inxt2[i];
      i2 = _MMG5_inxt2[i1];
      
      mesh->na++;
      pa = &mesh->edge[mesh->na];
      pa->a = pt->v[i1];
      pa->b = pt->v[i2];
      
      /* Case of an external boundary edge */
      if ( !iel ) {
        pa->tag = pt->tag[i];
        pa->ref = pt->edg[i];
      }
      /* Case of an internal boundary edge */
      else {
        pa->ref = mesh->info.iso ? MG_ISO : 0;
      }
    }
  }
  
  /* Set point tags */
  for (k=1; k<=mesh->na; k++) {
    pa = &mesh->edge[k];
    p0 = &mesh->point[pa->a];
    p0->tag |= MG_BDY;

    p0 = &mesh->point[pa->b];
    p0->tag |= MG_BDY;
  }
  
  return(1);
}





