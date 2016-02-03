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
 * \file mmg3d/hash_3d.c
 * \brief Functions for hash tables management and tetrahedra packing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmg3d.h"

#define KC    13

extern char  ddb;

/**
 * \param mesh pointer toward the mesh structure.
 *
 * tetra packing.
 *
 */
static void
_MMG5_paktet(MMG5_pMesh mesh) {
  MMG5_pTetra   pt,pt1;
  int      k;

  k = 1;
  do {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) {
      pt1 = &mesh->tetra[mesh->ne];
      assert(MG_EOK(pt1));
      memcpy(pt,pt1,sizeof(MMG5_Tetra));
      _MMG3D_delElt(mesh,mesh->ne);
    }
  }
  while ( ++k < mesh->ne );

  /* Recreate nil chain */
  mesh->nenil = mesh->ne + 1;

  for(k=mesh->nenil; k<=mesh->nemax-1; k++){
    mesh->tetra[k].v[3] = k+1;
  }
}

static int _MMG5_hashFace(MMG5_pMesh mesh,_MMG5_Hash *hash,int ia,int ib,int ic,int k) {
  _MMG5_hedge     *ph;
  int        key,mins,maxs,sum,j;

  mins = MG_MIN(ia,MG_MIN(ib,ic));
  maxs = MG_MAX(ia,MG_MAX(ib,ic));

  /* compute key */
  sum = ia + ib + ic;
  key = (_MMG5_KA*mins + _MMG5_KB*maxs) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a ) {
    if ( ph->a == mins && ph->b == maxs && ph->s == sum )
      return(ph->k);
    else {
      while ( ph->nxt && ph->nxt < hash->max ) {
        ph = &hash->item[ph->nxt];
        if ( ph->a == mins && ph->b == maxs && ph->s == sum )  return(ph->k);
      }
    }
    ph->nxt = hash->nxt;
    ph      = &hash->item[hash->nxt];
    ph->a   = mins;
    ph->b   = maxs;
    ph->s   = sum;
    ph->k   = k;
    hash->nxt = ph->nxt;
    ph->nxt = 0;

    if ( hash->nxt >= hash->max ) {
      _MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,0.2,_MMG5_hedge,"face",return(0));
      for (j=hash->nxt; j<hash->max; j++)  hash->item[j].nxt = j+1;
    }
    return(1);
  }

  /* insert new face */
  ph->a = mins;
  ph->b = maxs;
  ph->s = sum;
  ph->k = k;
  ph->nxt = 0;

  return(1);
}

/** return index of triangle ia ib ic */
static int _MMG5_hashGetFace(_MMG5_Hash *hash,int ia,int ib,int ic) {
  _MMG5_hedge  *ph;
  int     key,mins,maxs,sum;

  if ( !hash->item )  return(0);
  mins = MG_MIN(ia,MG_MIN(ib,ic));
  maxs = MG_MAX(ia,MG_MAX(ib,ic));

  /* compute key */
  sum = ia + ib + ic;
  key = (_MMG5_KA*mins + _MMG5_KB*maxs) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a ) {
    if ( ph->a == mins && ph->b == maxs && ph->s == sum )
      return(ph->k);
    else {
      while ( ph->nxt ) {
        ph = &hash->item[ph->nxt];
        if ( ph->a == mins && ph->b == maxs && ph->s == sum )  return(ph->k);
      }
    }
  }

  return(0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param pack we pack the mesh at function begining if \f$pack=1\f$.
 * \return 0 if failed, 1 otherwise.
 *
 * Create table of adjacency. Set pack variable to 0 for a compact
 * mesh and to 1 for a mesh that need to be packed.
 *
 */
int MMG3D_hashTetra(MMG5_pMesh mesh, int pack) {
  MMG5_pTetra    pt,pt1;
  int            k,kk,pp,l,ll,mins,mins1,maxs,maxs1,sum,sum1,iadr;
  int           *hcode,*link,hsize,inival;
  unsigned char  i,ii,i1,i2,i3;
  unsigned int   key;

  /* default */
  if ( mesh->adja ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: no re-build of adjacencies of mesh. ");
      fprintf(stdout,"mesh->adja must be freed to enforce analysis.\n");
    }
    return(1);
  }

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING STRUCTURE\n");

  /* packing : if not hash does not work */
  if ( pack )  _MMG5_paktet(mesh);

  /* memory alloc */
  _MMG5_ADD_MEM(mesh,(4*mesh->nemax+5)*sizeof(int),"adjacency table",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->adja,4*mesh->nemax+5,int);
  _MMG5_SAFE_CALLOC(hcode,mesh->ne+5,int);

  link  = mesh->adja;
  hsize = mesh->ne;

  /* init */
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 1: init\n");
  inival = 2147483647;
  iadr   = 0;
  for (k=0; k<=mesh->ne; k++)
    hcode[k] = -inival;

  /* hash tetras */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<4; i++) {
      i1 = _MMG5_idir[i][0];
      i2 = _MMG5_idir[i][1];
      i3 = _MMG5_idir[i][2];
      mins = MG_MIN(pt->v[i1],MG_MIN(pt->v[i2],pt->v[i3]));
      maxs = MG_MAX(pt->v[i1],MG_MAX(pt->v[i2],pt->v[i3]));

      /* compute key and insert */
      sum = pt->v[i1] + pt->v[i2] + pt->v[i3];
      key = _MMG5_KA*mins + _MMG5_KB*maxs + KC*sum;
      key = key % hsize + 1;
      iadr++;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 2: adjacencies\n");
  for (l=iadr; l>0; l--) {
    if ( link[l] >= 0 )  continue;

    /* current element */
    k = (l-1) / 4 + 1;
    i = (l-1) % 4;
    i1 = _MMG5_idir[i][0];
    i2 = _MMG5_idir[i][1];
    i3 = _MMG5_idir[i][2];
    pt = &mesh->tetra[k];
    mins = MG_MIN(pt->v[i1],MG_MIN(pt->v[i2],pt->v[i3]));
    maxs = MG_MAX(pt->v[i1],MG_MAX(pt->v[i2],pt->v[i3]));
    sum  = pt->v[i1] + pt->v[i2] + pt->v[i3];

    /* accross link */
    ll      = -link[l];
    pp      = 0;
    link[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 4 + 1;
      ii = (ll-1) % 4;
      i1 = _MMG5_idir[ii][0];
      i2 = _MMG5_idir[ii][1];
      i3 = _MMG5_idir[ii][2];
      pt1  = &mesh->tetra[kk];
      sum1 = pt1->v[i1] + pt1->v[i2] + pt1->v[i3];
      if ( sum1 == sum ) {
        mins1 = MG_MIN(pt1->v[i1],MG_MIN(pt1->v[i2],pt1->v[i3]));
        maxs1 = MG_MAX(pt1->v[i1],MG_MAX(pt1->v[i2],pt1->v[i3]));

        /* adjacent found */
        if ( mins1 == mins && maxs1 == maxs ) {
          if ( pp != 0 )  link[pp] = link[ll];
          link[l]  = 4*kk + ii;
          link[ll] = 4*k + i;
          break;
        }
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  _MMG5_SAFE_FREE(hcode);
  return(1);
}

/**
 * \param mesh pointer towar the mesh structure.
 * \param hash edges hash table.
 * \return 1 if success, 0 if failed.
 *
 * Set non-manifold tag at extremities of a non-manifold edge.  Check if a
 * non-manifold edge is at the interface of several distinct domains (thus we
 * can't travel from a domain through another one by adjacency) and in this
 * case, mark the edge and its extremities as required. Free the edge hash table
 * \a hash if success.
 *
 * \warning if fail, the edge hash table \a hash is not freed.
 *
 */
static inline
int _MMG5_setEdgeNmTag(MMG5_pMesh mesh, _MMG5_Hash *hash) {
  MMG5_pTetra         ptet;
  MMG5_pxTetra        pxt;
  MMG5_pTria          pt;
  _MMG5_Hash          hashF;
  _MMG5_hedge         *ph;
  int                 *adja,adj,pradj,piv,ilist;
  int                 k,i,i1,i2,ia,ib,l,it1,it2, nr;
  int                 v0,v1,v2,ipa,ipb,ok, kel,count,start;
  unsigned int        key;
  char                isbdy,iface;

  ok = 0;
  nr = 0;

  /* First: seek edges at the interface of two distinct domains and mark it as
   * required */
  for (k=1; k<=mesh->nt; k++) {
    pt  = &mesh->tria[k];

    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_NOM ) {
        i1 = _MMG5_inxt2[i];
        i2 = _MMG5_iprv2[i];

        /* compute key */
        ia  = MG_MIN(pt->v[i1],pt->v[i2]);
        ib  = MG_MAX(pt->v[i1],pt->v[i2]);
        key = (_MMG5_KA*ia + _MMG5_KB*ib) % hash->siz;
        ph  = &hash->item[key];

        assert(ph->a);
        while ( ph->a ) {
          if ( ph->a == ia && ph->b == ib ) break;
          assert(ph->nxt);
          ph = &hash->item[ph->nxt];
        }
        /* Set edge tag and point tags to MG_REQ if the non-manifold edge shared
         * separated domains */
        if ( ph->s > 3 ) {
          if ( !ok ) {

            /* If not already done, hash all the boundary tetra faces. */
            if ( ! _MMG5_hashNew(mesh,&hashF,0.51*mesh->nt,1.51*mesh->nt) )
              return(0);

            for ( kel=1; kel<=mesh->ne; ++kel ) {
              ptet = &mesh->tetra[kel];
              if ( !MG_EOK(ptet) || !ptet->xt ) continue;
              for ( iface=0; iface<4; ++iface) {
                v0 = ptet->v[_MMG5_idir[iface][0]];
                v1 = ptet->v[_MMG5_idir[iface][1]];
                v2 = ptet->v[_MMG5_idir[iface][2]];
                if ( !_MMG5_hashFace(mesh,&hashF,v0,v1,v2,kel) ) return(0);
              }
            }
            ok = 1;
          }
          /* Find the starting tetra, and the edge [ia,ib] in it. */
          start = _MMG5_hashGetFace(&hashF,pt->v[i],ia,ib);
          assert(start);
          ptet = &mesh->tetra[start];
          for (l=0; l<6; ++l) {
            ipa = _MMG5_iare[l][0];
            ipb = _MMG5_iare[l][1];
            if ( (ptet->v[ipa] == ia && ptet->v[ipb] == ib) ||
                 (ptet->v[ipa] == ib && ptet->v[ipb] == ia))  break;
          }
          assert(l<6);


          count = 0;

          /* Travel through the edge's shell and count the number of founded
           * boundaries. */
          ilist = 0;

          it1 = 0;
          it2 = 0;

          pradj = start;
          adja = &mesh->adja[4*(start-1)+1];
          adj = adja[_MMG5_ifar[l][0]] / 4;
          piv = ptet->v[_MMG5_ifar[l][1]];

          pxt = &mesh->xtetra[ptet->xt];

          iface = _MMG5_ifar[l][1];
          isbdy = pxt->ftag[iface];

          if ( isbdy ) {
            it1 = 4*start + iface;
            ++count;
          }

          while ( adj && (adj != start) ) {
            pradj = adj;
            /* travel through new tetra */
            if ( _MMG5_coquilTravel(mesh,ia,ib,&adj,&piv,&iface,&l) ) {
              if ( it1 == 0 ) {
                it1 = 4*pradj+iface;
                ++count;
              }
              else {
                it2 = 4*pradj+iface;
                ++count;
              }
            }

            /* overflow */
            if ( ++ilist > MMG3D_LMAX-2 ) {
              fprintf(stdout,"  ## Warning: problem in surface remesh process.");
              fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
                      _MMG3D_indPt(mesh,ia),_MMG3D_indPt(mesh,ib));
              fprintf(stdout,"  ##          Try to modify the hausdorff number,");
              fprintf(stdout," or/and the maximum mesh.\n");
              return(0);
            }
          }

          /* At this point, the first travel, in one direction, of the shell is
             complete. Now, analyze why the travel ended. */
          if ( adj == start ) {
            if ( (!it1 || !it2) || (it1 == it2) ) {
              _MMG5_coquilFaceErrorMessage(mesh, it1/4, it2/4);
            }
          }
          else {
            /* A boundary has been detected : slightly different configuration */
            assert(!adj);
            adj = pradj;
            ilist = 0;

            /* Start back everything from this tetra adj */
            ptet = &mesh->tetra[adj];

            assert(ptet->xt);
            pxt = &mesh->xtetra[ptet->xt];
            adja = &mesh->adja[4*(adj-1)+1];
            if ( ptet->v[ _MMG5_ifar[l][0] ] == piv ) {
              iface = _MMG5_ifar[l][1];
            }
            else {
              iface = _MMG5_ifar[l][0];
            }
            isbdy = pxt->ftag[iface];
            assert( isbdy );
            it1 = 4*pradj + iface;
            count = 1;

            while ( adj ) {
              pradj = adj;
              _MMG5_openCoquilTravel( mesh, ia, ib, &adj, &piv, &iface, &l );

              /* overflow */
              if ( ++ilist > MMG3D_LMAX-2 ) {
                fprintf(stdout,"  ## Warning: problem in surface remesh process.");
                fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
                        _MMG3D_indPt(mesh,ia),_MMG3D_indPt(mesh,ib));
                fprintf(stdout,"  ##          Try to modify the hausdorff number,");
                fprintf(stdout," or/and the maximum mesh.\n");
                return(0);
              }
              if ( ptet->xt ) {
                pxt = &mesh->xtetra[ptet->xt];
                if ( pxt->ftag[iface] ) ++count;
              }
              ptet = &mesh->tetra[adj];
            }

            assert(!adj);
            it2 = 4*pradj + iface;

            if ( (!it1 || !it2) || (it1 == it2) ) {
              _MMG5_coquilFaceErrorMessage(mesh, it1/4, it2/4);
              return(0);
            }
          }

          /* If ph->s do not match the number of encountred boundaries we have
             separated domains. */
          if ( count != ph->s ) {
            if ( !(pt->tag[i] & MG_REQ) ) {
              pt->tag[i] |= MG_REQ;
              ++nr;
            }
            mesh->point[pt->v[_MMG5_inxt2[i]]].tag |= MG_REQ;
            mesh->point[pt->v[_MMG5_iprv2[i]]].tag |= MG_REQ;
          }

          /* Work done for this edge: reset ph->s/ */
          ph->s = 0;
        }
      }
    }
  }
  if ( mesh->info.ddebug || abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"     %d required edges added\n",nr);

  /* Free the edge hash table */
  _MMG5_DEL_MEM(mesh,hash->item,(hash->max+1)*sizeof(_MMG5_hedge));
  return(1);
}

/**
 * \param mesh pointer towar the mesh structure.
 *
 * Seek the non-required non-manifold points and try to analyse whether they are
 * corner or required.
 *
 * \remark We don't know how to travel through the shell of a non-manifold point
 * by triangle adjacency. Thus the work done here can't be performed in the \ref
 * _MMG5_singul function.
 */
static inline
void _MMG5_setVertexNmTag(MMG5_pMesh mesh) {
  MMG5_pTetra         ptet;
  MMG5_pPoint         ppt;
  int                 k,i;
  int                 nc, nre, ng, nrp;

  /* Second: seek the non-required non-manifold points and try to analyse
   * whether they are corner or required. */
  nc = nre = 0;
  ++mesh->base;
  for (k=1; k<=mesh->ne; ++k) {
    ptet = &mesh->tetra[k];
    if ( !MG_EOK(ptet) ) continue;

    for ( i=0; i<4; ++i ) {
      ppt = &mesh->point[ptet->v[i]];
      if ( (!MG_VOK(ppt)) || (ppt->flag==mesh->base)  ) continue;
      ppt->flag = mesh->base;

      if ( (!(ppt->tag & MG_NOM)) || (ppt->tag & MG_REQ) ) continue;

      if ( !_MMG5_boulernm(mesh, k, i, &ng, &nrp) ) continue;
      if ( (ng+nrp) > 2 ) {
        ppt->tag |= MG_CRN + MG_REQ;
        nre++;
        nc++;
      }
      else if ( (ng == 1) && (nrp == 1) ) {
        ppt->tag |= MG_REQ;
        nre++;
      }
      else if ( ng == 1 && !nrp ){
        ppt->tag |= MG_CRN + MG_REQ;
        nre++;
        nc++;
      }
      else if ( ng == 1 && !nrp ){
        ppt->tag |= MG_CRN + MG_REQ;
        nre++;
        nc++;
      }
    }
  }

  if ( mesh->info.ddebug || abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"     %d corner and %d required vertices added\n",nc,nre);

}

/**
 * \param mesh pointer towar the mesh structure.
 * \param hash edges hash table.
 * \return 1 if success, 0 if failed.
 *
 * Set tags to non-manifold edges and vertices. Not done before because we need
 * the \ref MMG5_xTetra table.
 *
 * \warning if fail, the edge hash table \a hash is not freed.
 *
 */
int _MMG5_setNmTag(MMG5_pMesh mesh, _MMG5_Hash *hash) {

  /* First: seek edges at the interface of two distinct domains and mark it as
   * required */
  if ( !_MMG5_setEdgeNmTag(mesh,hash) ) return(0);

  /* Second: seek the non-required non-manifold points and try to analyse
   * whether they are corner or required. */
  _MMG5_setVertexNmTag(mesh);

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash Edges hash table.
 * \return 1 if success, 0 if failed.
 *
 * Create surface adjacency table. Allocate the edge hash table \a hash but
 * don't free it.
 *
 */
int _MMG3D_hashTria(MMG5_pMesh mesh, _MMG5_Hash *hash) {

  _MMG5_ADD_MEM(mesh,(3*mesh->nt+4)*sizeof(int),"surfacic adjacency table",return(0));
  _MMG5_SAFE_CALLOC(mesh->adjt,3*mesh->nt+4,int);

  return( _MMG5_mmgHashTria(mesh, mesh->adjt, hash, mesh->info.iso) );
}



/** remove edge from hash table */
int _MMG5_hashPop(_MMG5_Hash *hash,int a,int b) {
  _MMG5_hedge  *ph,*php;
  int          key,ia,ib,iph,iphp;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (_MMG5_KA*ia + _MMG5_KB*ib) % hash->siz;
  ph  = &hash->item[key];

  if ( !ph->a ) return(0);
  else if ( ph->a == ia && ph->b == ib ) {
    if ( !ph->nxt ) {
      memset(ph,0,sizeof(_MMG5_hedge));
      return(1);
    }
    else {
      iph = ph->nxt;
      php = ph;
      ph  = &hash->item[ph->nxt];
      memcpy(php,ph,sizeof(_MMG5_hedge));
      memset(ph,0,sizeof(_MMG5_hedge));
      ph->nxt   = hash->nxt;
      hash->nxt = iph;
      return(1);
    }
  }
  while ( ph->nxt ) {
    php = ph;
    ph  = &hash->item[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      if ( !ph->nxt ) {
        memset(ph,0,sizeof(_MMG5_hedge));
        ph->nxt   = hash->nxt;
        hash->nxt = php->nxt;
        php->nxt  = 0;
      }
      else {
        iph  = ph->nxt;
        iphp = php->nxt;
        php->nxt = iph;
        memset(ph,0,sizeof(_MMG5_hedge));
        ph->nxt   = hash->nxt;
        hash->nxt = iphp;
      }
      return(1);
    }
  }
  return(0);
}


/** set tag to edge on geometry */
int _MMG5_hTag(MMG5_HGeom *hash,int a,int b,int ref,char tag) {
  MMG5_hgeom  *ph;
  int     key,ia,ib;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (_MMG5_KA*ia + _MMG5_KB*ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( !ph->a )
    return(0);
  else if ( ph->a == ia && ph->b == ib ) {
    ph->tag |= tag;
    ph->ref  = ref;
    return(1);
  }
  while ( ph->nxt ) {
    ph = &hash->geom[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      ph->tag |= tag;
      ph->ref  = ref;
      return(1);
    }
  }
  return(0);
}

/** remove edge from hash table */
int _MMG5_hPop(MMG5_HGeom *hash,int a,int b,int *ref,char *tag) {
  MMG5_hgeom  *ph,*php;
  int     key,ia,ib,iph,iphp;

  *ref = 0;
  *tag = 0;
  if ( !hash->siz )  return(0);

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (_MMG5_KA*ia + _MMG5_KB*ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( !ph->a )  return(0);
  else if ( ph->a == ia && ph->b == ib ) {
    *ref = ph->ref;
    *tag = ph->tag;
    if ( !ph->nxt ) {
      memset(ph,0,sizeof(MMG5_hgeom));
    }
    else {
      iph = ph->nxt;
      php = ph;
      ph  = &hash->geom[ph->nxt];
      memcpy(php,ph,sizeof(MMG5_hgeom));
      memset(ph,0,sizeof(MMG5_hgeom));
      ph->nxt   = hash->nxt;
      hash->nxt = iph;
    }
    return(1);
  }
  while ( ph->nxt ) {
    php = ph;
    ph  = &hash->geom[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      *ref = ph->ref;
      *tag = ph->tag;
      if ( !ph->nxt ) {
        memset(ph,0,sizeof(MMG5_hgeom));
        ph->nxt   = hash->nxt;
        hash->nxt = php->nxt;
        php->nxt  = 0;
      }
      else {
        iph  = ph->nxt;
        iphp = php->nxt;
        php->nxt = iph;
        memset(ph,0,sizeof(MMG5_hgeom));
        ph->nxt   = hash->nxt;
        hash->nxt = iphp;
      }
      return(1);
    }
  }
  return(0);
}

/** get ref and tag to edge on geometry */
int _MMG5_hGet(MMG5_HGeom *hash,int a,int b,int *ref,char *tag) {
  MMG5_hgeom  *ph;
  int     key,ia,ib;

  *tag = 0;
  *ref = 0;
  if ( !hash->siz )  return(0);
  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (_MMG5_KA*ia + _MMG5_KB*ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( !ph->a )  return(0);
  else if ( ph->a == ia && ph->b == ib ) {
    *ref = ph->ref;
    *tag = ph->tag;
    return(1);
  }
  while ( ph->nxt ) {
    ph = &hash->geom[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      *ref = ph->ref;
      *tag = ph->tag;
      return(1);
    }
  }
  return(0);
}

/** store edge on geometry */
void _MMG5_hEdge(MMG5_pMesh mesh,int a,int b,int ref,char tag) {
  MMG5_hgeom  *ph;
  int     key,ia,ib,j;

  if ( !mesh->htab.siz )  return;
  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (_MMG5_KA*ia + _MMG5_KB*ib) % mesh->htab.siz;
  ph  = &mesh->htab.geom[key];

  if ( ph->a == ia && ph->b == ib )
    return;
  else if ( ph->a ) {
    while ( ph->nxt ) {
      ph = &mesh->htab.geom[ph->nxt];
      if ( ph->a == ia && ph->b == ib )  return;
    }
    ph->nxt = mesh->htab.nxt;
    ph      = &mesh->htab.geom[mesh->htab.nxt];
    ph->a   = ia;   ph->b   = ib;
    ph->ref = ref;  ph->tag = tag;
    mesh->htab.nxt = ph->nxt;
    ph->nxt = 0;
    if ( mesh->htab.nxt >= mesh->htab.max ) {
      if ( mesh->info.ddebug )
        fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",mesh->htab.max);
      _MMG5_TAB_RECALLOC(mesh,mesh->htab.geom,mesh->htab.max,0.2,MMG5_hgeom,
                         "larger htab table",
                         printf("  Exit program.\n");
                         exit(EXIT_FAILURE));
      for (j=mesh->htab.nxt; j<mesh->htab.max; j++)  mesh->htab.geom[j].nxt = j+1;
    }
    return;
  }
  /* insert new edge */
  ph->a   = ia;   ph->b   = ib;
  ph->ref = ref;  ph->tag = tag;
  ph->nxt = 0;
  return;
}

/** to store edge on geometry */
int _MMG5_hNew(MMG5_HGeom *hash,int hsiz,int hmax,int secure) {
  int   k;

  /* adjust hash table params */
  hash->geom = (MMG5_hgeom*)calloc(hmax+2,sizeof(MMG5_hgeom));
  if ( !hash->geom ) {
    perror("  ## Memory problem: calloc");
    if ( !secure )  return(0);
    else  exit(EXIT_FAILURE);
  }
  hash->siz  = hsiz;
  hash->max  = hmax + 1;
  hash->nxt  = hsiz;
  for (k=hsiz; k<hash->max; k++)
    hash->geom[k].nxt = k+1;
  return 1;
}

/**
 * \param mesh pointer toward he mesh structure.
 * \return 0 if failed, 1 otherwise
 *
 * Build hashtable for initial mesh edges.
 *
 */
int _MMG5_hGeom(MMG5_pMesh mesh) {
  MMG5_pTria   pt;
  MMG5_pEdge   pa;
  _MMG5_Hash   hash;
  int         *adja,k,kk,edg,ier;
  char         i,i1,i2,tag;

  /* if edges exist in mesh, hash special edges from existing field */
  if ( mesh->na ) {
    if ( !mesh->htab.geom ) {
      mesh->namax = MG_MAX(1.5*mesh->na,_MMG5_NAMAX);
      _MMG5_ADD_MEM(mesh,(3*mesh->namax+2)*sizeof(MMG5_hgeom),"htab",return(0));
      _MMG5_hNew(&mesh->htab,mesh->na,3*mesh->namax,1);
    }
    else {
      if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
        fprintf(stdout,"  ## Warning: no re-hash of edges of mesh. ");
        fprintf(stdout,"mesh->htab.geom must be freed to enforce analysis.\n");
      }
      _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(MMG5_Edge));
      mesh->na   = 0;
      return(1);
    }

    /* store initial edges */
    for (k=1; k<=mesh->na; k++) {
      pa = &mesh->edge[k];
      _MMG5_hEdge(mesh,pa->a,pa->b,pa->ref,pa->tag);
    }

    /* now check triangles */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      for (i=0; i<3; i++) {
        i1 = _MMG5_inxt2[i];
        i2 = _MMG5_iprv2[i];
        /* transfer non manifold tag to edges */
        if ( pt->tag[i] & MG_NOM )
          _MMG5_hTag(&mesh->htab,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]);

        _MMG5_hGet(&mesh->htab,pt->v[i1],pt->v[i2],&edg,&tag);
        pt->edg[i]  = edg;

        /* If we use the nosurf option and the edge is required, we don't want
         * to detect it as an edge whose tag has been modified for the option */
        if ( mesh->info.nosurf && (tag & MG_REQ) )
          pt->tag[i] &= ~MG_CRN;

        /* Mark edges as boundary edges */
        pt->tag[i] |= (tag | MG_BDY);
        _MMG5_hTag(&mesh->htab,pt->v[i1],pt->v[i2],edg,pt->tag[i]);
      }
    }
    _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(MMG5_Edge));
    mesh->na   = 0;
  }
  /* else, infer special edges from information carried by triangles */
  else {
    if ( !mesh->adjt ) {
      ier = _MMG3D_hashTria(mesh,&hash);
      _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
      if ( !ier ) return(0);
    }

    for (k=1; k<=mesh->nt; k++) {
      pt   = &mesh->tria[k];
      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        i1  = _MMG5_inxt2[i];
        i2  = _MMG5_iprv2[i];
        kk  = adja[i] / 3;
        if ( !kk || pt->tag[i] & MG_NOM )
          mesh->na++;
        else if ( (k < kk) && ( pt->edg[i] || pt->tag[i] ) )  mesh->na++;
      }
    }

    if ( mesh->htab.geom )
      _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));

    mesh->namax = MG_MAX(1.5*mesh->na,_MMG5_NAMAX);
    _MMG5_ADD_MEM(mesh,(3*mesh->namax+2)*sizeof(MMG5_hgeom),"htab",return(0));
    _MMG5_hNew(&mesh->htab,mesh->na,3*mesh->namax,1);
    mesh->na = 0;

    /* build hash for edges */
    for (k=1; k<=mesh->nt; k++) {
      pt   = &mesh->tria[k];
      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        i1  = _MMG5_inxt2[i];
        i2  = _MMG5_iprv2[i];
        kk  = adja[i] / 3;
        if ( !kk || pt->tag[i] & MG_NOM ) {
          if ( pt->tag[i] & MG_NOM ) {
            if ( mesh->info.iso )
              pt->edg[i] = ( pt->edg[i] != 0 ) ?  -abs(pt->edg[i]) : MG_ISO;
          }
          _MMG5_hEdge(mesh,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]);
        }
        else if ( k < kk && ( pt->edg[i] || pt->tag[i] ) )
          _MMG5_hEdge(mesh,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]);
      }
    }
    /* now check triangles */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      for (i=0; i<3; i++) {
        i1 = _MMG5_inxt2[i];
        i2 = _MMG5_iprv2[i];
        _MMG5_hGet(&mesh->htab,pt->v[i1],pt->v[i2],&edg,&tag);
        pt->edg[i]  = edg;
        pt->tag[i] |= tag;
      }
    }
  }
  return(1);
}

/**
 * \brief Check the matching between actual and given number of faces in the
 * mesh.
 * \param mesh pointer toward the mesh structure.
 * \return 1 if the number of counted faces match the number of given one,
 * 0 otherwise.
 *
 * Check the matching between actual and given number of faces in the
 * mesh: Count the number of faces in mesh and compare this number to
 * the number of given triangles. Delete the given triangles if they
 * didn't correspond to the founded faces and fill mesh->nt with the
 * real number of faces.
 *
 */
int _MMG5_chkNumberOfTri(MMG5_pMesh mesh) {
  MMG5_pTetra    pt,pt1;
  int      *adja,adj,k,i,nttmp;

  nttmp = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( !adj || ( pt->ref > pt1->ref) )
        nttmp++;
    }
  }
  if ( mesh->nt == nttmp ) return(1);
  else if ( mesh->nt ){
    if ( !mesh->info.iso && (mesh->info.imprim > 0 || mesh->info.ddebug) )
      fprintf(stdout,"  ## WARNING: INITIAL TRIANGLES ARE DELETED.\n");
    if ( !mesh->info.iso && (mesh->info.imprim > 3 || mesh->info.ddebug) ) {
      fprintf(stdout,"  Not enough or too much triangles for geometry (maybe");
      fprintf(stdout," you have 2 domains but only boundary/interface triangles).\n");
      fprintf(stdout," %d given triangles and %d counted triangles.\n",mesh->nt,nttmp);
    }
    _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(MMG5_Tria));
  }
  mesh->nt = nttmp;
  return(0);
}

/**
 * \param mesh pointer to the mesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Identify boundary triangles.
 *
 * \remark mesh->xtetra is not allocated when \ref _MMG5_bdryTria is called by
 * \ref _MMG3D_analys but it is allocated at save time.
 *
 */
int _MMG5_bdryTria(MMG5_pMesh mesh) {
  MMG5_pTetra    pt,pt1;
  MMG5_pTria     ptt;
  MMG5_pPoint    ppt;
  MMG5_pxTetra   pxt;
  int      *adja,adj,k;
  char      i;

  /* create triangles */
  _MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(MMG5_Tria),"triangles",return(0));
  _MMG5_SAFE_CALLOC(mesh->tria,mesh->nt+1,MMG5_Tria);

  mesh->nt = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    pxt = 0;
    if ( pt->xt )  pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( adj && ( pt->ref <= pt1->ref) )  continue;
      mesh->nt++;
      ptt = &mesh->tria[mesh->nt];
      ptt->v[0] = pt->v[_MMG5_idir[i][0]];
      ptt->v[1] = pt->v[_MMG5_idir[i][1]];
      ptt->v[2] = pt->v[_MMG5_idir[i][2]];
      if ( !adj ) {
        if ( pxt ) {
          if ( pxt->tag[_MMG5_iarf[i][0]] )  ptt->tag[0] = pxt->tag[_MMG5_iarf[i][0]];
          if ( pxt->tag[_MMG5_iarf[i][1]] )  ptt->tag[1] = pxt->tag[_MMG5_iarf[i][1]];
          if ( pxt->tag[_MMG5_iarf[i][2]] )  ptt->tag[2] = pxt->tag[_MMG5_iarf[i][2]];
          /* useful only when saving mesh */
          ptt->ref = pxt->ref[i];
        }
      }
      else {
        if ( pxt ) {
          if ( pxt->tag[_MMG5_iarf[i][0]] )  ptt->tag[0] = pxt->tag[_MMG5_iarf[i][0]];
          if ( pxt->tag[_MMG5_iarf[i][1]] )  ptt->tag[1] = pxt->tag[_MMG5_iarf[i][1]];
          if ( pxt->tag[_MMG5_iarf[i][2]] )  ptt->tag[2] = pxt->tag[_MMG5_iarf[i][2]];
          /* useful only when saving mesh */
        }
        ptt->ref = mesh->info.iso ? MG_ISO : 0;
      }
    }
  }

  /* set point tag */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ptt->v[i]];
      ppt->tag |= MG_BDY;
    }
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 if success.
 *
 * Set the triangles references to the tetrahedra faces and edges.
 *
 */
int _MMG5_bdrySet(MMG5_pMesh mesh) {
  MMG5_pTetra   pt,pt1;
  MMG5_pTria    ptt;
  MMG5_pxTetra  pxt;
  _MMG5_Hash     hash;
  int      *adja,adj,k,kt,ia,ib,ic,j,na;
  char     i,tag;

  if ( !mesh->nt )  return(1);

  if ( mesh->xtetra ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: no re-build of boundary tetras. ");
      fprintf(stdout,"mesh->xtetra must be freed to enforce analysis.\n");
    }
    return(1);
  }

  if ( ! _MMG5_hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) ) return(0);
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    _MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k);
  }
  na = 0;

  mesh->xt     = 0;
  mesh->xtmax  = mesh->ntmax + 2*na;

  _MMG5_ADD_MEM(mesh,(mesh->xtmax+1)*sizeof(MMG5_xTetra),"boundary tetrahedra",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->xtetra,mesh->xtmax+1,MMG5_xTetra);

  /* assign references to tetras faces */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( !adj || (adj > 0 && pt->ref != pt1->ref) ) {
        ia = pt->v[_MMG5_idir[i][0]];
        ib = pt->v[_MMG5_idir[i][1]];
        ic = pt->v[_MMG5_idir[i][2]];
        kt = _MMG5_hashGetFace(&hash,ia,ib,ic);
        assert(kt);
        if ( !pt->xt ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               printf("  Exit program.\n");
                               exit(EXIT_FAILURE));
          }
          pt->xt = mesh->xt;
        }
        ptt = &mesh->tria[kt];
        pxt = &mesh->xtetra[mesh->xt];
        pxt->ref[i]   = ptt->ref;
        pxt->ftag[i] |= MG_BDY;
        pxt->ftag[i] |= (ptt->tag[0] & ptt->tag[1] & ptt->tag[2]);
      }
    }
  }


  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( !pt->xt )  continue;
    pxt = &mesh->xtetra[pt->xt];
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      /* Set flag to know if tetra has the same orientation than the triangle */
      if ( adj && pt->ref < pt1->ref )  MG_CLR(pxt->ori,i);
      else  MG_SET(pxt->ori,i);
      /* Set edge tag */
      if ( pxt->ftag[i] ) {
        if ( adj && (pt->ref <= pt1->ref || (pt->ref == MG_PLUS)) ) {
          continue;
        }
        else {
          ia = pt->v[_MMG5_idir[i][0]];
          ib = pt->v[_MMG5_idir[i][1]];
          ic = pt->v[_MMG5_idir[i][2]];
          kt = _MMG5_hashGetFace(&hash,ia,ib,ic);
          ptt = &mesh->tria[kt];
          for (j=0; j<3; j++) {
            tag = pxt->ftag[i] | ptt->tag[j];
            if ( tag )
              _MMG5_settag(mesh,k,_MMG5_iarf[i][j],tag,ptt->edg[j]);
          }
        }
      }
    }
  }
  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  return(1);
}

/** Update tag and refs of tetra edges.
    If tetra is required, set the faces/edges to required */
int _MMG5_bdryUpdate(MMG5_pMesh mesh) {
  MMG5_pTetra   pt;
  MMG5_pTria    ptt;
  MMG5_pxTetra  pxt;
  _MMG5_Hash     hash;
  int      k,kt,ia,ib,ic,j;
  char     i,tag;

  if ( !mesh->nt )  return(1);
  if ( !_MMG5_hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) )  return(0);
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    _MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k);
  }

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( pt->tag & MG_REQ ) {
      mesh->point[mesh->tetra[k].v[0]].tag |= MG_REQ;
      mesh->point[mesh->tetra[k].v[1]].tag |= MG_REQ;
      mesh->point[mesh->tetra[k].v[2]].tag |= MG_REQ;
      mesh->point[mesh->tetra[k].v[3]].tag |= MG_REQ;
      _MMG5_settag(mesh,k,0,MG_REQ,0);
      _MMG5_settag(mesh,k,1,MG_REQ,0);
      _MMG5_settag(mesh,k,2,MG_REQ,0);
      _MMG5_settag(mesh,k,3,MG_REQ,0);
      _MMG5_settag(mesh,k,4,MG_REQ,0);
      _MMG5_settag(mesh,k,5,MG_REQ,0);
    }

    if ( !pt->xt )  continue;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++) {
      /* Set edge tag */
      if ( ! MG_GET(pxt->ori,i) ) continue;
      if ( pxt->ftag[i] & MG_BDY ) {
        ia = pt->v[_MMG5_idir[i][0]];
        ib = pt->v[_MMG5_idir[i][1]];
        ic = pt->v[_MMG5_idir[i][2]];
        kt = _MMG5_hashGetFace(&hash,ia,ib,ic);
        assert(kt);
        ptt = &mesh->tria[kt];
        if ( pt->tag & MG_REQ ) {
          pxt->ftag[i] |= MG_REQ;
          ptt->tag[0]   = MG_REQ;
          ptt->tag[1]   = MG_REQ;
          ptt->tag[2]   = MG_REQ;
        }
        for ( j=0; j<3; j++ ) {
          tag = ptt->tag[j];
          if ( tag || ptt->edg[j] )
            _MMG5_settag(mesh,k,_MMG5_iarf[i][j],tag,ptt->edg[j]);
        }
      }
    }
  }
  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Make orientation of triangles compatible with tetra faces.
 *
 */
int _MMG5_bdryPerm(MMG5_pMesh mesh) {
  MMG5_pTetra   pt,pt1;
  MMG5_pTria    ptt;
  MMG5_pPoint   ppt;
  _MMG5_Hash    hash;
  int     *adja,adj,k,kt,ia,ib,ic,nf;
  char     i;

  assert(mesh->nt);

  /* store triangles temporarily */
  if ( !_MMG5_hashNew(mesh,&hash,MG_MAX(0.51*mesh->nt,100),MG_MAX(1.51*mesh->nt,300)) )
    return(0);
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !_MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k) )  return(0);
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ptt->v[i]];
      if ( !mesh->info.iso ) ppt->tag |= MG_BDY;
    }
  }

  /* check orientation */
  nf = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( adj && (pt->ref <= pt1->ref || pt->ref == MG_PLUS) )
        continue;
      ia = pt->v[_MMG5_idir[i][0]];
      ib = pt->v[_MMG5_idir[i][1]];
      ic = pt->v[_MMG5_idir[i][2]];
      kt = _MMG5_hashGetFace(&hash,ia,ib,ic);
      if ( !kt ) {
        fprintf(stdout,"%s:%d: Error: function _MMG5_hashGetFace return 0.\n",__FILE__,__LINE__);
        fprintf(stdout," Maybe you have non-boundary triangles.");
        fprintf(stdout," Check triangle of vertices %d %d %d.\n",ia,ib,ic);
        exit(EXIT_FAILURE);
      }

      /* check orientation */
      ptt = &mesh->tria[kt];
      if ( ptt->v[0] == ia && ptt->v[1] == ib && ptt->v[2] == ic )
        continue;
      else {
        ptt->v[0] = ia;
        ptt->v[1] = ib;
        ptt->v[2] = ic;
        nf++;
      }
    }
  }
  if ( mesh->info.ddebug && nf > 0 )
    fprintf(stdout,"  ## %d faces reoriented\n",nf);

  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  return(1);
}
