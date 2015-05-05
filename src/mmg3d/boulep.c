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
 * \file mmg3d/boulep.c
 * \brief Functions for ball of points computation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmg3d.h"

extern MMG5_Info  info;

/** Return volumic ball (i.e. filled with tetrahedra) of point ip in tetra start.
    Results are stored under the form 4*kel + jel , kel = number of the tetra, jel = local
    index of p within kel */
int _MMG5_boulevolp (MMG5_pMesh mesh, int start, int ip, int * list){
  MMG5_pTetra  pt,pt1;
  int    *adja,nump,ilist,base,cur,k,k1;
  char    j,l,i;

  base = ++mesh->base;
  pt   = &mesh->tetra[start];
  nump = pt->v[ip];
  ilist = 0;

  /* Store initial tetrahedron */
  pt->flag = base;
  list[ilist] = 4*start + ip;
  ilist++;

  /* Explore list and travel by adjacency through elements sharing p */
  cur = 0;
  while ( cur < ilist ) {
    k = list[cur] / 4;
    i = list[cur] % 4; // index of point p in tetra k
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = _MMG5_inxt3[i];
      k1 = adja[i] / 4;
      if ( !k1 )  continue;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;
      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);
      /* overflow */
      if ( ilist > _MMG5_LMAX-3 )  return(0);
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }
  return(ilist);
}

/** Define normal and tangent vectors at a non manifold point (ip in start, supported by
    face iface), enumerating its (outer)surfacic ball ; return sng = whether point is singular
    or not */
int _MMG5_boulenm(MMG5_pMesh mesh,int start,int ip,int iface,
                  double n[3],double t[3]) {
  MMG5_pTetra   pt;
  MMG5_pPoint   p0,p1,ppt;
  double   dd,nt[3],l0,l1;
  int      base,nump,nr,nnm,k,piv,na,nb,adj,nvstart,fstart,aux,ip0,ip1;
  int     *adja;
  char     iopp,ipiv,indb,inda,i,ipa,ipb,isface,tag;
  char     indedg[4][4] = { {-1,0,1,2}, {0,-1,3,4}, {1,3,-1,5}, {2,4,5,-1} };

  base = ++mesh->base;
  nr  = nnm = 0;
  ip0 = ip1 = 0;

  memset(n,0.0,3*sizeof(double));
  memset(t,0.0,3*sizeof(double));

  pt   = &mesh->tetra[start];
  nump = pt->v[ip];
  k    = start;

  na   = pt->v[ip];
  nb   = pt->v[_MMG5_idir[iface][_MMG5_inxt2[_MMG5_idirinv[iface][ip]]]];
  piv  = pt->v[_MMG5_idir[iface][_MMG5_iprv2[_MMG5_idirinv[iface][ip]]]];

  iopp   = iface;
  fstart = 4*k+iopp;
  do {
    /* computation of normal and tangent at nump */
    if ( _MMG5_norface(mesh,k,iopp,nt) ) {
      n[0] += nt[0];
      n[1] += nt[1];
      n[2] += nt[2];
    }

    if ( pt->xt ) {
      for ( inda=0; inda<4; inda++ ){
        if ( pt->v[inda]==na ) break;
      }
      for ( indb=0; indb<4; indb++ ){
        if ( pt->v[indb]==nb ) break;
      }
      assert( (inda < 4) && (indb < 4));
      tag = mesh->xtetra[pt->xt].tag[indedg[inda][indb]];
    }

    else  tag = 0;

    if ( MG_EDG(tag) && !(tag & MG_NOM) )
      nr++;
    else if ( tag & MG_NOM ) {
      nnm++;
      if ( !ip0 )
        ip0 = nb;
      else
        ip1 = nb;
    }

    /* A boundary face has been hit : change travel edge */
    aux     = nb;
    nb      = piv;
    piv     = aux;
    nvstart = k;
    adj     = k;

    /* Now unfold shell of edge (na,nb) starting from k (included) */
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      if ( pt->flag != base ) {
        for (i=0; i<4; i++)
          if ( pt->v[i] == nump )  break;
        assert(i<4);
        pt->flag = base;
      }

      /* identification of edge number in tetra k */
      for (i=0; i<6; i++) {
        ipa = _MMG5_iare[i][0];
        ipb = _MMG5_iare[i][1];
        if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
             (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
        adj = adja[ _MMG5_ifar[i][0] ] / 4;
        ipiv = _MMG5_ifar[i][1];
        iopp = _MMG5_ifar[i][0];
        piv = pt->v[ipiv];
      }
      else {
        adj = adja[ _MMG5_ifar[i][1] ] / 4;
        ipiv = _MMG5_ifar[i][0];
        iopp = _MMG5_ifar[i][1];
        piv = pt->v[ipiv];
      }
      isface = (adja[iopp] == 0);
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  while ( 4*k+iopp != fstart );

  if ( (nr > 0 && nnm > 0) || nnm != 2 )  return(0);

  dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( dd > _MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    n[0] *= dd;
    n[1] *= dd;
    n[2] *= dd;
  }
  assert( ip0 && ip1 );
  if ( ip0 == ip1 )  return(0);

  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];
  ppt = &mesh->point[nump];

  l0 = (ppt->c[0] - p0->c[0])*(ppt->c[0] - p0->c[0]) \
    + (ppt->c[1] - p0->c[1])*(ppt->c[1] - p0->c[1]) + (ppt->c[2] - p0->c[2])*(ppt->c[2] - p0->c[2]);
  l1 = (ppt->c[0] - p1->c[0])*(ppt->c[0] - p1->c[0]) \
    + (ppt->c[1] - p1->c[1])*(ppt->c[1] - p1->c[1]) + (ppt->c[2] - p1->c[2])*(ppt->c[2] - p1->c[2]);
  l0 = sqrt(l0);
  l1 = sqrt(l1);

  if ( (l0 < _MMG5_EPSD2) || (l1 < _MMG5_EPSD2) ) {
    t[0] = p1->c[0] - p0->c[0];
    t[1] = p1->c[1] - p0->c[1];
    t[2] = p1->c[2] - p0->c[2];
  }
  else if ( l0 < l1 ) {
    dd = l0 / l1;
    t[0] = dd*(p1->c[0] - ppt->c[0]) + ppt->c[0] - p0->c[0];
    t[1] = dd*(p1->c[1] - ppt->c[1]) + ppt->c[1] - p0->c[1];
    t[2] = dd*(p1->c[2] - ppt->c[2]) + ppt->c[2] - p0->c[2];
  }
  else {
    dd = l1 / l0;
    t[0] = dd*(p0->c[0] - ppt->c[0]) + ppt->c[0] - p1->c[0];
    t[1] = dd*(p0->c[1] - ppt->c[1]) + ppt->c[1] - p1->c[1];
    t[2] = dd*(p0->c[2] - ppt->c[2]) + ppt->c[2] - p1->c[2];
  }
  dd = t[0]*n[0] + t[1]*n[1] + t[2]*n[2];
  t[0] -= dd*n[0];
  t[1] -= dd*n[1];
  t[2] -= dd*n[2];

  dd = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
  if ( dd > _MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    t[0] *= dd;
    t[1] *= dd;
    t[2] *= dd;
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param start index of the starting tetra.
 * \param ip index in \a start of the looked point.
 * \param iface index in \a start of the starting face.
 * \param listv pointer toward the computed volumic ball.
 * \param ilistv pointer toward the computed volumic ball size.
 * \param lists pointer toward the computed surfacic ball.
 * \param ilists pointer toward the computed surfacic ball size.
 * \param isnm is the looked point \a ip non-manifold?
 * \return -1 if fail, 1 otherwise.
 *
 * Compute the volumic ball of a SURFACE point \a p, as well as its surfacic
 * ball, starting from tetra \a start, with point \a ip, and face \a if in tetra
 * volumic ball.
 * \a listv[k] = 4*number of tet + index of point surfacic ball.
 * \a lists[k] = 4*number of tet + index of FAC.
 *
 */
int _MMG5_boulesurfvolp(MMG5_pMesh mesh,int start,int ip,int iface,
                        int *listv,int *ilistv,int *lists,int*ilists, int isnm)
{
  MMG5_pTetra  pt,pt1;
  MMG5_pxTetra pxt;
  int base,nump,k,k1,*adja,piv,na,nb,adj,cur,nvstart,fstart,aux;
  char iopp,ipiv,i,j,l,ipa,ipb,isface;

  base = ++mesh->base;
  *ilists = 0;
  *ilistv = 0;

  pt = &mesh->tetra[start];
  nump = pt->v[ip];
  k = start;

  na = pt->v[ip];
  nb = pt->v[_MMG5_idir[iface][_MMG5_inxt2[_MMG5_idirinv[iface][ip]]]];
  piv = pt->v[_MMG5_idir[iface][_MMG5_iprv2[_MMG5_idirinv[iface][ip]]]];

  iopp = iface;
  fstart = 4*k+iopp;

  do {
    /* A boundary face has been hit : change travel edge */
    lists[(*ilists)] = 4*k+iopp;
    (*ilists)++;
    if ( *ilists >= _MMG5_LMAX ) {
      fprintf(stdout,"  ## Warning: problem in surface remesh process.");
      fprintf(stdout," Surface ball of point %d contains too many elts.\n",
              _MMG5_indPt(mesh,nump));
      fprintf(stdout,"  ##          Try to modify the hausdorff number,");
      fprintf(stdout," or/and the maximum mesh.\n");
      return(-1);
    }

    aux = nb;
    nb = piv;
    piv = aux;
    nvstart = k;
    adj = k;

    /* Now unfold shell of edge (na,nb) starting from k (included)*/
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      if ( pt->flag != base ) {
        for (i=0; i<4; i++)
          if ( pt->v[i] == nump )  break;
        assert(i<4);
        listv[(*ilistv)] = 4*k+i;
        (*ilistv)++;
        pt->flag = base;
      }

      /* identification of edge number in tetra k */
      for (i=0; i<6; i++) {
        ipa = _MMG5_iare[i][0];
        ipb = _MMG5_iare[i][1];
        if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
             (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
        iopp = _MMG5_ifar[i][0];
        ipiv = _MMG5_ifar[i][1];
        adj = adja[ iopp ] / 4;
        piv = pt->v[ipiv];
      }
      else {
        ipiv = _MMG5_ifar[i][0];
        iopp = _MMG5_ifar[i][1];
        adj = adja[ iopp ] / 4;
        piv = pt->v[ipiv];
      }
      if ( isnm ) {
        isface = (adja[iopp] == 0);
      }
      else {
        isface = 0;
        if(pt->xt){
          pxt = &mesh->xtetra[pt->xt];
          isface = (MG_BDY & pxt->ftag[iopp]);
        }
      }
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  while ( 4*k+iopp != fstart );

  /* Now, surfacic ball is complete ; finish travel of volumic ball */
  cur = 0;  // Check numerotation
  while ( cur < (*ilistv) ) {
    k = listv[cur]/4;
    i = listv[cur]%4; // index of point p in tetra k
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = _MMG5_inxt3[i];
      k1 = adja[i]/4;
      if ( !k1 )  continue;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;

      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);

      /* overflow */
      if ( *ilistv > _MMG5_LMAX-3 ) {
        fprintf(stdout,"  ## Warning: problem in remesh process.");
        fprintf(stdout," Volumic ball of point %d contains too many elts.\n",
                _MMG5_indPt(mesh,nump));
        fprintf(stdout,"  ##          Try to modify the hausdorff number,");
        fprintf(stdout," or/and the maximum mesh.\n");
        return(-1);
      }
      listv[(*ilistv)] = 4*k1+j;
      (*ilistv)++;
    }
    cur++;
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param start index of the starting tetrahedron.
 * \param ip index of the looked ridge point.
 * \param iface index in \a start of the starting face.
 * \param il1 pointer toward the first ball size.
 * \param l1 pointer toward the first computed ball (associated to \a n_1's
 * side).
 * \param il2 pointer toward the second ball size.
 * \param l2 pointer toward the second computed ball (associated to \a n_2's
 * side).
 * \param ip1 index of the first extremity of the ridge.
 * \param ip2 index of the second extremity of the ridge.
 * \return 0 if fail, 1 otherwise.
 *
 * Computation of the two surface balls of a ridge point: the list \a l1 is
 * associated to the normal of face \a iface. \a ip0 and \a ip1 are the indices
 * of the 2 ending point of the ridge. Both lists are returned enumerated in
 * direct order.
 *
 */
int _MMG5_bouletrid(MMG5_pMesh mesh,int start,int iface,int ip,int *il1,int *l1,
                    int *il2,int *l2,int *ip0,int *ip1)
{
  MMG5_pTetra          pt;
  MMG5_pxTetra         pxt;
  MMG5_pPoint          ppt;
  int                  k,*adja,*ilist1,*ilist2,*list1,*list2,aux;
  int                  lists[_MMG5_LMAX+2], ilists;
  int                  idp,na, nb, base, iopp, ipiv, piv, fstart, nvstart, adj;
  int                  i,ifac,idx,idx2,idx_tmp,i1,ipa,ipb, isface;
  double               *n1,*n2,nt[3],ps1,ps2;

  pt = &mesh->tetra[start];
  if ( !MG_EOK(pt) )  return(0);

  assert(pt->xt);
  pxt = &mesh->xtetra[pt->xt];
  // We must call this function on a well orientated boundary face (to build the
  // direct surfacic ball).
  assert(pxt->ftag[iface] & MG_BDY);
  assert( MG_GET(pxt->ori,iface) );

  idp = pt->v[ip];
  k = start;

  ppt = &mesh->point[idp];
  assert( ppt->tag & MG_GEO );

  na  = pt->v[ip];
  nb  = pt->v[_MMG5_idir[iface][_MMG5_inxt2[_MMG5_idirinv[iface][ip]]]];
  piv = pt->v[_MMG5_idir[iface][_MMG5_iprv2[_MMG5_idirinv[iface][ip]]]];

  iopp = iface;
  fstart = 4*k+iopp;

  base = ++mesh->base;

  /* Set pointers on lists il1 and il2 to have il1 associated to the normal of
     the face iface.*/
  _MMG5_norpts(mesh, pt->v[_MMG5_idir[iface][0]],pt->v[_MMG5_idir[iface][1]],
               pt->v[_MMG5_idir[iface][2]],nt);

  n1 = &(mesh->xpoint[ppt->xp].n1[0]);
  n2 = &(mesh->xpoint[ppt->xp].n2[0]);
  ps1 = n1[0]*nt[0] + n1[1]*nt[1] + n1[2]*nt[2];
  ps2 = n2[0]*nt[0] + n2[1]*nt[1] + n2[2]*nt[2];

  if ( fabs(ps1) < fabs(ps2) ) {
    list1  = l2;
    list2  = l1;
    ilist1 = il2;
    ilist2 = il1;
  }
  else {
    list1  = l1;
    list2  = l2;
    ilist1 = il1;
    ilist2 = il2;
  }
  *ilist1 = 0;
  *ilist2 = 0;

  /* First: fill the surfacic ball. */
  ilists = 0;
  do {
    /* A boundary face has been hit : change travel edge */
    lists[ilists] = 4*k+iopp;
    ilists++;
    if ( ilists >= _MMG5_LMAX ) {
      fprintf(stdout,"  ## Warning: problem in surface remesh process.");
      fprintf(stdout," Surface ball of point %d contains too many elts.\n",
              _MMG5_indPt(mesh,idp));
      fprintf(stdout,"  ##          Try to modify the hausdorff number,");
      fprintf(stdout," or/and the maximum mesh.\n");
      return(-1);
    }

    aux = nb;
    nb = piv;
    piv = aux;
    nvstart = k;
    adj = k;

    /* Now unfold shell of edge (na,nb) starting from k (included) */
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      pt->flag = base;

      /* identification of edge number in tetra k */
      for (i=0; i<6; i++) {
        ipa = _MMG5_iare[i][0];
        ipb = _MMG5_iare[i][1];
        if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
             (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
        iopp = _MMG5_ifar[i][0];
        ipiv = _MMG5_ifar[i][1];
        adj = adja[ iopp ] / 4;
        piv = pt->v[ipiv];
      }
      else {
        ipiv = _MMG5_ifar[i][0];
        iopp = _MMG5_ifar[i][1];
        adj = adja[ iopp ] / 4;
        piv = pt->v[ipiv];
      }
      isface = 0;
      if(pt->xt){
        pxt = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  // Remark: here the test k!=start is a security bound: theorically it is
  // useless but in case of bad edge tag, it ensure that the loop is not
  // infinite.
  while ( 4*k+iopp != fstart );

  /* Second: travel through the surface ball until meeting a ridge. */
  for (idx=0; idx!=ilists; ++idx) {
    k    = lists[idx]/4;
    ifac = lists[idx]%4;
    pt   = &mesh->tetra[k];
    assert(pt->xt);
    pxt = &mesh->xtetra[pt->xt];

    for ( i=0; i<3; ++i ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    i1 = _MMG5_inxt2[i];
    if ( pxt->tag[_MMG5_iarf[ifac][i1]]  & MG_GEO ) break;
  }
  assert(idx < ilists);
  *ip0 = pt->v[_MMG5_idir[ifac][_MMG5_iprv2[i]]];

  /* Start from the hit boundary, until another boundary is hit and complete the
   * second ball */
  idx = (idx+1)%ilists;
  for (idx2=idx; idx2!=ilists+idx; ++idx2) {
    idx_tmp = idx2%ilists;
    k       = lists[idx_tmp]/4;
    ifac    = lists[idx_tmp]%4;
    pt      = &mesh->tetra[k];
    assert(pt->xt);
    pxt     = &mesh->xtetra[pt->xt];

    if ( (*ilist2) > _MMG5_LMAX-2 )  return(0);
    list2[(*ilist2)] = 4*k+ifac;
    (*ilist2)++;

    for ( i=0; i<3; ++i ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    i1 = _MMG5_inxt2[i];
    if ( pxt->tag[_MMG5_iarf[ifac][i1]]  & MG_GEO ) break;
  }
  assert(idx2 != ilists+idx);
  *ip1 = pt->v[_MMG5_idir[ifac][_MMG5_iprv2[i]]];

  /* Start again from the newly hit boundary, until another boundary is hit and
   * complete the first ball */
  idx = (idx2+1)%ilists;
  for (idx2=idx; idx2 != idx+ilists; ++idx2) {
    idx_tmp = idx2%ilists;
    k       = lists[idx_tmp]/4;
    ifac    = lists[idx_tmp]%4;
    pt      = &mesh->tetra[k];
    assert(pt->xt);
    pxt     = &mesh->xtetra[pt->xt];

    if ( (*ilist1) > _MMG5_LMAX-2 )  return(0);
    list1[(*ilist1)] = 4*k+ifac;
    (*ilist1)++;

    for ( i=0; i<3; ++i ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    i1 = _MMG5_inxt2[i];
    if ( pxt->tag[_MMG5_iarf[ifac][i1]]  & MG_GEO ) break;
  }
  assert(idx2 != ilists+idx);
  assert(*ip0 == pt->v[_MMG5_idir[ifac][_MMG5_iprv2[i]]]);

  return(1);
}

/** Get tag of edge ia in tetra start by travelling its shell until meeting a boundary face */
static inline int
_MMG5_gettag(MMG5_pMesh mesh,int start,int ia,int *tag,int *edg) {
  MMG5_pTetra        pt;
  MMG5_pxTetra       pxt;
  int           na,nb,*adja,adj,piv;
  unsigned char i,ipa,ipb;

  if ( start < 1 )  return(0);
  pt = &mesh->tetra[start];
  if ( !MG_EOK(pt) ) return(0);

  na   = pt->v[ _MMG5_iare[ia][0] ];
  nb   = pt->v[ _MMG5_iare[ia][1] ];

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[_MMG5_ifar[ia][0]] / 4;
  piv = pt->v[_MMG5_ifar[ia][1]];

  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    if ( (pxt->ftag[_MMG5_ifar[ia][0]] & MG_BDY) || (pxt->ftag[_MMG5_ifar[ia][1]] & MG_BDY) ) {
      *tag = pxt->tag[ia];
      *edg = pxt->edg[ia];
      return(1);
    }
  }

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = _MMG5_iare[i][0];
      ipb = _MMG5_iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      if ( (pxt->ftag[_MMG5_ifar[i][0]] & MG_BDY) || (pxt->ftag[_MMG5_ifar[i][1]] & MG_BDY) ) {
        *tag = pxt->tag[i];
        *edg = pxt->edg[i];
        return(1);
      }
    }

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ _MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ _MMG5_ifar[i][1] ];
    }
    else {
      adj = adja[ _MMG5_ifar[i][1] ] /4;
      piv = pt->v[ _MMG5_ifar[i][0] ];
    }
  }
  return(0);
}

/** Set tag and edg of edge ia (if need be) in tetra start by travelling its shell */
inline int
_MMG5_settag(MMG5_pMesh mesh,int start,int ia,int tag,int edg) {
  MMG5_pTetra        pt;
  MMG5_pxTetra       pxt;
  int           na,nb,*adja,adj,piv;
  unsigned char i,ipa,ipb;

  assert( start >= 1 );
  pt = &mesh->tetra[start];
  assert ( MG_EOK(pt) );

  na   = pt->v[ _MMG5_iare[ia][0] ];
  nb   = pt->v[ _MMG5_iare[ia][1] ];

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[_MMG5_ifar[ia][0]] / 4;
  piv = pt->v[_MMG5_ifar[ia][1]];

  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    if ( (pxt->ftag[_MMG5_ifar[ia][0]] & MG_BDY) ||
         (pxt->ftag[_MMG5_ifar[ia][1]] & MG_BDY) ) {
      pxt->tag[ia] |= tag;
      pxt->edg[ia]  = MG_MAX(pxt->edg[ia],edg);
    }
  }
  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = _MMG5_iare[i][0];
      ipb = _MMG5_iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      if ( (pxt->ftag[_MMG5_ifar[i][0]] & MG_BDY) ||
           (pxt->ftag[_MMG5_ifar[i][1]] & MG_BDY) ) {
        pxt->tag[i] |= tag;
        pxt->edg[i]  = MG_MAX(pxt->edg[i],edg);
      }
    }
    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ _MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ _MMG5_ifar[i][1] ];
    }
    else {
      adj = adja[ _MMG5_ifar[i][1] ] /4;
      piv = pt->v[ _MMG5_ifar[i][0] ];
    }
  }

  /* If all shell has been travelled, stop, else, travel it the other sense */
  if ( adj == start )  return(1);
  assert(!adj);

  pt = &mesh->tetra[start];
  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[_MMG5_ifar[ia][1]] / 4;
  piv = pt->v[_MMG5_ifar[ia][0]];

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = _MMG5_iare[i][0];
      ipb = _MMG5_iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      if ( (pxt->ftag[_MMG5_ifar[i][0]] & MG_BDY) ||
           (pxt->ftag[_MMG5_ifar[i][1]] & MG_BDY) ) {
        pxt->tag[i] |= tag;
        pxt->edg[i]  = MG_MAX(pxt->edg[i],edg);
      }
    }
    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ _MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ _MMG5_ifar[i][1] ];
    }
    else {
      adj = adja[ _MMG5_ifar[i][1] ] /4;
      piv = pt->v[ _MMG5_ifar[i][0] ];
    }
  }
  return(1);
}

/** Find all tets sharing edge ia of tetra start
    return 2*ilist if shell is closed, 2*ilist +1 otherwise
    return 0 if one of the tet of the shell is required */
int _MMG5_coquil(MMG5_pMesh mesh,int start,int ia,int * list) {
  MMG5_pTetra  pt;
  int     ilist,*adja,piv,adj,na,nb,ipa,ipb;
  char    i;

  assert ( start >= 1 );
  pt = &mesh->tetra[start];
  assert ( MG_EOK(pt) );

  na   = pt->v[ _MMG5_iare[ia][0] ];
  nb   = pt->v[ _MMG5_iare[ia][1] ];
  ilist = 0;
  list[ilist] = 6*start+ia;
  ilist++;

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[_MMG5_ifar[ia][0]] / 4; // start travelling by face (ia,0)
  piv = pt->v[_MMG5_ifar[ia][1]];

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    if ( pt->tag & MG_REQ )  return(0);
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = _MMG5_iare[i][0];
      ipb = _MMG5_iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    list[ilist] = 6*adj +i;
    ilist++;
    /* overflow */
    if ( ilist > _MMG5_LMAX-3 ) {
      fprintf(stdout,"  ## Warning: problem in remesh process.");
      fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
              _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
      fprintf(stdout,"  ##          Try to modify the hausdorff number,");
      fprintf(stdout," or/and the maximum mesh.\n");
      return(-1);
    }

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ _MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ _MMG5_ifar[i][1] ];
    }
    else {
      assert(pt->v[ _MMG5_ifar[i][1] ] == piv );
      adj = adja[ _MMG5_ifar[i][1] ] /4;
      piv = pt->v[ _MMG5_ifar[i][0] ];
    }
  }

  /* At this point, the first travel, in one direction, of the shell is
     complete. Now, analyze why the travel ended. */
  if ( adj == start )  return(2*ilist);
  assert(!adj); // a boundary has been detected

  adj = list[ilist-1] / 6;
  i   = list[ilist-1] % 6;
  ilist = 0;

  /* Start back everything from this tetra adj */
  list[ilist] = 6*adj + i;
  ilist++;
  /* overflow */
  if ( ilist > _MMG5_LMAX-3 ) {
    fprintf(stdout,"  ## Warning: problem in remesh process.");
    fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
            _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
    fprintf(stdout,"  ##          Try to modify the hausdorff number,");
    fprintf(stdout," or/and the maximum mesh.\n");
    return(-1);
  }

  adja = &mesh->adja[4*(adj-1)+1];
  if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
    adj = adja[ _MMG5_ifar[i][0] ] / 4;
    piv = pt->v[ _MMG5_ifar[i][1] ];
  }
  else {
    adj = adja[ _MMG5_ifar[i][1] ] /4;
    piv = pt->v[ _MMG5_ifar[i][0] ];
  }

  while ( adj ) {
    pt = &mesh->tetra[adj];
    if ( pt->tag & MG_REQ )  return(0);
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = _MMG5_iare[i][0];
      ipb = _MMG5_iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    list[ilist] = 6*adj +i;
    ilist++;
    /* overflow */
    if ( ilist > _MMG5_LMAX-2 ) {
      fprintf(stdout,"  ## Warning: problem in surface remesh process.");
      fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
              _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
      fprintf(stdout,"  ##          Try to modify the hausdorff number,");
      fprintf(stdout," or/and the maximum mesh.\n");
      return(-1);
    }

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ _MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ _MMG5_ifar[i][1] ];
    }
    else {
      adj = adja[ _MMG5_ifar[i][1] ] /4;
      piv = pt->v[ _MMG5_ifar[i][0] ];
    }
  }
  assert(!adj);
  return( 2*ilist+1 );
}

/** Identify whether edge ia in start is a boundary edge by unfolding its shell */
int _MMG5_srcbdy(MMG5_pMesh mesh,int start,int ia) {
  MMG5_pTetra      pt;
  MMG5_pxTetra     pxt;
  int         na,nb,adj,piv,*adja;
  char        ipa,ipb,iadj,i;

  pt = &mesh->tetra[start];
  na = pt->v[_MMG5_iare[ia][0]];
  nb = pt->v[_MMG5_iare[ia][1]];

  adja = &mesh->adja[4*(start-1)+1];
  iadj = _MMG5_ifar[ia][0];

  if(pt->xt){
    pxt = &mesh->xtetra[pt->xt];
    if( pxt->ftag[iadj] & MG_BDY )
      return(1);
  }

  adj = adja[iadj] / 4;
  piv = pt->v[_MMG5_ifar[ia][1]];

  while( adj && ( adj != start ) ) {
    pt = &mesh->tetra[adj];

    /* identification of edge number in tetra adj */
    for(i=0; i<6; i++) {
      ipa = _MMG5_iare[i][0];
      ipb = _MMG5_iare[i][1];
      if( ( pt->v[ipa] == na && pt->v[ipb] == nb ) ||
          ( pt->v[ipa] == nb && pt->v[ipb] == na ))
        break;

    }
    assert(i<6);

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
      iadj = _MMG5_ifar[i][0];
      adj = adja[ iadj ] / 4;
      piv = pt->v[ _MMG5_ifar[i][1] ];
    }
    else {
      iadj = _MMG5_ifar[i][1];
      adj = adja[ iadj ] /4;
      piv = pt->v[ _MMG5_ifar[i][0] ];
    }

    if(pt->xt){
      pxt = &mesh->xtetra[pt->xt];
      if( pxt->ftag[iadj] & MG_BDY )
        return(1);
    }
  }

  return(0);
}

/** print an error message if _MMG5_coquilFace detect a boundary topology problem */
static inline void
_MMG5_errorMessage(MMG5_pMesh mesh, int k1, int k2) {
  MMG5_pPoint ppt;
  MMG5_pTetra pt;
  int    np, ne, k, kel1, kel2;

  np = ne = kel1 = kel2 = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) )  ppt->tmp = ++np;
  }
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( MG_EOK(pt) ) {
      ne++;
      if ( k == k1 ) kel1 = ne;
      if ( k == k2 ) kel2 = ne;
    }
  }

  fprintf(stdout,"  ## Error: problem in surface remesh process");
  fprintf(stdout," (potential creation of a lonely boundary face):\n");

  if ( kel1 != 0 ) {
    pt = &mesh->tetra[k1];
    fprintf(stdout,"            look at elt %d:",kel1);
    fprintf(stdout," %d %d %d %d.\n", mesh->point[pt->v[0]].tmp,
            mesh->point[pt->v[1]].tmp,mesh->point[pt->v[2]].tmp,
            mesh->point[pt->v[3]].tmp);
    fprintf(stdout,"adjacent tetras %d %d %d %d\n",(&mesh->adja[3*(kel1-1)+1])[0],
            (&mesh->adja[3*(kel1-1)+1])[1],(&mesh->adja[3*(kel1-1)+1])[2],
            (&mesh->adja[3*(kel1-1)+1])[3]);
    fprintf(stdout,"vertex required? %d %d %d %d\n",mesh->point[pt->v[0]].tag & MG_REQ,
            mesh->point[pt->v[1]].tag & MG_REQ,
            mesh->point[pt->v[2]].tag & MG_REQ,mesh->point[pt->v[3]].tag & MG_REQ);
  } else if ( kel2 != 0 ) {
    fprintf(stdout,"            look at elt %d:",kel2);
    mesh->tetra[kel2].ref=5;
    fprintf(stdout," %d %d %d %d.\n", mesh->point[pt->v[0]].tmp,
            mesh->point[pt->v[1]].tmp,mesh->point[pt->v[2]].tmp,
            mesh->point[pt->v[3]].tmp);
  }
  fprintf(stdout,"  ##        Try to modify the hausdorff number,");
  fprintf(stdout," the maximum mesh size or/and the value of angle detection.\n");
  fprintf(stdout," You can also try to run with -noswap option but probably");
  fprintf(stdout," the final mesh will have poor quality.\n");
}

/** Find all tets sharing edge ia of tetra start, and stores boundary faces when
    met it1 & it2 = 6*iel + iface, iel = index of tetra, iface = index of face
    in tetra return 2*ilist if shell is closed, 2*ilist +1 otherwise */
int _MMG5_coquilface(MMG5_pMesh mesh,int start,int ia,int *list,int *it1,int *it2) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  int     *adja,piv,adj,na,nb,ipa,ipb,ilist,pradj;
  char     i,iface,isbdy;

  pt = &mesh->tetra[start];

  na   = pt->v[ _MMG5_iare[ia][0] ];
  nb   = pt->v[ _MMG5_iare[ia][1] ];

  ilist = 0;
  list[ilist] = 6*start+ia;
  ilist++;

  *it1 = 0;
  *it2 = 0;

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[_MMG5_ifar[ia][0]] / 4;
  piv = pt->v[_MMG5_ifar[ia][1]];

  pxt = &mesh->xtetra[pt->xt];

  iface = _MMG5_ifar[ia][1];
  isbdy = pxt->ftag[iface];
  if ( isbdy )
    *it1 = 4*start + iface;

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    pxt = 0;
    if ( pt->xt )
      pxt = &mesh->xtetra[pt->xt];

    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = _MMG5_iare[i][0];
      ipb = _MMG5_iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    list[ilist] = 6*adj +i;
    ilist++;
    /* overflow */
    if ( ilist > _MMG5_LMAX-2 ) {
      fprintf(stdout,"  ## Warning: problem in surface remesh process.");
      fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
              _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
      fprintf(stdout,"  ##          Try to modify the hausdorff number,");
      fprintf(stdout," or/and the maximum mesh.\n");
      return(-1);
    }

    /* set new tetra for travel */
    pradj = adj;
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ _MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ _MMG5_ifar[i][1] ];
      iface = _MMG5_ifar[i][1];
    }
    else {
      assert(pt->v[ _MMG5_ifar[i][1] ] == piv );
      adj = adja[ _MMG5_ifar[i][1] ] /4;
      piv = pt->v[ _MMG5_ifar[i][0] ];
      iface = _MMG5_ifar[i][0];
    }
    isbdy = pt->xt ? pxt->ftag[iface] : 0;

    if ( isbdy ) {
      if ( *it1 == 0 )
        *it1 = 4*pradj+iface;
      else {
        if ( *it2 ) {
          // Algiane: (assert commente) si on a 3 tetras de refs differentes dans la _MMG5_coquille??
          printf("  ## Warning: you have more than 2 boundaries in the shell of your edge.\n");
          printf("  Problem may occur during remesh process.\n");
        }
        //assert( *it2 == 0 );
        *it2 = 4*pradj+iface;
      }
    }
  }

  /* At this point, the first travel, in one direction, of the shell is complete. Now, analyze why
     the travel ended. */
  if ( adj == start ) {
    if ( (!(*it1) || !(*it2)) || ((*it1) == (*it2)) ) {
      _MMG5_errorMessage(mesh, (*it1)/4, (*it2)/4);
      return(-1);
    }
    return(2*ilist);
  }

  /* A boundary has been detected : slightly different configuration */
  assert(!adj);
  adj = list[ilist-1] / 6;
  i   = list[ilist-1] % 6;
  ilist = 0;

  /* Start back everything from this tetra adj */
  pradj = adj;
  /* overflow */
  if ( ilist > _MMG5_LMAX-2 ) {
    fprintf(stdout,"  ## Warning: problem in surface remesh process.");
    fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
            _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
    fprintf(stdout,"  ##          Try to modify the hausdorff number,");
    fprintf(stdout," or/and the maximum mesh.\n");
    return(-1);
  }

  pxt = 0;
  if ( pt->xt )
    pxt = &mesh->xtetra[pt->xt];

  adja = &mesh->adja[4*(adj-1)+1];
  if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
    iface = _MMG5_ifar[i][1];
  }
  else {
    iface = _MMG5_ifar[i][0];
  }
  isbdy = pt->xt ? pxt->ftag[iface] : 0;
  if ( isbdy ) {
    *it1 = 4*pradj + iface;
  }

  while ( adj ) {
    pt = &mesh->tetra[adj];

    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = _MMG5_iare[i][0];
      ipb = _MMG5_iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na) )  break;
    }
    assert(i<6);
    list[ilist] = 6*adj +i;
    ilist++;
    /* overflow */
    if ( ilist > _MMG5_LMAX-2 ) {
      fprintf(stdout,"  ## Warning: problem in surface remesh process.");
      fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
              _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
      fprintf(stdout,"  ##          Try to modify the hausdorff number,");
      fprintf(stdout," or/and the maximum mesh.\n");
      return(-1);
    }

    /* set new tetra for travel */
    pradj = adj;
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ _MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ _MMG5_ifar[i][1] ];
      iface = _MMG5_ifar[i][0];
    }
    else {
      adj = adja[ _MMG5_ifar[i][1] ] /4;
      piv = pt->v[ _MMG5_ifar[i][0] ];
      iface = _MMG5_ifar[i][1];
    }
  }

  assert(!adj);
  *it2 = 4*pradj + iface;

  if ( (!(*it1) || !(*it2)) || ((*it1) == (*it2)) ) {
    _MMG5_errorMessage(mesh, (*it1)/4, (*it2)/4);
    return(-1);
  }
  return ( 2*ilist+1 );
}

