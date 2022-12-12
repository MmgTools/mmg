/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
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
 * \file mmgs/boulep_s.c
 * \brief Functions for ball of points computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param start index of triangle to start.
 * \param ip index of point for wich we compute the ball.
 * \param ia index of edge in tria \a start from which we start to compute
 * the ball (so edge \a ia has \a ip as extremity).
 * \param list pointer toward the computed ball of \a ip.
 * \param ishell number of tria in the shell of \a ia (to fill)
 *
 * \return the size of the computed ball (ilist) or -ilist if fail.
 *
 * Find all triangles sharing \a ip, starting from \a start and storing first
 * the \a ishell trias in the shell of \a ia. Do not stop when crossing
 * ridge. Works for a non-manifold edge and if we cross nm edges within the ball of \a ip.
 *
 */
int MMGS_boulet(MMG5_pMesh mesh,int start,int ip,int ia,int *list, int *ishell) {
  MMG5_pTria    pt,pt1;
  MMG5_pPoint   ppt;
  int           *adja,cur,k,ilist,nump,k1,j,base;
  char          i,l;

  /* Check that ip is at extremity of edge ia */
  assert ( ip != ia );

  /* Functions related to collapse operators (boulechknm, colver...) implicitely
   * state that iq == MMG5_inxt2[ip] which implies ia == MMG5_iprv2[ip]. */
  assert ( MMG5_iprv2[ip] == ia );

  /* Initialization */
  base = ++mesh->base;
  pt   = &mesh->tria[start];

  nump = pt->v[ip];
  ppt  = &mesh->point[nump];

  /* store initial tria */
  pt->flag = base;
  list[0]  = 3*start + ip;
  *ishell  = 1;

  /** First stage: Fill the shell of ia */
  cur = mesh->adja[ 3*(start-1)+ 1 +ia ];
  while ( cur && (cur != 3*start+ia) ) {
    k = cur/3;
    i = cur%3;

    /* Store tria k and look for nump index */
    mesh->tria[k].flag = base;

    j = MMG5_inxt2[i];
    if ( mesh->tria[k].v[j] != nump ) {
      j = MMG5_inxt2[j];
      assert ( mesh->tria[k].v[j] == nump );
    }

    list[(*ishell)++] = 3*k + j;
    cur = mesh->adja[ 3*(k-1)+ 1 +i ];
  }

#ifndef NDEBUG
  if ( pt->tag[ia] & MG_NOM ) {
      assert ( (*ishell) >=3 );
  }
  else {
    assert ( (*ishell) < 3 );
  }
#endif

  /** Second stage: From each tria of the shell, explore list and travel by
   * adjacency through elements sharing nump */
  cur = 0;
  ilist = *ishell;
  while ( cur < ilist ) {
    k = list[cur] / 3;
    i = list[cur] % 3; // index of point nump in tria k

    adja = &mesh->adja[ 3*(k-1)+ 1 ];

    /* Pile up the 2 adjacent of k */
    for (l=0; l<2; l++) {
      i  = MMG5_inxt2[i];
      k1 = adja[i];
      if ( !k1 ) continue;

      j   = k1%3;
      k1 /= 3;

      pt1 = &mesh->tria[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;

      /* Find the local index of nump */
      j = MMG5_inxt2[j];
      if ( pt1->v[MMG5_inxt2[j]] == nump ) {
        j = MMG5_inxt2[j];
      }
      assert ( pt1->v[j] == nump );

      /* Check for overflow */
      if ( ilist > MMGS_LMAX-2 )  return -ilist;

      /* Store adjacent */
      list[ilist] = 3*k1 + j;
      ++ilist;
    }
    ++cur;
  }

  return ilist;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param start index of triangle to start.
 * \param ip index of point for wich we compute the ball.
 * \param list pointer toward the computed ball of \a ip.
 * \return the size of the computed ball or 0 if fail.
 *
 * Find all triangles sharing \a ip with ip a manifold point, \f$list[0] =\f$ \a
 * start do not stop when crossing ridge.
 *
 */
int MMGS_bouletmani(MMG5_pMesh mesh,int start,int ip,int *list) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  int           *adja,k,ilist;
  char          i,i1,i2;

  pt = &mesh->tria[start];

  ppt = &mesh->point[pt->v[ip]];
  ilist = 0;

  /* store neighbors */
  k = start;
  i = ip;
  do {
    if ( ilist > MMGS_LMAX-2 )  return 0;
    list[ilist] = 3*k + i;
    ++ilist;

    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = MMG5_inxt2[i];
  }
  while ( k && k != start );
  if ( k > 0 )  return ilist;

  if ( ppt->tag & MG_NOM )
    return 0;

  /* check if boundary hit */
  k = start;
  i = ip;
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i2 = MMG5_iprv2[i];
    k  = adja[i2] / 3;
    if ( k == 0 )  break;
    i  = adja[i2] % 3;
    i  = MMG5_iprv2[i];

    if ( ilist > MMGS_LMAX-2 )  return 0;
    list[ilist] = 3*k + i;
    ilist++;
  }
  while ( k );

  return ilist;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param start index of tria to start to compute the ball.
 * \param ip index of point in tria \a start for which we want to compute
 * the ball.
 * \param ie index of edge in tria \a start from which we start to compute
 * the ball (so edge \a ie has \a ip as extremity).
 * \param list pointer toward the computed ball of point.
 * \param ishell number of tria in the shell of \a ie (to fill)
 *
 * \return -ilist (number of tria in ball) if buffer overflow, ilist if the
 * collapse is ok, 0 if the collapse lead to a non manifold situation.
 *
 * Find all triangles sharing \a ip, starting from \a start and storing first
 * the \a ishell tria in the shell of \a ie. Do not stop when crossing
 * ridge. Check whether resulting configuration is manifold.
 *
 */
int MMGS_boulechknm(MMG5_pMesh mesh,int start,int ip,int ie,int *list,int *ishell) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  int           *adja,k,k1,ilist,base,iel,listq[MMGS_LMAX+2],ilistq,numq,cur;
  char          i,j,iq,voy;

  pt = &mesh->tria[start];
  if ( !MG_EOK(pt) )  return 0;

#warning ajeter : old test :  Ca doit marcher non?
  if ( mesh->point[pt->v[ip]].tag & MG_NOM ) return 0;

  /** Stage 1: fill the ball of point and store it in \a list (such as
   * list[k]%3 == ip ) */
  ilist = MMGS_boulet(mesh,start,ip,ie,list,ishell);

  if ( ilist < 1 ) return ilist;

  /** Stage 2: check if a collapse may lead to a non-convex situation */
  /* 1. Flags initialization: set all the flags of the points of the ball to base */
  base = ++mesh->base;
  for ( k=(*ishell); k<ilist; ++k ) {
    iel = list[k]/3;
    i   = list[k]%3;
    pt  = &mesh->tria[iel];
    mesh->point[pt->v[MMG5_inxt2[i]]].s = base;
    mesh->point[pt->v[MMG5_iprv2[i]]].s = base;
  }

  /* 2. Flags initialization: reset the fkags of the points of the triangles of
   * the shell ip-iq to 0. */
  for ( k=0; k<(*ishell); ++k ) {
    iel = list[k]/3;
    i   = list[k]%3;
    pt  = &mesh->tria[iel];
    mesh->point[pt->v[MMG5_inxt2[i]]].s = 0;
    mesh->point[pt->v[MMG5_iprv2[i]]].s = 0;
  }

  /* 3. Unfold the ball af iq and check that it doesn't contains common tria with
   * the ball of ip (except for the ip-iq shell) */
  pt       = &mesh->tria[start];
  j        = MMG5_iprv2[ip];
  assert ( j == ie && "Input edge index for the collapse don't match with the computed one." );
  iq       = MMG5_inxt2[ip];
  numq     = pt->v[iq];

  listq[0] = 3*start + iq;
  pt->flag = base;

  cur    = 0;
  ilistq = 1;
  while ( cur < ilistq ) {
    k = listq[cur] / 3;
    i = listq[cur] % 3; // index of point iq in tria k

    adja = &mesh->adja[ 3*(k-1)+ 1 ];

    /* Pile up the 2 adjacent of k */
    for (voy=0; voy<2; voy++) {
      i  = MMG5_inxt2[i];
      k1 = adja[i];
      if ( !k1 ) continue;

      j   = k1%3;
      k1 /= 3;

      pt = &mesh->tria[k1];
      if ( pt->flag == base )  continue;
      pt->flag = base;

      /* Find the local index of numq */
      j = MMG5_inxt2[j];
      if ( pt->v[MMG5_inxt2[j]] == numq ) {
        j = MMG5_inxt2[j];
      }
      assert ( pt->v[j] == numq );

      ppt = &mesh->point[pt->v[MMG5_inxt2[j]]];
      if ( ppt->s == base ) return 0;

      ppt = &mesh->point[pt->v[MMG5_iprv2[j]]];
      if ( ppt->s == base ) return 0;

      /* Check for overflow */
      if ( ilistq > MMGS_LMAX-2 )  return -ilist;

      /* Store adjacent */
      listq[ilistq] = 3*k1 + j;
      ++ilistq;
    }
    ++cur;
  }

  return ilist;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param start index of the starting triangle.
 * \param ip index of the looked ridge point.
 * \param il1 pointer toward the first ball size.
 * \param l1 pointer toward the first computed ball (associated to \a n1's
 * side).
 * \param il2 pointer toward the second ball size.
 * \param l2 pointer toward the second computed ball (associated to \a n2's
 * side).
 * \param ip0 index of the first extremity of the ridge.
 * \param ip1 index of the second extremity of the ridge.
 * \return 0 if fail, 1 otherwise.
 *
 * Computation of the two balls of a ridge point: the list \a l1 is associated
 * to normal \a n1's side. \a ip0 and \a ip1 are the indices of the 2 ending
 * point of the ridge. Both lists are returned enumerated in direct order.
 *
 * \warning can't be called from a non-manifold point
 */
int bouletrid(MMG5_pMesh mesh,int start,int ip,int *il1,int *l1,int *il2,int *l2,int *ip0,int *ip1) {
  MMG5_pTria           pt;
  MMG5_pPoint          ppt;
  int                  idp,k,kold,*adja,iel,*ilist1,*ilist2,*list1,*list2,aux;
  unsigned char        i,iold,i1,i2,ipn;
  double               *n1,*n2,nt[3],ps1,ps2;

  pt = &mesh->tria[start];
  if ( !MG_EOK(pt) )  return 0;

  idp = pt->v[ip];
  ppt = &mesh->point[idp];

  assert( ppt->tag & MG_GEO );
  assert( !(ppt->tag & MG_NOM) );

  /* set pointers: first manifold is on side of triangle */
  if ( !MMG5_nortri(mesh,pt,nt) )  return 0;

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

  /* First ball, first side (via i1)*/
  k = start;
  i = ip;
  do {
    pt   = &mesh->tria[k];
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];
    kold = k;
    iold = i;
    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_inxt2[i];
  }
  // Remark: here the test k!=start is a security bound: theorically it is
  // useless but in case of bad edge tag, it ensure that the loop is not
  // infinite.
  while ( k && !(pt->tag[i1] & MG_GEO) && k != start );
  *ip0 = pt->v[i2];

  /* Store the needed elements to start back in the new area,
     and complete first ball, second side (via i2) */
  iel = k;
  ipn = i;
  k  = kold;
  i  = iold;
  do {
    pt   = &mesh->tria[k];
    adja = &mesh->adja[3*(k-1)+1];
    if ( (*ilist1) > MMGS_LMAX-2 )  return 0;
    list1[(*ilist1)] = 3*k+i;
    (*ilist1)++;
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];
    k = adja[i2] / 3;
    i = adja[i2] % 3;
    i = MMG5_iprv2[i];
  }
  while ( k && !(MG_GEO & pt->tag[i2]) );
  *ip1 = pt->v[i1];

  /* Invert the order in list1, so that it is enumerated in DIRECT order */
  for (k=0; k<(*ilist1) / 2;k++) {
    aux = list1[k];
    list1[k] = list1[*ilist1-1-k];
    list1[*ilist1-1-k] = aux;
  }
  /* At this point, either something has been stored in iel or an open boundary has been hit */
  *ilist2 = 0;
  if ( !iel )  return 1;

  /* Else, start back from the hit boundary, until another boundary is hit */
  k  = iel;
  i  = ipn;
  do {
    pt   = &mesh->tria[k];
    adja = &mesh->adja[3*(k-1)+1];
    if ( *ilist2 > MMGS_LMAX-2 )  return 0;
    list2[*ilist2] = 3*k+i;
    (*ilist2)++;
    i1 = MMG5_inxt2[i];
    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_inxt2[i];
  }
  while ( k && !(MG_GEO & pt->tag[i1]) );

  if ( !(MG_GEO & pt->tag[i1]) )
    return 0;
  else
    return 1;
}
