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
 * \file mmg3d/swap.c
 * \brief Functions for swapping process over boundary.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg3d.h"

extern char ddb;

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param list pointer toward the shell of the edge.
 * \param ilist pointer toward the size of the shell of the edge.
 * \param it1 first element of the open shell.
 * \param it2 last element of the open shell.
 * \return 0 if fail, 1 otherwise.
 *
 * Check whether edge whose shell is provided should be swapped for
 * geometric approximation purposes (the 2 surface triangles are also
 * provided).
 *
 */
int _MMG5_chkswpbdy(MMG5_pMesh mesh, MMG5_pSol met, int *list,int ilist,int it1,int it2) {
  MMG5_pTetra   pt,pt0;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1,ppt0;
  MMG5_Tria     tt1,tt2;
  double        b0[3],b1[3],v[3],c[3],ux,uy,uz,ps,disnat,dischg;
  double        cal1,cal2,calnat,calchg,calold,calnew,caltmp,hausd;
  int           iel,iel1,iel2,np,nq,na1,na2,k,nminus,nplus;
  char          ifa1,ifa2,ia,ip,iq,ia1,ia2,j,isshell;
  MMG5_pPar     par;

  iel = list[0] / 6;
  ia  = list[0] % 6;
  pt  = &mesh->tetra[iel];
  pt0 = &mesh->tetra[0];
  ppt0= &mesh->point[0];

  np = pt->v[_MMG5_iare[ia][0]];
  nq = pt->v[_MMG5_iare[ia][1]];

  /* No swap of geometric edge */
  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    if ( (pxt->edg[ia]>0) || MG_EDG(pxt->tag[ia]) || (pxt->tag[ia] & MG_REQ) ||
         (pxt->tag[ia] & MG_NOM) )  return(0);
  }

  /* No swap when either internal or external component has only 1 element */
  //Algiane: pourquoi on ne check pas ca en multi-domaines?
  if ( mesh->info.iso ) {
    nminus = nplus = 0;
    for (k=0; k<ilist; k++) {
      iel = list[k] / 6;
      pt = &mesh->tetra[iel];
      if ( pt->ref == MG_MINUS )
        nminus++;
      else
        nplus++;
    }
    if ( nplus == 1 || nminus == 1 )  return(0);
  }
  iel1 = it1 / 4;
  ifa1 = it1 % 4;
  iel2 = it2 / 4;
  ifa2 = it2 % 4;
  _MMG5_tet2tri(mesh,iel1,ifa1,&tt1);
  _MMG5_tet2tri(mesh,iel2,ifa2,&tt2);

  for (ia1=0; ia1<3; ia1++) {
    if ( (tt1.v[ia1] != np) && (tt1.v[ia1] != nq) )  break;
  }
  assert( ia1 < 3 );
  assert( (tt1.v[_MMG5_inxt2[ia1]] == np && tt1.v[_MMG5_iprv2[ia1]] == nq) ||
          (tt1.v[_MMG5_inxt2[ia1]] == nq && tt1.v[_MMG5_iprv2[ia1]] == np) );
  na1 = tt1.v[ia1];

  for (ia2=0; ia2<3; ia2++) {
    if ( (tt2.v[ia2] != np) && (tt2.v[ia2] != nq) )  break;
  }

  assert ( ia2 < 3 );
  assert ( (tt2.v[_MMG5_inxt2[ia2]] == np && tt2.v[_MMG5_iprv2[ia2]] == nq) ||
           (tt2.v[_MMG5_inxt2[ia2]] == nq && tt2.v[_MMG5_iprv2[ia2]] == np) );
  na2 = tt2.v[ia2];

  /* Check non convexity (temporarily use b0,b1)*/
  _MMG5_norpts(mesh,tt1.v[ia1],tt1.v[_MMG5_inxt2[ia1]],tt2.v[ia2],b0);
  _MMG5_norpts(mesh,tt2.v[ia2],tt2.v[_MMG5_inxt2[ia2]],tt1.v[ia1],b1);
  ps = b0[0]*b1[0] + b0[1]*b1[1] + b0[2]*b1[2];
  if ( ps < _MMG5_ANGEDG ) return(0);

  /* Compare contributions to Hausdorff distance in both configurations */
  _MMG5_norface(mesh,iel1,ifa1,v);

  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];

  /* local hausdorff for triangles */
  hausd = mesh->info.hausd;
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    if ( (par->elt == MMG5_Triangle) &&
         ((tt1.ref == par->ref ) || (tt2.ref == par->ref)) ) {
      hausd = par->hausd;
    }
  }

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  _MMG5_BezierEdge(mesh,np,nq,b0,b1,0,v);
  c[0] = b0[0] - (p0->c[0] + _MMG5_ATHIRD*ux);
  c[1] = b0[1] - (p0->c[1] + _MMG5_ATHIRD*uy);
  c[2] = b0[2] - (p0->c[2] + _MMG5_ATHIRD*uz);

  disnat = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];

  c[0] = b1[0] - (p1->c[0] - _MMG5_ATHIRD*ux);
  c[1] = b1[1] - (p1->c[1] - _MMG5_ATHIRD*uy);
  c[2] = b1[2] - (p1->c[2] - _MMG5_ATHIRD*uz);

  disnat = MG_MAX(disnat, c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  disnat = MG_MAX(disnat,hausd * hausd);

  p0 = &mesh->point[na1];
  p1 = &mesh->point[na2];
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  _MMG5_BezierEdge(mesh,na1,na2,b0,b1,0,v);
  c[0] = b0[0] - (p0->c[0] + _MMG5_ATHIRD*ux);
  c[1] = b0[1] - (p0->c[1] + _MMG5_ATHIRD*uy);
  c[2] = b0[2] - (p0->c[2] + _MMG5_ATHIRD*uz);

  dischg = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];

  c[0] = b1[0] - (p1->c[0] - _MMG5_ATHIRD*ux);
  c[1] = b1[1] - (p1->c[1] - _MMG5_ATHIRD*uy);
  c[2] = b1[2] - (p1->c[2] - _MMG5_ATHIRD*uz);

  dischg = MG_MAX(dischg,c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  dischg = MG_MAX(dischg,hausd * hausd);

  if ( dischg > disnat )   return(0);
  if ( met->m ) {
    cal1 = _MMG5_caltri(mesh,met,&tt1);
    cal2 = _MMG5_caltri(mesh,met,&tt2);
  }
  else { // with -A option we don't have the metric here so we always use the
         // iso func.
    cal1 = _MMG5_caltri_iso(mesh,met,&tt1);
    cal2 = _MMG5_caltri_iso(mesh,met,&tt2);
  }

  calnat = MG_MIN(cal1,cal2);
  for (j=0; j<3; j++) {
    if ( tt1.v[j] == nq )  tt1.v[j] = na2;
    if ( tt2.v[j] == np )  tt2.v[j] = na1;
  }

  if ( met->m ) {
    cal1 = _MMG5_caltri(mesh,met,&tt1);
    cal2 = _MMG5_caltri(mesh,met,&tt2);
  }
  else { // with -A option we don't have the metric here so we always use the
         // iso func.
    cal1 = _MMG5_caltri_iso(mesh,met,&tt1);
    cal2 = _MMG5_caltri_iso(mesh,met,&tt2);
  }

  calchg = MG_MIN(cal1,cal2);
  if ( calchg < 1.01 * calnat )  return(0);

  /* Check mechanical validity of forthcoming operations */
  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];
  ppt0->c[0] = 0.5*(p0->c[0] + p1->c[0]);
  ppt0->c[1] = 0.5*(p0->c[1] + p1->c[1]);
  ppt0->c[2] = 0.5*(p0->c[2] + p1->c[2]);

  /* Check validity of insertion of midpoint on edge (pq), then collapse of m on a1 */
  calold = calnew = DBL_MAX;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 6;
    pt  = &mesh->tetra[iel];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    calold = MG_MIN(calold, pt->qual);

    ia1 = ia2 = ip = iq = -1;
    for (j=0; j< 4; j++) {
      if (pt->v[j] == np)  ip = j;
      else if (pt->v[j] == nq) iq = j;
      else if ( ia1 < 0 ) ia1 = j;
      else ia2 = j;
    }
    assert((ip >= 0) && (iq >= 0) && (ia1 >= 0) && (ia2 >= 0));
    isshell = (pt->v[ia1] == na1 || pt->v[ia2] == na1);

    /* 2 elts resulting from split and collapse */
    pt0->v[ip] = 0;
    if ( _MMG5_orcal(mesh,met,0) < _MMG5_NULKAL )  return(0);
    if ( !isshell ) {
      pt0->v[ip] = na1;
      caltmp = _MMG5_orcal(mesh,met,0);
      calnew = MG_MIN(calnew,caltmp);
    }
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[iq] = 0;
    if ( _MMG5_orcal(mesh,met,0) < _MMG5_NULKAL )  return(0);

    if ( !isshell ) {
      pt0->v[iq] = na1;
      caltmp = _MMG5_orcal(mesh,met,0);
      calnew = MG_MIN(calnew,caltmp);
    }
  }
  if ( calold < _MMG5_NULKAL && calnew <= calold )  return(0);
  else if ( calnew < 0.3 * calold )  return(0);

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure
 * \param met pointer toward the solution structure
 * \param list pointer toward the shell of the edge
 * \param ret dobble of the number of tetrahedra in the shell
 * \param it1 boundary face carrying the beforehand tested terminal
 * point for collapse
 * \return -1 if lack of memory, 0 if fail to swap, 1 otherwise
 *
 * Swap boundary edge whose shell is provided.
 *
 */
int _MMG5_swpbdy(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ret,int it1,_MMG5_pBucket bucket) {
  MMG5_pTetra   pt,pt1;
  MMG5_pPoint   p0,p1;
  int           iel,iel1,ilist,np,nq,nm;
  double        c[3],*mp, *mq, *mm;
  char          ia,iface1,j,ipa,im;
  int           ier;
#ifndef NDEBUG
  int           na;
#endif

  iel = list[0] / 6;
  ia  = list[0] % 6;
  pt  = &mesh->tetra[iel];

  np = pt->v[_MMG5_iare[ia][0]];
  nq = pt->v[_MMG5_iare[ia][1]];
#ifndef NDEBUG
  na = 0;
#endif

  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];

  /* search for na = the point on quadrangle surfacic configuration on which collapse
     validity has been checked in _MMG5_chkswpbdy */
  iel1 = it1 / 4;
  iface1 = it1 % 4;
  pt1 = &mesh->tetra[iel1];

  for (j=0; j<3;j++) {
    ipa = _MMG5_idir[iface1][j];
    if ( (pt1->v[ipa] != np)&&(pt1->v[ipa] != nq) ) {
#ifndef NDEBUG
      na = pt1->v[ipa];
#endif
      break;
    }
  }
  assert(na);

  /* Create midpoint m on edge (pq), then split edge */
  c[0] = 0.5*( p0->c[0] + p1->c[0]);
  c[1] = 0.5*( p0->c[1] + p1->c[1]);
  c[2] = 0.5*( p0->c[2] + p1->c[2]);
  nm = _MMG5_newPt(mesh,c,MG_BDY);
  if ( !nm ) {
    if ( bucket ) {
      _MMG5_POINT_AND_BUCKET_REALLOC(mesh,met,nm,mesh->gap,
                                     printf("  ## Error: unable to allocate a new point\n");
                                     _MMG5_INCREASE_MEM_MESSAGE();
                                     return(-1)
                                     ,c,MG_BDY);
    }
    else{
      _MMG5_POINT_REALLOC(mesh,met,np,mesh->gap,
                          printf("  ## Error: unable to allocate a new point\n");
                          _MMG5_INCREASE_MEM_MESSAGE();
                          return(-1)
                          ,c,MG_BDY);
    }
  }
  if ( met->m ) {
    mp = &met->m[met->size*np];
    mq = &met->m[met->size*nq];
    mm = &met->m[met->size*nm];
    if ( !_MMG5_intmetvol(mq,mp,mm,0.5) )  return(0);
  }

  ier = _MMG5_split1b(mesh,met,list,ret,nm,0);
  /* pointer adress may change if we need to realloc memory during split */
  pt  = &mesh->tetra[iel];
  pt1 = &mesh->tetra[iel1];

  if ( ier < 0 ) {
    fprintf(stdout,"  ## Warning: unable to swap boundary edge.\n");
    return(-1);
  }
  else if ( !ier )  return(0);

  /* Collapse m on na after taking (new) ball of m */
  memset(list,0,(_MMG5_LMAX+2)*sizeof(int));
  for (j=0; j<3; j++) {
    im = _MMG5_idir[iface1][j];
    if ( pt1->v[im] == nm )  break;
  }
  if ( pt1->v[im] != nm ){
    _MMG5_delPt(mesh,nm);
    fprintf(stdout,"%s:%d: Warning pt1->v[im] != nm\n",__FILE__,__LINE__);
    return(0);
  }
  ilist = _MMG5_boulevolp(mesh,iel1,im,list);

  assert(list[0]/4 == iel1);
  assert(pt1->v[ipa] == na);

  ier = _MMG5_colver(mesh,met,list,ilist,ipa);
  if ( ier < 0 ) {
    fprintf(stdout,"  ## Warning: unable to swap boundary edge.\n");
    return(-1);
  }
  else if ( ier ) {
    _MMG5_delPt(mesh,ier);
    ier = 1;
  }

  return(ier);
}
