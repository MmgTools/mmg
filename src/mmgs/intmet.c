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
 * \file mmgs/intmet.c
 * \brief Metric interpolations.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

extern char ddb;



/* Anisotropic metric interpolation between two points p1 and p2 such that edge
   0 = (p1p2) is ridge.  v is a direction vector, aimed at pointing towards
   direction of n1 at interpolated point */
int intridmet(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,double s,double v[3],double mr[6]) {
  MMG5_pTria     pt;
  MMG5_pxPoint     go1,go2;
  MMG5_pPoint    p1,p2;
  double   *m1,*m2,*n11,*n12,*n21,*n22,ps11,ps12,dd,hn1,hn2;
  int       ip1,ip2;
  char      i1,i2;

  pt  = &mesh->tria[k];
  i1  = _MMG5_inxt2[i];
  i2  = _MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];
  m1  = &met->m[6*ip1];
  m2  = &met->m[6*ip2];

  /* Case when both endpoints are singular */
  if ( MS_SIN(p1->tag) && MS_SIN(p2->tag) ) {
    dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      mr[0] = s < 0.5 ? m1[0] : m2[0];
      mr[1] = s < 0.5 ? m1[0] : m2[0];
      mr[2] = s < 0.5 ? m1[0] : m2[0];
    }
    else {
      mr[0] = m1[0]*m2[0] / dd;
      mr[1] = m1[0]*m2[0] / dd;
      mr[2] = m1[0]*m2[0] / dd;
    }
  }
  /* vertex p1 is singular, p2 is regular */
  else if ( MS_SIN(p1->tag) && (!MS_SIN(p2->tag)) ) {
    go2 = &mesh->xpoint[p2->xp];
    n21 = &go2->n1[0];
    n22 = &go2->n2[0];

    /* Interpolation of the eigenvalue associated to tangent vector */
    dd = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      mr[0] = s < 0.5 ? m1[0] : m2[0];
    }
    else {
      mr[0] = m1[0]*m2[0] / dd;
    }
    /* Interpolation of the two other eigenvalues */
    dd  = (1-s)*sqrt(m2[1]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      hn1 = s < 0.5 ? m1[0] : m2[1];
    }
    else {
      hn1 = m1[0]*m2[1] / dd;
    }
    dd = (1-s)*sqrt(m2[2]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      hn2 = s < 0.5 ? m1[0] : m2[2];
    }
    else {
      hn2 = m1[0]*m2[2] / dd;
    }

    /* Decision of the ordering of hn1 and hn2 */
    ps11 = n21[0]*v[0] + n21[1]*v[1] + n21[2]*v[2];
    ps12 = n22[0]*v[0] + n22[1]*v[1] + n22[2]*v[2];
    if ( fabs(ps11) > fabs(ps12) ) {
      mr[1] = hn1;
      mr[2] = hn2;
    }
    else {
      mr[1] = hn2;
      mr[2] = hn1;
    }
  }
  /* vertex p2 is singular, p1 is regular */
  else if ( MS_SIN(p2->tag) && (!MS_SIN(p1->tag)) ) {
    go1 = &mesh->xpoint[p2->xp];
    n11 = &go1->n1[0];
    n12 = &go1->n2[0];

    /* Interpolation of the eigenvalue associated to tangent vector */
    dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      mr[0] = s < 0.5 ? m1[0] : m2[0];
    }
    else {
      mr[0] = m1[0]*m2[0] / dd;
    }
    /* Interpolation of the two other eigenvalues */
    dd = (1-s)*sqrt(m2[0]) + s*sqrt(m1[1]);
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      hn1 = s < 0.5 ? m1[1] : m2[0];
    }
    else {
      hn1 = m1[1]*m2[0] / dd;
    }
    dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[2]);
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      hn2 = s < 0.5 ? m1[2] : m2[0];
    }
    else {
      hn2 = m1[2]*m2[0] / dd;
    }

    /* Decision of the ordering of hn1 and hn2 */
    ps11 = n11[0]*v[0] + n11[1]*v[1] + n11[2]*v[2];
    ps12 = n12[0]*v[0] + n12[1]*v[1] + n12[2]*v[2];
    if ( fabs(ps11) > fabs(ps12) ) {
      mr[1] = hn1;
      mr[2] = hn2;
    }
    else {
      mr[1] = hn2;
      mr[2] = hn1;
    }
  }
  /* p1,p2 : nonsingular vertices */
  else {
    go1 = &mesh->xpoint[p1->xp];
    go2 = &mesh->xpoint[p2->xp];

    /* Interpolation of the eigenvalue associated to tangent vector */
    dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      mr[0] = s < 0.5 ? m1[0] : m2[0];
    }
    else {
      mr[0] = m1[0]*m2[0] / dd;
    }

    /* Pairing of normal vectors at p1 and p2 */
    n11 = &go1->n1[0];
    n12 = &go1->n2[0];
    n21 = &go2->n1[0];
    n22 = &go2->n2[0];
    ps11 = n11[0]*n21[0] + n11[1]*n21[1] + n11[2]*n21[2];
    ps12 = n11[0]*n22[0] + n11[1]*n22[1] + n11[2]*n22[2];
    if ( fabs(ps11) > fabs(ps12) ) {   //n11 and n21 go together
      dd  = (1-s)*sqrt(m2[1]) + s*sqrt(m1[1]);
      dd *= dd;
      if ( dd < _MMG5_EPSD ) {
        hn1 = s < 0.5 ? m1[1] : m2[1];
      }
      else {
        hn1 = m1[1]*m2[1] / dd;
      }
      dd = (1-s)*sqrt(m2[2]) + s*sqrt(m1[2]);
      dd *= dd;
      if ( dd < _MMG5_EPSD ) {
        hn2 = s < 0.5 ? m1[2] : m2[2];
      }
      else {
        hn2 = m1[2]*m2[2] / dd;
      }
    }
    else {
      dd  = (1-s)*sqrt(m2[2]) + s*sqrt(m1[1]);
      dd *= dd;
      if ( dd < _MMG5_EPSD ) {
        hn1 = s < 0.5 ? m1[1] : m2[2];
      }
      else {
        hn1 = m1[1]*m2[2] / dd;
      }
      dd  = (1-s)*sqrt(m2[1]) + s*sqrt(m1[2]);
      dd *= dd;
      if ( dd < _MMG5_EPSD ) {
        hn2 = s < 0.5 ? m1[2] : m2[1];
      }
      else {
        hn2 = m1[2]*m2[1] / dd;
      }
    }

    /* Now, hn1 is the eigenvalue associated to the direction at interpolated point,
       closest to n11 (hn2 -> n12) ; one may need a different orientation, and put eigenvalue of
       direction closest to v (= interpolated normal) first */
    ps11 = n11[0]*v[0] + n11[1]*v[1] + n11[2]*v[2];
    ps12 = n12[0]*v[0] + n12[1]*v[1] + n12[2]*v[2];
    if ( fabs(ps11) > fabs(ps12) ) {
      mr[1] = hn1;
      mr[2] = hn2;
    }
    else {
      mr[1] = hn2;
      mr[2] = hn1;
    }
  }
  mr[3] = 0.0;
  mr[4] = 0.0;
  mr[5] = 0.0;

  return(1);
}

/* Metric interpolation between points p1 and p2, in tria it at parameter 0 <= s0 <= 1 from p1
   result is stored in mr. edge p1p2 must not be a ridge */
int intregmet(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,double s,double mr[6]) {
  MMG5_pTria     pt;

  pt  = &mesh->tria[k];
  return(_MMG5_interpreg_ani(mesh,met,pt,i,s,mr));

  return(1);
}

/* linear interpolation of sizemap along edge i of tria k */
void intmet_iso(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s) {
  MMG5_pTria  pt;
  int    ip1,ip2;
  char   i1,i2;

  pt  = &mesh->tria[k];
  i1  = _MMG5_inxt2[i];
  i2  = _MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  met->m[ip] = s * (met->m[ip1] + met->m[ip2]);
}

void intmet_ani(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  MMG5_pxPoint    go;
  double  *m;

  pt = &mesh->tria[k];
  m  = &met->m[6*ip];
  if ( pt->tag[i] & MG_GEO ) {
    ppt = &mesh->point[ip];
    assert(ppt->xp);
    go = &mesh->xpoint[ppt->xp];
    intridmet(mesh,met,k,i,s,go->n1,m);
  }
  else {
    intregmet(mesh,met,k,i,s,m);
  }
}

/* Interpolation of size features for a 3*3 metric tensor field, between points np and nq,
   at parameter s from np */
int intmet33(MMG5_pMesh mesh,MMG5_pSol met,int np,int nq,int ip,double s) {
  int     order;
  double  *m,*n,*mr,lambda[3],vp[3][3],mu[3],is[6],isnis[6],mt[9],P[9],dd;
  char    i;

  m  = &met->m[6*np];
  n  = &met->m[6*nq];
  mr = &met->m[6*ip];

  /* Compute inverse of square root of matrix M : is = P*diag(1/sqrt(lambda))*{^t}P */
  order = _MMG5_eigenv(1,m,lambda,vp);
  if ( !order ) return(0);

  for (i=0; i<3; i++) {
    if ( lambda[i] < _MMG5_EPSD ) return(0);
    lambda[i] = sqrt(lambda[i]);
    lambda[i] = 1.0 / lambda[i];
  }

  is[0] = lambda[0]*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0] + lambda[1]*vp[2][0]*vp[2][0];
  is[1] = lambda[0]*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1] + lambda[1]*vp[2][0]*vp[2][1];
  is[2] = lambda[0]*vp[0][0]*vp[0][2] + lambda[1]*vp[1][0]*vp[1][2] + lambda[1]*vp[2][0]*vp[2][2];
  is[3] = lambda[0]*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1] + lambda[1]*vp[2][1]*vp[2][1];
  is[4] = lambda[0]*vp[0][1]*vp[0][2] + lambda[1]*vp[1][1]*vp[1][2] + lambda[1]*vp[2][1]*vp[2][2];
  is[5] = lambda[0]*vp[0][2]*vp[0][2] + lambda[1]*vp[1][2]*vp[1][2] + lambda[1]*vp[2][2]*vp[2][2];

  mt[0] = n[0]*is[0] + n[1]*is[1] + n[2]*is[2];
  mt[1] = n[0]*is[1] + n[1]*is[3] + n[2]*is[4];
  mt[2] = n[0]*is[2] + n[1]*is[4] + n[2]*is[5];
  mt[3] = n[1]*is[0] + n[3]*is[1] + n[4]*is[2];
  mt[4] = n[1]*is[1] + n[3]*is[3] + n[4]*is[4];
  mt[5] = n[1]*is[2] + n[3]*is[4] + n[4]*is[5];
  mt[6] = n[2]*is[0] + n[4]*is[1] + n[5]*is[2];
  mt[7] = n[2]*is[1] + n[4]*is[3] + n[5]*is[4];
  mt[8] = n[2]*is[2] + n[4]*is[4] + n[5]*is[5];

  isnis[0] = is[0]*mt[0] + is[1]*mt[3] + is[2]*mt[6];
  isnis[1] = is[0]*mt[1] + is[1]*mt[4] + is[2]*mt[7];
  isnis[2] = is[0]*mt[2] + is[1]*mt[5] + is[2]*mt[8];
  isnis[3] = is[1]*mt[1] + is[3]*mt[4] + is[4]*mt[7];
  isnis[4] = is[1]*mt[2] + is[3]*mt[5] + is[4]*mt[8];
  isnis[5] = is[2]*mt[2] + is[4]*mt[5] + is[5]*mt[8];

  order = _MMG5_eigenv(1,isnis,lambda,vp);
  if ( !order ) return(0);

  /* P = is * (vp) */
  P[0] = is[0]*vp[0][0] + is[1]*vp[0][1] + is[2]*vp[0][2];
  P[1] = is[0]*vp[1][0] + is[1]*vp[1][1] + is[2]*vp[1][2];
  P[2] = is[0]*vp[2][0] + is[1]*vp[2][1] + is[2]*vp[2][2];
  P[3] = is[1]*vp[0][0] + is[3]*vp[0][1] + is[4]*vp[0][2];
  P[4] = is[1]*vp[1][0] + is[3]*vp[1][1] + is[4]*vp[1][2];
  P[5] = is[1]*vp[2][0] + is[3]*vp[2][1] + is[4]*vp[2][2];
  P[6] = is[2]*vp[0][0] + is[4]*vp[0][1] + is[5]*vp[0][2];
  P[7] = is[2]*vp[1][0] + is[4]*vp[1][1] + is[5]*vp[1][2];
  P[8] = is[2]*vp[2][0] + is[4]*vp[2][1] + is[5]*vp[2][2];

  /* At this point, theory states that ^tPMP = I, {^t}PNP=\Lambda */
  /* Linear interpolation between sizes */
  for(i=0; i<3; i++) {
    if ( lambda[i] < 0.0 ) return(0);
    dd = s*sqrt(lambda[i]) + (1.0-s);
    dd = dd*dd;
    if ( dd < _MMG5_EPSD )  return(0);
    mu[i] = lambda[i]/dd;
  }

  if ( !invmatg(P,mt) )  return(0);

  /* Resulting matrix = ^tP^{-1} diag(mu) P^{-1} */
  mr[0] = mu[0]*mt[0]*mt[0] + mu[1]*mt[3]*mt[3] + mu[2]*mt[6]*mt[6];
  mr[1] = mu[0]*mt[0]*mt[1] + mu[1]*mt[3]*mt[4] + mu[2]*mt[6]*mt[7];
  mr[2] = mu[0]*mt[0]*mt[2] + mu[1]*mt[3]*mt[5] + mu[2]*mt[6]*mt[8];
  mr[3] = mu[0]*mt[1]*mt[1] + mu[1]*mt[4]*mt[4] + mu[2]*mt[7]*mt[7];
  mr[4] = mu[0]*mt[1]*mt[2] + mu[1]*mt[4]*mt[5] + mu[2]*mt[7]*mt[8];
  mr[5] = mu[0]*mt[2]*mt[2] + mu[1]*mt[5]*mt[5] + mu[2]*mt[8]*mt[8];

  return(1);
}
