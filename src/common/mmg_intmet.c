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
 * \file common/mmg_intmet.c
 * \brief Functions to compute metric interpolation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg.h"


/* Compute the interpolated (2 x 2) metric from metrics m and n, at parameter s :
   mr = (1-s)*m +s*n, both metrics being expressed in the simultaneous reduction basis:
   linear interpolation of sizes */
static int MMG5_intmet22(double *m,double *n,double *mr,double s) {
  double  det,imn[4],dd,sqDelta,trimn,lambda[2],vp0[2],vp1[2],dm[2],dn[2],vnorm,d0,d1,ip[4];

  /* Compute imn = M^{-1}N */
  det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < _MMG5_EPS*_MMG5_EPS ) {
    printf("BEWARE : function intmet : null metric det : %E \n",det);
    return(0);
  }
  det = 1.0 / det;

  imn[0] = det * ( m[2]*n[0] - m[1]*n[1]);
  imn[1] = det * ( m[2]*n[1] - m[1]*n[2]);
  imn[2] = det * (-m[1]*n[0] + m[0]*n[1]);
  imn[3] = det * (-m[1]*n[1] + m[0]*n[2]);
  dd = imn[0] - imn[3];
  sqDelta = sqrt(fabs(dd*dd + 4.0*imn[1]*imn[2]));
  trimn = imn[0] + imn[3];

  lambda[0] = 0.5 * (trimn - sqDelta);
  if ( lambda[0] < 0.0 ) {
    printf("Les valeurs propres : %f \n",lambda[0]);
    return(0);
  }

  /* First case : matrices m and n are homothetic = n = lambda0*m */
  if ( sqDelta < _MMG5_EPS ) {
    dd  = (1.0-s)*sqrt(lambda[0]) + s;
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      mr[0] = m[0];
      mr[1] = m[1];
      mr[2] = m[2];
      return(1);
    }
    dd = lambda[0] / dd;
    mr[0] = dd * m[0];
    mr[1] = dd * m[1];
    mr[2] = dd * m[2];
    return(1);
  }

  /* Second case : both eigenvalues of imn are distinct ; theory says qf associated to m and n
     are diagonalizable in basis (vp0, vp1) - the coreduction basis */
  else {
    lambda[1] = 0.5 * (trimn + sqDelta);
    assert(lambda[1] >= 0.0);

    vp0[0] = imn[1];
    vp0[1] = (lambda[0] - imn[0]);
    vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);
    if ( vnorm < _MMG5_EPS ) {
      vp0[0] = (lambda[0] - imn[3]);
      vp0[1] = imn[2];
      vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);
    }
    vnorm   = 1.0 / vnorm;
    vp0[0] *= vnorm;
    vp0[1] *= vnorm;

    vp1[0] = imn[1];
    vp1[1] = (lambda[1] - imn[0]);
    vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);
    if ( vnorm < _MMG5_EPS ) {
      vp1[0] = (lambda[1] - imn[3]);
      vp1[1] = imn[2];
      vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);
    }
    vnorm   = 1.0 / vnorm;
    vp1[0] *= vnorm;
    vp1[1] *= vnorm;

    /* Compute diagonal values in simultaneous reduction basis */
    dm[0] = m[0]*vp0[0]*vp0[0] + 2.0*m[1]*vp0[0]*vp0[1] + m[2]*vp0[1]*vp0[1];
    dm[1] = m[0]*vp1[0]*vp1[0] + 2.0*m[1]*vp1[0]*vp1[1] + m[2]*vp1[1]*vp1[1];
    dn[0] = n[0]*vp0[0]*vp0[0] + 2.0*n[1]*vp0[0]*vp0[1] + n[2]*vp0[1]*vp0[1];
    dn[1] = n[0]*vp1[0]*vp1[0] + 2.0*n[1]*vp1[0]*vp1[1] + n[2]*vp1[1]*vp1[1];

    /* Diagonal values of the interpolated metric */
    dd  = (1.0-s)*sqrt(dn[0]) + s*sqrt(dm[0]);
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      d0 = s < 0.5 ? dm[0] : dn[0];
    }
    else {
      d0 = dm[0]*dn[0] / dd;
    }
    dd = (1.0-s)*sqrt(dn[1]) + s*sqrt(dm[1]);
    dd *= dd;
    if ( dd < _MMG5_EPSD ) {
      d1 = s < 0.5 ? dm[1] : dn[1];
    }
    else{
      d1 = dm[1]*dn[1] / dd;
    }

    /* Intersected metric = tP^-1 diag(d0,d1)P^-1, P = (vp0, vp1) stored in columns */
    det = vp0[0]*vp1[1] - vp0[1]*vp1[0];
    if ( fabs(det) < _MMG5_EPS )  return(0);
    det = 1.0 / det;

    ip[0] =  vp1[1]*det;
    ip[1] = -vp1[0]*det;
    ip[2] = -vp0[1]*det;
    ip[3] =  vp0[0]*det;

    mr[0] = d0*ip[0]*ip[0] + d1*ip[2]*ip[2];
    mr[1] = d0*ip[0]*ip[1] + d1*ip[2]*ip[3];
    mr[2] = d0*ip[1]*ip[1] + d1*ip[3]*ip[3];
  }

  return(1);
}


/**
 * \param ma pointer on a metric
 * \param mb pointer on a metric
 * \param mp pointer on the computed interpolated metric
 * \param t interpolation parameter (comprise between 0 and 1)
 *
 *
 * Linear interpolation of isotropic sizemap along an edge
 *
 */
int _MMG5_interp_iso(double *ma,double *mb,double *mp,double t) {

  *mp = (1.0-t)*(*ma) + t*(*mb);

  return(1);

}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param pt pointer toward the triangle structure.
 * \param i char : edge of the triangle pt
 * \param s interpolated parameter (comprise between 0 and 1)
 * \param mr computed interpolated metric
 *
 * Metric interpolation between points p1 and p2, in tria \a pt at parameter 0
 * <= \a s <= 1 from p1 result is stored in \a mr. edge p1p2 must not be a ridge
 */
int _MMG5_interpreg_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria pt,char i,double s,double mr[6]) {
  MMG5_pPoint    p1,p2;
  _MMG5_Bezier   b;
  double         b1[3],b2[3],bn[3],c[3],nt[3],cold[3],n[3],nold[3],mold[6],m1[6],m2[6];
  double        *n1,*n2,step,u,r[3][3],mt1[3],mt2[3],dd;
  int            ip1,ip2,nstep,l;
  char           i1,i2;

  /* Number of steps for parallel transport */
  nstep = 4;
  _MMG5_nortri(mesh,pt,nt);
  i1  = _MMG5_inxt2[i];
  i2  = _MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];

  if ( !_MMG5_bezierCP(mesh,pt,&b,1) )  return(0);

  n1 = &b.n[i1][0];
  n2 = &b.n[i2][0];
  memcpy(bn,&b.n[i+3][0],3*sizeof(double));
  memcpy(b1,&b.b[2*i+3][0],3*sizeof(double));
  memcpy(b2,&b.b[2*i+4][0],3*sizeof(double));

  /* Parallel transport of metric at p1 to point p(s) */
  step = s / nstep;
  cold[0] = p1->c[0];
  cold[1] = p1->c[1];
  cold[2] = p1->c[2];

  nold[0] = n1[0];
  nold[1] = n1[1];
  nold[2] = n1[2];

  if ( MG_SIN(p1->tag) || (p1->tag & MG_NOM)) {
    memcpy(m1,&met->m[6*ip1],6*sizeof(double));
  }
  else {
    if ( MG_GEO & p1->tag ) {
      if ( !_MMG5_buildridmetnor(mesh,met,pt->v[i1],nt,m1) )  return(0);
    }
    else {
      memcpy(m1,&met->m[6*ip1],6*sizeof(double));
    }
    memcpy(mold,m1,6*sizeof(double));

    /* Go from point (l-1)step, to point l step */
    for (l=1; l<=nstep; l++) {
      u    = l*step;
      c[0] = p1->c[0]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[0]\
        + 3.0*u*u*(1.0-u)*b2[0] + u*u*u*p2->c[0];
      c[1] = p1->c[1]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[1]\
        + 3.0*u*u*(1.0-u)*b2[1] + u*u*u*p2->c[1];
      c[2] = p1->c[2]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[2]\
        + 3.0*u*u*(1.0-u)*b2[2] + u*u*u*p2->c[2];

      n[0] = (1.0-u)*(1.0-u)*n1[0] + 2.0*u*(1.0-u)*bn[0] + u*u*n2[0];
      n[1] = (1.0-u)*(1.0-u)*n1[1] + 2.0*u*(1.0-u)*bn[1] + u*u*n2[1];
      n[2] = (1.0-u)*(1.0-u)*n1[2] + 2.0*u*(1.0-u)*bn[2] + u*u*n2[2];
      dd   = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
      if ( dd < _MMG5_EPSD )  return(0);
      dd = 1.0 / sqrt(dd);
      n[0] *= dd;
      n[1] *= dd;
      n[2] *= dd;

      if ( !_MMG5_paratmet(cold,nold,mold,c,n,m1) )  return(0);

      memcpy(cold,c,3*sizeof(double));
      memcpy(nold,n,3*sizeof(double));
      memcpy(mold,m1,6*sizeof(double));
    }
  }

  /* Parallel transport of metric at p2 to point p(s) */
  step = (1.0-s) / nstep;
  cold[0] = p2->c[0];
  cold[1] = p2->c[1];
  cold[2] = p2->c[2];

  nold[0] = n2[0];
  nold[1] = n2[1];
  nold[2] = n2[2];

  if ( MG_SIN(p2->tag) || (p2->tag & MG_NOM)) {
    memcpy(m2,&met->m[6*ip2],6*sizeof(double));

    /* In this pathological case, n is empty */
    if ( MG_SIN(p1->tag) || (p1->tag & MG_NOM))
      memcpy(n,n2,3*sizeof(double));
  }
  else {
    if ( p2->tag & MG_GEO ) {
      if ( !_MMG5_buildridmetnor(mesh,met,pt->v[i2],nt,m2))  return(0);
    }
    else {
      memcpy(m2,&met->m[6*ip2],6*sizeof(double));
    }
    memcpy(mold,m2,6*sizeof(double));

    /* Go from point (l-1)step, to point l step */
    for (l=1; l<=nstep; l++) {
      u    = 1.0 - l*step;
      c[0] = p1->c[0]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[0]\
        + 3.0*u*u*(1.0-u)*b2[0] + u*u*u*p2->c[0];
      c[1] = p1->c[1]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[1]\
        + 3.0*u*u*(1.0-u)*b2[1] + u*u*u*p2->c[1];
      c[2] = p1->c[2]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[2]\
        + 3.0*u*u*(1.0-u)*b2[2] + u*u*u*p2->c[2];

      n[0] = (1.0-u)*(1.0-u)*n1[0] + 2.0*u*(1.0-u)*bn[0] + u*u*n2[0];
      n[1] = (1.0-u)*(1.0-u)*n1[1] + 2.0*u*(1.0-u)*bn[1] + u*u*n2[1];
      n[2] = (1.0-u)*(1.0-u)*n1[2] + 2.0*u*(1.0-u)*bn[2] + u*u*n2[2];
      dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
      if ( dd < _MMG5_EPSD )  return(0);
      dd = 1.0 / sqrt(dd);
      n[0] *= dd;
      n[1] *= dd;
      n[2] *= dd;

      if ( !_MMG5_paratmet(cold,nold,mold,c,n,m2) )  return(0);

      memcpy(cold,c,3*sizeof(double));
      memcpy(nold,n,3*sizeof(double));
      memcpy(mold,m2,6*sizeof(double));
    }
  }
  /* At this point, c is point p(s), n is the normal at p(s), m1 and m2 are the 3*3
     transported metric tensors from p1 and p2 to p(s) */

  /* Rotate both matrices to the tangent plane */
  if ( !_MMG5_rotmatrix(n,r) )  return(0);
  _MMG5_rmtr(r,m1,mold);
  mt1[0] = mold[0];
  mt1[1] = mold[1];
  mt1[2] = mold[3];

  _MMG5_rmtr(r,m2,mold);
  mt2[0] = mold[0];
  mt2[1] = mold[1];
  mt2[2] = mold[3];

  /* Interpolate both metrics expressed in the same tangent plane : cold is reused */
  if ( !MMG5_intmet22(mt1,mt2,cold,s) ) {
    printf("Impossible interpolation between points : %d %d\n",pt->v[i1],pt->v[i2]);
    printf("m1 : %E %E %E %E %E %E \n",m1[0],m1[1],m1[2],m1[3],m1[4],m1[5]);
    printf("m2 : %E %E %E %E %E %E \n",m2[0],m2[1],m2[2],m2[3],m2[4],m2[5]);
    printf("mt1 : %E %E %E et det %E \n",mt1[0],mt1[1],mt1[2],mt1[0]*mt1[2]-mt1[1]*mt1[1]);
    printf("mt2 : %E %E %E et det : %E \n",mt2[0],mt2[1],mt2[2],mt2[0]*mt2[2]-mt2[1]*mt2[1]);
    exit(0);
    return(0);
  }

  /* End rotating back tangent metric into canonical basis of R^3 : mr = {^t}R*cold*R
     mt1 and mt2 serve for nothing ; let them be the lines of cold * R  */
  mt1[0] = cold[0]*r[0][0] + cold[1]*r[1][0];  mt1[1] = cold[0]*r[0][1] + cold[1]*r[1][1];   mt1[2] = cold[0]*r[0][2] + cold[1]*r[1][2] ;
  mt2[0] = cold[1]*r[0][0] + cold[2]*r[1][0];  mt2[1] = cold[1]*r[0][1] + cold[2]*r[1][1];   mt2[2] = cold[1]*r[0][2] + cold[2]*r[1][2] ;

  mr[0] = r[0][0] * mt1[0] + r[1][0] * mt2[0];
  mr[1] = r[0][0] * mt1[1] + r[1][0] * mt2[1];
  mr[2] = r[0][0] * mt1[2] + r[1][0] * mt2[2];
  mr[3] = r[0][1] * mt1[1] + r[1][1] * mt2[1];
  mr[4] = r[0][1] * mt1[2] + r[1][1] * mt2[2];
  mr[5] = r[0][2] * mt1[2] + r[1][2] * mt2[2];

//todo : warning false
//reappliquer la rotation aux deux metriques initiales
//extraire les vp : 2 devraient etre egales et la troisieme differente
//prendre le min des vp

  return(1);

}
