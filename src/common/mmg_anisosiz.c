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
 * \file common/mmg_anisosiz.c
 * \brief Fonctions for anisotropic size map computation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgs.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param np0 index of edge's extremity.
 * \param np1 index of edge's extremity.
 * \param isedg 1 if the edge is a ridge, 0 otherwise.
 * \return length of edge according to the prescribed metric.
 *
 * Compute length of edge \f$[i0;i1]\f$ according to the prescribed aniso.
 * metric.
 *
 */
double _MMG5_lenedg_ani(MMG5_pMesh mesh,MMG5_pSol met,int np0,int np1,char isedg) {
  MMG5_pPoint   p0,p1;
  double        gammaprim0[3],gammaprim1[3],t[3],*n1,*n2,ux,uy,uz,ps1,ps2,l0,l1;
  double        *m0,*m1,met0[6],met1[6];

  p0 = &mesh->point[np0];
  p1 = &mesh->point[np1];

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  /* computation of the two tangent vectors to the underlying curve of [i0i1] */
  if ( MS_SIN(p0->tag) ) {
    gammaprim0[0] = ux;
    gammaprim0[1] = uy;
    gammaprim0[2] = uz;
  }
  else if ( isedg ) {
    memcpy(t,p0->n,3*sizeof(double));
    ps1 = ux*t[0] + uy*t[1] + uz*t[2];
    gammaprim0[0] = ps1*t[0];
    gammaprim0[1] = ps1*t[1];
    gammaprim0[2] = ps1*t[2];
  }
  else {
    if ( MG_GEO & p0->tag ) {
      //assert(p0->xp);
      n1 = &mesh->xpoint[p0->xp].n1[0];
      n2 = &mesh->xpoint[p0->xp].n2[0];
      ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
      ps2 = ux*n2[0] + uy*n2[1] + uz*n2[2];

      if ( fabs(ps2) < fabs(ps1) ) {
        n1  = &mesh->xpoint[p0->xp].n2[0];
        ps1 = ps2;
      }
    }
    else if ( MG_REF & p0->tag ) {
      n1  = &mesh->xpoint[p0->xp].n1[0];
      ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
    }
    else {
      n1  = &(p0->n[0]);
      ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
    }
    gammaprim0[0] = ux - ps1*n1[0];
    gammaprim0[1] = uy - ps1*n1[1];
    gammaprim0[2] = uz - ps1*n1[2];
  }

  if ( MS_SIN(p1->tag) ) {
    gammaprim1[0] = -ux;
    gammaprim1[1] = -uy;
    gammaprim1[2] = -uz;
  }
  else if ( isedg ) {
    memcpy(t,p1->n,3*sizeof(double));
    ps1 = -ux*t[0] - uy*t[1] - uz*t[2];
    gammaprim1[0] = ps1*t[0];
    gammaprim1[1] = ps1*t[1];
    gammaprim1[2] = ps1*t[2];
  }
  else {
    if ( MG_GEO & p1->tag ) {
      n1 = &mesh->xpoint[p1->xp].n1[0];
      n2 = &mesh->xpoint[p1->xp].n2[0];
      ps1 = -ux*n1[0] - uy*n1[1] - uz*n1[2];
      ps2 = -ux*n2[0] - uy*n2[1] - uz*n2[2];

      if ( fabs(ps2) < fabs(ps1) ) {
        n1  = &mesh->xpoint[p1->xp].n2[0];
        ps1 = ps2;
      }
    }
    else if ( MG_REF & p1->tag ) {
      n1  = &mesh->xpoint[p1->xp].n1[0];
      ps1 = - ux*n1[0] - uy*n1[1] - uz*n1[2];
    }
    else {
      n1  = &(p1->n[0]);
      ps1 = -ux*n1[0] - uy*n1[1] - uz*n1[2];
    }
    gammaprim1[0] = - ux - ps1*n1[0];
    gammaprim1[1] = - uy - ps1*n1[1];
    gammaprim1[2] = - uz - ps1*n1[2];
  }

  /* Set metrics */
  if ( MS_SIN(p0->tag) ) {
    m0 = &met->m[6*np0];
  }
  else if ( MG_GEO & p0->tag ) {
    if ( !_MMG5_buildridmet(mesh,met,np0,ux,uy,uz,met0) )  return(-1.0);
    m0 = met0;
  }
  else {
    m0 = &met->m[6*np0];
  }

  if ( MS_SIN(p1->tag) ) {
    m1 = &met->m[6*np1];
  }
  else if ( MG_GEO & p1->tag ) {
    if ( !_MMG5_buildridmet(mesh,met,np1,ux,uy,uz,met1) )  return(-1.0);
    m1 = met1;
  }
  else {
    m1 = &met->m[6*np1];
  }

  /* computation of the length of the two tangent vectors in their respective tangent plane */
  l0 = m0[0]*gammaprim0[0]*gammaprim0[0] + m0[3]*gammaprim0[1]*gammaprim0[1] \
    + m0[5]*gammaprim0[2]*gammaprim0[2] \
    + 2.0*m0[1]*gammaprim0[0]*gammaprim0[1]  + 2.0*m0[2]*gammaprim0[0]*gammaprim0[2] \
    + 2.0*m0[4]*gammaprim0[1]*gammaprim0[2];

  l1 = m1[0]*gammaprim1[0]*gammaprim1[0] + m1[3]*gammaprim1[1]*gammaprim1[1] \
    + m1[5]*gammaprim1[2]*gammaprim1[2] \
    +2.0*m1[1]*gammaprim1[0]*gammaprim1[1]  + 2.0*m1[2]*gammaprim1[0]*gammaprim1[2] \
    + 2.0*m1[4]*gammaprim1[1]*gammaprim1[2];

  if(l0 < 0) {
    printf("neg %e\n",l0);
    l0 =1;
  }
  if(l1 < 0) {
    printf("neg1 %e\n",l1);
    l1 = 1;
  }
  l0 = 0.5*(sqrt(l0) + sqrt(l1));
  return(l0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the meric structure.
 * \param ptt pointer toward the triangle structure.
 * \return The computed area.
 *
 * Compute the area of the surface triangle \a ptt with respect to
 * the anisotropic metric \a met.
 *
 */
double _MMG5_surftri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) {
  MMG5_pPoint    p[3];
  _MMG5_Bezier    b;
  int       np[3];
  double    surf,ux,uy,uz,dens,m[3][6],J[3][2],mJ[3][2],tJmJ[2][2];
  char      i,i1,i2;

  surf = 0.0;

  for (i=0; i<3; i++) {
    np[i] = ptt->v[i];
    p[i]  = &mesh->point[np[i]];
  }
  if ( !_MMG5_bezierCP(mesh,ptt,&b,1) ) return(0.0);

  /* Set metric tensors at vertices of tria iel */
  for(i=0; i<3; i++) {
    i1 = _MMG5_inxt2[i];
    i2 = _MMG5_iprv2[i];
    ux = 0.5*(p[i1]->c[0]+p[i2]->c[0]) - p[i]->c[0];
    uy = 0.5*(p[i1]->c[1]+p[i2]->c[1]) - p[i]->c[1];
    uz = 0.5*(p[i1]->c[2]+p[i2]->c[2]) - p[i]->c[2];

    if ( MS_SIN(p[i]->tag) ) {
      memcpy(&m[i][0],&met->m[6*np[i]],6*sizeof(double));
    }
    else if ( p[i]->tag & MG_GEO ) {
      if ( !_MMG5_buildridmet(mesh,met,np[i],ux,uy,uz,&m[i][0]) )  return(0.0);
    }
    else {
      memcpy(&m[i][0],&met->m[6*np[i]],6*sizeof(double));
    }
  }

  /* Compute density integrand of volume at the 3 vertices of T */
  for (i=0; i<3; i++) {
    if ( i == 0 ) {
      J[0][0] = 3.0*( b.b[7][0] - b.b[0][0] ) ; J[0][1] = 3.0*( b.b[6][0] - b.b[0][0] );
      J[1][0] = 3.0*( b.b[7][1] - b.b[0][1] ) ; J[1][1] = 3.0*( b.b[6][1] - b.b[0][1] );
      J[2][0] = 3.0*( b.b[7][2] - b.b[0][2] ) ; J[2][1] = 3.0*( b.b[6][2] - b.b[0][2] );
    }
    else if ( i == 1 ) {
      J[0][0] = 3.0*( b.b[1][0] - b.b[8][0] ) ; J[0][1] = 3.0*( b.b[3][0] - b.b[8][0] );
      J[1][0] = 3.0*( b.b[1][1] - b.b[8][1] ) ; J[1][1] = 3.0*( b.b[3][1] - b.b[8][1] );
      J[2][0] = 3.0*( b.b[1][2] - b.b[8][2] ) ; J[2][1] = 3.0*( b.b[3][2] - b.b[8][2] );
    }
    else {
      J[0][0] = 3.0*( b.b[4][0] - b.b[5][0] ) ; J[0][1] = 3.0*( b.b[2][0] - b.b[5][0] );
      J[1][0] = 3.0*( b.b[4][1] - b.b[5][1] ) ; J[1][1] = 3.0*( b.b[2][1] - b.b[5][1] );
      J[2][0] = 3.0*( b.b[4][2] - b.b[5][2] ) ; J[2][1] = 3.0*( b.b[2][2] - b.b[5][2] );
    }

    mJ[0][0] = m[i][0]*J[0][0] + m[i][1]*J[1][0] + m[i][2]*J[2][0];
    mJ[1][0] = m[i][1]*J[0][0] + m[i][3]*J[1][0] + m[i][4]*J[2][0];
    mJ[2][0] = m[i][2]*J[0][0] + m[i][4]*J[1][0] + m[i][5]*J[2][0];

    mJ[0][1] = m[i][0]*J[0][1] + m[i][1]*J[1][1] + m[i][2]*J[2][1];
    mJ[1][1] = m[i][1]*J[0][1] + m[i][3]*J[1][1] + m[i][4]*J[2][1];
    mJ[2][1] = m[i][2]*J[0][1] + m[i][4]*J[1][1] + m[i][5]*J[2][1];

    /* dens = sqrt(tJacsigma * M * Jacsigma )*/
    tJmJ[0][0] = J[0][0]*mJ[0][0] + J[1][0]*mJ[1][0] + J[2][0]*mJ[2][0];
    tJmJ[0][1] = J[0][0]*mJ[0][1] + J[1][0]*mJ[1][1] + J[2][0]*mJ[2][1];
    tJmJ[1][0] = J[0][1]*mJ[0][0] + J[1][1]*mJ[1][0] + J[2][1]*mJ[2][0];
    tJmJ[1][1] = J[0][1]*mJ[0][1] + J[1][1]*mJ[1][1] + J[2][1]*mJ[2][1];

    dens = tJmJ[0][0]*tJmJ[1][1] - tJmJ[1][0]*tJmJ[0][1];
    if ( dens < 0.0 ) {
      //fprintf(stdout,"  ## Density should be positive : %E for elt %d %d %d \n",dens,ptt->v[0],ptt->v[1],ptt->v[2]);
      //return(0.0);
    }
    surf += sqrt(fabs(dens));
  }

  surf *= _MMG5_ATHIRD;
  return(surf);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ismet has the user provided a metric?
 *
 * Search for points with unintialized metric and define anisotropic size at
 * this points.
 *
 */
void _MMG5_defUninitSize(MMG5_pMesh mesh,MMG5_pSol met,char ismet)
{
  MMG5_pPoint   ppt;
  double        *m,*n,r[3][3],isqhmax;
  int           k;

  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) || ppt->flag == 1 )  continue;

    m = &met->m[6*k];
    if(ismet) {
      ppt->flag = 1;
      continue;
    }
    memset(m,0,6*sizeof(double));
    if ( MS_SIN(ppt->tag) ) {
      m[0] = m[3] = m[5] = isqhmax;
    }
    else if ( ppt->tag & MG_GEO ) {
      /* We store the size in the tangent dir in m[0], in the n1 dir in m[1] and
       * in the n2 dir in m[2]. */
      m[0] = m[1] = m[2] = isqhmax;
    }
    else {
      n = ppt->tag & MG_REF ? &mesh->xpoint[ppt->xp].n1[0] : ppt->n;
      _MMG5_rotmatrix(n,r);
      m[0] = isqhmax*(r[0][0]*r[0][0]+r[1][0]*r[1][0]+r[2][0]*r[2][0]);
      m[1] = isqhmax*(r[0][0]*r[0][1]+r[1][0]*r[1][1]+r[2][0]*r[2][1]);
      m[2] = isqhmax*(r[0][0]*r[0][2]+r[1][0]*r[1][2]+r[2][0]*r[2][2]);
      m[3] = isqhmax*(r[0][1]*r[0][1]+r[1][1]*r[1][1]+r[2][1]*r[2][1]);
      m[4] = isqhmax*(r[0][1]*r[0][2]+r[1][1]*r[1][2]+r[2][1]*r[2][2]);
      m[5] = isqhmax*(r[0][2]*r[0][2]+r[1][2]*r[1][2]+r[2][2]*r[2][2]);
    }
    ppt->flag = 1;
  }
}

/**
 * \param k index of the tetrahedra from which we come.
 * \param p0 pointer toward the point on which we want to def the metric.
 * \param i0 pointer toward the local index of the point in tria.
 * \param b control polygon of triangle.
 * \param r rotation matrix.
 * \param c physical coordinates of the curve edge mid-point.
 * \param lispoi list of incident vertices to p0
 * \param tAA matrix to fill
 * \param tAb second member
 *
 * Fill matrice tAA and second member tAb with A=(\sum X_{P_i}^2 \sum
 * Y_{P_i}^2 \sum X_{P_i}Y_{P_i}) and b=\sum Z_{P_i} with P_i the physical
 * points at edge [i0;i1] extremities and middle.
 * Compute the physical coor \a c of the curve edge's
 * mid-point for a regular or reference point.
 *
 */
void _MMG5_fillDefmetregSys( int k, MMG5_pPoint p0, int i0, _MMG5_Bezier b,
                             double r[3][3], double c[3],
                             double *lispoi,
                             double tAA[6], double tAb[3] )
{
  double b0[3],b1[3],d[3];
  int    j;

  for(j=0; j<10; j++){
    c[0] = b.b[j][0] - p0->c[0];
    c[1] = b.b[j][1] - p0->c[1];
    c[2] = b.b[j][2] - p0->c[2];

    b.b[j][0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
    b.b[j][1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
    b.b[j][2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];
  }

/* Mid-point along edge [i0;i1] and endpoint in the rotated frame */
  if(i0 == 0){
    memcpy(b0,&(b.b[7][0]),3*sizeof(double));
    memcpy(b1,&(b.b[8][0]),3*sizeof(double));
  }
  else if(i0 == 1){
    memcpy(b0,&(b.b[3][0]),3*sizeof(double));
    memcpy(b1,&(b.b[4][0]),3*sizeof(double));
  }
  else{
    memcpy(b0,&(b.b[5][0]),3*sizeof(double));
    memcpy(b1,&(b.b[6][0]),3*sizeof(double));
  }

/* At this point, the two control points b0 and b1 are expressed in the
 * rotated frame. We compute the physical coor of the curve edge's
 * mid-point. */
  c[0] = 3.0/8.0*b0[0] + 3.0/8.0*b1[0] + 1.0/8.0*lispoi[3*k+1];
  c[1] = 3.0/8.0*b0[1] + 3.0/8.0*b1[1] + 1.0/8.0*lispoi[3*k+2];
  c[2] = 3.0/8.0*b0[2] + 3.0/8.0*b1[2] + 1.0/8.0*lispoi[3*k+3];

/* Fill matrice tAA and second member tAb with A=(\sum X_{P_i}^2 \sum
 * Y_{P_i}^2 \sum X_{P_i}Y_{P_i}) and b=\sum Z_{P_i} with P_i the physical
 * points at edge [i0;i1] extremities and middle. */
  tAA[0] += c[0]*c[0]*c[0]*c[0];
  tAA[1] += c[0]*c[0]*c[1]*c[1];
  tAA[2] += c[0]*c[0]*c[0]*c[1];
  tAA[3] += c[1]*c[1]*c[1]*c[1];
  tAA[4] += c[0]*c[1]*c[1]*c[1];
  tAA[5] += c[0]*c[0]*c[1]*c[1];

  tAb[0] += c[0]*c[0]*c[2];
  tAb[1] += c[1]*c[1]*c[2];
  tAb[2] += c[0]*c[1]*c[2];

  tAA[0] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1];
  tAA[1] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2];
  tAA[2] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2];
  tAA[3] += lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2];
  tAA[4] += lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2];
  tAA[5] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2];

  tAb[0] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+3];
  tAb[1] += lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+3];
  tAb[2] += lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+3];

/* Mid-point along median edge (coor of mid-point in parametric space : (1/4
 * 1/4)) and endpoint in the rotated frame (coor of end-point in parametric
 * space (1/2 1/2)). */
  if(i0 == 0){
    c[0] = A64TH*(b.b[1][0] + b.b[2][0] + 3.0*(b.b[3][0] + b.b[4][0])) \
      + 3.0*A16TH*(b.b[6][0] + b.b[7][0] + b.b[9][0]) + A32TH*(b.b[5][0] + b.b[8][0]);
    c[1] = A64TH*(b.b[1][1] + b.b[2][1] + 3.0*(b.b[3][1] + b.b[4][1])) \
      + 3.0*A16TH*(b.b[6][1] + b.b[7][1] + b.b[9][1]) + A32TH*(b.b[5][1] + b.b[8][1]);
    c[2] = A64TH*(b.b[1][2] + b.b[2][2] + 3.0*(b.b[3][2] + b.b[4][2])) \
      + 3.0*A16TH*(b.b[6][2] + b.b[7][2] + b.b[9][2]) + A32TH*(b.b[5][2] + b.b[8][2]);

    d[0] = 0.125*b.b[1][0] + 0.375*(b.b[3][0] + b.b[4][0]) + 0.125*b.b[2][0];
    d[1] = 0.125*b.b[1][1] + 0.375*(b.b[3][1] + b.b[4][1]) + 0.125*b.b[2][1];
    d[2] = 0.125*b.b[1][2] + 0.375*(b.b[3][2] + b.b[4][2]) + 0.125*b.b[2][2];
  }
  else if(i0 == 1){
    c[0] = A64TH*(b.b[0][0] + b.b[2][0] + 3.0*(b.b[5][0] + b.b[6][0])) \
      + 3.0*A16TH*(b.b[3][0] + b.b[8][0] + b.b[9][0]) + A32TH*(b.b[4][0] + b.b[7][0]);
    c[1] = A64TH*(b.b[0][1] + b.b[2][1] + 3.0*(b.b[5][1] + b.b[6][1])) \
      + 3.0*A16TH*(b.b[3][1] + b.b[8][1] + b.b[9][1]) + A32TH*(b.b[4][1] + b.b[7][1]);
    c[2] = A64TH*(b.b[0][2] + b.b[2][2] + 3.0*(b.b[5][2] + b.b[6][2])) \
      + 3.0*A16TH*(b.b[3][2] + b.b[8][2] + b.b[9][2]) + A32TH*(b.b[4][2] + b.b[7][2]);

    d[0] = 0.125*b.b[2][0] + 0.375*(b.b[5][0] + b.b[6][0]) + 0.125*b.b[0][0];
    d[1] = 0.125*b.b[2][1] + 0.375*(b.b[5][1] + b.b[6][1]) + 0.125*b.b[0][1];
    d[2] = 0.125*b.b[2][2] + 0.375*(b.b[5][2] + b.b[6][2]) + 0.125*b.b[0][2];
  }
  else{
    c[0] = A64TH*(b.b[0][0] + b.b[1][0] + 3.0*(b.b[7][0] + b.b[8][0])) \
      + 3.0*A16TH*(b.b[4][0] + b.b[5][0] + b.b[9][0]) + A32TH*(b.b[3][0] + b.b[6][0]);
    c[1] = A64TH*(b.b[0][1] + b.b[1][1] + 3.0*(b.b[7][1] + b.b[8][1])) \
      + 3.0*A16TH*(b.b[4][1] + b.b[5][1] + b.b[9][1]) + A32TH*(b.b[3][1] + b.b[6][1]);
    c[2] = A64TH*(b.b[0][2] + b.b[1][2] + 3.0*(b.b[7][2] + b.b[8][2])) \
      + 3.0*A16TH*(b.b[4][2] + b.b[5][2] + b.b[9][2]) + A32TH*(b.b[3][2] + b.b[6][2]);

    d[0] = 0.125*b.b[0][0] + 0.375*(b.b[7][0] + b.b[8][0]) + 0.125*b.b[1][0];
    d[1] = 0.125*b.b[0][1] + 0.375*(b.b[7][1] + b.b[8][1]) + 0.125*b.b[1][1];
    d[2] = 0.125*b.b[0][2] + 0.375*(b.b[7][2] + b.b[8][2]) + 0.125*b.b[1][2];
  }

/* Fill matrice tAA and second member tAb*/
  tAA[0] += c[0]*c[0]*c[0]*c[0];
  tAA[1] += c[0]*c[0]*c[1]*c[1];
  tAA[2] += c[0]*c[0]*c[0]*c[1];
  tAA[3] += c[1]*c[1]*c[1]*c[1];
  tAA[4] += c[0]*c[1]*c[1]*c[1];
  tAA[5] += c[0]*c[0]*c[1]*c[1];

  tAb[0] += c[0]*c[0]*c[2];
  tAb[1] += c[1]*c[1]*c[2];
  tAb[2] += c[0]*c[1]*c[2];

  tAA[0] += d[0]*d[0]*d[0]*d[0];
  tAA[1] += d[0]*d[0]*d[1]*d[1];
  tAA[2] += d[0]*d[0]*d[0]*d[1];
  tAA[3] += d[1]*d[1]*d[1]*d[1];
  tAA[4] += d[0]*d[1]*d[1]*d[1];
  tAA[5] += d[0]*d[0]*d[1]*d[1];

  tAb[0] += d[0]*d[0]*d[2];
  tAb[1] += d[1]*d[1]*d[2];
  tAb[2] += d[0]*d[1]*d[2];

  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param r pointer toward the rotation matrix.
 * \param c physical coordinates of the curve edge mid-point.
 * \param tAA matrix of the system to solve.
 * \param tAb second member.
 * \param m pointer toward the metric.
 * \param isqhmax maximum size for edge.
 * \param isqhmin minimum size for edge.
 * \param hausd hausdorff value at point.
 * \return 1 if success, 0 if fail.
 *
 * Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) for a regular
 * point.
 *
 */
int _MMG5_solveDefmetregSys( MMG5_pMesh mesh, double r[3][3], double c[3],
                             double tAA[6], double tAb[3], double *m,
                             double isqhmin, double isqhmax, double hausd)
{
  double intm[3], kappa[2], vp[2][2], b0[3], b1[3], b2[3];

  memset(intm,0.0,3*sizeof(double));

  /* case planar surface : tAb = 0 => no curvature */
  /* isotropic metric with hmax size*/

  if((tAb[0]*tAb[0] + tAb[1]*tAb[1] + tAb[2]*tAb[2]) < _MMG5_EPSD) {
    m[0] = isqhmax;
    m[1] = 0;
    m[2] = 0;
    m[3] = isqhmax;
    m[4] = 0;
    m[5] = isqhmax;
    return(1);
  }

  /* solve now (a b c) = tAA^{-1} * tAb */
  if ( !_MMG5_sys33sym(tAA,tAb,c) ) {
    printf(" La matrice %f %f %f %f %f %f \n",tAA[0],tAA[1],tAA[2],tAA[3],tAA[4],tAA[5]);
    return(0);
  }
  intm[0] = 2.0*c[0];
  intm[1] = c[2];
  intm[2] = 2.0*c[1];

  /* At this point, intm stands for the integral matrix of Taubin's approach : vp[0] and vp[1]
     are the two pr. directions of curvature, and the two curvatures can be inferred from lambdas*/
  _MMG5_eigensym(intm,kappa,vp);

  /* Truncation of eigenvalues */
  kappa[0] = 2.0/9.0 * fabs(kappa[0])/hausd;
  kappa[0] = MG_MIN(kappa[0],isqhmin);
  kappa[0] = MG_MAX(kappa[0],isqhmax);

  kappa[1] = 2.0/9.0 * fabs(kappa[1])/hausd;
  kappa[1] = MG_MIN(kappa[1],isqhmin);
  kappa[1] = MG_MAX(kappa[1],isqhmax);

  /* Send back the metric to the canonical basis of tangent plane :
     diag(lambda) = {^t}vp * M * vp, M = vp * diag(lambda) * {^t}vp */
  intm[0] = kappa[0]*vp[0][0]*vp[0][0] + kappa[1]*vp[1][0]*vp[1][0];
  intm[1] = kappa[0]*vp[0][0]*vp[0][1] + kappa[1]*vp[1][0]*vp[1][1];
  intm[2] = kappa[0]*vp[0][1]*vp[0][1] + kappa[1]*vp[1][1]*vp[1][1];

  /* At this point, intm (with 0 replaced by isqhmax in the z direction) is the
     desired metric, except it is expressed in the rotated bc, that is intm = R
     * metric in bc * ^t R, so metric in bc = ^tR*intm*R */

  /* b0, b1 and b2 are the lines of matrix intm*R  */
  // intm = intm[0]  intm[1]    0
  //        intm[1]  intm[2]    0
  //           0       0     isqhmax

  b0[0] = intm[0]*r[0][0] + intm[1]*r[1][0];
  b0[1] = intm[0]*r[0][1] + intm[1]*r[1][1];
  b0[2] = intm[0]*r[0][2] + intm[1]*r[1][2];
  b1[0] = intm[1]*r[0][0] + intm[2]*r[1][0];
  b1[1] = intm[1]*r[0][1] + intm[2]*r[1][1];
  b1[2] = intm[1]*r[0][2] + intm[2]*r[1][2];
  b2[0] = isqhmax*r[2][0];
  b2[1] = isqhmax*r[2][1];
  b2[2] = isqhmax*r[2][2];

  m[0] = r[0][0] * b0[0] + r[1][0] * b1[0] + r[2][0] * b2[0];
  m[1] = r[0][0] * b0[1] + r[1][0] * b1[1] + r[2][0] * b2[1];
  m[2] = r[0][0] * b0[2] + r[1][0] * b1[2] + r[2][0] * b2[2];

  m[3] = r[0][1] * b0[1] + r[1][1] * b1[1] + r[2][1] * b2[1];
  m[4] = r[0][1] * b0[2] + r[1][1] * b1[2] + r[2][1] * b2[2];

  m[5] = r[0][2] * b0[2] + r[1][2] * b1[2] + r[2][2] * b2[2];

  /* Security check : normal in the kernel */
  /*if((fabs(p0->m[0]*n[0] + p0->m[1]*n[1] + p0->m[2]*n[2] ) > 1.0e-4)){
    printf("VALEUR ETRANGE... %f \n",fabs(p0->m[0]*n[0] + p0->m[1]*n[1] + p0->m[2]*n[2] ));
    }
    if((fabs(p0->m[1]*n[0] + p0->m[3]*n[1] + p0->m[4]*n[2] ) > 1.0e-4)){
    printf("VALEUR ETRANGE... %f \n",fabs(p0->m[1]*n[0] + p0->m[3]*n[1] + p0->m[4]*n[2] ));
    }

    if((fabs(p0->m[2]*n[0] + p0->m[4]*n[1] + p0->m[5]*n[2] ) > 1.0e-4)){
    printf("VALEUR ETRANGE... %f \n",fabs(p0->m[2]*n[0] + p0->m[4]*n[1] + p0->m[5]*n[2]));
    } */

  /*if(ddb) {
    printf("La matrice %f %f %f\n",p0->m[0],p0->m[1],p0->m[2]);
    printf("            %f %f %f\n",p0->m[1],p0->m[3],p0->m[4]);
    printf("            %f %f %f\n",p0->m[2],p0->m[4],p0->m[5]);

    }*/
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param p0 pointer toward the point on which we want to define the metric.
 * \param ipref table containing the indices of the edge extremities.
 * \param r pointer toward the rotation matrix.
 * \param c physical coordinates of the curve edge mid-point.
 * \param tAA matrix of the system to solve.
 * \param tAb second member.
 * \param m pointer toward the metric.
 * \param isqhmax maximum size for edge.
 * \param isqhmin minimum size for edge.
 * \param hausd hausdorff value at point.
 * \return 1 if success, 0 if fail.
 *
 * Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) for a ref
 * point.
 *
 */
int _MMG5_solveDefmetrefSys( MMG5_pMesh mesh, MMG5_pPoint p0, int ipref[2],
                             double r[3][3], double c[3],
                             double tAA[6], double tAb[3], double *m,
                             double isqhmin, double isqhmax, double hausd)
{
  MMG5_pPoint  p1;
  double       intm[3], kappa[2], vp[2][2], b0[3], b1[3], b2[3], kappacur;
  double       gammasec[3],tau[2], ux, uy, uz, ps1, l, ll, *t, *t1;
  int          i;

  memset(intm,0.0,3*sizeof(double));

  /* case planar surface : tAb = 0 => no curvature */
  /* isotropic metric with hmax size*/
  if((tAb[0]*tAb[0] + tAb[1]*tAb[1] + tAb[2]*tAb[2]) < _MMG5_EPSD) {
    m[0] = isqhmax;
    m[1] = 0;
    m[2] = 0;
    m[3] = isqhmax;
    m[4] = 0;
    m[5] = isqhmax;
    return(1);
  }

  /* solve now (a b c) = tAA^{-1} * tAb */
  if ( !_MMG5_sys33sym(tAA,tAb,c) )  return(0);

  intm[0] = 2.0*c[0];
  intm[1] = c[2];
  intm[2] = 2.0*c[1];

  /* At this point, intm stands for the integral matrix of Taubin's approach :
     vp[0] and vp[1] are the two pr. directions of curvature, and the two
     curvatures can be inferred from lambdas*/
  _MMG5_eigensym(intm,kappa,vp);

  /* Truncation of eigenvalues */
  kappa[0] = 2.0/9.0 * fabs(kappa[0])/mesh->info.hausd;
  kappa[0] = MG_MIN(kappa[0],isqhmin);
  kappa[0] = MG_MAX(kappa[0],isqhmax);

  kappa[1] = 2.0/9.0 * fabs(kappa[1])/mesh->info.hausd;
  kappa[1] = MG_MIN(kappa[1],isqhmin);
  kappa[1] = MG_MAX(kappa[1],isqhmax);

  /* Send back the metric to the canonical basis of tangent plane :
     diag(lambda) = {^t}vp * M * vp, M = vp * diag(lambda) * {^t}vp */
  intm[0] = kappa[0]*vp[0][0]*vp[0][0] + kappa[1]*vp[1][0]*vp[1][0];
  intm[1] = kappa[0]*vp[0][0]*vp[0][1] + kappa[1]*vp[1][0]*vp[1][1];
  intm[2] = kappa[0]*vp[0][1]*vp[0][1] + kappa[1]*vp[1][1]*vp[1][1];

  /* Now express metric with respect to the approx of the underlying ref
   * curve */
  t = &p0->n[0];
  kappacur = 0.0;

  for (i=0; i<2; i++) {
    p1 = &mesh->point[ipref[i]];
    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    ps1 =  ux*t[0] + uy*t[1] + uz*t[2];
    c[0] = _MMG5_ATHIRD*ps1*t[0];
    c[1] = _MMG5_ATHIRD*ps1*t[1];
    c[2] = _MMG5_ATHIRD*ps1*t[2];

    b0[0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
    b0[1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
    b0[2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];

    if ( (MG_CRN & p1->tag) || (MG_NOM & p1->tag) ) {
      c[0] = p1->c[0] - _MMG5_ATHIRD*ux;
      c[1] = p1->c[1] - _MMG5_ATHIRD*uy;
      c[2] = p1->c[2] - _MMG5_ATHIRD*uz;
    }
    else {
      assert(MG_REF & p1->tag);
      t1 = &(p1->n[0]);
      ps1 =  -(ux*t1[0] + uy*t1[1] + uz*t1[2]);
      c[0] = p1->c[0] + _MMG5_ATHIRD*ps1*t1[0];
      c[1] = p1->c[1] + _MMG5_ATHIRD*ps1*t1[1];
      c[2] = p1->c[2] + _MMG5_ATHIRD*ps1*t1[2];
    }
    c[0] -= p0->c[0];
    c[1] -= p0->c[1];
    c[2] -= p0->c[2];

    b1[0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
    b1[1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
    b1[2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];

    /* Everything is expressed in the rotated frame */
    tau[0] = 3.0*b0[0];
    tau[1] = 3.0*b0[1];
    ll = tau[0]*tau[0] + tau[1]*tau[1];
    if ( ll < _MMG5_EPSD ) {
      kappacur = isqhmax;
      continue;
    }
    l = 1.0 / sqrt(ll);
    tau[0] *= l;
    tau[1] *= l;

    gammasec[0] = -12.0*b0[0] + 6.0*b1[0];
    gammasec[1] = -12.0*b0[1] + 6.0*b1[1];
    gammasec[2] = -12.0*b0[2] + 6.0*b1[2];

    ps1 = tau[0]*gammasec[0] + tau[1]*gammasec[1];
    c[0] = gammasec[0] - ps1*tau[0];
    c[1] = gammasec[1] - ps1*tau[1];
    c[2] = gammasec[2];

    // p.s. with normal at p0
    kappacur = MG_MAX(kappacur,MG_MAX(0.0,1.0/ll*fabs(c[2])));
  }

  /* Rotation of tangent vector : tau is reused */
  c[0] =  r[0][0]*t[0] + r[0][1]*t[1] + r[0][2]*t[2];
  c[1] =  r[1][0]*t[0] + r[1][1]*t[1] + r[1][2]*t[2];
  c[2] =  r[2][0]*t[0] + r[2][1]*t[1] + r[2][2]*t[2];
  memcpy(tau,&c[0],2*sizeof(double));

  /* Truncation of curvature */
  kappacur = 1.0/8.0 * kappacur/mesh->info.hausd;
  kappacur = MG_MIN(kappacur,isqhmin);
  kappacur = MG_MAX(kappacur,isqhmax);

  /* The associated matrix in basis (rt, orth rt) */
  c[0] = kappacur*tau[0]*tau[0] + isqhmax*tau[1]*tau[1];
  c[1] = (kappacur - isqhmax)*tau[0]*tau[1];
  c[2] = kappacur*tau[1]*tau[1] + isqhmax*tau[0]*tau[0];

  /* Reuse b0 for commodity */
  _MMG5_intmetsavedir(mesh,c,intm,b0);
  memcpy(intm,b0,3*sizeof(double));

  /* At this point, intm (with 0 replaced by isqhmax in the z direction) is the
     desired metric, except it is expressed in the rotated bc, that is intm = R
     * metric in bc * ^t R, so metric in bc = ^tR*intm*R */
  // intm = intm[0]  intm[1]    0
  //        intm[1]  intm[2]    0
  //           0       0     isqhmax

  /* b0 and b1 serve now for nothing : let them be the lines of matrix intm*R */
  b0[0] = intm[0]*r[0][0] + intm[1]*r[1][0];
  b0[1] = intm[0]*r[0][1] + intm[1]*r[1][1];
  b0[2] = intm[0]*r[0][2] + intm[1]*r[1][2];
  b1[0] = intm[1]*r[0][0] + intm[2]*r[1][0];
  b1[1] = intm[1]*r[0][1] + intm[2]*r[1][1];
  b1[2] = intm[1]*r[0][2] + intm[2]*r[1][2];
  b2[0] = isqhmax*r[2][0];
  b2[1] = isqhmax*r[2][1];
  b2[2] = isqhmax*r[2][2];

  m[0] = r[0][0] * b0[0] + r[1][0] * b1[0] + r[2][0] * b2[0];
  m[1] = r[0][0] * b0[1] + r[1][0] * b1[1] + r[2][0] * b2[1];
  m[2] = r[0][0] * b0[2] + r[1][0] * b1[2] + r[2][0] * b2[2];

  m[3] = r[0][1] * b0[1] + r[1][1] * b1[1] + r[2][1] * b2[1];
  m[4] = r[0][1] * b0[2] + r[1][1] * b1[2] + r[2][1] * b2[2];

  m[5] = r[0][2] * b0[2] + r[1][2] * b1[2] + r[2][2] * b2[2];

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param p0 pointer toward the point at which we define the metric.
 * \param idp global index of the point at which we define the metric.
 * \param iprid pointer toward the two extremities of the ridge.
 * \param isqhmin minimum edge size.
 * \param isqhmax maximum edge size.
 * \return the computed ridge size in the tangent direction.
 *
 * Compute the specific size of a ridge in the direction of the tangent of the
 * ridge.
 *
 **/
double _MMG5_ridSizeInTangentDir(MMG5_pMesh mesh, MMG5_pPoint p0, int idp,
                                 int* iprid, double isqhmin,double isqhmax)
{
  int    i;
  double n0[3],tau[3],gammasec[3],c[3],ps,ll,l,m;
  double b0[3],b1[3],kappacur;

  m = isqhmax;
  for (i=0; i<2; i++) {
    kappacur = 0.0;
    // Remark: bezierEdge don't use n0 in case of a ridge so it's ok to call it
    // with an undefined n0.
    _MMG5_bezierEdge(mesh,idp,iprid[i],b0,b1,1,n0);

    /* tau is the bezier curve derivative in p0 (parametric coor t=0) */
    tau[0] = 3.0*(b0[0] - p0->c[0]);
    tau[1] = 3.0*(b0[1] - p0->c[1]);
    tau[2] = 3.0*(b0[2] - p0->c[2]);
    ll = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];
    if ( ll < _MMG5_EPSD )  continue;
    l = 1.0 / sqrt(ll);
    tau[0] *= l;
    tau[1] *= l;
    tau[2] *= l;

    gammasec[0] = 6.0*p0->c[0] -12.0*b0[0] + 6.0*b1[0];
    gammasec[1] = 6.0*p0->c[1] -12.0*b0[1] + 6.0*b1[1];
    gammasec[2] = 6.0*p0->c[2] -12.0*b0[2] + 6.0*b1[2];

    ps = tau[0]*gammasec[0] + tau[1]*gammasec[1] + tau[2]*gammasec[2];
    c[0] = gammasec[0] - ps*tau[0];
    c[1] = gammasec[1] - ps*tau[1];
    c[2] = gammasec[2] - ps*tau[2];

    kappacur = MG_MAX(0.0,1.0/ll*sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]));
    kappacur = 1.0/8.0*kappacur/mesh->info.hausd;
    kappacur = MG_MIN(kappacur,isqhmin);
    kappacur = MG_MAX(kappacur,isqhmax);
    m = MG_MAX(m,kappacur);
  }
  return(m);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param i0 local index in the face of the point on which we want to compute
 * the metric
 * \param bcu pointer toward the barycentric coordinates of vector \a u in the
 * looked face.
 * \param b bezier control polygon for the looked face.
 * \param isqhmin minimum edge size.
 * \param isqhmax maximum edge size.
 * \return the computed ridge size in first or second normal direction
 * (depending of i0).
 *
 * Compute the specific size of a ridge in the direction of the normal of the
 * looked face.
 *
 **/
double _MMG5_ridSizeInNormalDir(MMG5_pMesh mesh,int i0,double* bcu,
                                _MMG5_Bezier *b,int isqhmin,int isqhmax)
{
  double lambda[2],Jacb[3][2],Hb[3][3],tau[3],ll,l,gammasec[3],c[3];
  double ps,kappacur;

  if ( i0 == 0 ) { // w = 1, u,v = 0
    lambda[0] = bcu[1];
    lambda[1] = bcu[2];

    Jacb[0][0] = 3.0*(b->b[7][0]-b->b[0][0]);
    Jacb[1][0] = 3.0*(b->b[7][1]-b->b[0][1]);
    Jacb[2][0] = 3.0*(b->b[7][2]-b->b[0][2]);

    Jacb[0][1] = 3.0*(b->b[6][0]-b->b[0][0]);
    Jacb[1][1] = 3.0*(b->b[6][1]-b->b[0][1]);
    Jacb[2][1] = 3.0*(b->b[6][2]-b->b[0][2]);

    /* Hb[i] = hessian matrix of i-th component of b at point p0 */
    Hb[0][0] = 6.0*(b->b[0][0] - 2.0*b->b[7][0] + b->b[8][0]);
    Hb[1][0] = 6.0*(b->b[0][1] - 2.0*b->b[7][1] + b->b[8][1]);
    Hb[2][0] = 6.0*(b->b[0][2] - 2.0*b->b[7][2] + b->b[8][2]);

    Hb[0][1] = 6.0*(b->b[0][0] - b->b[7][0] - b->b[6][0] + b->b[9][0]);
    Hb[1][1] = 6.0*(b->b[0][1] - b->b[7][1] - b->b[6][1] + b->b[9][1]);
    Hb[2][1] = 6.0*(b->b[0][2] - b->b[7][2] - b->b[6][2] + b->b[9][2]);

    Hb[0][2] = 6.0*(b->b[0][0] - 2.0*b->b[6][0] + b->b[5][0]);
    Hb[1][2] = 6.0*(b->b[0][1] - 2.0*b->b[6][1] + b->b[5][1]);
    Hb[2][2] = 6.0*(b->b[0][2] - 2.0*b->b[6][2] + b->b[5][2]);
  }
  else if ( i0 == 1 ) {  //w = v = 0, u = 1;
    lambda[0] = bcu[0];
    lambda[1] = bcu[1];

    Jacb[0][0] = 3.0*(b->b[1][0]-b->b[8][0]);
    Jacb[1][0] = 3.0*(b->b[1][1]-b->b[8][1]);
    Jacb[2][0] = 3.0*(b->b[1][2]-b->b[8][2]);

    Jacb[0][1] = 3.0*(b->b[3][0]-b->b[8][0]);
    Jacb[1][1] = 3.0*(b->b[3][1]-b->b[8][1]);
    Jacb[2][1] = 3.0*(b->b[3][2]-b->b[8][2]);

    Hb[0][0] = 6.0*(b->b[1][0] - 2.0*b->b[8][0] + b->b[7][0]);
    Hb[1][0] = 6.0*(b->b[1][1] - 2.0*b->b[8][1] + b->b[7][1]);
    Hb[2][0] = 6.0*(b->b[1][2] - 2.0*b->b[8][2] + b->b[7][2]);

    Hb[0][1] = 6.0*(b->b[7][0] - b->b[8][0] - b->b[9][0] + b->b[3][0]);
    Hb[1][1] = 6.0*(b->b[7][1] - b->b[8][1] - b->b[9][1] + b->b[3][1]);
    Hb[2][1] = 6.0*(b->b[7][2] - b->b[8][2] - b->b[9][2] + b->b[3][2]);

    Hb[0][2] = 6.0*(b->b[4][0] - 2.0*b->b[9][0] + b->b[7][0]);
    Hb[1][2] = 6.0*(b->b[4][1] - 2.0*b->b[9][1] + b->b[7][1]);
    Hb[2][2] = 6.0*(b->b[4][2] - 2.0*b->b[9][2] + b->b[7][2]);
  }
  else {   //w =u = 0, v =1
    lambda[0] = bcu[2];
    lambda[1] = bcu[0];

    Jacb[0][0] = 3.0*(b->b[4][0]-b->b[5][0]);
    Jacb[1][0] = 3.0*(b->b[4][1]-b->b[5][1]);
    Jacb[2][0] = 3.0*(b->b[4][2]-b->b[5][2]);

    Jacb[0][1] = 3.0*(b->b[2][0]-b->b[5][0]);
    Jacb[1][1] = 3.0*(b->b[2][1]-b->b[5][1]);
    Jacb[2][1] = 3.0*(b->b[2][2]-b->b[5][2]);

    Hb[0][0] = 6.0*(b->b[3][0] - 2.0*b->b[9][0] + b->b[6][0]);
    Hb[1][0] = 6.0*(b->b[3][1] - 2.0*b->b[9][1] + b->b[6][1]);
    Hb[2][0] = 6.0*(b->b[3][2] - 2.0*b->b[9][2] + b->b[6][2]);

    Hb[0][1] = 6.0*(b->b[4][0] - b->b[5][0] - b->b[9][0] + b->b[6][0]);
    Hb[1][1] = 6.0*(b->b[4][1] - b->b[5][1] - b->b[9][1] + b->b[6][1]);
    Hb[2][1] = 6.0*(b->b[4][2] - b->b[5][2] - b->b[9][2] + b->b[6][2]);

    Hb[0][2] = 6.0*(b->b[2][0] - 2.0*b->b[5][0] + b->b[6][0]);
    Hb[1][2] = 6.0*(b->b[2][1] - 2.0*b->b[5][1] + b->b[6][1]);
    Hb[2][2] = 6.0*(b->b[2][2] - 2.0*b->b[5][2] + b->b[6][2]);
  }

  /* tau = jacb *(lambda0,lambda1)*/
  tau[0] = Jacb[0][0]*lambda[0] + Jacb[0][1]*lambda[1];
  tau[1] = Jacb[1][0]*lambda[0] + Jacb[1][1]*lambda[1];
  tau[2] = Jacb[2][0]*lambda[0] + Jacb[2][1]*lambda[1];
  ll = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];
  if ( ll < _MMG5_EPSD )  return(0);

  l = 1.0 / sqrt(ll);
  tau[0] *= l;
  tau[1] *= l;
  tau[2] *= l;

  gammasec[0] = Hb[0][0]*lambda[0]*lambda[0] + 2.0*Hb[0][1]*lambda[0]*lambda[1] + Hb[0][2]*lambda[1]*lambda[1];
  gammasec[1] = Hb[1][0]*lambda[0]*lambda[0] + 2.0*Hb[1][1]*lambda[0]*lambda[1] + Hb[1][2]*lambda[1]*lambda[1];
  gammasec[2] = Hb[2][0]*lambda[0]*lambda[0] + 2.0*Hb[2][1]*lambda[0]*lambda[1] + Hb[2][2]*lambda[1]*lambda[1];

  ps = tau[0]*gammasec[0] + tau[1]*gammasec[1] + tau[2]*gammasec[2];
  c[0] = gammasec[0] - ps*tau[0];
  c[1] = gammasec[1] - ps*tau[1];
  c[2] = gammasec[2] - ps*tau[2];

  kappacur = MG_MAX(0.0,1.0/ll*sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]));
  kappacur = 1.0/8.0 * kappacur/mesh->info.hausd;
  kappacur = MG_MIN(kappacur,isqhmin);
  kappacur = MG_MAX(kappacur,isqhmax);

  return(kappacur);
}
