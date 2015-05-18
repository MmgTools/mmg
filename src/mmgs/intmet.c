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



/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ip1 global index of ridge extremity.
 * \param ip2 global index of ridge extremity.
 * \param s interpolation parameter (between 0 and 1).
 * \param v normal at the point at which we want to compute the metric.
 * \param mr computed anisotropic size.
 * \return 1 if success, 0 otherwise.
 *
 * Anisotropic metric interpolation between two points \f$p_1\f$ and \f$p_2\f$
 * such that \f$edge_0 = (p_1p_2)\f$ is ridge. \a v is a direction vector, aimed
 * at pointing towards direction of n1 at interpolated point.
 *
 */
int _MMG5_intridmet(MMG5_pMesh mesh,MMG5_pSol met,int ip1, int ip2,double s,
                    double v[3],double mr[6]) {
  MMG5_pxPoint   go1,go2;
  MMG5_pPoint    p1,p2;
  double         *m1,*m2,*n11,*n12,*n21,*n22,ps11,ps12,dd,hn1,hn2;

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

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k element index.
 * \param i local index of edge in \a k.
 * \param ip global index of the new point in which we want to compute the metric.
 * \param s interpolation parameter (between 0 and 1).
 * \return 0 if fail, 1 otherwise.
 *
 * Interpolation of anisotropic sizemap at parameter \a s along edge \a i of elt
 * \a k.
 *
 */
void intmet_ani(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  MMG5_pxPoint  go;
  double        *m;
  int           ip1, ip2, i1, i2;

  pt  = &mesh->tria[k];
  i1  = _MMG5_inxt2[i];
  i2  = _MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];

  m  = &met->m[6*ip];
  if ( pt->tag[i] & MG_GEO ) {
    ppt = &mesh->point[ip];
    assert(ppt->xp);
    go = &mesh->xpoint[ppt->xp];
    _MMG5_intridmet(mesh,met,ip1,ip2,s,go->n1,m);
  }
  else {
    intregmet(mesh,met,k,i,s,m);
  }
}
