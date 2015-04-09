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
 * \file common/mmg_tools.c
 * \brief Fonctions for anisotropic size map computation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param np0 index of edge's extremity.
 * \param ux distance \f$[p0;p1]\f$ along x axis.
 * \param uy distance \f$[p0;p1]\f$ along y axis.
 * \param uz distance \f$[p0;p1]\f$ along z axis.
 * \param mr computed metric tensor.
 * \return 1.
 *
 * Build metric tensor at ridge point p0, when computations with respect to p1
 * are to be held.
 *
 */
int _MMG5_buildridmet(MMG5_pMesh mesh,MMG5_pSol met,int np0,
                      double ux,double uy,double uz,double mr[6]) {
  MMG5_pPoint p0;
  MMG5_pxPoint  go;
  double ps1,ps2,*n1,*n2,*t,*m,dv,u[3],r[3][3];

  p0 = &mesh->point[np0];
  if ( !(MG_GEO & p0->tag) )  return(0);
  m = &met->m[6*(np0)+1];
  t = &p0->n[0];
  go = &mesh->xpoint[p0->ig];

  /* Decide between the two possible configurations */
  n1 = &go->n1[0];
  n2 = &go->n2[0];

  ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
  ps2 = ux*n2[0] + uy*n2[1] + uz*n2[2];

  if ( fabs(ps2)<fabs(ps1) ) {
    n1 = &go->n2[0];
    dv = m[2];
  }
  else{
    dv = m[1];
  }

  u[0] = n1[1]*t[2] - n1[2]*t[1];
  u[1] = n1[2]*t[0] - n1[0]*t[2];
  u[2] = n1[0]*t[1] - n1[1]*t[0];

  /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(m[0],dv,0)*/
  r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n1[0];
  r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n1[1];
  r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n1[2];

  mr[0] = m[0]*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1];
  mr[1] = m[0]*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1];
  mr[2] = m[0]*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1];
  mr[3] = m[0]*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1];
  mr[4] = m[0]*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1];
  mr[5] = m[0]*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1];

  return(1);
}

/**
 * \param mesh pointer toward the mesh stucture.
 * \param ip1 first point of face.
 * \param ip2 second point of face.
 * \param ip3 third point of face.
 * \param n pointer to store the computed normal.
 * \return 1
 *
 * Compute face normal.
 *
 */
inline int _MMG5_norpts(MMG5_pMesh mesh,int ip1,int ip2, int ip3,double *n) {
  MMG5_pPoint   p1,p2,p3;
  double   dd,abx,aby,abz,acx,acy,acz,det;

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  p3 = &mesh->point[ip3];

  /* area */
  abx = p2->c[0] - p1->c[0];
  aby = p2->c[1] - p1->c[1];
  abz = p2->c[2] - p1->c[2];

  acx = p3->c[0] - p1->c[0];
  acy = p3->c[1] - p1->c[1];
  acz = p3->c[2] - p1->c[2];

  n[0] = aby*acz - abz*acy;
  n[1] = abz*acx - abx*acz;
  n[2] = abx*acy - aby*acx;
  det  = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];

  if ( det < _MMG5_EPSD2 )  return(0);

  dd = 1.0 / sqrt(det);
  n[0] *= dd;
  n[1] *= dd;
  n[2] *= dd;

  return(1);
}

/**
 * \param mesh pointer toward the mesh stucture.
 * \param pt pointer toward the triangle structure.
 * \param n pointer to store the computed normal.
 * \return 1
 *
 * Compute face normal.
 *
 */
inline int _MMG5_nortri(MMG5_pMesh mesh,MMG5_pTria pt,double *n) {
  double   *a,*b,*c,dd,abx,aby,abz,acx,acy,acz,det;

  a = mesh->point[pt->v[0]].c;
  b = mesh->point[pt->v[1]].c;
  c = mesh->point[pt->v[2]].c;

  /* area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];

  n[0] = aby*acz - abz*acy;
  n[1] = abz*acx - abx*acz;
  n[2] = abx*acy - aby*acx;
  det  = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( det < _MMG5_EPSD2 )  return(0);

  dd = 1.0 / sqrt(det);
  n[0] *= dd;
  n[1] *= dd;
  n[2] *= dd;
  return(1);
}
