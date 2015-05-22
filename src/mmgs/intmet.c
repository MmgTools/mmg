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
 * \param k element index.
 * \param i local index of edge in \a k.
 * \param s interpolation parameter.
 * \param mr computed metric.
 * \return call to _MMG5_interpreg_ani (thus, 0 if fail, 1 otherwise).
 *
 * Metric interpolation on edge \a i in elt \a it at
 * parameter \f$ 0 <= s0 <= 1 \f$ from \a p1 result is stored in \a mr. edge
 * \f$ p_1-p_2 \f$ must not be a ridge.
 *
 * */
int intregmet(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,double s,double mr[6]) {
  MMG5_pTria     pt;

  pt  = &mesh->tria[k];
  return(_MMG5_interpreg_ani(mesh,met,pt,i,s,mr));
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
