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
 * \file mmg3d/intmet.c
 * \brief Metric interpolations.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmg3d.h"

//* ** */
 /* * \param mesh pointer toward the mesh structure. */
 /* * \param met pointer toward the metric structure. */
 /* * \param k element index. */
 /* * \param i local index of edge in \a k. */
 /* * \param ip global index of the new point in which we want to compute the metric. */
 /* * \param s interpolation parameter (between 0 and 1). */
 /* * \return 0 if fail, 1 otherwise. */
 /* * */
 /* * Interpolation of anisotropic sizemap at parameter \a s along edge \a i of elt */
 /* * \a k. */
 /* * */
 /* *\/ */
/* void _MMG5_intmet_ani(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip, */
/*                       double s) { */
/*   MMG5_pTetra   pt; */
/*   MMG5_pPoint   ppt; */
/*   MMG5_pxPoint  pxp; */
/*   double        *m; */

/*   pt = &mesh->tetra[k]; */
/*   m  = &met->m[6*ip]; */
/*   if ( pt->tag[i] & MG_GEO ) { */
/*     ppt = &mesh->point[ip]; */
/*     assert(ppt->xp); */
/*     pxp = &mesh->xpoint[ppt->xp]; */
/*     _MMG5_intridmet(mesh,met,k,i,s,pxp->n1,m); */
/*   } */
/*   else { */
/*     _MMG5_intregmet(mesh,met,k,i,s,m); */
/*   } */
/* } */

/**
 * \param ma pointer on a metric
 * \param mb pointer on a metric
 * \param mp pointer on the computed interpolated metric
 * \param t interpolation parameter (comprise between 0 and 1)
 * \return 1 if success, 0 if fail.
 *
 * Linear interpolation of anisotropic sizemap along an internal edge
 *
 */
int _MMG5_intmetvol_ani(double *ma,double *mb,double *mp,double t) {
  double        dma[6],dmb[6],mai[6],mbi[6],mi[6];
  int           i;

  for (i=0; i<6; i++) {
    dma[i] = ma[i];
    dmb[i] = mb[i];
  }
  if ( !_MMG5_invmat(dma,mai) || !_MMG5_invmat(dmb,mbi) ) {
    fprintf(stderr,"  ## INTERP INVALID METRIC.\n");
    return(0);
  }
  for (i=0; i<6; i++)
    mi[i] = (1.0-t)*mai[i] + t*mbi[i];

  if ( !_MMG5_invmat(mi,mai) ) {
    fprintf(stderr,"  ## INTERP INVALID METRIC.\n");
    return(0);
  }

  for (i=0; i<6; i++)  mp[i] = mai[i];
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of the tetra.
 * \param ip index of the point on which we compute the metric.
 * \param cb barycentric coordinates of \a ip in \a k.
 * \return 1.
 *
 * Linear interpolation of isotropic sizemap in a tetra given the barycentric
 * coordinates of the new point in \a k.
 *
 */
int _MMG5_interp4bar_iso(MMG5_pMesh mesh, MMG5_pSol met, int k, int ip,
                         double cb[4]) {
  MMG5_pTetra pt;

  pt = &mesh->tetra[k];

  met->m[ip] = cb[0]*met->m[pt->v[0]]+cb[1]*met->m[pt->v[1]] +
    cb[2]*met->m[pt->v[2]]+cb[3]*met->m[pt->v[3]];

  return(1);

}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of the tetra.
 * \param ip index of the point on which we compute the metric.
 * \param cb barycentric coordinates of \a ip in \a k.
 * \return 1 if success, 0 if fail.
 *
 * Linear interpolation of anisotropic sizemap in a tetra given the barycentric
 * coordinates of the new point in \a k.
 *
 */
#warning add test for the ridge point metric
int _MMG5_interp4bar_ani(MMG5_pMesh mesh, MMG5_pSol met, int k, int ip,
                         double cb[4]) {
  MMG5_pTetra   pt;
  double        dm0[6],dm1[6],dm2[6],dm3[6];
  double        m0i[6],m1i[6],m2i[6],m3i[6],mi[6];
  int           i;

  pt = &mesh->tetra[k];

  for (i=0; i<6; i++) {
    dm0[i] = met->m[met->size*pt->v[0]+i];
    dm1[i] = met->m[met->size*pt->v[1]+i];
    dm2[i] = met->m[met->size*pt->v[2]+i];
    dm3[i] = met->m[met->size*pt->v[3]+i];
  }
  if ( !_MMG5_invmat(dm0,m0i) || !_MMG5_invmat(dm1,m1i) ||
       !_MMG5_invmat(dm2,m2i) || !_MMG5_invmat(dm3,m3i) ) {
    fprintf(stderr,"  ## INTERP INVALID METRIC.\n");
    return(0);
  }
  for (i=0; i<6; i++)
    mi[i] = cb[0]*m0i[i] + cb[1]*m1i[i] + cb[2]*m2i[i] + cb[3]*m3i[i];

  if ( !_MMG5_invmat(mi,m0i) ) {
    fprintf(stderr,"  ## INTERP INVALID METRIC.\n");
    return(0);
  }

  for (i=0; i<6; i++)  met->m[met->size*ip] = m0i[i];

  return 1;
}
