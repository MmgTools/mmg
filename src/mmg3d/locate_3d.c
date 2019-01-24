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
 * \file locate_3d.c
 * \brief research of the element to which a point belongs
 * \author L. Cirrottola
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 */

#include "mmg3d.h"


/**
 * \param mesh pointer to the mesh structure
 * \param pt pointer to the current tetra
 * \param coord pointer to the point coordinates
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 *
 * \return 0 if fail, 1 if success
 *
 *  Compute the barycentric coordinates of a given point in a given tetrahedron.
 *
 */
int MMG3D_compute_baryCoord( MMG5_pMesh mesh, MMG5_pTetra pt,
                    double *coord, MMG5_baryCoord *barycoord ) {
  double *c0,*c1,*c2,*c3,vol;

  vol = MMG5_orvol( mesh->point, pt->v );

  c0 = mesh->point[pt->v[0]].c;
  c1 = mesh->point[pt->v[1]].c;
  c2 = mesh->point[pt->v[2]].c;
  c3 = mesh->point[pt->v[3]].c;

  barycoord[0].val = MMG5_det4pt( coord, c1,    c2,    c3    )/vol;
  barycoord[1].val = MMG5_det4pt( c0,    coord, c2,    c3    )/vol;
  barycoord[2].val = MMG5_det4pt( c0,    c1,    coord, c3    )/vol;
  barycoord[3].val = MMG5_det4pt( c0,    c1,    c2,    coord )/vol;

  barycoord[0].idx = 0;
  barycoord[1].idx = 1;
  barycoord[2].idx = 2;
  barycoord[3].idx = 3;

  return 1;
}

/**
 * \param a pointer to point barycentric coordinates
 * \param b pointer to point barycentric coordinates
 *
 * \return -1 if (a < b), +1 if (a > b), 0 if equal
 *
 *  Compare the barycentric coordinates of a given point in a given tetrahedron.
 *
 */
int MMG3D_compare_baryCoord( const void *a,const void *b ) {
  MMG5_baryCoord *coord_a;
  MMG5_baryCoord *coord_b;

  coord_a = (MMG5_baryCoord *)a;
  coord_b = (MMG5_baryCoord *)b;

  if( coord_a->val > coord_b-> val ) return 1;
  if( coord_a->val < coord_b-> val ) return -1;

  return 0;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 *
 * \return ie if positive, index of the target element; if negative, index of
 * the closest element; 0 if not found
 *
 *  Locate a point in a background mesh by traveling the elements adjacency.
 *
 */
int MMG3D_locatePoint( MMG5_pMesh mesh, MMG5_pPoint ppt ) {
  MMG5_pTetra    ptr,pt1;
  MMG5_baryCoord barycoord[4];
  int            *adja,iel,ip,idxTet,step,closestTet;
  double         vol,eps,closestDist;

  eps = MMG5_EPS;

  idxTet = 1;

  step = 0;
  ++mesh->base;
  while(step <= mesh->ne) {
    step++;

    /** Get tetra */
    ptr = &mesh->tetra[idxTet];
    if ( !MG_EOK(ptr) ) continue;

    adja = &mesh->adja[4*(idxTet-1)+1];
    vol = MMG5_orvol( mesh->point, ptr->v );

    /** Mark tetra */
    ptr->flag = mesh->base;

    /** Get barycentric coordinates and sort them in ascending order */
    MMG3D_compute_baryCoord(mesh, ptr, ppt->c, barycoord);
    qsort(barycoord,4,sizeof(MMG5_baryCoord),MMG3D_compare_baryCoord);

    /** Exit if inside the element */
    if( barycoord[0].val > -eps ) break;

    /** Compute new direction */
    for( ip=0; ip<4; ip++ ) {
      iel = adja[barycoord[ip].idx]/4;

      /* Skip if on boundary */
      if (!iel) continue;

      /* Skip if already marked */
      pt1 = &mesh->tetra[iel];
      if(pt1->flag == mesh->base) continue;

      /* Get next otherwise */
      idxTet = iel;
      break;
    }

    /** Stuck: Start exhaustive research */
    if (ip == 4) step = mesh->ne+1;

  }

  /** Boundary hit or cyclic path: Perform exhaustive research */
  if( step == (mesh->ne+1) ) {
    fprintf(stderr,"\n  ## Warning: %s: Cannot locate point, performing"
            " exhaustive research.\n",__func__);

    closestTet = 0;
    closestDist = 1.0e10;
    for( idxTet=1; idxTet<mesh->ne+1; idxTet++ ) {

      /** Get tetra */
      ptr = &mesh->tetra[idxTet];

      if ( !MG_EOK(ptr) ) continue;

      adja = &mesh->adja[4*(idxTet-1)+1];
      vol = MMG5_orvol( mesh->point, ptr->v );

      /** Mark tetra */
      ptr->flag = mesh->base;

      /** Get barycentric coordinates and sort them in ascending order */
      MMG3D_compute_baryCoord(mesh, ptr, ppt->c, barycoord);
      qsort(barycoord,4,sizeof(MMG5_baryCoord),MMG3D_compare_baryCoord);

      /** Exit if inside the element */
      if( barycoord[0].val > -eps ) break;

      /** Save element index (with negative sign) if it is the closest one */
      if( fabs(barycoord[0].val)*vol < closestDist ) {
        closestDist = fabs(barycoord[0].val)*vol;
        closestTet = -idxTet;
      }

    }

    /** Element not found: Return the closest one with negative sign (if found) */
    if ( idxTet == mesh->ne+1 ) {
      fprintf(stderr,"\n  ## Warning: %s: Point not located, smallest external"
              " volume %e.",__func__,closestDist);
      return closestTet;
    }

  }
  return idxTet;
}
