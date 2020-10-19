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
 * \file mmg2d/scalem_2d.c
 * \brief Scale and unscale mesh and solution
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "mmg2d.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a solution structure.
 * \param sol pointer toward a solution structure (level-set or displacement).
 * \return 1 if success, 0 if fail (computed bounding box too small
 * or one af the anisotropic input metric is not valid).
 *
 * Scale the mesh and the size informations between 0 and 1.
 * Compute a default value for the hmin/hmax parameters if needed.
 * Truncate the metric sizes between hmin/hmax
 *
 */
int MMG2D_scaleMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  double         dd,isqhmin,isqhmax;
  double         *m;
  double         lambda[2],v[2][2];
  int            i,k,iadr;
  int8_t         sethmin,sethmax;
  static int8_t  mmgWarn0=0, mmgWarn1=0;

  if ( !MMG5_scale_meshAndSol(mesh,met,sol,&dd,&sethmin,&sethmax) ) {
    return 0;
  }

  if ( (!met) || (met && !met->np) ) {
    return 1;
  }

  /* metric truncature and normalization and default values for hmin/hmax if not
   * provided by the user ( 0.1 \times the minimum of the metric sizes for hmin
   * and 10 \times the max of the metric sizes for hmax ). */
  switch (met->size) {
  case 1:
    if ( !MMG5_scale_scalarMetric ( mesh, met, dd, sethmin, sethmax ) ) {
      return 0;
    }
    break;

  case 3:
    dd = 1.0 / (dd*dd);
    /* Normalization */
    for (k=1; k<=mesh->np; k++) {
      iadr = k*met->size;
      for (i=0; i<met->size; i++)  met->m[iadr+i] *= dd;
    }

    /* compute hmin and hmax parameters if not provided by the user */
    if ( !sethmin ) {
      mesh->info.hmin = FLT_MAX;
    }
    if ( !sethmax ) {
      mesh->info.hmax = 0.;
    }

    for (k=1; k<=mesh->np; k++)  {
      iadr = k*met->size;
      m    = &met->m[iadr];

      /* Check the input metric */
      if ( !MMG5_eigensym(m,lambda,v) ) {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  ## Error: %s: at least 1 wrong metric.\n",
                  __func__);
        }
        return 0;
      }
      for (i=0; i<2; i++) {
        if(lambda[i]<=0) {
          if ( !mmgWarn1 ) {
            mmgWarn1 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 wrong metric"
                    " (eigenvalue : %e %e -- det %e\n",__func__,lambda[0],
                    lambda[1],m[0]*m[2]-m[1]*m[1]);
          }
          return 0;
        }
        if ( !sethmin ) {
          mesh->info.hmin = MG_MIN(mesh->info.hmin,1./sqrt(lambda[i]));
        }

        if ( !sethmax ) {
          mesh->info.hmax = MG_MAX(mesh->info.hmax,1./sqrt(lambda[i]));
        }
      }
    }

    /* Check the compatibility between the user settings and the automatically
     * computed values */
    MMG5_check_hminhmax(mesh,sethmin,sethmax);

    /* Truncature */
    isqhmin  = 1.0 / (mesh->info.hmin*mesh->info.hmin);
    isqhmax  = 1.0 / (mesh->info.hmax*mesh->info.hmax);
    for (k=1; k<=mesh->np; k++) {
      iadr = k*met->size;

      m    = &met->m[iadr];
      if ( !MMG5_eigensym(m,lambda,v) ) {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  ## Error: %s: at least 1 wrong metric.\n",
                  __func__);
        }
        return 0;
      }
      for (i=0; i<2; i++) {
        if(lambda[i]<=0) {
          if ( !mmgWarn1 ) {
            mmgWarn1 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 wrong metric"
                    " (eigenvalue : %e %e -- det %e\n",__func__,lambda[0],
                    lambda[1],m[0]*m[2]-m[1]*m[1]);
          }
          return 0;
        }
        lambda[i]=MG_MIN(isqhmin,lambda[i]);
        lambda[i]=MG_MAX(isqhmax,lambda[i]);
      }
      m[0] = v[0][0]*v[0][0]*lambda[0] + v[1][0]*v[1][0]*lambda[1];
      m[1] = v[0][0]*v[0][1]*lambda[0] + v[1][0]*v[1][1]*lambda[1];
      m[2] = v[0][1]*v[0][1]*lambda[0] + v[1][1]*v[1][1]*lambda[1];
    }
    break;
  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected metric size (%d)\n",__func__,met->size);
    break;
  }

  return 1;
}
