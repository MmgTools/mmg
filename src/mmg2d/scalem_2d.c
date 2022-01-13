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
  double         dd;

  if ( !MMG5_scale_meshAndSol ( mesh, met, sol, &dd ) ) {
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
    if ( !MMG5_scale_scalarMetric ( mesh, met, dd ) ) {
      return 0;
    }
    break;

  case 3:
    /* Set function pointer before it is used */
    MMG5_solTruncature_ani = MMG2D_solTruncature_ani;
    if ( !MMG5_scale_tensorMetric ( mesh, met, dd ) ) {
      return 0;
    }

    break;
  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected metric size (%d)\n",__func__,met->size);
    break;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the solution structure.
 *
 * \return 0 if fail, 1 if succeed.
 *
 * Compute hmin and hmax values from unflagged points (if not setted by the
 * user), check the values and truncate the 2D metric.
 *
 */
int MMG2D_solTruncature_ani(MMG5_pMesh mesh, MMG5_pSol met) {
  MMG5_pPoint ppt;
  int         k,iadr;
  double      isqhmin, isqhmax;
  double      lambda[2],vp[2][2];

 /* Security check: if hmin (resp. hmax) is not setted, it means that sethmin
   * (resp. sethmax) is not setted too */
  if ( !MMG5_check_setted_hminhmax(mesh) ) {
    return 0;
  }

  /* If not provided by the user, compute hmin/hmax from the metric computed by
   * the DoSol function. */
  isqhmin = 0.;
  isqhmax = FLT_MAX;
  if ( (!mesh->info.sethmin) || (!mesh->info.sethmax) ) {
    for (k=1; k<=mesh->np; k++)  {
      ppt = &mesh->point[k];
      if ( (!MG_VOK(ppt)) || ppt->flag ) continue;
      iadr = met->size*k;

      MMG5_eigensym(met->m+iadr,lambda,vp);

      assert (lambda[0] > 0. && lambda[1] > 0. && "Negative eigenvalue");

      isqhmin = MG_MAX(isqhmin,lambda[0]);
      isqhmin = MG_MAX(isqhmin,lambda[1]);

      isqhmax = MG_MIN(isqhmax,lambda[0]);
      isqhmax = MG_MIN(isqhmax,lambda[1]);
    }
  }

  if ( !mesh->info.sethmin ) {
    mesh->info.hmin = 1./sqrt(isqhmin);
  }

  if ( !mesh->info.sethmax ) {
    mesh->info.hmax = 1./sqrt(isqhmax);
  }

  MMG5_check_hminhmax(mesh, mesh->info.sethmin, mesh->info.sethmax);

  /* vertex size */
  isqhmin = 1./(mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1./(mesh->info.hmax*mesh->info.hmax);

  for (k=1; k<=mesh->np; k++) {
    iadr = 3*k;

    ppt = &mesh->point[k];
    if ( ppt->flag ) {
      met->m[iadr] = met->m[iadr+2] = isqhmax;
      met->m[iadr+1] = 0;
    }
    else {
      MMG5_eigensym(met->m+iadr,lambda,vp);

      lambda[0]=MG_MAX(isqhmax,MG_MIN(isqhmin,lambda[0]));
      lambda[1]=MG_MAX(isqhmax,MG_MIN(isqhmin,lambda[1]));

      met->m[iadr]   = vp[0][0]*vp[0][0]*lambda[0] + vp[1][0]*vp[1][0]*lambda[1];
      met->m[iadr+1] = vp[0][0]*vp[0][1]*lambda[0] + vp[1][0]*vp[1][1]*lambda[1];
      met->m[iadr+2] = vp[0][1]*vp[0][1]*lambda[0] + vp[1][1]*vp[1][1]*lambda[1];
    }
  }

  if ( mesh->info.ddebug ) {
    /* print unscaled values for debug purpose */
    fprintf(stdout,"     After truncature computation:   hmin %lf (user setted %d)\n"
            "                                     hmax %lf (user setted %d)\n",
            mesh->info.delta * mesh->info.hmin,mesh->info.sethmin,
            mesh->info.delta * mesh->info.hmax,mesh->info.sethmax);
  }
  return 1;
}
