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
 * \file common/scalem.c
 * \brief Scale and unscale mesh and solution.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 if fail (computed bounding box too small).
 *
 * Compute the mesh bounding box and fill the \a min, \a max and \a delta fields
 * of the \a MMG5_info structure.
 *
 */
int MMG5_boundingBox(MMG5_pMesh mesh) {
  MMG5_pPoint    ppt;
  int            k,i;
  double         dd;

  /* compute bounding box */
  for (i=0; i<mesh->dim; i++) {
    mesh->info.min[i] =  DBL_MAX;
    mesh->info.max[i] = -DBL_MAX;
  }
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    for (i=0; i<mesh->dim; i++) {
      if ( ppt->c[i] > mesh->info.max[i] )  mesh->info.max[i] = ppt->c[i];
      if ( ppt->c[i] < mesh->info.min[i] )  mesh->info.min[i] = ppt->c[i];
    }
    ppt->tmp = 0;
  }
  mesh->info.delta = 0.0;
  for (i=0; i<mesh->dim; i++) {
    dd = mesh->info.max[i] - mesh->info.min[i];
    if ( dd > mesh->info.delta )  mesh->info.delta = dd;
  }
  if ( mesh->info.delta < MMG5_EPSD ) {
    fprintf(stderr,"\n  ## Error: %s: unable to scale mesh:"
            " Check that your mesh contains non-zero points and "
            "valid elements.\n",__func__);
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sethmin 1 if hmin is setted by the user.
 * \param sethmax 1 if hmax is setted by the user.
 *
 * Check the compatibility between the automatically computed hmin/hmax values
 * and the user settings.
 *
 */
void MMG5_check_hminhmax(MMG5_pMesh mesh, int8_t sethmin, int8_t sethmax) {

  if ( !sethmin ) {
    mesh->info.hmin *=.1;
    /* Check that user has not given a hmax value lower that the founded
     * hmin. */
    if ( mesh->info.hmin > mesh->info.hmax ) {
      mesh->info.hmin = 0.1*mesh->info.hmax;
    }
  }
  if ( !sethmax ) {
    mesh->info.hmax *=10.;
    /* Check that user has not given a hmin value bigger that the founded
     * hmax. */
    if ( mesh->info.hmax < mesh->info.hmin ) {
      mesh->info.hmax = 10.*mesh->info.hmin;
    }
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param dd scaling value.
 * \param sethmin 1 if hmin must not be automatically computed
 * \param sethmax 1 if hmin must not be automatically computed
 *
 * \return 1 if success, 0 if fail.
 *
 * Scale and truncate by hmin/hmax the scalar metric stored in met.  If
 * hmin/hmax are not provided by the user, it is automatically computed from the
 * metric.
 *
 */
int MMG5_scale_scalarMetric(MMG5_pMesh mesh, MMG5_pSol met, double dd,
                            int8_t sethmin, int8_t sethmax) {
  int    k;
  static int8_t mmgWarn0 = 0;

  for (k=1; k<=mesh->np; k++)  {
    if( !MG_VOK( &mesh->point[k] ) ) continue;
    /* Check the metric */
    if ( met->m[k] <= 0 ) {
      if ( !mmgWarn0 ) {
        mmgWarn0 = 1;
        fprintf(stderr,"\n  ## Error: %s: at least 1 wrong metric.\n",
                __func__);
        return 0;
      }
    }
    /* normalization */
    met->m[k] *= dd;
  }

  /* compute hmin and hmax parameters if not provided by the user */
  if ( !sethmin ) {
    mesh->info.hmin = FLT_MAX;
    for (k=1; k<=mesh->np; k++)  {
      mesh->info.hmin = MG_MIN(mesh->info.hmin,met->m[k]);
    }
  }
  if ( !sethmax ) {
    mesh->info.hmax = 0.;
    for (k=1; k<=mesh->np; k++)  {
      mesh->info.hmax = MG_MAX(mesh->info.hmax,met->m[k]);
    }
  }

  /* Check the compatibility between the user settings and the automatically
   * computed values */
  MMG5_check_hminhmax(mesh,sethmin,sethmax);

  /* Truncature */
  for (k=1; k<=mesh->np; k++)  {
    met->m[k]=MG_MAX(mesh->info.hmin,met->m[k]);
    met->m[k]=MG_MIN(mesh->info.hmax,met->m[k]);
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a metric
 * \param sol pointer toward a solution structure (level-set or displacement).
 * \param dd pointer toward the scaling value (to fill)
 * \param sethmin setted to 1 if hmin must not be computed from the metric.
 * \param sethmax setted to 1 if hmax must not be computed from the metric.
 *
 * \return 1 if success, 0 if fail.
 *
 * Scale the mesh and the size informations between 0 and 1.
 * Compute a default value for the hmin/hmax parameters if needed.
 *
 */
int MMG5_scale_meshAndSol(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,double *dd,
                          int8_t* sethmin,int8_t *sethmax ) {
  MMG5_pPoint    ppt;
  MMG5_pPar      par;
  int            k,i;
  int8_t         hsizOrOptim;

  /* sol is a level-set or a displacement so it cannot be an aniso metric */
  if ( sol ) { assert ( sol->type == MMG5_Scalar || sol->type == MMG5_Vector ); }

  /* met is a metric so it cannot be a vector */
  if ( met ) { assert ( met->type != MMG5_Vector ); }

  /* if we are not in iso or lagrangian mode, check that sol isn't allocated */
  if ( (!mesh->info.iso) && mesh->info.lag < 0 ) {
    assert ( !(sol && sol->m) );
  }

  /* compute bounding box */
  if ( ! MMG5_boundingBox(mesh) ) return 0;

  /* normalize coordinates */
  *dd = 1.0 / mesh->info.delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    for (i=0 ; i<mesh->dim ; i++)
      ppt->c[i] = (*dd) * (ppt->c[i] - mesh->info.min[i]);
  }

  mesh->info.hausd *= (*dd);
  mesh->info.ls    *= (*dd);
  mesh->info.hsiz  *= (*dd);

  /* normalize local parameters */
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    par->hmin  *= (*dd);
    par->hmax  *= (*dd);
    par->hausd *= (*dd);
  }

  /* Check if hmin/hmax have been provided by the user and scale it if yes */
  *sethmin = 0;
  *sethmax = 0;

  if ( mesh->info.hmin > 0. ) {
    mesh->info.hmin  *= (*dd);
    *sethmin = 1;
  }
  if ( mesh->info.hmax > 0. ) {
    mesh->info.hmax  *= (*dd);
    *sethmax = 1;
  }

  hsizOrOptim = ( mesh->info.hsiz > 0. || mesh->info.optim )? 1 : 0;

  /* if we are not in optim or hsiz mode and if we don't have a metric, compute
   * default hmin/hmax */
  if ( (!hsizOrOptim) && (!(met && met->np)) ) {
    /* Set default values to hmin/hmax from the bounding box if not provided by
     * the user */
    if ( !MMG5_Set_defaultTruncatureSizes(mesh,*sethmin,*sethmax) ) {
      fprintf(stderr,"\n  ## Error: %s: Exit program.\n",__func__);
      return 0;
    }
    *sethmin = 1;
    *sethmax = 1;
  }

  if ( sol && sol->np ) {
    for ( k=sol->size; k<sol->size*(mesh->np+1); k++ ) {
      sol->m[k]   *= (*dd);
    }
  }

  return 1;

}
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param sol pointer toward a solution structure (level-set or displacement).
 *
 * \return 1 if success, 0 if fail (computed bounding box too small
 * or one af the anisotropic input metric is not valid).
 *
 * Scale the mesh and the size informations between 0 and 1.
 * Compute a default value for the hmin/hmax parameters if needed.
 *
 */
int MMG5_scaleMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  double         dd,d1;
  int            k,i;
  double         *m;
  double         lambda[3],v[3][3],isqhmin,isqhmax;
  int8_t         sethmin,sethmax;

  if ( !MMG5_scale_meshAndSol(mesh,met,sol,&dd,&sethmin,&sethmax) ) {
    return 0;
  }

  if ( (!met) || (met && !met->np) || (!met->m) ) {
    return 1;
  }

  switch ( met->size ) {
  case 1:
    if ( !MMG5_scale_scalarMetric ( mesh, met, dd, sethmin, sethmax ) ) {
      return 0;
    }
    break;

  case 6:
      d1 = 1.0 / (dd*dd);
      /* Normalization */
      for (k=1; k<=mesh->np; k++) {

        if( !MG_VOK( &mesh->point[k] ) ) continue;

        for ( i=0; i<met->size; ++i ) {
          met->m[6*k+i] *= d1;
        }
      }

      /* compute hmin and hmax parameters if not provided by the user and check
       * the input metric */
      if ( !sethmin ) {
        mesh->info.hmin = FLT_MAX;
      }

      if ( !sethmax ) {
        mesh->info.hmax = 0.;
      }

      for (k=1; k<=mesh->np; k++) {
        if( !MG_VOK( &mesh->point[k] ) ) continue;

        m    = &met->m[6*k];

        /* Check the input metric */
        if ( !MMG5_eigenv3d(1,m,lambda,v) ) {
          fprintf(stderr,"\n  ## Error: %s: unable to diagonalize at least"
                  " 1 metric (point %d).\n",__func__,k);
          return 0;
        }
        for (i=0; i<3; i++) {
          if(lambda[i]<=0) {
            fprintf(stderr,"\n  ## Error: %s: at least 1 wrong metric "
                    "(point %d -> eigenvalues : %e %e %e).\n"
                    "            metric tensor: %e %e %e %e %e %e.\n",
                    __func__,k,lambda[0],lambda[1],lambda[2],
                              m[0],m[1],m[2],m[3],m[4],m[5]);
            return 0;
          }
          if ( !sethmin )
            mesh->info.hmin = MG_MIN(mesh->info.hmin,1./sqrt(lambda[i]));
          if ( !sethmax )
            mesh->info.hmax = MG_MAX(mesh->info.hmax,1./sqrt(lambda[i]));
        }
      }

      /* Check the compatibility between the user settings and the automatically
       * computed values */
      MMG5_check_hminhmax(mesh,sethmin,sethmax);

      /* Truncature */
      isqhmin  = 1.0 / (mesh->info.hmin*mesh->info.hmin);
      isqhmax  = 1.0 / (mesh->info.hmax*mesh->info.hmax);
      for (k=1; k<=mesh->np; k++) {

        if( !MG_VOK( &mesh->point[k] ) ) continue;

        m    = &met->m[6*k];

        if ( !MMG5_eigenv3d(1,m,lambda,v) ) {
          fprintf(stderr,"\n  ## Error: %s: unable to diagonalize at least"
                  " 1 metric (point %d).\n",__func__,k);
          return 0;
        }
        for (i=0; i<3; i++) {
          if(lambda[i]<=0) {
            fprintf(stderr,"\n  ## Error: %s: at least 1 wrong metric "
                    "(point %d -> eigenvalues : %e %e %e).\n"
                    "            metric tensor: %e %e %e %e %e %e.\n",
                    __func__,k,lambda[0],lambda[1],lambda[2],
                              m[0],m[1],m[2],m[3],m[4],m[5]);
            return 0;
          }
          lambda[i]=MG_MIN(isqhmin,lambda[i]);
          lambda[i]=MG_MAX(isqhmax,lambda[i]);
        }
        m[0] = v[0][0]*v[0][0]*lambda[0] + v[1][0]*v[1][0]*lambda[1] + v[2][0]*v[2][0]*lambda[2];
        m[1] = v[0][0]*v[0][1]*lambda[0] + v[1][0]*v[1][1]*lambda[1] + v[2][0]*v[2][1]*lambda[2];
        m[2] = v[0][0]*v[0][2]*lambda[0] + v[1][0]*v[1][2]*lambda[1] + v[2][0]*v[2][2]*lambda[2];
        m[3] = v[0][1]*v[0][1]*lambda[0] + v[1][1]*v[1][1]*lambda[1] + v[2][1]*v[2][1]*lambda[2];
        m[4] = v[0][1]*v[0][2]*lambda[0] + v[1][1]*v[1][2]*lambda[1] + v[2][1]*v[2][2]*lambda[2];
        m[5] = v[0][2]*v[0][2]*lambda[0] + v[1][2]*v[1][2]*lambda[1] + v[2][2]*v[2][2]*lambda[2];
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
 * \param met pointer toward a metric.
 * \param sol pointer toward a solution structure (level-set or displacement).
 *
 * \return 1.
 *
 * Unscale the mesh and the size informations to their initial sizes.
 *
 */
int MMG5_unscaleMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  MMG5_pPoint     ppt;
  double          dd;
  int             k,i,iadr;
  MMG5_pPar       par;

  /* de-normalize coordinates */
  dd = mesh->info.delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    for ( i=0; i<mesh->dim; ++i ) {
      ppt->c[i] = ppt->c[i] * dd + mesh->info.min[i];
    }
  }

  /* unscale paramter values */
  if ( !mesh->info.sethmin ) {
    mesh->info.hmin = MMG5_NONSET_HMIN;
  }
  else {
    mesh->info.hmin  *= dd;
  }

  if ( !mesh->info.sethmax ) {
    mesh->info.hmax = MMG5_NONSET_HMAX;
  }
  else {
    mesh->info.hmax  *= dd;
  }
  mesh->info.hausd *= dd;
  mesh->info.ls    *= dd;
  mesh->info.hsiz  *= dd;

  /* normalize local parameters */
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    par->hmin  *= dd;
    par->hmax  *= dd;
    par->hausd *= dd;
  }

  /* de-normalize level-set or displacement */
  if ( sol && sol->np && sol->m ) {
    assert ( sol->size != MMG5_Tensor );
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      iadr = k*sol->size;
      for (i=0; i<sol->size; i++)  sol->m[iadr+i] *= dd;
    }
  }

  /* reset the scaling data to ensure that if we try to unscale again, we will
   * do nothing */
  mesh->info.delta = 1.;
  mesh->info.min[0]= 0.;
  mesh->info.min[1]= 0.;
  mesh->info.min[2]= 0.;

  /* de-normalize metric */
  if ( !(met && met->np && met->m) )  return 1;

  /* unscale sizes */
  switch (met->type) {
  case 1:
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      met->m[k] *= dd;
    }
    break;
  case 3:
    dd = 1.0 / (dd*dd);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      for (i=0; i<met->size; i++)  met->m[met->size*k+i] *= dd;
    }
    break;
  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected metric size (%d)\n",__func__,met->size);
    break;
  }

  return 1;
}
