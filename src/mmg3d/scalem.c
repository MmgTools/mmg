/* =============================================================================
**  This file is part of the Mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
**
**  Mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with Mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the Mmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file mmg3d/scalem.c
 * \brief Scale and unscale mesh and solution.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmg3d.h"

int _MMG5_scaleMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSingul sing) {
  MMG5_pPoint    ppt;
#ifdef SINGUL
  MMG5_psPoint   ppts;
  double    deltb,delta[3];
#endif
  double    dd;
  int       i,k;
  MMG5_pPar      par;

  /* compute bounding box */
  for (i=0; i<3; i++) {
    mesh->info.min[i] =  DBL_MAX;
    mesh->info.max[i] = -DBL_MAX;
  }
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    for (i=0; i<3; i++) {
      if ( ppt->c[i] > mesh->info.max[i] )  mesh->info.max[i] = ppt->c[i];
      if ( ppt->c[i] < mesh->info.min[i] )  mesh->info.min[i] = ppt->c[i];
    }
    ppt->tmp = 0;
  }
  mesh->info.delta = 0.0;
  for (i=0; i<3; i++) {
    dd = mesh->info.max[i] - mesh->info.min[i];
    if ( dd > mesh->info.delta )  mesh->info.delta = dd;
  }
  if ( mesh->info.delta < _MMG5_EPSD ) {
    fprintf(stdout,"  ## Unable to scale mesh.\n");
    return(0);
  }

  /* normalize coordinates */
  dd = 1.0 / mesh->info.delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->c[0] = dd * (ppt->c[0] - mesh->info.min[0]);
    ppt->c[1] = dd * (ppt->c[1] - mesh->info.min[1]);
    ppt->c[2] = dd * (ppt->c[2] - mesh->info.min[2]);
  }

  /* normalize values */
  mesh->info.hmin  *= dd;
  mesh->info.hmax  *= dd;
  mesh->info.hausd *= dd;

  /* normalize sizes */
  if ( met->size == 1 && met->m ) {
    for (k=1; k<=mesh->np; k++)
      met->m[k] *= dd;
  }

#ifdef SINGUL
  /* 2nd mesh (sing) is quarter sized */
  /* Get the size of sing in every direction */
  if ( mesh->info.sing && sing->ns ) {
    deltb = 0.0;

    for (i=0; i<mesh->dim; i++) {
      delta[i] = fabs(sing->max[i]-sing->min[i]);
      if ( delta[i] > deltb )  deltb = delta[i];   // deltb = max dimension
    }
    if ( deltb < _MMG5_EPSD ) {
      fprintf(stdout,"  ## Unable to scale mesh\n");
      return(0);
    }

    /* centering */
    dd = 1.0 / deltb;
    for (k=1; k<=sing->ns; k++) {
      ppts = &sing->point[k];
      for (i=0; i<mesh->dim; i++) {
        ppts->c[i] = MMG5_SIZE * (dd * (ppts->c[i]-sing->min[i])) +
          0.5 * (1.0 - MMG5_SIZE*dd*delta[i]);
      }
    }
  }
#endif
  /* normalize local parameters */
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    par->hausd *= dd;
  }

  return(1);
}

int _MMG5_unscaleMesh(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pPoint     ppt;
  double     dd;
  int        k;
  MMG5_pPar       par;

  /* de-normalize coordinates */
  dd = mesh->info.delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->c[0] = ppt->c[0] * dd + mesh->info.min[0];
    ppt->c[1] = ppt->c[1] * dd + mesh->info.min[1];
    ppt->c[2] = ppt->c[2] * dd + mesh->info.min[2];
  }

  /* unscale sizes */
  if(met->m){
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) )	met->m[k] *= dd;
    }
  }

  /* unscale paramter values */
  mesh->info.hmin  *= dd;
  mesh->info.hmax  *= dd;
  mesh->info.hausd *= dd;

  /* normalize local parameters */
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    par->hausd *= dd;
  }

  return(1);
}
