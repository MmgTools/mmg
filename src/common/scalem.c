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
 * \file common/scalem.c
 * \brief Scale and unscale mesh and solution.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 if fail (computed bounding box too small).
 *
 * Compute the mesh bounding box and fill the \a min, \a max and \a delta fields
 * of the \ref _MMG5_info structure.
 *
 */
int _MMG5_boundingBox(MMG5_pMesh mesh) {
  MMG5_pPoint    ppt;
  int            k,i;
  double         dd;

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

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric or solution structure.
 * \return 1 if success, 0 if fail (computed bounding box too small).
 *
 * Scale the mesh and the size informations between 0 and 1.
 *
 */
int _MMG5_scaleMesh(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pPoint    ppt;
  double         dd,d1;
  int            k;
  MMG5_pPar      par;


  /* compute bounding box */
  if ( ! _MMG5_boundingBox(mesh) ) return(0);

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
  if ( mesh->info.hmin > 0. )  mesh->info.hmin  *= dd;
  else  mesh->info.hmin = 0.01;

  if ( mesh->info.hmax > 0. )  mesh->info.hmax  *= dd;
  else mesh->info.hmax  = 1.;

  mesh->info.hausd *= dd;

  /* normalize sizes */
  if ( met->m ) {
    if ( met->size == 1 ) {
      for (k=1; k<=mesh->np; k++)    met->m[k] *= dd;
    }
    else if ( met->size==3 ){
//CECILE : a verifier si size==3 ou non (pour disp)
      d1 = 1.0 / (dd*dd);
      for (k=1; k<=6*mesh->np; k++)  met->m[k] *= d1;
    }
  }

  /* normalize local parameters */
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    par->hmin  *= dd;
    par->hmax  *= dd;
    par->hausd *= dd;
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric or solution structure.
 * \return 1.
 *
 * Unscale the mesh and the size informations to their initial sizes.
 *
 */
int _MMG5_unscaleMesh(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pPoint     ppt;
  double          dd;
  int             k,i;
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
  if ( met->m ) {
//CECILE : a mettre d'equerre : pour disp ete aniso
    if ( met->size == 6 ) {
      dd = 1.0 / (dd*dd);
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) )  continue;
        for (i=0; i<6; i++)  met->m[6*(k)+1+i] *= dd;
      }
    }
    else {
      for (k=1; k<=mesh->np ; k++) {
        ppt = &mesh->point[k];
        if ( MG_VOK(ppt) )  met->m[k] *= dd;
      }
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
