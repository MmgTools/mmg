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
#include "mmg2d.h"


int MMG2_scaleMesh(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria     pt;
  //Displ     pd;
  MMG5_pPoint    ppt;
  MMG5_Info     *info;
  double    dd;
  int       i,k,iadr,sethmin,sethmax;

  // pd  = mesh->disp;

  /* compute bounding box */
  info = &mesh->info;
  for (i=0; i<2; i++) {
    info->min[i] =  DBL_MAX;
    info->max[i] = -DBL_MAX;
  }
  // #warning no normalization
  //return(1);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !M_VOK(ppt) ) continue;
    for (i=0; i<2; i++) {
      if ( ppt->c[i] > info->max[i] )  info->max[i] = ppt->c[i];
      if ( ppt->c[i] < info->min[i] )  info->min[i] = ppt->c[i];
    }
  }
  info->delta = info->max[0]-info->min[0];
  dd = info->max[1]-info->min[1];
  if ( dd > info->delta )
    info->delta = dd;
  if ( info->delta < EPS30 ) {
    fprintf(stdout,"  ## Unable to scale mesh.\n");
    return(0);
  }

  /* normalize coordinates */
  dd = PRECI / info->delta;


  sethmin = 0;
  sethmax = 0;
  if ( mesh->info.hmin > 0. ) {
    mesh->info.hmin  *= dd;
    sethmin = 1;
  }
  else
    mesh->info.hmin  = 0.01;

  if ( mesh->info.hmax > 0. ) {
    mesh->info.hmax  *= dd;
    sethmax = 1;
  }
  else
    mesh->info.hmax  = 1.;
  if ( mesh->info.hmax < mesh->info.hmin ) {
    if ( sethmin && sethmax ) {
      fprintf(stdout,"  ## ERROR: MISMATCH PARAMETERS:"
              "MINIMAL MESH SIZE LARGER THAN MAXIMAL ONE.\n");
      fprintf(stdout,"  Exit program.\n");
      exit(EXIT_FAILURE);
    }
    else if ( sethmin )
      mesh->info.hmax = 100. * mesh->info.hmin;
    else
      mesh->info.hmin = 0.01 * mesh->info.hmax;
  }

  mesh->info.hausd *= dd;

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !M_VOK(ppt) )  continue;
    ppt->c[0] = dd * (ppt->c[0] - info->min[0]);
    ppt->c[1] = dd * (ppt->c[1] - info->min[1]);
  }

  /* normalize metric */
  if ( !sol->np )  return(1);
  switch (sol->size) {
  case 1:
    for (k=1; k<=mesh->np; k++)  {
      sol->m[k] *= dd;
      if(sethmin) sol->m[k]=MG_MAX(mesh->info.hmin,sol->m[k]);
      if(sethmax) sol->m[k]=MG_MIN(mesh->info.hmax,sol->m[k]);
    }
    break;

  case 3:
    if(sethmin || sethmax) {
      printf("warning imposed hmin and/or hmax ignored\n");
    }
    dd = 1.0 / (dd*dd);
    for (k=1; k<=mesh->np; k++) {
      iadr = (k-1)*sol->size + 1;
      for (i=0; i<sol->size; i++)  sol->m[iadr+i] *= dd;
      if(sethmin || sethmax) {
//#warning todo : metric troncature
      }
    }
    break;
  }

  /* compute quality */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    pt->qual = MMG2_caltri_in(mesh,sol,pt);
  }

  return(1);
}


int MMG2_unscaleMesh(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pPoint     ppt;
  MMG5_Info      *info;
  double     dd;
  int        i,k,iadr;

  info = &mesh->info;
  //return(1);

  /* de-normalize coordinates */
  dd = info->delta / (double)PRECI;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !M_VOK(ppt) )  continue;
    ppt->c[0] = ppt->c[0] * dd + info->min[0];
    ppt->c[1] = ppt->c[1] * dd + info->min[1];
  }

  /* de-normalize metric */
  if ( !sol->np )  return(1);
  switch (sol->size) {
  case 1:
    for (k=1; k<=mesh->np; k++)  sol->m[k] *= dd;
    break;

  case 3:
    dd = 1.0 / (dd*dd);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !M_VOK(ppt) )  continue;
      iadr = (k-1)*sol->size + 1;
      for (i=0; i<sol->size; i++)  sol->m[iadr+i] *= dd;
    }
    break;
  }

  return(1);
}

