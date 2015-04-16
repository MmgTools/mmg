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
 * \file mmg3d/mmg3d3.c
 * \brief Lagragian meshing.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmg3d.h"

extern char  ddb;

/** Check if moving mesh with disp for a fraction t yields a valid mesh */
int _MMG5_chkmovmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t) {
  MMG5_pTetra  pt;
  MMG5_pPoint  ppt;
  double       *v,c[4][3],tau,vol;
  int          k,np;
  char         i,j;
  
  /* Pseudo time-step = fraction of disp to perform */
  tau = (double)t / _MMG5_SHORTMAX;
  
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    
    for (i=0; i<4; i++) {
      np = pt->v[i];
      ppt = &mesh->point[np];
      v = &disp->m[3*(np-1)+1];
      for (j=0; j<3; j++)
         c[i][j] = ppt->c[j]+tau*v[j];
    }
    
    vol = _MMG5_det4pt(c[0],c[1],c[2],c[3]);
    
    if ( vol < _MMG5_NULKAL )        // Other criteria : eg. a rate of degradation, etc... ?
      return(0);
  }

  return(1);
}

/** Return the largest fraction t that makes the motion along disp valid */
short _MMG5_dikomv(MMG5_pMesh mesh,MMG5_pSol disp) {
  int     it,maxit;
  short   t,tmin,tmax;
  char    ier;
 
  maxit = 200;
  it = 0;
  
  tmin = 0;
  tmax = _MMG5_SHORTMAX;
  
  /* If full displacement can be achieved */
  if ( _MMG5_chkmovmesh(mesh,disp,tmax) )
    return(tmax);

  /* Else, find the largest displacement by dichotomy */
  while( tmin != tmax && it < maxit ) {
    t = (tmin+tmax)/2;
  
    /* Case that tmax = tmin +1 : check move with tmax */
    if ( t == tmin ) {
      ier = _MMG5_chkmovmesh(mesh,disp,tmax);
      if ( ier )
        return(tmax);
      else
        return(tmin);
    }
    
    /* General case: check move with t */
    ier = _MMG5_chkmovmesh(mesh,disp,t);
    if ( ier )
      tmin = t;
    else
      tmax = t;
  
    it++;
  }
  
  return(tmin);
}

/** Perform mesh motion along disp, for a fraction t, and the corresponding updates */
int _MMG5_dispmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t) {
  MMG5_pPoint ppt;
  double      *v,tau,ctau;
  int         k;
  char        i;
  
  tau = (double)t /_MMG5_SHORTMAX;
  ctau = 1.0 - tau;
  
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    
    if ( !MG_VOK(ppt) ) continue;
    v = &disp->m[3*(k-1)+1];
    
    for (i=0; i<3; i++) {
      ppt->c[i] = ppt->c[i] + tau*v[i];
      v[i] *= ctau;
    }
  }
  
  return(1);
}

/** Lagrangian node displacement and meshing */
int _MMG5_mmg3d3(MMG5_pMesh mesh,MMG5_pSol disp) {
  short   t;
  char    ier;
  
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** LAGRANGIAN MOTION\n");
  
  t = _MMG5_dikomv(mesh,disp);
  printf("Greatest possible displacement : %d \n",t);
  
  ier = _MMG5_dispmesh(mesh,disp,t);
  if ( !ier ) {
    fprintf(stdout,"  ** Impossible motion\n");
    return(0);
  }
  
  /* Clean memory (but not pointer) */
  /* Doing this, memcur of mesh is decreased by size of displacement... ????? */
  _MMG5_DEL_MEM(mesh,disp->m,(disp->size*disp->npmax+1)*sizeof(double));
  
  /* Generates errors saying that the pointer being freed was not allocated */
  //_MMG5_DEL_MEM(mesh,disp,sizeof(MMG5_pSol));

  return(1);
}