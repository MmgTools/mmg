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
 * \file mmg2d/API_functionsf.c
 * \brief Fortran API functions for MMG2D library.
 * \author Cecile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 07 2015
 * \copyright GNU Lesser General Public License.
 * \note Please, refer to the \ref mmg2d/libmmg2d.h file for functions
 * documentation.
 *
 * Define the Fortran API functions for MMG2D library: adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */
#include "mmg2d.h"
/**
 * See \ref MMG2_Init_mesh function in mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2_INIT_MESH, mmg2_init_mesh,(MMG5_pMesh *mesh, MMG5_pSol *sol
               ),(mesh,sol
                 )) {
  MMG2_Init_mesh(mesh,sol);
  return;
}
/**
 * See \ref MMG5_Set_iparameter function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_SET_IPARAMETER,mmg5_set_iparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *iparam, int *val, int* retval),
             (mesh,sol,iparam,val,retval)){
  *retval = MMG5_Set_iparameter(*mesh,*sol,*iparam,*val);
  return;
}
/**
 * See \ref MMG5_Set_dparameter function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_SET_DPARAMETER,mmg5_set_dparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *dparam, double *val, int* retval),
             (mesh,sol,dparam,val,retval)){
  *retval = MMG5_Set_dparameter(*mesh,*sol,*dparam,*val);
  return;
}
/**
 * See \ref MMG5_Set_meshSize function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_SET_MESHSIZE,mmg5_set_meshsize,
             (MMG5_pMesh *mesh, int *np, int *nt, int *na, int *retval),
             (mesh,np,nt,na,retval)) {
  *retval = MMG5_Set_meshSize(*mesh,*np,*nt,*na);
  return;
}
/**
 * See \ref MMG5_Set_solSize function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_SET_SOLSIZE,mmg5_set_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity,
              int* np, int* typSol, int* retval),
             (mesh, sol, typEntity, np, typSol, retval)) {
  *retval = MMG5_Set_solSize(*mesh,*sol,*typEntity,*np,*typSol);
  return;
}
/**
 * See \ref MMG5_Set_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_VERTEX,mmg5_set_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, int* ref,
              int* pos, int* retval),
             (mesh,c0,c1,ref,pos,retval)) {

  *retval = MMG5_Set_vertex(*mesh,*c0,*c1,*ref,*pos);
  return;
}
/**
 * See \ref MMG5_Set_triangle function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_SET_TRIANGLE,mmg5_set_triangle,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref,int* pos,
              int* retval),
             (mesh,v0,v1,v2,ref,pos,retval)) {
  *retval = MMG5_Set_triangle(*mesh, *v0, *v1, *v2, *ref, *pos);
  return;
}
/**
 * See \ref MMG5_Set_edge function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_SET_EDGE,mmg5_set_edge,
             (MMG5_pMesh *mesh, int *v0, int *v1, int *ref, int *pos, int* retval),
             (mesh,v0,v1,ref,pos,retval)){
  *retval = MMG5_Set_edge(*mesh,*v0,*v1,*ref,*pos);
  return;
}
/**
 * See \ref MMG5_Get_meshSize function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_GET_MESHSIZE,mmg5_get_meshsize,
             (MMG5_pMesh *mesh, int* np, int* nt, int* na, int* retval),
             (mesh,np,nt, na,retval)) {

  *retval = MMG5_Get_meshSize(*mesh,np,nt,na);
  return;
}
/**
 * See \ref MMG5_Set_scalarSol function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_SET_SCALARSOL,mmg5_set_scalarsol,
             (MMG5_pSol *met, double *s, int *pos, int* retval),
             (met,s,pos,retval)) {
  *retval = MMG5_Set_scalarSol(*met,*s,*pos);
  return;
}
/**
 * See \ref MMG5_Set_tensorSol function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_SET_TENSORSOL,mmg5_set_tensorsol,
             (MMG5_pSol *met, double *s, int *pos, int* retval),
             (met,s,pos,retval)) {
  *retval = MMG5_Set_tensorSol(*met,s,*pos);
  return;
}
/**
 * See \ref MMG5_Chk_meshData function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_CHK_MESHDATA,mmg5_chk_meshdata,
             (MMG5_pMesh *mesh,MMG5_pSol *met, int* retval),
             (mesh,met,retval)) {
  *retval = MMG5_Chk_meshData(*mesh,*met);
  return;
}
/**
 * See \ref MMG5_Free_structures function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_FREE_STRUCTURES,mmg5_free_structures,(MMG5_pMesh *mesh,MMG5_pSol *met
               ),(mesh,met
                 )){
  MMG5_Free_structures(*mesh,*met);
  return;
}


/**
 * See \ref MMG2_loadMesh function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2_LOADMESH,mmg5_loadmesh,(MMG5_pMesh *mesh,char* meshin,int* strlen,int* retval),(mesh, meshin,strlen,retval)){

  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG2_loadMesh(*mesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG5_loadSol function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_LOADSOL,mmg5_loadSol,(MMG5_pSol *met,char *meshin,int* strlen,int *npmax,int *msh,int* retval),(met,meshin,strlen,npmax,msh,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG2_loadSol(*met,tmp,*npmax,*msh);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG2_saveSol function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2_SAVESOL,mmg2_saveSol,(MMG5_pMesh *mesh,MMG5_pSol *met,char *meshin,int *strlen,int* retval),(mesh,met,meshin,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG2_saveSol(*mesh,*met,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}
