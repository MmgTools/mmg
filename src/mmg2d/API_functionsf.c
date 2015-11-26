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
 * See \ref MMG2D_Init_mesh function in mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_INIT_MESH, mmg2d_init_mesh,(MMG5_pMesh *mesh, MMG5_pSol *sol),
             (mesh,sol)) {
  MMG2D_Init_mesh(mesh,sol);
  return;
}
/**
 * See \ref MMG2D_Init_fileNames function in mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_INIT_FILENAMES, mmg2d_init_filenames,(MMG5_pMesh *mesh, MMG5_pSol *sol),
             (mesh,sol)) {
  MMG2D_Init_fileNames(*mesh,*sol);
  return;
}
/**
 * See \ref MMG2D_Init_parameters function in mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_INIT_PARAMETERS, mmg2d_init_parameters,(MMG5_pMesh *mesh),
             (mesh)) {
  MMG2D_Init_parameters(*mesh);
  return;
}
/**
 * See \ref MMG2D_Set_inputMeshName function in \ref common/libmmgcommon.h file.
 */
FORTRAN_NAME(MMG2D_SET_INPUTMESHNAME, mmg2d_set_inputmeshname,
             (MMG5_pMesh *mesh, char* meshin, int *strlen, int* retval),
             (mesh,meshin,strlen,retval)) {
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMG2D_Set_inputMeshName(*mesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG2D_Set_inputSolName function in \ref common/libmmgcommon.h file.
 */
FORTRAN_NAME(MMG2D_SET_INPUTSOLNAME, mmg2d_set_inputsolname,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solin, int* strlen, int* retval),
             (mesh,sol,solin,strlen,retval)) {

  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,solin,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMG2D_Set_inputSolName(*mesh,*sol,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG2D_Set_outputMeshName function in mmgs/libmmgs.h or
 * mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG2D_SET_OUTPUTMESHNAME,mmg2d_set_outputmeshname,
             (MMG5_pMesh *mesh, char* meshout, int* strlen,int* retval),
             (mesh,meshout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshout,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMG2D_Set_outputMeshName(*mesh, tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG2D_Set_outputSolName function in \ref common/libmmgcommon.h file.
 */
FORTRAN_NAME(MMG2D_SET_OUTPUTSOLNAME,mmg2d_set_outputsolname,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solout,int* strlen, int* retval),
             (mesh,sol,solout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,solout,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMG2D_Set_outputSolName(*mesh,*sol,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG2D_Set_iparameter function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SET_IPARAMETER,mmg2d_set_iparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *iparam, int *val, int* retval),
             (mesh,sol,iparam,val,retval)){
  *retval = MMG2D_Set_iparameter(*mesh,*sol,*iparam,*val);
  return;
}
/**
 * See \ref MMG2D_Set_dparameter function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SET_DPARAMETER,mmg2d_set_dparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *dparam, double *val, int* retval),
             (mesh,sol,dparam,val,retval)){
  *retval = MMG2D_Set_dparameter(*mesh,*sol,*dparam,*val);
  return;
}
/**
 * See \ref MMG2D_Set_meshSize function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SET_MESHSIZE,mmg2d_set_meshsize,
             (MMG5_pMesh *mesh, int *np, int *nt, int *na, int *retval),
             (mesh,np,nt,na,retval)) {
  *retval = MMG2D_Set_meshSize(*mesh,*np,*nt,*na);
  return;
}
/**
 * See \ref MMG2D_Set_solSize function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SET_SOLSIZE,mmg2d_set_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity,
              int* np, int* typSol, int* retval),
             (mesh, sol, typEntity, np, typSol, retval)) {
  *retval = MMG2D_Set_solSize(*mesh,*sol,*typEntity,*np,*typSol);
  return;
}
/**
 * See \ref MMG2D_Set_vertex function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SET_VERTEX,mmg2d_set_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, int* ref,
              int* pos, int* retval),
             (mesh,c0,c1,ref,pos,retval)) {

  *retval = MMG2D_Set_vertex(*mesh,*c0,*c1,*ref,*pos);
  return;
}
/**
 * See \ref MMG2D_Get_vertex function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_GET_VERTEX,mmg2d_get_vertex,
             (MMG5_pMesh *mesh, int *num,double* c0, double* c1, int* ref,
              int* isCorner, int* isRequired, int* retval),
             (mesh,num,c0,c1,ref,isCorner,isRequired, retval)) {
  *retval = MMG2D_Get_vertex(*mesh,*num,c0,c1,ref,isCorner,isRequired);
  return;
}
/**
 * See \ref MMG2D_Set_triangle function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SET_TRIANGLE,mmg2d_set_triangle,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref,int* pos,
              int* retval),
             (mesh,v0,v1,v2,ref,pos,retval)) {
  *retval = MMG2D_Set_triangle(*mesh, *v0, *v1, *v2, *ref, *pos);
  return;
}
/**
 * See \ref MMG2D_Set_edge function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SET_EDGE,mmg2d_set_edge,
             (MMG5_pMesh *mesh, int *v0, int *v1, int *ref, int *pos, int* retval),
             (mesh,v0,v1,ref,pos,retval)){
  *retval = MMG2D_Set_edge(*mesh,*v0,*v1,*ref,*pos);
  return;
}
/**
 * See \ref MMG2D_Get_meshSize function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_GET_MESHSIZE,mmg2d_get_meshsize,
             (MMG5_pMesh *mesh, int* np, int* nt, int* na, int* retval),
             (mesh,np,nt, na,retval)) {

  *retval = MMG2D_Get_meshSize(*mesh,np,nt,na);
  return;
}
/**
 * See \ref MMG2D_Set_scalarSol function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SET_SCALARSOL,mmg2d_set_scalarsol,
             (MMG5_pSol *met, double *s, int *pos, int* retval),
             (met,s,pos,retval)) {
  *retval = MMG2D_Set_scalarSol(*met,*s,*pos);
  return;
}
/**
 * See \ref MMG2D_Set_tensorSol function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SET_TENSORSOL,mmg2d_set_tensorsol,
             (MMG5_pSol *met, double *s, int *pos, int* retval),
             (met,s,pos,retval)) {
  *retval = MMG2D_Set_tensorSol(*met,s,*pos);
  return;
}
/**
 * See \ref MMG2D_Chk_meshData function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_CHK_MESHDATA,mmg2d_chk_meshdata,
             (MMG5_pMesh *mesh,MMG5_pSol *met, int* retval),
             (mesh,met,retval)) {
  *retval = MMG2D_Chk_meshData(*mesh,*met);
  return;
}
/**
 * See \ref MMG2D_Free_structures function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_FREE_STRUCTURES,mmg2d_free_structures,(MMG5_pMesh *mesh,MMG5_pSol *met
               ),(mesh,met
                 )){
  MMG2D_Free_structures(*mesh,*met);
  return;
}
/**
 * See \ref MMG2D_Free_names function in \ref common/libmmgcommon.h file.
 */
FORTRAN_NAME(MMG2D_FREE_NAMES,mmg2d_free_names,(MMG5_pMesh *mesh,MMG5_pSol *met ),
             (mesh,met)){
  MMG2D_Free_names(*mesh,*met);

  return;
}


/**
 * See \ref MMG2D_loadMesh function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_LOADMESH,mmg2d_loadmesh,(MMG5_pMesh *mesh,char* meshin,int* strlen,int* retval),(mesh, meshin,strlen,retval)){

  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG2D_loadMesh(*mesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG2D_loadSol function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_LOADSOL,mmg2d_loadsol,(MMG5_pMesh *mesh,MMG5_pSol *met,char *meshin,int* strlen,int *msh,int* retval),(mesh,met,meshin,strlen,msh,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG2D_loadSol(*mesh,*met,tmp,*msh);
  _MMG5_SAFE_FREE(tmp);

  return;
}


/**
 * See \ref MMG2D_saveSol function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SAVESOL,mmg2d_savesol,(MMG5_pMesh *mesh,MMG5_pSol *met,char *meshin,int *strlen,int* retval),(mesh,met,meshin,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG2D_saveSol(*mesh,*met,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}
