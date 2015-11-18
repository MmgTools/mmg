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
 * \file mmgs/API_functionsf.c
 * \brief Fortran API functions for MMGS library.
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \note Please, refer to the \ref mmgs/libmmgs.h file for functions
 * documentation.
 *
 * Define the Fortran API functions for MMGS library: adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */

#include "mmgs.h"

/**
 * See \ref MMGS_Init_mesh function in common/libmmg.h file.
 */
FORTRAN_NAME(MMGS_INIT_MESH, mmgs_init_mesh,(MMG5_pMesh *mesh, MMG5_pSol *sol,
                                             MMG5_pSol *dummy),
             (mesh,sol,dummy) ){

MMGS_Init_mesh(mesh,sol);

return;
}

/**
 * See \ref MMGS_Init_parameters function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_INIT_PARAMETERS,mmgs_init_parameters,(MMG5_pMesh *mesh),(mesh)) {
  MMGS_Init_parameters(*mesh);
  return;
}
/**
 * See \ref MMGS_Set_inputMeshName function in \ref common/libmmg.h file.
 */
FORTRAN_NAME(MMGS_SET_INPUTMESHNAME, mmgs_set_inputmeshname,
             (MMG5_pMesh *mesh, char* meshin, int *strlen, int* retval),
             (mesh,meshin,strlen,retval)) {
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMGS_Set_inputMeshName(*mesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMGS_Set_inputSolName function in \ref common/libmmg.h file.
 */
FORTRAN_NAME(MMGS_SET_INPUTSOLNAME, mmgs_set_inputsolname,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solin, int* strlen, int* retval),
             (mesh,sol,solin,strlen,retval)) {

  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,solin,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMGS_Set_inputSolName(*mesh,*sol,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMGS_Set_outputMeshName function in mmgs/libmmgs.h or
 * mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_OUTPUTMESHNAME,mmgs_set_outputmeshname,
             (MMG5_pMesh *mesh, char* meshout, int* strlen,int* retval),
             (mesh,meshout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshout,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMGS_Set_outputMeshName(*mesh, tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMGS_Set_outputSolName function in \ref common/libmmg.h file.
 */
FORTRAN_NAME(MMGS_SET_OUTPUTSOLNAME,mmgs_set_outputsolname,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solout,int* strlen, int* retval),
             (mesh,sol,solout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,solout,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMGS_Set_outputSolName(*mesh,*sol,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMGS_Set_meshSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_MESHSIZE,mmgs_set_meshsize,
             (MMG5_pMesh *mesh, int *np, int *nt, int *na, int *retval),
             (mesh,np,nt,na,retval)) {
  *retval = MMGS_Set_meshSize(*mesh,*np,*nt,*na);
  return;
}

/**
 * See \ref MMGS_Get_solSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_GET_SOLSIZE,mmgs_get_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity, int* np, int* typSol, int* retval),
             (mesh,sol,typEntity,np,typSol,retval)) {

  *retval = MMGS_Get_solSize(*mesh,*sol,typEntity,np,typSol);
  return;
}

/**
 * See \ref MMGS_Get_meshSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_GET_MESHSIZE,mmgs_get_meshsize,
             (MMG5_pMesh *mesh, int* np, int* nt, int* na, int* retval),
             (mesh,np,nt, na,retval)) {

  *retval = MMGS_Get_meshSize(*mesh,np,nt,na);
  return;
}

/**
 * See \ref MMGS_Set_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_VERTEX,mmgs_set_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
              int* pos, int* retval),
             (mesh,c0,c1,c2,ref,pos,retval)) {

  *retval = MMGS_Set_vertex(*mesh,*c0,*c1,*c2,*ref,*pos);
  return;
}

/**
 * See \ref MMGS_Get_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_GET_VERTEX,mmgs_get_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
              int* isCorner, int* isRequired, int* retval),
             (mesh,c0,c1,c2,ref,isCorner,isRequired, retval)) {
  *retval = MMGS_Get_vertex(*mesh,c0,c1,c2,ref,isCorner,isRequired);
  return;
}
/**
 * See \ref MMGS_Set_triangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_TRIANGLE,mmgs_set_triangle,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref,int* pos,
              int* retval),
             (mesh,v0,v1,v2,ref,pos,retval)) {
  *retval = MMGS_Set_triangle(*mesh, *v0, *v1, *v2, *ref, *pos);
  return;
}

/**
 * See \ref MMGS_Get_triangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_GET_TRIANGLE,mmgs_get_triangle,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref
              ,int* isRequired, int* retval),
             (mesh,v0,v1,v2,ref,isRequired,retval)) {
  *retval = MMGS_Get_triangle(*mesh,v0,v1,v2,ref,isRequired);
  return;
}

/**
 * See \ref MMGS_Set_edge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_EDGE,mmgs_set_edge,
             (MMG5_pMesh *mesh, int *v0, int *v1, int *ref, int *pos, int* retval),
             (mesh,v0,v1,ref,pos,retval)){
  *retval = MMGS_Set_edge(*mesh,*v0,*v1,*ref,*pos);
  return;
}

/**
 * See \ref MMGS_Get_edge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_GET_EDGE,mmgs_get_edge,(MMG5_pMesh *mesh, int* e0, int* e1, int* ref
                                          ,int* isRidge, int* isRequired, int* retval),
             (mesh,e0,e1,ref,isRidge,isRequired,retval)) {
  *retval = MMGS_Get_edge(*mesh,e0,e1,ref,isRidge,isRequired);
  return;
}

/**
 * See \ref MMGS_Set_corner function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_CORNER,mmgs_set_corner,(MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval =  MMGS_Set_corner(*mesh,*k);
  return;
}

/**
 * See \ref MMGS_Set_requiredVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_REQUIREDVERTEX,mmgs_set_requiredvertex,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval =  MMGS_Set_requiredVertex(*mesh,*k);
  return;
}
/**
 * See \ref MMGS_Set_requiredTriangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_REQUIREDTRIANGLE,mmgs_set_requiredtriangle,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMGS_Set_requiredTriangle(*mesh, *k);
  return;
}

/**
 * See \ref MMGS_Set_ridge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_RIDGE,mmgs_set_ridge,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMGS_Set_ridge(*mesh,*k);
  return;
}

/**
 * See \ref MMGS_Set_requiredEdge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_REQUIREDEDGE,mmgs_set_requirededge,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMGS_Set_requiredEdge(*mesh,*k);
  return;
}

/**
 * See \ref MMGS_Set_scalarSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_SCALARSOL,mmgs_set_scalarsol,
             (MMG5_pSol *met, double *s, int *pos, int* retval),
             (met,s,pos,retval)) {
  *retval = MMGS_Set_scalarSol(*met,*s,*pos);
  return;
}

/**
 * See \ref MMGS_Get_scalarSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_GET_SCALARSOL,mmgs_get_scalarsol,
             (MMG5_pSol *met, double* s, int* retval),
             (met,s,retval)) {
  *retval = MMGS_Get_scalarSol(*met,s);
  return;
}

/**
 * See \ref MMGS_Set_vectorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_VECTORSOL,mmgs_set_vectorsol,
             (MMG5_pSol *met, double *vx, double *vy, double *vz,
              int *pos, int* retval),
             (met,vx,vy,vz,pos,retval)) {
  *retval = MMGS_Set_vectorSol(*met,*vx,*vy,*vz,*pos);
  return;
}

/**
 * See \ref MMGS_Get_vectorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_GET_VECTORSOL,mmgs_get_vectorsol,
             (MMG5_pSol *met, double* vx,double *vy, double *vz, int* retval),
             (met,vx,vy,vz,retval)) {
  *retval = MMGS_Get_vectorSol(*met,vx,vy,vz);
  return;
}

/**
 * See \ref MMGS_Set_tensorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_TENSORSOL,mmgs_set_tensorsol,
             (MMG5_pSol *met, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int *pos, int* retval),
             (met,m11,m12,m13,m22,m23,m33,pos,retval)) {
  *retval = MMGS_Set_tensorSol(*met,*m11,*m12,*m13,*m22,*m23,*m33,*pos);
  return;
}

/**
 * See \ref MMGS_Get_tensorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_GET_TENSORSOL,mmgs_get_tensorsol,
             (MMG5_pSol *met, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int* retval),
             (met,m11,m12,m13,m22,m23,m33,retval)) {
  *retval = MMGS_Get_tensorSol(*met,m11,m12,m13,m22,m23,m33);
  return;
}

/**
 * See \ref MMGS_Chk_meshData function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_CHK_MESHDATA,mmgs_chk_meshdata,
             (MMG5_pMesh *mesh,MMG5_pSol *met, int* retval),
             (mesh,met,retval)) {
  *retval = MMGS_Chk_meshData(*mesh,*met);
  return;
}

/**
 * See \ref MMGS_Set_iparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_IPARAMETER,mmgs_set_iparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *iparam, int *val, int* retval),
             (mesh,sol,iparam,val,retval)){
  *retval = MMGS_Set_iparameter(*mesh,*sol,*iparam,*val);
  return;
}

/**
 * See \ref MMGS_Get_iparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_GET_IPARAMETER,mmgs_get_iparameter,
             (MMG5_pMesh *mesh, int *iparam, int* retval),
             (mesh,iparam,retval)){
  *retval = MMGS_Get_iparameter(*mesh,*iparam);
  return;
}

/**
 * See \ref MMGS_Set_dparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_DPARAMETER,mmgs_set_dparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *dparam, double *val, int* retval),
             (mesh,sol,dparam,val,retval)){
  *retval = MMGS_Set_dparameter(*mesh,*sol,*dparam,*val);
  return;
}

/**
 * See \ref MMGS_Set_localParameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SET_LOCALPARAMETER,mmgs_set_localparameter,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, int *typ, int *ref,
              double *hmin, double *hmax,double *hausd, int* retval),
             (mesh,sol,typ,ref,hmin,hmax,hausd,retval)){
  *retval = MMGS_Set_localParameter(*mesh,*sol,*typ,*ref,*hmin,*hmax,*hausd);
  return;
}


/**
 * See \ref MMGS_Free_structures function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_FREE_STRUCTURES,mmgs_free_structures,(MMG5_pMesh *mesh,
                                                        MMG5_pSol *met,
                                                        MMG5_pSol *dummy),
             (mesh,met,dummy) ){
  MMGS_Free_structures(*mesh,*met);
  return;
}

/**
 * See \ref MMGS_Free_names function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_FREE_NAMES,mmgs_free_names,(MMG5_pMesh *mesh,MMG5_pSol *met,
                                              MMG5_pSol *dummy),
             (mesh,met,dummy))
{
  MMGS_Free_names(*mesh,*met);
  return;
}

/**
 * See \ref MMGS_loadMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_LOADMESH,mmgs_loadmesh,(MMG5_pMesh *mesh,int* retval),(mesh, retval)){

  *retval = MMGS_loadMesh(*mesh);

  return;
}

/**
 * See \ref MMGS_loadMet function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_LOADMET,mmgs_loadmet,(MMG5_pMesh *mesh,MMG5_pSol *met,int* retval),(mesh,met,retval)){

  *retval = MMGS_loadMet(*mesh,*met);

  return;
}

/**
 * See \ref MMGS_saveMet function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMGS_SAVEMET,mmgs_savemet,(MMG5_pMesh *mesh,MMG5_pSol *met,int* retval),(mesh,met,retval)){

  *retval = MMGS_saveMet(*mesh,*met);

  return;
}
