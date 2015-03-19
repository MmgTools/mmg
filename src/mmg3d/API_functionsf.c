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
 * \file mmg3d/API_functionsf.c
 * \brief Fortran API functions for MMG3D library.
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \note Please, refer to the \ref mmg3d/libmmg3d.h file for functions
 * documentation.
 *
 * Define the Fortran API functions for MMG3D library: adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */

#include "mmg3d.h"

/**
 * \def FORTRAN_NAME(nu,nl,pl,pc)
 * \brief Adds function definitions.
 * \param nu function name in upper case.
 * \param nl function name in lower case.
 * \param pl type of arguments.
 * \param pc name of arguments.
 * \note Macro coming from Scotch library.
 *
 * Adds function definitions with upcase, underscore and double
 * underscore to match any fortran compiler.
 *
 */
#define FORTRAN_NAME(nu,nl,pl,pc)               \
    void nu pl;                                 \
    void nl pl                                  \
    { nu pc; }                                  \
    void nl##_ pl                               \
    { nu pc; }                                  \
    void nl##__ pl                              \
    { nu pc; }                                  \
    void nu pl

/**
 * See \ref MMG5_Init_mesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_INIT_MESH, mmg5_init_mesh,(MMG5_pMesh *mesh, MMG5_pSol *sol
#ifdef SINGUL
                                             , MMG5_pSingul *sing
#endif
                 ),(mesh,sol
#ifdef SINGUL
                    ,sing
#endif
                     )) {
#ifdef SINGUL
    MMG5_Init_mesh(mesh,sol,sing);
#else
    MMG5_Init_mesh(mesh,sol);
#endif

    return;
}

/**
 * See \ref MMG5_Init_parameters function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_INIT_PARAMETERS,mmg5_init_parameters,(MMG5_pMesh *mesh),(mesh)) {
    MMG5_Init_parameters(*mesh);
    return;
}

/**
 * See \ref MMG5_Init_fileNames function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(INIT_FILENAMES,init_filenames,(MMG5_pMesh *mesh,MMG5_pSol *sol
#ifdef SINGUL
                                            ,MMG5_pSingul *sing
#endif
                 ),(mesh,sol
#ifdef SINGUL
                    ,sing
#endif
                     )) {
#ifdef SINGUL
    MMG5_Init_fileNames(*mesh,*sol,*sing);
#else
    MMG5_Init_fileNames(*mesh,*sol);
#endif
    return;
}

/**
 * See \ref MMG5_Set_inputMeshName function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_INPUTMESHNAME, mmg5_set_inputmeshname,
             (MMG5_pMesh *mesh, char* meshin, int *strlen, int* retval),
             (mesh,meshin,strlen,retval)) {
    char *tmp = NULL;

    tmp = (char*)malloc((*strlen+1)*sizeof(char));
    strncpy(tmp,meshin,*strlen);
    tmp[*strlen] = '\0';
    *retval = MMG5_Set_inputMeshName(*mesh,tmp);
    _MMG5_SAFE_FREE(tmp);

    return;
}

/**
 * See \ref MMG5_Set_inputSolName function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_INPUTSOLNAME, mmg5_set_inputsolname,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solin, int* strlen, int* retval),
             (mesh,sol,solin,strlen,retval)) {

    char *tmp = NULL;

    tmp = (char*)malloc((*strlen+1)*sizeof(char));
    strncpy(tmp,solin,*strlen);
    tmp[*strlen] = '\0';
    *retval = MMG5_Set_inputSolName(*mesh,*sol,tmp);
    _MMG5_SAFE_FREE(tmp);

    return;
}

/**
 * See \ref MMG5_Set_outputMeshName function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_OUTPUTMESHNAME,mmg5_set_outputmeshname,
             (MMG5_pMesh *mesh, char* meshout, int* strlen,int* retval),
             (mesh,meshout,strlen,retval)){
    char *tmp = NULL;

    tmp = (char*)malloc((*strlen+1)*sizeof(char));
    strncpy(tmp,meshout,*strlen);
    tmp[*strlen] = '\0';
    *retval = MMG5_Set_outputMeshName(*mesh, tmp);
    _MMG5_SAFE_FREE(tmp);

    return;
}

/**
 * See \ref MMG5_Set_outputSolName function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_OUTPUTSOLNAME,mmg5_set_outputsolname,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solout,int* strlen, int* retval),
             (mesh,sol,solout,strlen,retval)){
    char *tmp = NULL;

    tmp = (char*)malloc((*strlen+1)*sizeof(char));
    strncpy(tmp,solout,*strlen);
    tmp[*strlen] = '\0';
    *retval = MMG5_Set_outputSolName(*mesh,*sol,tmp);
    _MMG5_SAFE_FREE(tmp);

    return;
}

#ifdef SINGUL
/**
 * See \ref MMG5_Set_inputSingulName function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_INPUTSINGULNAME,mmg5_set_inputsingulname,
             (MMG5_pMesh *mesh,MMG5_pSingul *sing, char* singin, int* strlen, int* retval),
             (mesh,sing,singin,strlen,retval)) {
    char *tmp = NULL;

    tmp = (char*)malloc((*strlen+1)*sizeof(char));
    strncpy(tmp,singin,*strlen);
    tmp[*strlen] = '\0';
    *retval = MMG5_Set_inputSingulName(*mesh,*sing,tmp);
    _MMG5_SAFE_FREE(tmp);

    return;
}
#endif

/**
 * See \ref MMG5_Set_solSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SOLSIZE,mmg5_set_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity,
              int* np, int* typSol, int* retval),
             (mesh, sol, typEntity, np, typSol, retval)) {
    *retval = MMG5_Set_solSize(*mesh,*sol,*typEntity,*np,*typSol);
    return;
}

/**
 * See \ref MMG5_Set_meshSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_MESHSIZE,mmg5_set_meshsize,
             (MMG5_pMesh *mesh, int *np, int *ne, int *nt, int *na, int *retval),
             (mesh,np,ne,nt,na,retval)) {
    *retval = MMG5_Set_meshSize(*mesh,*np,*ne,*nt,*na);
    return;
}

#ifdef SINGUL
/**
 * See \ref MMG5_Set_singulSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SINGULSIZE,mmg5_set_singulsze,
             (MMG5_pMesh *mesh,MMG5_pSingul *sing, int *np, int *na, int *retval),
             (mesh,sing,np,na,retval)) {
    *retval = MMG5_Set_singulSize(*mesh,*sing, *np, *na);
    return;
}
#endif

/**
 * See \ref MMG5_Get_solSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_SOLSIZE,mmg5_get_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity, int* np, int* typSol, int* retval),
             (mesh,sol,typEntity,np,typSol,retval)) {

    *retval = MMG5_Get_solSize(*mesh,*sol,typEntity,np,typSol);
    return;
}

/**
 * See \ref MMG5_Get_meshSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_MESHSIZE,mmg5_get_meshsize,
             (MMG5_pMesh *mesh, int* np, int* ne, int* nt, int* na, int* retval),
             (mesh,np,ne,nt, na,retval)) {

    *retval = MMG5_Get_meshSize(*mesh,np,ne,nt,na);
    return;
}

/**
 * See \ref MMG5_Set_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_VERTEX,mmg5_set_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
              int* pos, int* retval),
             (mesh,c0,c1,c2,ref,pos,retval)) {

    *retval = MMG5_Set_vertex(*mesh,*c0,*c1,*c2,*ref,*pos);
    return;
}

/**
 * See \ref MMG5_Get_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_VERTEX,mmg5_get_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
              int* isCorner, int* isRequired, int* retval),
             (mesh,c0,c1,c2,ref,isCorner,isRequired, retval)) {
    *retval = MMG5_Get_vertex(*mesh,c0,c1,c2,ref,isCorner,isRequired);
    return;
}

/**
 * See \ref MMG5_Set_tetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_TETRAHEDRON,mmg5_set_tetrahedron,
             (MMG5_pMesh *mesh, int *v0, int *v1, int *v2, int *v3, int *ref,
              int *pos, int* retval),
             (mesh,v0,v1,v2,v3,ref,pos,retval)){
    *retval = MMG5_Set_tetrahedron(*mesh,*v0,*v1,*v2,*v3,*ref,*pos);
    return;
}

/**
 * See \ref MMG5_Get_tetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_TETRAHEDRON,mmg5_get_tetrahedron,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* v3,
              int* ref, int* isRequired, int* retval),
             (mesh,v0,v1,v2,v3,ref,isRequired,retval)) {
    *retval = MMG5_Get_tetrahedron(*mesh,v0,v1,v2,v3,ref,isRequired);
    return;
}

/**
 * See \ref MMG5_Set_triangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_TRIANGLE,mmg5_set_triangle,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref,int* pos,
              int* retval),
             (mesh,v0,v1,v2,ref,pos,retval)) {
    *retval = MMG5_Set_triangle(*mesh, *v0, *v1, *v2, *ref, *pos);
    return;
}

/**
 * See \ref MMG5_Get_triangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_TRIANGLE,mmg5_get_triangle,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref
              ,int* isRequired, int* retval),
             (mesh,v0,v1,v2,ref,isRequired,retval)) {
    *retval = MMG5_Get_triangle(*mesh,v0,v1,v2,ref,isRequired);
    return;
}

/**
 * See \ref MMG5_Set_edge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_EDGE,mmg5_set_edge,
             (MMG5_pMesh *mesh, int *v0, int *v1, int *ref, int *pos, int* retval),
             (mesh,v0,v1,ref,pos,retval)){
    *retval = MMG5_Set_edge(*mesh,*v0,*v1,*ref,*pos);
    return;
}

/**
 * See \ref MMG5_Get_edge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_EDGE,mmg5_get_edge,(MMG5_pMesh *mesh, int* e0, int* e1, int* ref
                                          ,int* isRidge, int* isRequired, int* retval),
             (mesh,e0,e1,ref,isRidge,isRequired,retval)) {
    *retval = MMG5_Get_edge(*mesh,e0,e1,ref,isRidge,isRequired);
    return;
}

/**
 * See \ref MMG5_Set_corner function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_CORNER,mmg5_set_corner,(MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
    *retval =  MMG5_Set_corner(*mesh,*k);
    return;
}

/**
 * See \ref MMG5_Set_requiredVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_REQUIREDVERTEX,mmg5_set_requiredvertex,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
    *retval =  MMG5_Set_requiredVertex(*mesh,*k);
    return;
}

/**
 * See \ref MMG5_Set_requiredTetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_REQUIREDTETRAHEDRON,mmg5_set_requiredtetrahedron,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
    *retval = MMG5_Set_requiredTetrahedron(*mesh,*k);
    return;
}

/**
 * See \ref MMG5_Set_requiredTriangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_REQUIREDTRIANGLE,mmg5_set_requiredtriangle,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
    *retval = MMG5_Set_requiredTriangle(*mesh, *k);
    return;
}

/**
 * See \ref MMG5_Set_ridge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_RIDGE,mmg5_set_ridge,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
    *retval = MMG5_Set_ridge(*mesh,*k);
    return;
}

/**
 * See \ref MMG5_Set_requiredEdge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_REQUIREDEDGE,mmg5_set_requirededge,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
    *retval = MMG5_Set_requiredEdge(*mesh,*k);
    return;
}

/**
 * See \ref MMG5_Set_scalarSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SCALARSOL,mmg5_set_scalarsol,
             (MMG5_pSol *met, double *s, int *pos, int* retval),
             (met,s,pos,retval)) {
    *retval = MMG5_Set_scalarSol(*met,*s,*pos);
    return;
}

/**
 * See \ref MMG5_Get_scalarSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_SCALARSOL,mmg5_get_scalarsol,
             (MMG5_pSol *met, double* s, int* retval),
             (met,s,retval)) {
    *retval = MMG5_Get_scalarSol(*met,s);
    return;
}

#ifdef SINGUL
/**
 * See \ref MMG5_Set_singulVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SINGULVERTEX,mmg5_set_singluvertex,
             (MMG5_pSingul *sing, double *c0, double *c1,
              double *c2, int *typ, int *pos, int* retval),
             (sing,c0,c1,c2,typ,pos,retval)) {
    *retval = MMG5_Set_singulVertex(*sing,*c0,*c1,*c2,*typ,*pos);
    return;
}

/**
 * See \ref MMG5_Get_singulVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SINGULEDGE,mmg5_set_singuledge,
             (MMG5_pSingul *sing, int *v0, int *v1, int *ref, int *pos, int* retval),
             (sing,v0,v1,ref,pos,retval)) {
    *retval = MMG5_Set_singulEdge(*sing,*v0,*v1,*ref,*pos);
    return;
}

/**
 * See \ref MMG5_Set_singulCorner function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SINGULCORNER,mmg5_set_singulcorner,
             (MMG5_pSingul *sing, int *k, int* retval),
             (sing,k,retval)) {
    *retval = MMG5_Set_singulCorner(*sing,*k);
    return;
}

/**
 * See \ref MMG5_Set_singulRequiredVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SINGULREQUIREDVERTEX,mmg5_set_singulrequiredvertex,
             (MMG5_pSingul *sing, int *k, int* retval),
             (sing,k,retval)) {
    *retval = MMG5_Set_singulRequiredVertex(*sing,*k);
    return;
}

/**
 * See \ref MMG5_Set_singulRidge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SINGULRIDGE,mmg5_set_singulridge,
             (MMG5_pSingul *sing, int *k, int* retval),
             (sing,k,retval)) {
    *retval = MMG5_Set_singulRidge(*sing,*k);
    return;
}

/**
 * See \ref MMG5_Set_singulRequiredEdge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SINGULREQUIREDEDGE,set_singulrequirededge,
             (MMG5_pSingul *sing, int *k, int* retval),
             (sing,k,retval)) {
    *retval =  MMG5_Set_singulRequiredEdge(*sing,*k);
    return;
}
#endif

/**
 * See \ref MMG5_Set_handGivenMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_HANDGIVENMESH,mmg5_set_handgivenmesh,
             (MMG5_pMesh *mesh),
             (mesh)) {
    MMG5_Set_handGivenMesh(*mesh);
    return;
}

/**
 * See \ref MMG5_Chk_meshData function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_CHK_MESHDATA,mmg5_chk_meshdata,
             (MMG5_pMesh *mesh,MMG5_pSol *met, int* retval),
             (mesh,met,retval)) {
    *retval = MMG5_Chk_meshData(*mesh,*met);
    return;
}

/**
 * See \ref MMG5_Set_iparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_IPARAMETER,mmg5_set_iparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *iparam, int *val, int* retval),
             (mesh,sol,iparam,val,retval)){
    *retval = MMG5_Set_iparameter(*mesh,*sol,*iparam,*val);
    return;
}

/**
 * See \ref MMG5_Set_dparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_DPARAMETER,mmg5_set_dparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *dparam, double *val, int* retval),
             (mesh,sol,dparam,val,retval)){
    *retval = MMG5_Set_dparameter(*mesh,*sol,*dparam,*val);
    return;
}

/**
 * See \ref MMG5_Set_localParameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_LOCALPARAMETER,mmg5_set_localparameter,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, int *typ, int *ref, double *val, int* retval),
             (mesh,sol,typ,ref,val,retval)){
    *retval = MMG5_Set_localParameter(*mesh,*sol,*typ,*ref,*val);
    return;
}


/**
 * See \ref MMG5_Free_names function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_FREE_NAMES,mmg5_free_names,(MMG5_pMesh *mesh,MMG5_pSol *met
#ifdef SINGUL
                                              ,MMG5_pSingul *singul
#endif
                 ),(mesh,met
#ifdef SINGUL
                    ,singul
#endif
                     )){
#ifdef SINGUL
    MMG5_Free_names(*mesh,*met,*singul);
#else
    MMG5_Free_names(*mesh,*met);
#endif
    return;
}

/**
 * See \ref MMG5_Free_structures function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_FREE_STRUCTURES,mmg5_free_structures,(MMG5_pMesh *mesh,MMG5_pSol *met
#ifdef SINGUL
                                                        ,MMG5_pSingul *singul
#endif
                 ),(mesh,met
#ifdef SINGUL
                    ,singul
#endif
                     )){
#ifdef SINGUL
    MMG5_Free_structures(*mesh,*met,*singul);
#else
    MMG5_Free_structures(*mesh,*met);
#endif
    return;
}

/**
 * See \ref MMG5_loadMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_LOADMESH,mmg5_loadmesh,(MMG5_pMesh *mesh,int* retval),(mesh, retval)){

    *retval = MMG5_loadMesh(*mesh);

    return;
}

/**
 * See \ref MMG5_loadMet function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_LOADMET,mmg5_loadmet,(MMG5_pMesh *mesh,MMG5_pSol *met,int* retval),(mesh,met,retval)){

    *retval = MMG5_loadMet(*mesh,*met);

    return;
}

/**
 * See \ref MMG5_saveMet function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SAVEMET,mmg5_savemet,(MMG5_pMesh *mesh,MMG5_pSol *met,int* retval),(mesh,met,retval)){

    *retval = MMG5_saveMet(*mesh,*met);

    return;
}
