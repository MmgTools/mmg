/* =============================================================================
**  This file is part of the MMG3D 5 software package for the tetrahedral
**  mesh modification.
**  Copyright (c) 2014 Inria / Université de Bordeaux, IMB / UPMC, LJLL.
**
**  MMG3D 5 is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  MMG3D 5 is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with MMG3D 5 (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the MMG3D 5 distribution only if you accept them.
** =============================================================================
*/

/**
 *
 * Written by Cecile Dobrzynski (IMB), Charles Dapogny,
 * Pascal Frey (LJLL) and Algiane Froehly
 * Copyright (c) 2004- IMB/LJLL.
 * All rights reserved.
 *
 * API functions for fortran call
 *
 * Integers parameters:
 *    MMG5_IPARAM_verbose           = [-10..10] , Tune level of verbosity;
 *    MMG5_IPARAM_mem               = [n/-1]    , Set maximal memory size to n Mbytes/keep the default value;
 *    MMG5_IPARAM_debug             = [1/0]     , Turn on/off debug mode;
 *    MMG5_IPARAM_angle             = [1/0]     , Turn on/off angle detection;
 *    MMG5_IPARAM_iso               = [1/0]     , Turn on/off levelset meshing;
 *    MMG5_IPARAM_noinsert          = [1/0]     , avoid/allow point insertion/deletion;
 *    MMG5_IPARAM_noswap            = [1/0]     , avoid/allow edge or face flipping;
 *    MMG5_IPARAM_nomove            = [1/0]     , avoid/allow point relocation;
 *    MMG5_IPARAM_numberOflocalParam= [n]       , number of local parameters;
 *    MMG5_IPARAM_renum             = [1/0]     , Turn on/off the renumbering using SCOTCH;
 *    MMG5_IPARAM_sing              = [1/0]     , Turn on/off the insertion of singularities
 *                                        (need to compile with -DSINGUL flag);
 *    MMG5_IPARAM_bucket            = [val]     , Specify the size of the bucket per dimension (Delaunay)
 *                                        (need to compile with PATTERN=NO);
 * Double parameters:
 *    MMG5_DPARAM_angleDetection   = [val]     , angle detection;
 *    MMG5_DPARAM_hmin             = [val]     , minimal mesh size;
 *    MMG5_DPARAM_hmax             = [val]     , maximal mesh size;
 *    MMG5_DPARAM_hausd = [val]     , control global Hausdorff distance
 *                                    (on all the boundary surfaces of the mesh);
 *    MMG5_DPARAM_hgrad            = [val]     , control gradation;
 *    MMG5_DPARAM_ls               = [val]     , level set value;
 **/

#include "mmg3d.h"

/** Macro from Scotch **/
#define FORTRAN_NAME(nu,nl,pl,pc) \
  void nu pl;			  \
  void nl pl			  \
  { nu pc; }			  \
  void nl##_ pl			  \
  { nu pc; }			  \
  void nl##__ pl		  \
  { nu pc; }			  \
  void nu pl

/** Allocate the mesh and sol structures and initialize it to default values */
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
  Init_mesh(mesh,sol,sing);
#else
  Init_mesh(mesh,sol);
#endif

  return;
}

/** Initialization of parameters stored in the Info structure */
FORTRAN_NAME(MMG5_INIT_PARAMETERS,mmg5_init_parameters,(pMesh *mesh),(mesh)) {
  Init_parameters(*mesh);
  return;
}

/** default values for file names */
FORTRAN_NAME(INIT_FILENAMES,init_filenames,(pMesh *mesh,pSol *sol
#ifdef SINGUL
					    ,pSingul *sing
#endif
					    ),(mesh,sol
#ifdef SINGUL
					       ,sing
#endif
					       )) {
#ifdef SINGUL
  Init_fileNames(*mesh,*sol,*sing);
#else
  Init_fileNames(*mesh,*sol);
#endif
  return;
}


/** Set the name of input mesh */
FORTRAN_NAME(MMG5_SET_INPUTMESHNAME, mmg5_set_inputmeshname,
	     (MMG5_pMesh *mesh, char* meshin, int *strlen, int* retval),
	     (mesh,meshin,strlen,retval)) {
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';
  *retval = Set_inputMeshName(*mesh,tmp);
  SAFE_FREE(tmp);

  return;
}

/** Set the name of input sol */
FORTRAN_NAME(MMG5_SET_INPUTSOLNAME, mmg5_set_inputsolname,
	     (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solin, int* strlen, int* retval),
	     (mesh,sol,solin,strlen,retval)) {

  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,solin,*strlen);
  tmp[*strlen] = '\0';
  *retval = Set_inputSolName(*mesh,*sol,tmp);
  SAFE_FREE(tmp);

  return;
}

/** Set the name of output mesh */
FORTRAN_NAME(MMG5_SET_OUTPUTMESHNAME,mmg5_set_outputmeshname,
	     (MMG5_pMesh *mesh, char* meshout, int* strlen,int* retval),
 	     (mesh,meshout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshout,*strlen);
  tmp[*strlen] = '\0';
  *retval = Set_outputMeshName(*mesh, tmp);
  SAFE_FREE(tmp);

  return;
}

/** Set the name of output sol */
FORTRAN_NAME(MMG5_SET_OUTPUTSOLNAME,mmg5_set_outputsolname,
	     (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solout,int* strlen, int* retval),
	     (mesh,sol,solout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,solout,*strlen);
  tmp[*strlen] = '\0';
  *retval = Set_outputSolName(*mesh,*sol,tmp);
  SAFE_FREE(tmp);

  return;
}

#ifdef SINGUL
/** Set the name of input singularities file */
FORTRAN_NAME(MMG5_SET_INPUTSINGULNAME,mmg5_set_inputsingulname,
	     (MMG5_pMesh *mesh,MMG5_pSingul *sing, char* singin, int* strlen, int* retval),
	     (mesh,sing,singin,strlen,retval)) {
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,singin,*strlen);
  tmp[*strlen] = '\0';
  *retval = Set_inputSingulName(*mesh,*sing,tmp);
  SAFE_FREE(tmp);

  return;
}
#endif


/** Set the solution number, dimension and type */
FORTRAN_NAME(MMG5_SET_SOLSIZE,mmg5_set_solsize,
	     (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity,
	      int* np, int* typSol, int* retval),
	     (mesh, sol, typEntity, np, typSol, retval)) {
  *retval = Set_solSize(*mesh,*sol,*typEntity,*np,*typSol);
  return;
}

/** Set the number of vertices, tetrahedra, triangles and edges of the mesh
    and allocate the associated tables. If call twice, reset the whole mesh to
    realloc it at the new size */
FORTRAN_NAME(MMG5_SET_MESHSIZE,mmg5_set_meshsize,
	     (MMG5_pMesh *mesh, int *np, int *ne, int *nt, int *na, int *retval),
	     (mesh,np,ne,nt,na,retval)) {
  *retval = Set_meshSize(*mesh,*np,*ne,*nt,*na);
  return;
}

#ifdef SINGUL
/** Set the number of singular vertices and edges and allocate
    the associated tables. */
FORTRAN_NAME(MMG5_SET_SINGULSIZE,mmg5_set_singulsze,
	     (MMG5_pMesh *mesh,MMG5_pSingul *sing, int *np, int *na, int *retval),
	     (mesh,sing,np,na,retval)) {
  *retval = Set_singulSize(*mesh,*sing, *np, *na);
  return;
}
#endif

/** Get the solution number, dimension and type */
FORTRAN_NAME(MMG5_GET_SOLSIZE,mmg5_get_solsize,
	     (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity, int* np, int* typSol, int* retval),
	     (mesh,sol,typEntity,np,typSol,retval)) {

  *retval = Get_solSize(*mesh,*sol,typEntity,np,typSol);
  return;
}

/** Get the number of vertices, tetrahedra, triangles and edges of the mesh. */
FORTRAN_NAME(MMG5_GET_MESHSIZE,mmg5_get_meshsize,
	     (MMG5_pMesh *mesh, int* np, int* ne, int* nt, int* na, int* retval),
	     (mesh,np,ne,nt, na,retval)) {

  *retval = Get_meshSize(*mesh,np,ne,nt,na);
  return;
}

/** Set vertex of coordinates c0,c1,c2 and reference ref at position pos
    in mesh structure */
FORTRAN_NAME(MMG5_SET_VERTEX,mmg5_set_vertex,
	     (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
	      int* pos, int* retval),
	     (mesh,c0,c1,c2,ref,pos,retval)) {

  *retval = Set_vertex(*mesh,*c0,*c1,*c2,*ref,*pos);
  return;
}

/** Get coordinates c0,c1,c2 and reference ref of next vertex of mesh  */
FORTRAN_NAME(MMG5_GET_VERTEX,mmg5_get_vertex,
	     (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
	      int* isCorner, int* isRequired, int* retval),
	     (mesh,c0,c1,c2,ref,isCorner,isRequired, retval)) {
  *retval = Get_vertex(*mesh,c0,c1,c2,ref,isCorner,isRequired);
  return;
}

/** Set tetrahedra of vertices v0,v1,v2,v3 and reference ref at position pos
    in mesh structure */
FORTRAN_NAME(MMG5_SET_TETRAHEDRA,mmg5_set_tetrahedra,
	     (MMG5_pMesh *mesh, int *v0, int *v1, int *v2, int *v3, int *ref,
	      int *pos, int* retval),
	     (mesh,v0,v1,v2,v3,ref,pos,retval)){
  *retval = Set_tetrahedra(*mesh,*v0,*v1,*v2,*v3,*ref,*pos);
  return;
}

/** Get vertices v0,v1,v2,v3 and reference ref of next tetra of mesh  */
FORTRAN_NAME(MMG5_GET_TETRAHEDRA,mmg5_get_tetrahedra,
	     (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* v3,
	      int* ref, int* isRequired, int* retval),
	     (mesh,v0,v1,v2,v3,ref,isRequired,retval)) {
  *retval = Get_tetrahedra(*mesh,v0,v1,v2,v3,ref,isRequired);
  return;
}

/** Set triangle of vertices v0,v1,v2 and reference ref at position pos
    in mesh structure */
FORTRAN_NAME(MMG5_SET_TRIANGLE,mmg5_set_triangle,
	     (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref,int* pos,
	      int* retval),
	     (mesh,v0,v1,v2,ref,pos,retval)) {
  *retval = Set_triangle(*mesh, *v0, *v1, *v2, *ref, *pos);
  return;
}

/** Get vertices v0,v1,v2 and reference ref of next triangle of mesh  */
FORTRAN_NAME(MMG5_GET_TRIANGLE,mmg5_get_triangle,
	     (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref
	      ,int* isRequired, int* retval),
	     (mesh,v0,v1,v2,ref,isRequired,retval)) {
  *retval = Get_triangle(*mesh,v0,v1,v2,ref,isRequired);
  return;
}

/** Set edges of extremities v0,v1 and reference ref at position pos
    in mesh structure */
FORTRAN_NAME(MMG5_SET_EDGES,mmg5_set_edges,
	     (MMG5_pMesh *mesh, int *v0, int *v1, int *ref, int *pos, int* retval),
	     (mesh,v0,v1,ref,pos,retval)){
  *retval = Set_edges(*mesh,*v0,*v1,*ref,*pos);
  return;
}

/** Get extremities e0,e1 and reference ref of next edge of mesh  */
FORTRAN_NAME(MMG5_GET_EDGE,mmg5_get_edge,(MMG5_pMesh *mesh, int* e0, int* e1, int* ref
				 ,int* isRidge, int* isRequired, int* retval),
	     (mesh,e0,e1,ref,isRidge,isRequired,retval)) {
  *retval = Get_edge(*mesh,e0,e1,ref,isRidge,isRequired);
  return;
}

/** Set corner at point k  */
FORTRAN_NAME(MMG5_SET_CORNER,mmg5_set_corner,(MMG5_pMesh *mesh, int *k, int* retval),
	     (mesh,k,retval)) {
  *retval =  Set_corner(*mesh,*k);
  return;
}

/** Set point k as required  */
FORTRAN_NAME(MMG5_SET_REQUIREDVERTEX,mmg5_set_requiredvertex,
	     (MMG5_pMesh *mesh, int *k, int* retval),
	     (mesh,k,retval)) {
  *retval =  Set_requiredVertex(*mesh,*k);
  return;
}

/** Set element k as required  */
FORTRAN_NAME(MMG5_SET_REQUIREDTETRAHEDRA,mmg5_set_requiredtetrahedra,
	     (MMG5_pMesh *mesh, int *k, int* retval),
	     (mesh,k,retval)) {
  *retval = Set_requiredTetrahedra(*mesh,*k);
  return;
}

/** Set triangle k as required  */
FORTRAN_NAME(MMG5_SET_REQUIREDTRIANGLE,mmg5_set_requiredtriangle,
	     (MMG5_pMesh *mesh, int *k, int* retval),
	     (mesh,k,retval)) {
  *retval = Set_requiredTriangle(*mesh, *k);
  return;
}

/** Set ridge at edge k  */
FORTRAN_NAME(MMG5_SET_RIDGE,mmg5_set_ridge,
	     (MMG5_pMesh *mesh, int *k, int* retval),
	     (mesh,k,retval)) {
  *retval = Set_ridge(*mesh,*k);
  return;
}

/** Set edge k as required  */
FORTRAN_NAME(MMG5_SET_REQUIREDEDGE,mmg5_set_requirededge,
	     (MMG5_pMesh *mesh, int *k, int* retval),
	     (mesh,k,retval)) {
  *retval = Set_requiredEdge(*mesh,*k);
  return;
}

/** Set scalar value s at position pos in solution structure */
FORTRAN_NAME(MMG5_SET_SCALARSOL,mmg5_set_scalarsol,
	     (MMG5_pSol *met, double *s, int *pos, int* retval),
	     (met,s,pos,retval)) {
  *retval = Set_scalarSol(*met,*s,*pos);
  return;
}

/** Get solution s of next vertex of mesh  */
FORTRAN_NAME(MMG5_GET_SCALARSOL,mmg5_get_scalarsol,
	     (MMG5_pSol *met, double* s, int* retval),
	     (met,s,retval)) {
  *retval =  Get_scalarSol(*met,s);
  return;
}

#ifdef SINGUL
/** Set singular point of coordinates c0,c1,c2 at position pos
    in the singularities structure */
FORTRAN_NAME(MMG5_SET_SINGULVERTEX,mmg5_set_singluvertex,
	     (MMG5_pSingul *sing, double *c0, double *c1,
	      double *c2, int *typ, int *pos, int* retval),
	     (sing,c0,c1,c2,typ,pos,retval)) {
  *retval = Set_singulVertex(*sing,*c0,*c1,*c2,*typ,*pos);
  return;
}


/** Set singular edge of extremities v0,v1 and reference ref
    at position pos in the singularities structure */
FORTRAN_NAME(MMG5_SET_SINGULEDGE,mmg5_set_singuledge,
	     (MMG5_pSingul *sing, int *v0, int *v1, int *ref, int *pos, int* retval),
	     (sing,v0,v1,ref,pos,retval)) {
  *retval =  Set_singulEdge(*sing,*v0,*v1,*ref,*pos);
  return;
}

/** Set corner at singular vertex k  */
FORTRAN_NAME(MMG5_SET_SINGULCORNER,mmg5_set_singulcorner,
	     (MMG5_pSingul *sing, int *k, int* retval),
	     (sing,k,retval)) {
  *retval = Set_singulCorner(*sing,*k);
  return;
}

/** Set required vertex at singular vertex k  */
FORTRAN_NAME(MMG5_SET_SINGULREQUIREDVERTEX,mmg5_set_singulrequiredvertex,
	     (MMG5_pSingul *sing, int *k, int* retval),
	     (sing,k,retval)) {
  *retval = Set_singulRequiredVertex(*sing,*k);
  return;
}

/** Set required edge at singular edge k  */
FORTRAN_NAME(MMG5_SET_SINGULRIDGE,mmg5_set_singulridge,
	     (MMG5_pSingul *sing, int *k, int* retval),
	     (sing,k,retval)) {
  *retval =  Set_singulRidge(*sing,*k);
  return;
}

/** Set required edge at singular edge k  */
FORTRAN_NAME(MMG5_SET_SINGULREQUIREDEDGE,set_singulrequirededge,
	     (MMG5_pSingul *sing, int *k, int* retval),
	     (sing,k,retval)) {
  *retval =  Set_singulRequiredEdge(*sing,*k);
  return;
}
#endif

/** To mark as ended a mesh given without using the API functions
    (for example, mesh given by mesh->point[i] = 0 ...)*/
FORTRAN_NAME(MMG5_SET_HANDGIVENMESH,mmg5_set_handgivenmesh,
	     (MMG5_pMesh *mesh),
	     (mesh)) {
  Set_handGivenMesh(*mesh);
  return;
}

/** Check if the number of given entities match with mesh and sol size (not mandatory)
    and check mesh datas  */
FORTRAN_NAME(MMG5_CHK_MESHDATA,mmg5_chk_meshdata,
	     (MMG5_pMesh *mesh,MMG5_pSol *met, int* retval),
	     (mesh,met,retval)) {
  *retval = Chk_meshData(*mesh,*met);
  return;
}

/** Set integer parameter iparam at value val */
FORTRAN_NAME(MMG5_SET_IPARAMETERS,mmg5_set_iparameters,
	     (MMG5_pMesh *mesh, MMG5_pSol *sol, int *iparam, int *val, int* retval),
	     (mesh,sol,iparam,val,retval)){
  *retval = Set_iparameters(*mesh,*sol,*iparam,*val);
  return;
}

/** Set double parameter dparam at value val */
FORTRAN_NAME(MMG5_SET_DPARAMETERS,mmg5_set_dparameters,
	     (MMG5_pMesh *mesh, MMG5_pSol *sol, int *dparam, double *val, int* retval),
	     (mesh,sol,dparam,val,retval)){
  *retval = Set_dparameters(*mesh,*sol,*dparam,*val);
  return;
}

/** Set local parameters: set the hausdorff value at val for all elements
    of type typ and reference ref  */
FORTRAN_NAME(MMG5_SET_LOCALPARAMETERS,mmg5_set_localparameters,
	     (MMG5_pMesh *mesh,MMG5_pSol *sol, int *typ, int *ref, double *val, int* retval),
	     (mesh,sol,typ,ref,val,retval)){
  *retval = Set_localParameters(*mesh,*sol,*typ,*ref,*val);
  return;
}


/** File names deallocations before return */
FORTRAN_NAME(MMG5_FREE_NAMES,mmg5_free_names,(pMesh *mesh,pSol *met
#ifdef SINGUL
				    ,pSingul *singul
#endif
				    ),(mesh,met
#ifdef SINGUL
				       ,singul
#endif
				       )){
#ifdef SINGUL
  Free_names(*mesh,*met,*singul);
#else
  Free_names(*mesh,*met);
#endif
  return;
}

/** Structure deallocations before return */
FORTRAN_NAME(MMG5_FREE_STRUCTURES,mmg5_free_structures,(pMesh *mesh,pSol *met
#ifdef SINGUL
					      ,pSingul *singul
#endif
					      ),(mesh,met
#ifdef SINGUL
						 ,singul
#endif
						 )){
#ifdef SINGUL
  Free_structures(*mesh,*met,*singul);
#else
  Free_structures(*mesh,*met);
#endif
  return;
}

FORTRAN_NAME(MMG5_LOADMESH,mmg5_loadmesh,(pMesh *mesh,int* retval),(mesh, retval)){

  *retval = loadMesh(*mesh);

  return;
}

FORTRAN_NAME(MMG5_LOADMET,mmg5_loadmet,(pMesh *mesh,pSol *met,int* retval),(mesh,met,retval)){

  *retval = loadMet(*mesh,*met);

  return;
}

FORTRAN_NAME(MMG5_SAVEMET,mmg5_savemet,(pMesh *mesh,pSol *met,int* retval),(mesh,met,retval)){

  *retval = saveMet(*mesh,*met);

  return;
}
