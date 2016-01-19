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
 * \file mmgs/libmmgs.h
 * \brief API headers for the mmgs library
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning To keep the genheader working, don't break line between the enum
 * name and the opening brace (it creates errors under windows)
 *
 */

#ifndef _MMGSLIB_H
#define _MMGSLIB_H

#include "mmgcommon.h"

/**
 * Maximum array size when storing adjacent points (or ball) of a vertex.
 */
#define MMGS_LMAX      1024

/**
 * \enum MMGS_Param
 * \brief Input parameters for mmg library.
 *
 * Input parameters for mmg library. Options prefixed by \a
 * MMGS_IPARAM asked for integers values ans options prefixed by \a
 * MMGS_DPARAM asked for real values.
 *
 */
enum MMGS_Param {
  MMGS_IPARAM_verbose,           /*!< [-10..10], Tune level of verbosity */
  MMGS_IPARAM_mem,               /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  MMGS_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  MMGS_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  MMGS_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
  MMGS_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  MMGS_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
  MMGS_IPARAM_nreg,              /*!< [0/1], Disabled/enabled normal regularization */
  MMGS_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
  MMGS_IPARAM_renum,             /*!< [1/0], Turn on/off point relocation with Scotch */
  MMGS_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
  MMGS_DPARAM_hmin,              /*!< [val], Minimal mesh size */
  MMGS_DPARAM_hmax,              /*!< [val], Maximal mesh size */
  MMGS_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  MMGS_DPARAM_hgrad,             /*!< [val], Control gradation */
  MMGS_PARAM_size,               /*!< [n], Number of parameters */
};

/*----------------------------- functions header -----------------------------*/
/* Initialization functions */
/* init structures */
/**
 * \param starter dummy argument used to initialize the variadic argument list
 * \param ... variadic arguments. For now, you need to call the \a
 * MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh,
 * MMG5_ARG_ppMet, your_metric,MMG5_ARG_end). Here, \a your_mesh is a pointer
 * toward \a MMG5_pMesh and \a your_metric a pointer toward \a MMG5_pSol.
 *
 * MMG structures allocation and initialization.
 *
 */
void MMGS_Init_mesh(enum MMG5_arg starter,...);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Initialize file names to their default values.
 *
 */
void  MMGS_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 */
void  MMGS_Init_parameters(MMG5_pMesh mesh);

/* init file names */
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * Set the name of input mesh.
 *
 */
int  MMGS_Set_inputMeshName(MMG5_pMesh mesh, char* meshin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * Set the name of output mesh file.
 *
 */
int  MMGS_Set_outputMeshName(MMG5_pMesh mesh, char* meshout);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * Set the name of input solution file.
 *
 */
int  MMGS_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solout name of the output solution file.
 * \return 0 if failed, 1 otherwise.
 *
 *  Set the name of output solution file.
 *
 */
int  MMGS_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solout);

/* init structure sizes */
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typEntity type of solutions entities (vertices, triangles...).
 * \param np number of solutions.
 * \param typSol type of solution (scalar, vectorial...).
 * \return 0 if failed, 1 otherwise.
 *
 * Set the solution number, dimension and type.
 *
 */
int  MMGS_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, int np, int typSol);
/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param nt number of triangles.
 * \param na number of edges.
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of vertices, triangles and edges of the
 * mesh and allocate the associated tables. If call twice, reset the
 * whole mesh to realloc it at the new size
 *
 */
int  MMGS_Set_meshSize(MMG5_pMesh mesh, int np, int nt, int na);

/* init structure datas */
/**
 * \param mesh pointer toward the mesh structure.
 * \param c0 coordinate of the point along the first dimension.
 * \param c1 coordinate of the point along the second dimension.
 * \param c2 coordinate of the point along the third dimension.
 * \param ref point reference.
 * \param pos position of the point in the mesh.
 * \return 1.
 *
 * Set vertex of coordinates \a c0, \a c1,\a c2 and reference \a ref
 * at position \a pos in mesh structure
 *
 */
int  MMGS_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                     double c2, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of triangle.
 * \param v1 second vertex of triangle.
 * \param v2 third vertex of triangle.
 * \param ref triangle reference.
 * \param pos triangle position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set triangle of vertices \a v0, \a v1, \a v2 and reference \a ref
 * at position \a pos in mesh structure.
 *
 */
int  MMGS_Set_triangle(MMG5_pMesh mesh, int v0, int v1,
                       int v2, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first extremity of the edge.
 * \param v1 second extremity of the edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set edges of extremities \a v0, \a v1 and reference \a ref at
 * position \a pos in mesh structure
 *
 */
int  MMGS_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set corner at point \a pos.
 *
 */
int  MMGS_Set_corner(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set point \a k as required.
 *
 */
int  MMGS_Set_requiredVertex(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Set triangle \a k as required.
 *
 */
int  MMGS_Set_requiredTriangle(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set ridge at edge \a k.
 *
 */
int  MMGS_Set_ridge(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set edge \a k as required.
 *
 */
int  MMGS_Set_requiredEdge(MMG5_pMesh mesh, int k);
/**
 * \param met pointer toward the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure
 *
 */
int  MMGS_Set_scalarSol(MMG5_pSol met, double s,int pos);
/**
 * \param met pointer toward the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set vectorial value \f$(v_x,v_y,v_z)\f$ at position \a pos in solution
 * structure.
 *
 */
int MMGS_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz, int pos);
/**
 * \param met pointer toward the sol structure.
 * \param m11 value of the tensorial solution at position (1,1) in the tensor.
 * \param m12 value of the tensorial solution at position (1,2) in the tensor.
 * \param m13 value of the tensorial solution at position (1,3) in the tensor.
 * \param m22 value of the tensorial solution at position (2,2) in the tensor.
 * \param m23 value of the tensorial solution at position (2,3) in the tensor.
 * \param m33 value of the tensorial solution at position (3,3) in the tensor.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in solution
 * structure.
 *
 */
int MMGS_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos);

/* check init */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Check if the number of given entities match with mesh and sol size
 * (not mandatory) and check mesh datas.
 *
 */
int MMGS_Chk_meshData(MMG5_pMesh mesh, MMG5_pSol met);

/** functions to set parameters */
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param iparam integer parameter to set (see \a MMGS_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 */
int  MMGS_Set_iparameter(MMG5_pMesh mesh,MMG5_pSol sol, int iparam, int val);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param dparam double parameter to set (see \a MMGS_Param structure).
 * \param val value of the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
 *
 */
int  MMGS_Set_dparameter(MMG5_pMesh mesh,MMG5_pSol sol, int dparam, double val);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typ type of entity (triangle, edge,...).
 * \param ref reference of the entity.
 * \param hmin minimal edge size.
 * \param hmax maximal edge size.
 * \param hausd value of the Hausdorff number.
 * \return 0 if failed, 1 otherwise.
 *
 * Set local parameters: set the hausdorff value at \a hausd, the minmal edge
 * size value at \a hmin and the maximal edge size value at \a hmax for all
 * elements of type \a typ and reference \a ref.
 *
 */
int  MMGS_Set_localParameter(MMG5_pMesh mesh, MMG5_pSol sol, int typ, int ref,
                             double hmin, double hmax, double hausd);

/** recover datas */
/**
 * \param mesh pointer toward the mesh structure.
 * \param np pointer toward the number of vertices.
 * \param nt pointer toward the number of triangles.
 * \param na pointer toward the number of edges.
 * \return 1.
 *
 * Get the number of vertices, triangles and edges of the mesh.
 *
 */
int  MMGS_Get_meshSize(MMG5_pMesh mesh, int* np, int* nt, int* na);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typEntity pointer toward the type of entities to which solutions are applied.
 * \param np pointer toward the number of solutions.
 * \param typSol pointer toward the type of the solutions (scalar, vectorial...)
 * \return 1.
 *
 * Get the solution number, dimension and type.
 *
 */
int  MMGS_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, int* np,
                      int* typSol);
/**
 * \param mesh pointer toward the mesh structure.
 * \param c0 pointer toward the coordinate of the point along the first dimension.
 * \param c1 pointer toward the coordinate of the point along the second dimension.
 * \param c2 pointer toward the coordinate of the point along the third dimension.
 * \param ref pointer to the point reference.
 * \param isCorner pointer toward the flag saying if point is corner.
 * \param isRequired pointer toward the flag saying if point is required.
 * \return 1.
 *
 * Get coordinates \a c0, \a c1,\a c2 and reference \a ref of next
 * vertex of mesh.
 *
 */
int  MMGS_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, int* ref,
                     int* isCorner, int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 pointer toward the first vertex of triangle.
 * \param v1 pointer toward the second vertex of triangle.
 * \param v2 pointer toward the third vertex of triangle.
 * \param ref pointer toward the triangle reference.
 * \param isRequired pointer toward the flag saying if triangle is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0,\a v1,\a v2 and reference \a ref of next
 * triangle of mesh.
 *
 */
int  MMGS_Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref,
                       int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param e0 pointer toward the first extremity of the edge.
 * \param e1 pointer toward the second  extremity of the edge.
 * \param ref pointer toward the edge reference.
 * \param isRidge pointer toward the flag saying if the edge is ridge.
 * \param isRequired pointer toward the flag saying if the edge is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get extremities \a e0, \a e1 and reference \a ref of next edge of mesh.
 *
 */
int  MMGS_Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired);
/**
 * \param met pointer toward the sol structure.
 * \param s pointer toward the scalar solution value.
 * \return 0 if failed, 1 otherwise.
 *
 * Get solution \a s of next vertex of mesh.
 *
 */
int  MMGS_Get_scalarSol(MMG5_pSol met, double* s);
/**
 * \param met pointer toward the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vectorial solution \f$(v_x,v_y,vz)\f$ of next vertex of mesh.
 *
 */
int MMGS_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz);
/**
 * \param met pointer toward the sol structure.
 * \param m11 pointer toward the position (1,1) in the solution tensor.
 * \param m12 pointer toward the position (1,2) in the solution tensor.
 * \param m13 pointer toward the position (1,3) in the solution tensor.
 * \param m22 pointer toward the position (2,2) in the solution tensor.
 * \param m23 pointer toward the position (2,3) in the solution tensor.
 * \param m33 pointer toward the position (3,3) in the solution tensor.
 * \return 0 if failed, 1 otherwise.
 *
 * Get tensorial solution of next vertex of mesh.
 *
 */
int MMGS_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33);
/**
 * \param mesh pointer toward the mesh structure.
 * \param iparam integer parameter to set (see \a MMGS_Param structure).
 * \return The value of integer parameter.
 *
 * Get the value of integer parameter \a iparam.
 *
 */
int MMGS_Get_iparameter(MMG5_pMesh mesh, int iparam);

/* input/output functions */
/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh data.
 *
 */
int  MMGS_loadMesh(MMG5_pMesh mesh, char* filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 */
int  MMGS_saveMesh(MMG5_pMesh mesh, char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Load metric field.
 *
 */
int  MMGS_loadSol(MMG5_pMesh mesh,MMG5_pSol met, char* filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric.
 *
 */
int  MMGS_saveSol(MMG5_pMesh mesh, MMG5_pSol met, char *filename);

/* deallocations */
/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments. For now, you need to call the \a
 * MMGS_Free_all function with the following arguments :
 * MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh,
 * MMG5_ARG_ppMet, your_metric,MMG5_ARG_end). Here, \a your_mesh is a pointer
 * toward \a MMG5_pMesh and \a your_metric a pointer toward \a MMG5_pSol.
 *
 * Deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void MMGS_Free_all(enum MMG5_arg starter,...);

/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments. For now, you need to call the \a
 * MMGS_Free_structures function with the following arguments :
 * MMGS_Free_structures(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh,
 * MMG5_ARG_ppMet, your_metric,MMG5_ARG_end). Here, \a your_mesh is a pointer
 * toward \a MMG5_pMesh and \a your_metric a pointer toward \a MMG5_pSol.
 *
 * Structure deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void MMGS_Free_structures(enum MMG5_arg starter,...);

/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments. For now, you need to call the \a
 * MMGS_Free_names function with the following arguments :
 * MMGS_Free_names(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh,
 * MMG5_ARG_ppMet, your_metric,MMG5_ARG_end). Here, \a your_mesh is a pointer
 * toward \a MMG5_pMesh and \a your_metric a pointer toward \a MMG5_pSol.
 *
 * Structure deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void MMGS_Free_names(enum MMG5_arg starter,...);

/* library */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol (metric) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the library.
 *
 */
int  MMGS_mmgslib(MMG5_pMesh mesh, MMG5_pSol met);

/** To associate function pointers without calling MMGS_mmg3dlib */
/**
 * \param mesh pointer toward the mesh structure (unused).
 *
 * Set pointer for MMGS_saveMesh function.
 *
 */
void  MMGS_Set_saveFunc(MMG5_pMesh mesh);
/**
 * \param mesh pointer toward the mesh structure (unused).
 * \param met pointer toward the sol structure (unused).
 *
 * Set function pointers for caltet, lenedg, defsiz and gradsiz.
 *
 */
void  MMGS_setfunc(MMG5_pMesh mesh,MMG5_pSol met);

/* Tools for the library */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * Set function pointers depending if case is iso or aniso.
 *
 */
void MMGS_setfunc(MMG5_pMesh mesh,MMG5_pSol met);
/**
 * \param prog pointer toward the program name.
 *
 * Print help for mmgs options.
 *
 */
void MMGS_usage(char *prog);
/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
void MMGS_defaultValues(MMG5_pMesh mesh);
/**
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 * \return 1.
 *
 * Store the info structure in the mesh structure.
 *
 */
int MMGS_stockOptions(MMG5_pMesh mesh, MMG5_Info *info);
/**
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 *
 * Recover the info structure stored in the mesh structure.
 *
 */
void MMGS_destockOptions(MMG5_pMesh mesh, MMG5_Info *info);

/**
 * \brief Return adjacent elements of a triangle.
 * \param mesh pointer toward the mesh structure.
 * \param kel triangle index.
 * \param listri pointer toward the table of the indices of the three adjacent
 * triangles of the elt \a kel (the index is 0 if there is no adjacent).
 * \return 1.
 *
 * Find the indices of the 3 adjacent elements of triangle \a
 * kel. \f$v_i = 0\f$ if the \f$i^{th}\f$ face has no adjacent element
 * (so we are on a boundary face).
 *
 */
int MMGS_Get_adjaTri(MMG5_pMesh mesh, int kel, int listri[3]);

/**
 * \brief Return adjacent elements of a triangle.
 * \param mesh pointer toward the mesh structure.
 * \param ip vertex index.
 * \param start index of a triangle holding \a ip.
 * \param lispoi pointer toward an array of size MMGS_LMAX that will contain
 * the indices of adjacent vertices to the vertex \a ip.
 * \return nbpoi the number of adjacent points if success, 0 if fail.
 *
 * Find the indices of the adjacent vertices of the vertex \a
 * ip of the triangle \a start.
 *
 */
extern
int MMGS_Get_adjaVerticesFast(MMG5_pMesh mesh, int ip,int start, int lispoi[MMGS_LMAX]);

#endif
