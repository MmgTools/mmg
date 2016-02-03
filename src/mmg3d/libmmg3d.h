/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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
 * \file mmg3d/libmmg3d.h
 * \brief API headers for the mmg3d library
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning To keep the genheader working, don't break line between the enum
 * name and the opening brace (it creates errors under windows)
 * \warning Use the MMG3D_ prefix: MMG5_ prefix will became obsolete...
 *
 * \include libexamples/mmg3d/adaptation_example0/example0_a/main.c
 * \include libexamples/mmg3d/example0/adaptation_example0_b/main.c
 * \include libexamples/mmg3d/adaptation_example0_fortran/example0_a/main.F90
 * \include libexamples/mmg3d/adaptation_example0_fortran/example0_b/main.F90
 * \include libexamples/mmg3d/adaptation_example1/main.c
 * \include libexamples/mmg3d/adaptation_example2/main.c
 * \include libexamples/mmg3d/IsosurfDiscretization_example0/main.c
 * \include libexamples/mmg3d/LagrangianMotion_example0/main.c
 */

#ifndef _MMG3DLIB_H
#define _MMG3DLIB_H

#include "mmgcommon.h"

/**
 * Maximum array size when storing adjacent points (or ball) of a vertex.
 */
#define MMG3D_LMAX      10240

/**
 * \enum MMG3D_Param
 * \brief Input parameters for mmg library.
 *
 * Input parameters for mmg library. Options prefixed by \a
 * MMG3D_IPARAM asked for integers values ans options prefixed by \a
 * MMG3D_DPARAM asked for real values.
 *
 */
enum MMG3D_Param {
  MMG3D_IPARAM_verbose,           /*!< [-10..10], Tune level of verbosity */
  MMG3D_IPARAM_mem,               /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  MMG3D_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  MMG3D_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  MMG3D_IPARAM_iso,               /*!< [1/0], Level-set meshing */
  MMG3D_IPARAM_lag,               /*!< [-1/0/1/2], Lagrangian option */
  MMG3D_IPARAM_optim,             /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  MMG3D_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
  MMG3D_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  MMG3D_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
  MMG3D_IPARAM_nosurf,            /*!< [1/0], Avoid/allow surface modifications */
  MMG3D_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
  MMG3D_IPARAM_renum,             /*!< [1/0], Turn on/off point relocation with Scotch */
  MMG3D_IPARAM_bucket,            /*!< [n], Specify the size of the bucket per dimension (DELAUNAY) */
  MMG3D_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
  MMG3D_DPARAM_hmin,              /*!< [val], Minimal mesh size */
  MMG3D_DPARAM_hmax,              /*!< [val], Maximal mesh size */
  MMG3D_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  MMG3D_DPARAM_hgrad,             /*!< [val], Control gradation */
  MMG3D_DPARAM_ls,                /*!< [val], Value of level-set (not use for now) */
  MMG3D_PARAM_size,               /*!< [n], Number of parameters */
};

/*----------------------------- functions header -----------------------------*/
/* Initialization functions */
/* init structures */
/**
 * \param starter dummy argument used to initialize the variadic argument list
 * \param ... variadic arguments that depend to the library function that you
 * want to call. For the MMG3D_mmg3dlib or the MMG3D_mmg3dls functions, you need
 * to call the \a MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppSol,
 * &your_metric,MMG5_ARG_end). For the MMG3D_mmg3dmov function, you must call
 * : MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppSol,
 * &your_metric,MMG5_ARG_ppDisp, &your_displacement,MMG5_ARG_end). Here,
 * \a your_mesh is a \a MMG5_pMesh, \a your_metric a \a MMG5_pSol and \a
 * your_displacement a \a MMG5_pSol.
 *
 * MMG structures allocation and initialization.
 *
 */
void MMG3D_Init_mesh(enum MMG5_arg starter,...);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Initialize file names to their default values.
 *
 */
void  MMG3D_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 */
void  MMG3D_Init_parameters(MMG5_pMesh mesh);

/* init file names */
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * Set the name of input mesh.
 *
 */
int  MMG3D_Set_inputMeshName(MMG5_pMesh mesh, char* meshin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * Set the name of output mesh file.
 *
 */
int  MMG3D_Set_outputMeshName(MMG5_pMesh mesh, char* meshout);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * Set the name of input solution file.
 *
 */
int  MMG3D_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solout name of the output solution file.
 * \return 0 if failed, 1 otherwise.
 *
 *  Set the name of output solution file.
 *
 */
int  MMG3D_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solout);

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
int  MMG3D_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, int np, int typSol);
/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param ne number of elements (tetrahedra).
 * \param nt number of triangles.
 * \param na number of edges.
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of vertices, tetrahedra, triangles and edges of the
 * mesh and allocate the associated tables. If call twice, reset the
 * whole mesh to realloc it at the new size
 *
 */
int  MMG3D_Set_meshSize(MMG5_pMesh mesh, int np, int ne, int nt, int na);

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
int  MMG3D_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                     double c2, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of tetrahedron.
 * \param v1 second vertex of tetrahedron.
 * \param v2 third vertex of tetrahedron.
 * \param v3 fourth vertex of tetrahedron.
 * \param ref tetrahedron reference.
 * \param pos tetrahedron position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set tetrahedra of vertices \a v0, \a v1,\a v2,\a v3 and reference
 * \a ref at position \a pos in mesh structure.
 *
 */
int  MMG3D_Set_tetrahedron(MMG5_pMesh mesh, int v0, int v1,
                          int v2, int v3, int ref, int pos);
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
int  MMG3D_Set_triangle(MMG5_pMesh mesh, int v0, int v1,
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
int  MMG3D_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set corner at point \a pos.
 *
 */
int  MMG3D_Set_corner(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set point \a k as required.
 *
 */
int  MMG3D_Set_requiredVertex(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k element index.
 * \return 1.
 *
 * Set element \a k as required.
 *
 */
int  MMG3D_Set_requiredTetrahedron(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Set triangle \a k as required.
 *
 */
int  MMG3D_Set_requiredTriangle(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set ridge at edge \a k.
 *
 */
int  MMG3D_Set_ridge(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set edge \a k as required.
 *
 */
int  MMG3D_Set_requiredEdge(MMG5_pMesh mesh, int k);
/**
 * \param met pointer toward the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure
 *
 */
int  MMG3D_Set_scalarSol(MMG5_pSol met, double s,int pos);
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
int MMG3D_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz, int pos);
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
int MMG3D_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * To mark as ended a mesh given without using the API functions
 * (for example, mesh given by mesh->point[i] = 0 ...). Not recommanded.
 *
 */
void MMG3D_Set_handGivenMesh(MMG5_pMesh mesh);

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
int MMG3D_Chk_meshData(MMG5_pMesh mesh, MMG5_pSol met);

/** functions to set parameters */
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param iparam integer parameter to set (see \a MMG3D_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 */
int  MMG3D_Set_iparameter(MMG5_pMesh mesh,MMG5_pSol sol, int iparam, int val);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param dparam double parameter to set (see \a MMG3D_Param structure).
 * \param val value of the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
 *
 */
int  MMG3D_Set_dparameter(MMG5_pMesh mesh,MMG5_pSol sol, int dparam, double val);
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
 * Set local parameters: set the hausdorff value at \a val for all
 * elements of type \a typ and reference \a ref.
 *
 */
int  MMG3D_Set_localParameter(MMG5_pMesh mesh, MMG5_pSol sol, int typ, int ref,
                              double hmin,double hmax,double hausd);

/** recover datas */
/**
 * \param mesh pointer toward the mesh structure.
 * \param np pointer toward the number of vertices.
 * \param ne pointer toward the number of elements (tetrahedra).
 * \param nt pointer toward the number of triangles.
 * \param na pointer toward the number of edges.
 * \return 1.
 *
 * Get the number of vertices, tetrahedra, triangles and edges of the mesh.
 *
 */
int  MMG3D_Get_meshSize(MMG5_pMesh mesh, int* np, int* ne, int* nt, int* na);
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
int  MMG3D_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, int* np,
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
int  MMG3D_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, int* ref,
                     int* isCorner, int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 pointer toward the first vertex of tetrahedron.
 * \param v1 pointer toward the second vertex of tetrahedron.
 * \param v2 pointer toward the third vertex of tetrahedron.
 * \param v3 pointer toward the fourth vertex of tetrahedron.
 * \param ref pointer toward the tetrahedron reference.
 * \param isRequired pointer toward the flag saying if tetrahedron is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0, \a v1, \a v2, \a v3 and reference \a ref of
 * next tetra of mesh.
 *
 */
int  MMG3D_Get_tetrahedron(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* v3,
                          int* ref, int* isRequired);
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
int  MMG3D_Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref,
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
int  MMG3D_Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired);
/**
 * \param met pointer toward the sol structure.
 * \param s pointer toward the scalar solution value.
 * \return 0 if failed, 1 otherwise.
 *
 * Get solution \a s of next vertex of mesh.
 *
 */
int  MMG3D_Get_scalarSol(MMG5_pSol met, double* s);
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
int MMG3D_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz);
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
int MMG3D_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33);
/**
 * \param mesh pointer toward the mesh structure.
 * \param iparam integer parameter to set (see \a MMG3D_Param structure).
 * \return The value of integer parameter.
 *
 * Get the value of integer parameter \a iparam.
 *
 */
int MMG3D_Get_iparameter(MMG5_pMesh mesh, int iparam);

/* input/output functions */
/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh data.
 *
 */
int MMG3D_loadMesh(MMG5_pMesh mesh,char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param filename pointer toward the name of file.

 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 */
int MMG3D_saveMesh(MMG5_pMesh mesh, char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load metric field.
 *
 */
int MMG3D_loadSol(MMG5_pMesh mesh,MMG5_pSol met, char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric.
 *
 */
int MMG3D_saveSol(MMG5_pMesh mesh,MMG5_pSol met, char *filename);

/* deallocations */
/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments that depend to the library function that you
 * have call. For the MMG3D_mmg3dlib or the MMG3D_mmg3dls functions, you need to
 * call the \a MMG3D_Free_all function with the following arguments :
 * MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh, MMG5_ARG_ppMet,
 * your_metric,MMG5_ARG_end). For the MMG3D_mmg3dmov function, you must call :
 * MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh, MMG5_ARG_ppMet,
 * your_metric,MMG5_ARG_ppDisp, your_displacement,MMG5_ARG_end). Here,
 * \a your_mesh is a pointer toward \a MMG5_pMesh, \a your_metric a pointer
 * toward \a MMG5_pSol and \a your_displacement a pointer toward \a MMG5_pSol.
 *
 * Deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void MMG3D_Free_all(enum MMG5_arg starter,...);

/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments that depend to the library function that you
 * have call. For the MMG3D_mmg3dlib or the MMG3D_mmg3dls functions, you need to
 * call the \a MMG3D_Free_structures function with the following arguments :
 * MMG3D_Free_structures(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh, MMG5_ARG_ppMet,
 * your_metric,MMG5_ARG_end). For the MMG3D_mmg3dmov function, you must call :
 * MMG3D_Free_structures(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh, MMG5_ARG_ppMet,
 * your_metric,MMG5_ARG_ppDisp, your_displacement,MMG5_ARG_end). Here,
 * \a your_mesh is a pointer toward \a MMG5_pMesh, \a your_metric a pointer
 * toward \a MMG5_pSol and \a your_displacement a pointer toward \a MMG5_pSol.
 *
 * Structure deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void MMG3D_Free_structures(enum MMG5_arg starter,...);

/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments that depend to the library function that you
 * have call. For the MMG3D_mmg3dlib or the MMG3D_mmg3dls functions, you need to
 * call the \a MMG3D_Free_names function with the following arguments :
 * MMG3D_Free_names(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh, MMG5_ARG_ppMet,
 * your_metric,MMG5_ARG_end). For the MMG3D_mmg3dmov function, you must call :
 * MMG3D_Free_names(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh, MMG5_ARG_ppMet,
 * your_metric,MMG5_ARG_ppDisp, your_displacement,MMG5_ARG_end). Here,
 * \a your_mesh is a pointer toward \a MMG5_pMesh, \a your_metric a pointer
 * toward \a MMG5_pSol and \a your_displacement a pointer toward \a MMG5_pSol.
 *
 * Structure deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void MMG3D_Free_names(enum MMG5_arg starter,...);

/* library */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol (metric) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the remesh library.
 *
 */
int  MMG3D_mmg3dlib(MMG5_pMesh mesh, MMG5_pSol met );

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol (level-set) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the level-set discretization library.
 *
 */
int  MMG3D_mmg3dls(MMG5_pMesh mesh, MMG5_pSol met );

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol (output metric) structure.
 * \param disp pointer toward a sol (displacement for the lagrangian motion
 * mode) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the rigidbody movement library.
 *
 */
int  MMG3D_mmg3dmov(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol disp );

/** Tools for the library */
/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
void MMG3D_defaultValues(MMG5_pMesh mesh);

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 1.
 *
 * Store command line arguments.
 *
 */
int  MMG3D_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 1.
 *
 * Read local parameters file. This file must have the same name as
 * the mesh with the \a .mmg3d extension or must be named \a
 * DEFAULT.mmg3d.
 *
 */
int  MMG3D_parsop(MMG5_pMesh mesh,MMG5_pSol met);
/**
 * \param prog pointer toward the program name.
 *
 * Print help for mmg3d options.
 *
 */
void  MMG3D_usage(char *prog);
/**
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 * \return 1.
 *
 * Store the info structure in the mesh structure.
 *
 */
int  MMG3D_stockOptions(MMG5_pMesh mesh, MMG5_Info *info);
/**
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 *
 * Recover the info structure stored in the mesh structure.
 *
 */
void  MMG3D_destockOptions(MMG5_pMesh mesh, MMG5_Info *info);

/** Checks */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param critmin minimum quality for elements.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param eltab table of invalid elements.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage
 * (before the _MMG5_defsiz call), 1 for special storage (after this call).
 *
 * Search invalid elements (in term of quality or edge length).
 *
 */
int MMG3D_mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol met,double critmin,
                    double lmin, double lmax, int *eltab,char metRidTyp);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param critmin minimum quality for elements.
 * \param eltab pointer toward the table of invalid elements.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage
 * (before the _MMG5_defsiz call), 1 for special storage (after this call).
 *
 * Store elements which have worse quality than \a critmin in \a eltab,
 * \a eltab is allocated and could contain \a mesh->ne elements.
 *
 */
void  MMG3D_searchqua(MMG5_pMesh mesh, MMG5_pSol met, double critmin, int *eltab,
                     char metRidTyp);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param eltab table of invalid elements.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage
 * (before the _MMG5_defsiz call), 1 for special storage (after this call).
 *
 * \return 1 if success, 0 otherwise.
 *
 * Store in \a eltab elements which have edge lengths shorter than \a lmin
 * or longer than \a lmax, \a eltab is allocated and could contain \a mesh->ne
 * elements.
 *
 */
int  MMG3D_searchlen(MMG5_pMesh mesh, MMG5_pSol met, double lmin, double lmax,
                    int *eltab,char  metRidTyp);

/** Utils */
/**
 * \brief Return adjacent elements of a tetrahedron.
 * \param mesh pointer toward the mesh structure.
 * \param kel tetrahedron index.
 * \param listet pointer toward the table of the 4 tetra adjacent to \a kel.
 * (the index is 0 if there is no adjacent)
 * \return 1.
 *
 * Find the indices of the 4 adjacent elements of tetrahedron \a
 * kel. \f$listet[i] = 0\f$ if the \f$i^{th}\f$ face has no adjacent element
 * (so we are on a boundary face).
 *
 */
int MMG3D_Get_adjaTet(MMG5_pMesh mesh,int kel, int listet[4]);
/**
 * \param ca pointer toward the coordinates of the first edge's extremity.
 * \param cb pointer toward the coordinates of the second edge's extremity.
 * \param ma pointer toward the metric associated to the first edge's extremity.
 * \param mb pointer toward the metric associated to the second edge's extremity.
 * \return edge length.
 *
 * Compute length of edge \f$[ca,cb]\f$ (with \a ca and \a cb
 * coordinates of edge extremities) according to the size
 * prescription.
 *
 */
double (*MMG3D_lenedgCoor)(double *ca,double *cb,double *sa,double *sb);

/**
 * \param mesh pointer toward the mesh structure.
 * \param pack we pack the mesh at function begining if \f$pack=1\f$.
 * \return 0 if failed, 1 otherwise.
 *
 * Create table of adjacency. Set pack variable to 0 for a compact
 * mesh and to 1 for a mesh that need to be packed.
 *
 */
int  MMG3D_hashTetra(MMG5_pMesh mesh, int pack);

/**
 * \param mesh pointer toward the mesh structure
 * \param met pointer toward the sol structure
 * \return 1 if success
 *
 * Compute isotropic size map according to the mean of the length of the edges
 * passing through a point.
 *
 */
int MMG3D_DoSol(MMG5_pMesh mesh,MMG5_pSol met);

/** To associate function pointers without calling MMG3D_mmg3dlib */
/**
 * \param mesh pointer toward the mesh structure (unused).
 * \param met pointer toward the sol structure (unused).
 *
 * Set function pointers for caltet, lenedg, lenedgCoor defsiz, gradsiz...
 * depending if the readed metric is anisotropic or isotropic
 *
 */
void  MMG3D_setfunc(MMG5_pMesh mesh,MMG5_pSol met);

/** Old API °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°*/

/**
 * \enum MMG5_Param
 * \brief Input parameters for mmg library.
 *
 * Input parameters for mmg library. Options prefixed by \a
 * MMG5_IPARAM asked for integers values ans options prefixed by \a
 * MMG5_DPARAM asked for real values.
 *
 */
enum MMG5_Param {
  MMG5_IPARAM_verbose,           /*!< [-10..10], Tune level of verbosity */
  MMG5_IPARAM_mem,               /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  MMG5_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  MMG5_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  MMG5_IPARAM_iso,               /*!< [1/0], Level-set meshing */
  MMG5_IPARAM_lag,               /*!< [-1/0/1/2], Lagrangian option */
  MMG5_IPARAM_optim,             /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  MMG5_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
  MMG5_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  MMG5_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
  MMG5_IPARAM_nosurf,            /*!< [1/0], Avoid/allow surface modifications */
  MMG5_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
  MMG5_IPARAM_renum,             /*!< [1/0], Turn on/off point relocation with Scotch */
  MMG5_IPARAM_bucket,            /*!< [n], Specify the size of the bucket per dimension (DELAUNAY) */
  MMG5_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
  MMG5_DPARAM_hmin,              /*!< [val], Minimal mesh size */
  MMG5_DPARAM_hmax,              /*!< [val], Maximal mesh size */
  MMG5_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  MMG5_DPARAM_hgrad,             /*!< [val], Control gradation */
  MMG5_DPARAM_ls,                /*!< [val], Value of level-set (not use for now) */
  MMG5_PARAM_size,               /*!< [n], Number of parameters */
};

/*----------------------------- functions header -----------------------------*/
/* Initialization functions */
/* init structures */
/**
 * \param mesh adress of a pointer toward a pointer toward the mesh structure.
 * \param sol adress of a pointer toward a sol structure (metric or level-set).
 * \param disp adress of a pointer toward a sol structure
 * (displacement for the lagrangian mode).
 *
 * Allocate the mesh and solution structures and initialize it to
 * their default values.
 *
 */
void  MMG5_Init_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol, MMG5_pSol *disp );

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
int  MMG5_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, int np, int typSol);
/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param ne number of elements (tetrahedra).
 * \param nt number of triangles.
 * \param na number of edges.
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of vertices, tetrahedra, triangles and edges of the
 * mesh and allocate the associated tables. If call twice, reset the
 * whole mesh to realloc it at the new size
 *
 */
int  MMG5_Set_meshSize(MMG5_pMesh mesh, int np, int ne, int nt, int na);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 */
void  MMG5_Init_parameters(MMG5_pMesh mesh);

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
int  MMG5_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                     double c2, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of tetrahedron.
 * \param v1 second vertex of tetrahedron.
 * \param v2 third vertex of tetrahedron.
 * \param v3 fourth vertex of tetrahedron.
 * \param ref tetrahedron reference.
 * \param pos tetrahedron position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set tetrahedra of vertices \a v0, \a v1,\a v2,\a v3 and reference
 * \a ref at position \a pos in mesh structure.
 *
 */
int  MMG5_Set_tetrahedron(MMG5_pMesh mesh, int v0, int v1,
                          int v2, int v3, int ref, int pos);
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
int  MMG5_Set_triangle(MMG5_pMesh mesh, int v0, int v1,
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
int  MMG5_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set corner at point \a pos.
 *
 */
int  MMG5_Set_corner(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set point \a k as required.
 *
 */
int  MMG5_Set_requiredVertex(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k element index.
 * \return 1.
 *
 * Set element \a k as required.
 *
 */
int  MMG5_Set_requiredTetrahedron(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Set triangle \a k as required.
 *
 */
int  MMG5_Set_requiredTriangle(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set ridge at edge \a k.
 *
 */
int  MMG5_Set_ridge(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set edge \a k as required.
 *
 */
int  MMG5_Set_requiredEdge(MMG5_pMesh mesh, int k);
/**
 * \param met pointer toward the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure
 *
 */
int  MMG5_Set_scalarSol(MMG5_pSol met, double s,int pos);
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
int MMG5_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz, int pos);
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
int MMG5_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * To mark as ended a mesh given without using the API functions
 * (for example, mesh given by mesh->point[i] = 0 ...). Not recommanded.
 *
 */
void MMG5_Set_handGivenMesh(MMG5_pMesh mesh);

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
int MMG5_Chk_meshData(MMG5_pMesh mesh, MMG5_pSol met);

/** functions to set parameters */
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param iparam integer parameter to set (see \a MMG5_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 */
int  MMG5_Set_iparameter(MMG5_pMesh mesh,MMG5_pSol sol, int iparam, int val);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param dparam double parameter to set (see \a MMG5_Param structure).
 * \param val value of the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
 *
 */
int  MMG5_Set_dparameter(MMG5_pMesh mesh,MMG5_pSol sol, int dparam, double val);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typ type of entity (triangle, edge,...).
 * \param ref reference of the entity.
 * \param val value of the Hausdorff number.
 * \return 0 if failed, 1 otherwise.
 *
 * Set local parameters: set the hausdorff value at \a val for all
 * elements of type \a typ and reference \a ref.
 *
 */
int  MMG5_Set_localParameter(MMG5_pMesh mesh, MMG5_pSol sol, int typ, int ref, double val);

/** recover datas */
/**
 * \param mesh pointer toward the mesh structure.
 * \param np pointer toward the number of vertices.
 * \param ne pointer toward the number of elements (tetrahedra).
 * \param nt pointer toward the number of triangles.
 * \param na pointer toward the number of edges.
 * \return 1.
 *
 * Get the number of vertices, tetrahedra, triangles and edges of the mesh.
 *
 */
int  MMG5_Get_meshSize(MMG5_pMesh mesh, int* np, int* ne, int* nt, int* na);
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
int  MMG5_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, int* np,
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
int  MMG5_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, int* ref,
                     int* isCorner, int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 pointer toward the first vertex of tetrahedron.
 * \param v1 pointer toward the second vertex of tetrahedron.
 * \param v2 pointer toward the third vertex of tetrahedron.
 * \param v3 pointer toward the fourth vertex of tetrahedron.
 * \param ref pointer toward the tetrahedron reference.
 * \param isRequired pointer toward the flag saying if tetrahedron is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0, \a v1, \a v2, \a v3 and reference \a ref of
 * next tetra of mesh.
 *
 */
int  MMG5_Get_tetrahedron(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* v3,
                          int* ref, int* isRequired);
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
int  MMG5_Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref,
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
int  MMG5_Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired);
/**
 * \param met pointer toward the sol structure.
 * \param s pointer toward the scalar solution value.
 * \return 0 if failed, 1 otherwise.
 *
 * Get solution \a s of next vertex of mesh.
 *
 */
int  MMG5_Get_scalarSol(MMG5_pSol met, double* s);
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
int MMG5_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz);
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
int MMG5_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33);
/**
 * \param mesh pointer toward the mesh structure.
 * \param iparam integer parameter to set (see \a MMG5_Param structure).
 * \return The value of integer parameter.
 *
 * Get the value of integer parameter \a iparam.
 *
 */
int MMG5_Get_iparameter(MMG5_pMesh mesh, int iparam);

/* input/output functions */
/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh data.
 *
 */
int  MMG5_loadMesh(MMG5_pMesh mesh);
/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 */
int  MMG5_saveMesh(MMG5_pMesh mesh);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Load metric field.
 *
 */
int  MMG5_loadMet(MMG5_pMesh mesh,MMG5_pSol met);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric.
 *
 */
int  MMG5_saveMet(MMG5_pMesh mesh, MMG5_pSol met);

/* deallocations */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure (metric or solution).
 * \param disp pointer toward a sol structure (displacement for the
 * lagrangian mode).
 *
 * Deallocations before return.
 *
 */
void MMG5_Free_all(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol disp );

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure (metric or solution).
 * \param disp pointer toward a sol structure (displacement).
 *
 * Structure deallocations before return.
 *
 */
void MMG5_Free_structures(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol disp );

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (metric or solution).
 * \param disp pointer toward a sol structure (displacement).
 *
 * File name deallocations before return.
 *
 */
void MMG5_Free_names(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol disp);

/* library */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol (metric or level-set) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the remesh library.
 *
 */
int  MMG5_mmg3dlib(MMG5_pMesh mesh, MMG5_pSol met );


/* Library tools*/
/**
 * \brief Return adjacent elements of a tetrahedron.
 * \param mesh pointer toward the mesh structure.
 * \param kel tetrahedron index.
 * \param v0 pointer toward the index of the adjacent element of \a kel through
 * its face number 0.
 * \param v1 pointer toward the index of the adjacent element of \a kel through
 * its face number 1.
 * \param v2 pointer toward the index of the adjacent element of \a kel through
 * its face number 2.
 * \param v3 pointer toward the index of the adjacent element of \a kel through
 * its face number 3.
 * \return 1.
 *
 * Find the indices of the 4 adjacent elements of tetrahedron \a
 * kel. \f$v_i = 0\f$ if the \f$i^{th}\f$ face has no adjacent element
 * (so we are on a boundary face).
 *
 */
int MMG5_Get_adjaTet(MMG5_pMesh mesh, int kel, int *v0, int *v1, int *v2, int *v3);
/**
 * \param prog pointer toward the program name.
 *
 * Print help for mmgs options.
 *
 */
void MMG5_usage(char *prog);
/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
void MMG5_defaultValues(MMG5_pMesh mesh);

int MMG5_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met);

int MMG5_parsop(MMG5_pMesh mesh,MMG5_pSol met);
/**
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 * \return 1.
 *
 * Store the info structure in the mesh structure.
 *
 */
int MMG5_stockOptions(MMG5_pMesh mesh, MMG5_Info *info);

int MMG5_mmg3dcheck(MMG5_pMesh,MMG5_pSol,double, double,double, int *,char);

void MMG5_searchqua(MMG5_pMesh,MMG5_pSol,double, int *, char);

int MMG5_searchlen(MMG5_pMesh,MMG5_pSol, double, double, int *,char);

/**
 * \brief Return adjacent elements of a triangle.
 * \param mesh pointer toward the mesh structure.
 * \param ip vertex index.
 * \param vtab pointer toward an array of size MMG3D_LMAX that will contain
 * the indices of adjacent vertices to the vertex \a k.
 * \return 1 if success.
 *
 * Find the indices of the adjacent vertices of the vertex \a
 * ip.
 *
 */
int MMG3D_Get_adjaVertices(MMG5_pMesh mesh, int ip, int vtab[MMG3D_LMAX]);

#endif
