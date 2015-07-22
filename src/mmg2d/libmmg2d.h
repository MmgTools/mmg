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
 * \file mmg2d/libmmg2d.h
 * \brief C API for MMG2D library.
 * \author Cecile Dobrzynski and Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 */

#ifndef _MMG2DLIB_H
#define _MMG2DLIB_H

#include "mmg.h"
/**
 * \def MMG5_SUCCESS
 *
 * Return value for success.
 *
 */
#define MMG5_SUCCESS       0
/**
 * \def MMG5_LOWFAILURE
 *
 * Return value if the remesh process failed but we can save a conform
 * mesh.
 *
 */
#define MMG5_LOWFAILURE    1
/**
 * \def MMG5_STRONGFAILURE
 *
 * Return value if the remesh process failed and the mesh is
 * non-conform.
 *
 */
#define MMG5_STRONGFAILURE 2
/**
 * \enum MMG5_type
 * \brief Type of solutions.
 */
enum MMG5_type
{
  MMG5_Notype, /*!< Undefined type (unusable) */
  MMG5_Scalar, /*!< Scalar solution */
  MMG5_Vector, /*!< Vectorial solution */
  MMG5_Tensor  /*!< Tensorial solution */
};
/**
 * \enum MMG5_entities
 * \brief Type of mesh entities to which solutions are applied.
 */
enum MMG5_entities
{
  MMG5_Noentity, /*!< Undefined type (unusable) */
  MMG5_Vertex, /*!< Vertex entity */
  MMG5_Triangle, /*!< Triangle entity */
};
/**
 * \enum MMG5_Param
 * \brief Input parameters for mmg library.
 *
 * Input parameters for mmg library. Options prefixed by \a
 * MMG5_IPARAM asked for integers values ans options prefixed by \a
 * MMG5_DPARAM asked for real values.
 *
 */
enum MMG5_Param
{
  MMG5_IPARAM_verbose,           /*!< [-10..10], Tune level of verbosity */
  MMG5_IPARAM_mem,               /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  MMG5_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  MMG5_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  MMG5_IPARAM_iso,               /*!< [1/0], Level-set meshing */
  MMG5_IPARAM_lag,               /*!< [-1/0/1/2], Lagrangian option */
  MMG5_IPARAM_msh,               /*!< [0/1/2], Read/write to gmsh visu if val=1 (out) if val=2 (in/out) */
  MMG5_IPARAM_numsubdomain,       /*!<only if no given triangle, save the subdomain nb (0==all subdomain) */
  MMG5_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
  MMG5_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  MMG5_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
  MMG5_IPARAM_nosurf,            /*!< [1/0], Avoid/allow surface modifications */
  MMG5_IPARAM_bucket,            /*!< [n], Specify the size of the bucket per dimension (DELAUNAY) */
  MMG5_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
  MMG5_DPARAM_hmin,              /*!< [val], Minimal mesh size */
  MMG5_DPARAM_hmax,              /*!< [val], Maximal mesh size */
  MMG5_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  MMG5_DPARAM_hgrad,             /*!< [val], Control gradation */
  MMG5_DPARAM_ls,                /*!< [val], Value of level-set (not use for now) */
};

/*----------------------------- functions header -----------------------------*/
/* Initialization functions */
/* init structures */
/**
 * \param mesh pointer toward a pointer toward the mesh structure.
 * \param sol pointer toward a pointer toward the sol structure.
 *
 * Allocate the mesh and solution structures and initialize it to
 * their default values.
 *
 */
void  MMG2_Init_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol);
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
int MMG5_Set_iparameter(MMG5_pMesh mesh, MMG5_pSol sol, int iparam, int val);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param dparam double parameter to set (see \a MMG5_Param structure).
 * \val value of the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
 *
 */
int MMG5_Set_dparameter(MMG5_pMesh mesh, MMG5_pSol sol, int dparam, double val);
/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param nt number of triangles.
 * \param na number of edges.
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of vertices, tetrahedra, triangles and edges of the
 * mesh and allocate the associated tables. If call twice, reset the
 * whole mesh to realloc it at the new size
 *
 */
int  MMG5_Set_meshSize(MMG5_pMesh mesh, int np, int nt, int na);
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
int MMG5_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity,
                     int np, int typSol);
/* init structure datas */
/**
 * \param mesh pointer toward the mesh structure.
 * \param c0 coordinate of the point along the first dimension.
 * \param c1 coordinate of the point along the second dimension.
 * \param ref point reference.
 * \param pos position of the point in the mesh.
 * \return 1.
 *
 * Set vertex of coordinates \a c0, \a c1 and reference \a ref
 * at position \a pos in mesh structure
 *
 */
int  MMG5_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                     int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param num integer
 * \param c0 pointer toward the coordinate of the point along the first dimension.
 * \param c1 pointer toward the coordinate of the point along the second dimension.
 * \param ref pointer to the point reference.
 * \param isCorner pointer toward the flag saying if point is corner.
 * \param isRequired pointer toward the flag saying if point is required.
 * \return 1.
 *
 * Get coordinates \a c0, \a c1 and reference \a ref of 
 * vertex num of mesh.
 *
 */
int  MMG5_Get_vertex(MMG5_pMesh mesh, int num,double* c0, double* c1, int* ref,
                     int* isCorner, int* isRequired);
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
int MMG5_Set_triangle(MMG5_pMesh mesh, int v0, int v1, 
                      int v2, int ref, int pos);

/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of edge.
 * \param v1 second vertex of edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set edge of vertices \a v0, \a v1 and reference \a ref
 * at position \a pos in mesh structure.
 *
 */
int MMG5_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref, int pos);

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
int  MMG5_Get_meshSize(MMG5_pMesh mesh, int* np, int* nt, int* na);
/**
 * \param met pointer toward the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure
 *
 */
int MMG5_Set_scalarSol(MMG5_pSol met, double s, int pos);

/**
 * \param met pointer toward the sol structure.
 * \param s solution symetric tensor value (s11 s12 s22)
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensor value \a s at position \a pos in solution structure
 *
 */
int MMG5_Set_tensorSol(MMG5_pSol met, double* s, int pos);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Check if the number of given entities match with mesh and sol size
 * (not mandatory) and check mesh datas.
 *
 */
int MMG5_Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * Deallocations before return.
 *
 */
void MMG5_Free_all(MMG5_pMesh mesh,MMG5_pSol met
  );
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * Structure deallocations before return.
 *
 */
void MMG5_Free_structures(MMG5_pMesh mesh,MMG5_pSol met
  );

int MMG2_loadMesh(MMG5_pMesh ,char *);
int MMG2_loadSol(MMG5_pSol ,char *,int,int);
int MMG2_loadVect(MMG5_pMesh ,char *);
int MMG2_saveMesh(MMG5_pMesh ,char *);
int MMG2_saveSol(MMG5_pMesh ,MMG5_pSol ,char *);
int MMG2_saveVect(MMG5_pMesh mesh,MMG5_pSol sol,char *filename,double lambda);

int MMG2_mmg2dlib(MMG5_pMesh mesh,MMG5_pSol sol,void (*titi)(int ,int,int,int,int));
void (*MMG2_callbackinsert) (int ,int ,int ,int, int);

#endif
