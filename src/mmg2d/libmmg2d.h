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
 * \brief API headers for the mmg2d library
 * \author Cecile Dobrzynski and Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning To keep the genheader working, don't break line between the enum
 * name and the opening brace (it creates errors under windows)
 */

#ifndef _MMG2DLIB_H
#define _MMG2DLIB_H

#include "mmgcommon.h"

/**
 * Maximum array size when storing adjacent points (or ball) of a vertex.
 */
#define MMG2D_LMAX   1024

/**
 * \enum MMG2D_Param
 * \brief Input parameters for mmg library.
 *
 * Input parameters for mmg library. Options prefixed by \a
 * MMG2D_IPARAM asked for integers values ans options prefixed by \a
 * MMG2D_DPARAM asked for real values.
 *
 */
enum MMG2D_Param {
  MMG2D_IPARAM_verbose,           /*!< [-10..10], Tune level of verbosity */
  MMG2D_IPARAM_mem,               /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  MMG2D_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  MMG2D_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  MMG2D_IPARAM_iso,               /*!< [1/0], Level-set meshing */
  MMG2D_IPARAM_lag,               /*!< [-1/0/1/2], Lagrangian option */
  MMG2D_IPARAM_msh,               /*!< [0/1/2], Read/write to gmsh visu if val=1 (out) if val=2 (in/out) */
  MMG2D_IPARAM_numsubdomain,       /*!<only if no given triangle, save the subdomain nb (0==all subdomain) */
  MMG2D_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
  MMG2D_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  MMG2D_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
  MMG2D_IPARAM_nosurf,            /*!< [1/0], Avoid/allow surface modifications */
  MMG2D_IPARAM_bucket,            /*!< [n], Specify the size of the bucket per dimension (DELAUNAY) */
  MMG2D_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
  MMG2D_DPARAM_hmin,              /*!< [val], Minimal mesh size */
  MMG2D_DPARAM_hmax,              /*!< [val], Maximal mesh size */
  MMG2D_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  MMG2D_DPARAM_hgrad,             /*!< [val], Control gradation */
  MMG2D_DPARAM_ls,                /*!< [val], Value of level-set (not use for now) */
};

/*----------------------------- functions header -----------------------------*/
/* Initialization functions */
/* init structures */

/**
 * \param starter dummy argument used to initialize the variadic argument list
 * \param ... variadic arguments. For now, you need to call the \a
 * MMG2D_Init_mesh function with the following arguments :
 * MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh,
 * MMG5_ARG_ppMet, your_metric,MMG5_ARG_end). Here, \a your_mesh is a pointer
 * toward \a MMG5_pMesh and \a your_metric a pointer toward \a MMG5_pSol.
 *
 * MMG structures allocation and initialization.
 *
 */
void MMG2D_Init_mesh(enum MMG5_arg starter,...);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Initialize file names to their default values.
 *
 */
void  MMG2D_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 */
void  MMG2D_Init_parameters(MMG5_pMesh mesh);

/* init file names */
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * Set the name of input mesh.
 *
 */
int  MMG2D_Set_inputMeshName(MMG5_pMesh mesh, char* meshin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * Set the name of output mesh file.
 *
 */
int  MMG2D_Set_outputMeshName(MMG5_pMesh mesh, char* meshout);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * Set the name of input solution file.
 *
 */
int  MMG2D_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solout name of the output solution file.
 * \return 0 if failed, 1 otherwise.
 *
 *  Set the name of output solution file.
 *
 */
int  MMG2D_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solout);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param iparam integer parameter to set (see \a MMG2D_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 */
int MMG2D_Set_iparameter(MMG5_pMesh mesh, MMG5_pSol sol, int iparam, int val);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param dparam double parameter to set (see \a MMG2D_Param structure).
 * \param val value of the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
 *
 */
int MMG2D_Set_dparameter(MMG5_pMesh mesh, MMG5_pSol sol, int dparam, double val);
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
int  MMG2D_Set_meshSize(MMG5_pMesh mesh, int np, int nt, int na);
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
int MMG2D_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity,
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
int  MMG2D_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
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
int  MMG2D_Get_vertex(MMG5_pMesh mesh, int num,double* c0, double* c1, int* ref,
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
int MMG2D_Set_triangle(MMG5_pMesh mesh, int v0, int v1,
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
int MMG2D_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref, int pos);

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
int  MMG2D_Get_meshSize(MMG5_pMesh mesh, int* np, int* nt, int* na);
/**
 * \param met pointer toward the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure
 *
 */
int MMG2D_Set_scalarSol(MMG5_pSol met, double s, int pos);

/**
 * \param met pointer toward the sol structure.
 * \param s solution symetric tensor value (s11 s12 s22)
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensor value \a s at position \a pos in solution structure
 *
 */
int MMG2D_Set_tensorSol(MMG5_pSol met, double* s, int pos);

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Check if the number of given entities match with mesh and sol size
 * (not mandatory) and check mesh datas.
 *
 */
int MMG2D_Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met);

/* deallocations */
/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments. For now, you need to call the \a
 * MMG2D_Free_all function with the following arguments :
 * MMG2D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh,
 * MMG5_ARG_ppMet, your_metric,MMG5_ARG_end). Here, \a your_mesh is a pointer
 * toward \a MMG5_pMesh and \a your_metric a pointer toward \a MMG5_pSol.
 *
 * Deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void MMG2D_Free_all(enum MMG5_arg starter,...);

/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments. For now, you need to call the \a
 * MMG2D_Free_structures function with the following arguments :
 * MMG2D_Free_structures(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh,
 * MMG5_ARG_ppMet, your_metric,MMG5_ARG_end). Here, \a your_mesh is a pointer
 * toward \a MMG5_pMesh and \a your_metric a pointer toward \a MMG5_pSol.
 *
 * Structure deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void MMG2D_Free_structures(enum MMG5_arg starter,...);

/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments. For now, you need to call the \a
 * MMG2D_Free_names function with the following arguments :
 * MMG2D_Free_names(MMG5_ARG_start,MMG5_ARG_ppMesh, your_mesh,
 * MMG5_ARG_ppMet, your_metric,MMG5_ARG_end). Here, \a your_mesh is a pointer
 * toward \a MMG5_pMesh and \a your_metric a pointer toward \a MMG5_pSol.
 *
 * Structure deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void MMG2D_Free_names(enum MMG5_arg starter,...);

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of the readed file.
 * \return 0 or -1 if fail, 1 otherwise
 *
 * Read mesh data.
 *
 */
int MMG2D_loadMesh(MMG5_pMesh mesh,char * filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure..
 * \param filename name of the solution file.
 * \param msh if 1, read the 2D solution saved in 3D (.sol files saved by gmsh)
 * \return 0 or -1 if fail, 1 otherwise.
 *
 * Load solution field.
 *
 */
int MMG2D_loadSol(MMG5_pMesh mesh,MMG5_pSol sol,char * filename,int msh);

int MMG2D_loadVect(MMG5_pMesh ,char *);

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of the readed file.
 * \return 0 or -1 if fail, 1 otherwise.
 *
 * Save mesh data.
 *
 */
int MMG2D_saveMesh(MMG5_pMesh ,char *);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure..
 * \param filename name of the solution file.
 * \return 0 or -1 if fail, 1 otherwise.
 *
 * Save metric field.
 *
 */
int MMG2D_saveSol(MMG5_pMesh  mesh,MMG5_pSol sol ,char *filename);
int MMG2D_saveVect(MMG5_pMesh mesh,MMG5_pSol sol,char *filename,double lambda);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward a sol structure (metric).
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 *
 * Main program for the mesh adaptation library .
 *
 */
int MMG2D_mmg2dlib(MMG5_pMesh mesh,MMG5_pSol sol);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward a sol structure (metric).
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 *
 * Main program for the mesh generation library .
 *
 */
int MMG2D_mmg2dmesh(MMG5_pMesh mesh,MMG5_pSol sol);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward a sol structure (metric).
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 *
 * Main program for the level-set discretization library .
 *
 */
int MMG2D_mmg2dls(MMG5_pMesh mesh,MMG5_pSol sol) ;
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward a sol structure (displacement).
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 *
 * Main program for the rigid body movement library .
 *
 */
int MMG2D_mmg2dmov(MMG5_pMesh mesh,MMG5_pSol sol);

/* Tools for the library */
void (*MMG2D_callbackinsert) (int ,int ,int ,int, int);

/**
 * \param type metric type
 *
 * Set function pointers for length, caltri, buckin... depending if case is iso or aniso
 *
 */
void MMG2D_setfunc(MMG5_pMesh mesh,MMG5_pSol met);

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
int MMG2D_Get_adjaTri(MMG5_pMesh mesh, int kel, int listri[3]);

/**
 * \brief Return adjacent elements of a triangle.
 * \param mesh pointer toward the mesh structure.
 * \param ip vertex index.
 * \param lispoi pointer toward an array of size MMG2D_LMAX that will contain
 * the indices of adjacent vertices to the vertex \a ip.
 * \return nbpoi the number of adjacent points if success, 0 if fail.
 *
 * Find the indices of the adjacent vertices of the vertex \a
 * ip.
 *
 */
extern
int MMG2D_Get_adjaVertices(MMG5_pMesh mesh, int ip, int lispoi[MMG2D_LMAX]);

/**
 * \brief Return adjacent elements of a triangle.
 * \param mesh pointer toward the mesh structure.
 * \param ip vertex index.
 * \param start index of a triangle holding \a ip.
 * \param lispoi pointer toward an array of size MMG2D_LMAX that will contain
 * the indices of adjacent vertices to the vertex \a ip.
 * \return nbpoi the number of adjacent points if success, 0 if fail.
 *
 * Find the indices of the adjacent vertices of the vertex \a
 * ip of the triangle \a start.
 *
 */
int MMG2D_Get_adjaVerticesFast(MMG5_pMesh mesh, int ip,int start,
                               int lispoi[MMG2D_LMAX]);

/**
 * \param mesh pointer toward the mesh structure
 *
 * Free the mesh elements (and the adjacency).
 *
 */
void MMG2D_Free_Triangles(MMG5_pMesh mesh);

/**
 * \param mesh pointer toward the mesh structure
 *
 * Free the mesh edges (and the associated xpoints).
 *
 */
void MMG2D_Free_Edges(MMG5_pMesh mesh);

/**
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward the solution structure
 *
 * Free the solution.
 *
 */
void MMG2D_Free_Solutions(MMG5_pMesh mesh,MMG5_pSol sol);

#endif
