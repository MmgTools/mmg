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
 * \file common/libmmg.h
 * \brief C API for MMG libraries.
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 */

#ifndef _MMGLIB_H
#define _MMGLIB_H

#include "chrono.h"

/**
 * \struct MMG5_Par
 * number) associated to a specific reference.
 *
 * Store the local values (minimal and maximal sizes and Hausdorff number)
 * associated to the given reference of an element of type \a elt (point,
 * edge... ).
 *
 */
typedef struct {
  double   hmin; /*!< minimal size for edges */
  double   hmax; /*!< maximal size for edges */
  double   hausd; /*!< Hausdorff value */
  int      ref; /*!< Reference value */
  char     elt; /*!< Element type */
} MMG5_Par;
typedef MMG5_Par * MMG5_pPar;

/**
 * \struct MMG5_Point
 * \brief Structure to store points of a MMG mesh.
 * \todo What to do with n[3], try to remove ig,s.
 */
typedef struct {
  double   c[3]; /*!< Coordinates of point */
  double   n[3]; /*!< Tangeant for mmgs */
  int      ref; /*!< Reference of point */
  int      xp; /*!< Surface point number */
  int      tmp; /*!< Index of point in the saved mesh (we don't count
                  the unused points)*/
  int      flag; /*!< Flag to know if we have already treated the point */
  char     tag; /*!< Contains binary flags : if \f$tag=23=16+4+2+1\f$, then
                  the point is \a MG_REF, \a MG_GEO, \a MG_REQ and \a MG_BDY */
  char     tagdel; /*!< Tag for delaunay */
  int      ig,s;
} MMG5_Point;
typedef MMG5_Point * MMG5_pPoint;

/**
 * \struct MMG5_xPoint
 * \brief Structure to store surface points of a MMG mesh.
 */
typedef struct {
  double   n1[3],n2[3]; /*!< Normals at boundary vertex;
                          n1!=n2 if the vertex belong to a ridge */
  double   t[3]; /*!< Tangeant at vertex */
} MMG5_xPoint;
typedef MMG5_xPoint * MMG5_pxPoint;

/**
 * \struct MMG5_Edge
 * \brief Structure to store edges of a MMG mesh.
 */
typedef struct {
  int      a,b; /*!< Extremities of the edge */
  int      ref; /*!< Reference of the edge */
  char     tag; /*!< Binary flags */
} MMG5_Edge;
typedef MMG5_Edge * MMG5_pEdge;

/**
 * \struct MMG5_Tria
 * \brief Structure to store triangles of a MMG mesh.
 * \todo try to remove cc.
 */
typedef struct {
  int      v[3]; /*!< Vertices of the triangle */
  int      ref; /*!< Reference of the triangle */
  int      base;
  int      cc;
  int      edg[3]; /*!< edg[i] contains the ref of the \f$i^{th}\f$ edge
                     of triangle */
  int      flag;
  char     tag[3]; /*!< tag[i] contains the tag associated to the
                     \f$i^{th}\f$ edge of triangle */
} MMG5_Tria;
typedef MMG5_Tria * MMG5_pTria;

/**
 * \struct MMG5_Tetra
 * \brief Structure to store tetrahedra of a MMG mesh.
 */
typedef struct {
  int      v[4]; /*!< Vertices of the tetrahedron */
  int      ref; /*!< Reference of the tetrahedron */
  int      base;
  int      mark; /*!< Used for delaunay */
  int      xt; /*!< Index of the surface \ref MMG5_xTetra associated to
                 the tetrahedron*/
  int      flag;
  char     tag;
  double   qual; /*!< Quality of the element */
} MMG5_Tetra;
typedef MMG5_Tetra * MMG5_pTetra;


/**
 * \struct MMG5_xTetra
 * \brief Structure to store the surface tetrahedra of a MMG mesh.
 */
typedef struct {
  int      ref[4]; /*!< ref[i] is the reference of the opposite triangle to the
                     \f$i^{th}\f$ vertex of the tetrahedron;*/
  int      edg[6]; /*!< edg[i] contains the reference of the
                     \f$i^{th}\f$ edge of the tetrahedron */
  char     ftag[4]; /*!< ftag[i] contains the tag associated to the
                      \f$i^{th}\f$ face of the tetrahedron */
  char     tag[6]; /*!< tag[i] contains the tag associated to the
                     \f$i^{th}\f$ edge of the tetrahedron */
  char     ori; /*!< Orientation of the triangles of the tetrahedron:
                  the $\f$i^{th}\f$ bit of ori is set to 0 when the
                  \f$i^{th}\f$ face is bad orientated */
} MMG5_xTetra;
typedef MMG5_xTetra * MMG5_pxTetra;

/**
 * \struct MMG5_Info
 * \brief Store input parameters of the run.
 */
typedef struct {
  double        dhd,hmin,hmax,hgrad,hausd,min[3],max[3],delta,ls;
  int           mem,npar,npari;
  char          nreg;
  int           renum;
  char          imprim,ddebug,badkal,iso,fem;
  unsigned char noinsert, noswap, nomove;
  int           bucket;
  MMG5_pPar     par;
} MMG5_Info;

/**
 * \struct MMG5_hgeom
 * \brief To store geometric edges.
 */
typedef struct {
  int   a; /*!< First extremity of edge */
  int   b;  /*!< Second extremity of edge */
  int   ref; /*!< Reference of edge */
  int   nxt; /*!< Next element of hash table */
  char  tag; /*!< tag of edge */
} MMG5_hgeom;

typedef struct {
  int         siz,max,nxt;
  MMG5_hgeom  *geom;
} MMG5_HGeom;

/**
 * \struct MMG5_Mesh
 * \brief MMG mesh structure.
 * \todo try to remove nc1;
 */
typedef struct {
  int       ver; /*!< Version of the mesh file */
  int       dim; /*!< Dimension of the mesh */
  int       type; /*!< Type of the mesh */
  long long memMax; /*!< Maximum memory available */
  long long memCur; /*!< Current memory used */
  double    gap; /*!< Gap for table reallocation */
  int       npi,nti,nai,nei,np,na,nt,ne,npmax,namax,ntmax,nemax,xpmax,xtmax;
  int       nc1;

  int       base; /*!< Used with \a flag to know if an entity has been
                    treated */
  int       mark; /*!< Flag for delaunay (to know if an entity has
                    been treated) */
  int       xp,xt; /*!< Number of surfaces points/triangles */
  int       npnil; /*!< Index of first unused point */
  int       nenil; /*!< Index of first unused element */
  int      *adja; /*!< Table of tetrahedron adjacency: if
                    \f$adjt[4*i+1+j]=4*k+l\f$ then the \f$i^{th}\f$ and
                    \f$k^th\f$ tetrahedra are adjacent and share their
                    faces \a j and \a l (resp.) */
  int      *adjt; /*!< Table of triangles adjacency: if
                    \f$adjt[3*i+1+j]=3*k+l\f$ then the \f$i^{th}\f$ and
                    \f$k^th\f$ triangles are adjacent and share their
                    edges \a j and \a l (resp.) */
  char     *namein; /*!< Input mesh name */
  char     *nameout; /*!< Output mesh name */

  MMG5_pPoint    point; /*!< Pointer toward the \ref MMG5_Point structure */
  MMG5_pxPoint   xpoint; /*!< Pointer toward the \ref MMG5_xPoint structure */
  MMG5_pTetra    tetra; /*!< Pointer toward the \ref MMG5_Tetra structure */
  MMG5_pxTetra   xtetra; /*!< Pointer toward the \ref MMG5_xTetra structure */
  MMG5_pTria     tria; /*!< Pointer toward the \ref MMG5_Tria structure */
  MMG5_pEdge     edge; /*!< Pointer toward the \ref MMG5_Edge structure */
  MMG5_HGeom     htab; /*!< \ref MMG5_HGeom structure */
  MMG5_Info      info; /*!< \ref MMG5_Info structure */
} MMG5_Mesh;
typedef MMG5_Mesh  * MMG5_pMesh;

/**
 * \struct MMG5_sol
 * \brief MMG Solution structure (for solution or metric).
 */
typedef struct {
  int       ver; /* Version of the solution file */
  int       dim; /* Dimension of the solution file*/
  int       np; /* Number of points of the solution */
  int       npmax; /* Maximum number of points */
  int       npi; /* Temporary number of points (internal use only) */
  int       size; /* Number of solutions per entity */
  int       type; /* Type of the solution (scalar, vectorial of tensorial) */
  double   *m; /*!< Solution values */
  char     *namein; /*!< Input solution file name */
  char     *nameout; /*!< Output solution file name */
} MMG5_Sol;
typedef MMG5_Sol * MMG5_pSol;

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
void  MMG5_Init_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Initialize file names to their default values.
 *
 */
void  MMG5_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 */
void  (*MMG5_Init_parameters)(MMG5_pMesh mesh);

/* init file names */
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * Set the name of input mesh.
 *
 */
int  MMG5_Set_inputMeshName(MMG5_pMesh mesh, char* meshin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * Set the name of output mesh file.
 *
 */
int  MMG5_Set_outputMeshName(MMG5_pMesh mesh, char* meshout);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * Set the name of input solution file.
 *
 */
int  MMG5_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solout name of the output solution file.
 * \return 0 if failed, 1 otherwise.
 *
 *  Set the name of output solution file.
 *
 */
int  MMG5_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solout);

/* deallocations */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * File name deallocations before return.
 *
 */
void MMG5_Free_names(MMG5_pMesh mesh, MMG5_pSol met);

#endif
