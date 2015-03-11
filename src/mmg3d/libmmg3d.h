/* =============================================================================
**  This file is part of the Mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
**
**  Mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with Mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the Mmg distribution only if you accept them.
** =============================================================================
*/

 /**
 * \file mmg3d/libmmg3d.h
 * \brief C API for MMG3D library.
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 */

#ifndef _MMG3DLIB_H
#define _MMG3DLIB_H

#include "chrono.h"

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
 * \def SIZE
 *
 * Size of the mesh of singularities inside the main mesh (SINGUL mode only).
 *
 */
#define MMG5_SIZE 0.75

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
 * \brief Input parameters for Mmg library.
 *
 * Input parameters for Mmg library. Options prefixed by \a
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
    MMG5_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
    MMG5_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
    MMG5_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
    MMG5_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
    MMG5_IPARAM_renum,             /*!< [1/0], Turn on/off point relocation with Scotch */
    MMG5_IPARAM_sing,              /*!< [1/0], Turn on/off the insertion of singularities */
    MMG5_IPARAM_bucket,            /*!< [n], Specify the size of the bucket per dimension (DELAUNAY) */
    MMG5_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
    MMG5_DPARAM_hmin,              /*!< [val], Minimal mesh size */
    MMG5_DPARAM_hmax,              /*!< [val], Maximal mesh size */
    MMG5_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
    MMG5_DPARAM_hgrad,             /*!< [val], Control gradation */
    MMG5_DPARAM_ls,                /*!< [val], Value of level-set (not use for now) */
    MMG5_PARAM_size,               /*!< [n], Number of parameters */
  };

/**
 * \struct MMG5_Par
 * \brief Local Hausdorff number associated to a specific reference.
 *
 * Store the local Hausdorff number associated to the given reference
 * of an element of type \a elt (point, edge... ).
 *
 */
typedef struct {
  double   hausd; /*!< Hausdorff value */
  int      ref; /*!< Reference value */
  char     elt; /*!< Element type */
} MMG5_Par;
typedef MMG5_Par * MMG5_pPar;

/**
 * \struct MMG5_Point
 * \brief Structure to store points of a MMG3D mesh.
 */
typedef struct {
  double   c[3]; /*!< Coordinates of point */
  int      ref; /*!< Reference of point */
  int      xp; /*!< Surface point number */
  int      tmp; /*!< Index of point in the saved mesh (we don't count
                   the unused points)*/
  int      flag; /*!< Flag to know if we have already treated the point */
  char     tag; /*!< Contains binary flags : if \f$tag=23=16+4+2+1\f$, then
                   the point is \a MG_REF, \a MG_GEO, \a MG_REQ and \a MG_BDY */
  char     tagdel; /*!< Tag for delaunay */
} MMG5_Point;
typedef MMG5_Point * MMG5_pPoint;

/**
 * \struct MMG5_xPoint
 * \brief Structure to store surface points of a MMG3D mesh.
 */
typedef struct {
  double   n1[3],n2[3]; /*!< Normals at boundary vertex;
                           n1!=n2 if the vertex belong to a ridge */
  double   t[3]; /*!< Tangeant at vertex */
} MMG5_xPoint;
typedef MMG5_xPoint * MMG5_pxPoint;

/**
 * \struct MMG5_Edge
 * \brief Structure to store edges of a MMG3D mesh.
 */
typedef struct {
  int      a,b; /*!< Extremities of the edge */
  int      ref; /*!< Reference of the edge */
  char     tag; /*!< Binary flags */
} MMG5_Edge;
typedef MMG5_Edge * MMG5_pEdge;

/**
 * \struct MMG5_Tria
 * \brief Structure to store triangles of a MMG3D mesh.
 */
typedef struct {
  int      v[3]; /*!< Vertices of the triangle */
  int      ref; /*!< Reference of the triangle */
  int      base;
  int      edg[3]; /*!< edg[i] contains the ref of the \f$i^{th}\f$ edge
                      of triangle */
  int      flag;
  char     tag[3]; /*!< tag[i] contains the tag associated to the
                      \f$i^{th}\f$ edge of triangle */
} MMG5_Tria;
typedef MMG5_Tria * MMG5_pTria;

/**
 * \struct MMG5_Tetra
 * \brief Structure to store tetrahedra of a MMG3D mesh.
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
 * \brief Structure to store the surface tetrahedra of a MMG3D mesh.
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
 * \struct MMG5_hgeom
 * \brief To store geometric edges.
 */
typedef struct {
  int   a,b,ref,nxt;
  char  tag;
} MMG5_hgeom;

typedef struct {
  int         siz,max,nxt;
  MMG5_hgeom  *geom;
} MMG5_HGeom;

/**
 * \struct MMG5_Info
 * \brief Store input parameters of the run.
 */
typedef struct {
  double        dhd,hmin,hmax,hgrad,hausd,min[3],max[3],delta,ls;
  int           mem,sing,npar,npari;
  int           renum;
  char          imprim,ddebug,badkal,iso,fem;
  unsigned char noinsert, noswap, nomove;
  int           bucket;
  MMG5_pPar     par;
} MMG5_Info;

/**
 * \struct MMG5_Mesh
 * \brief MMG3D mesh structure.
 */
typedef struct {
  int       ver; /*!< Version of the mesh file */
  int       dim; /*!< Dimension of the mesh */
  int       type; /*!< Type of the mesh */
  long long memMax; /*!< Maximum memory available */
  long long memCur; /*!< Current memory used */
  double    gap; /*!< Gap for table reallocation */
  int       npi,nti,nai,nei,np,na,nt,ne,npmax,namax,ntmax,nemax,xpmax,xtmax;
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
 * \brief MMG3D Solution structure (for solution or metric).
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

/**
 * \struct MMG5_sPoint
 * \brief Structure to store MMG3D singular point.
 * (only for singularities insertion).
 *
 * Structure to store MMG3D singular points (only used for singularities
 * insertion: \a SINGUL precompilator flag).
 *
 */
typedef struct {
  double         c[3]; /*!< Coordinates of the point */
  double         n[3]; /*!< Normal to the point */
  int            flag,tmp,tet;
  unsigned char  tag;
} MMG5_sPoint;
typedef MMG5_sPoint * MMG5_psPoint;

/**
 * \struct MMG5_Singul
 * \brief Structure to store the singularities of a mesh.
 * (only for singularities insertion).
 */
typedef struct {
  char     *namein; /*!< Name of the mesh containing the singularities */
  double   min[3]; /*!< Minimum coordinates for rescaling */
  double   max[3]; /*!< Maximum coordinates for rescaling */
  int      nsi;
  int      ns; /*!< Number of singular vertices */
  int      na; /*!< Number of singular edges */
  MMG5_psPoint  point; /*!< Pointer toward \ref MMG5_sPoint structure */
  MMG5_pEdge    edge; /*!< Pointer toward \ref MMG5_Edge structure */
} MMG5_Singul;
typedef MMG5_Singul * MMG5_pSingul;


/*----------------------------- functions header -----------------------------*/
/** Initialization functions */
/* init structures */
#ifndef SINGUL
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
 * singularities mode).
 *
 * Initialize file names to their default values.
 *
 */
void  MMG5_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);
#else
/**
 * \param mesh pointer toward a pointer toward the mesh structure.
 * \param sol pointer toward a pointer toward the sol structure.
 * \param sing pointer toward a pointer toward the sing structure
 * (only for insertion of singularities mode).
 *
 * Allocate the mesh and solution structures and initialize it to
 * their default values.
 *
 */
void  MMG5_Init_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol, MMG5_pSingul *sing);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param sing pointer toward the sing structure (only for insertion of
 * singularities mode).
 *
 * Initialize file names to their default values.
 *
 */
void  MMG5_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSingul sing);
#endif
/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 */
void  MMG5_Init_parameters(MMG5_pMesh mesh);

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
 * \param solout name of the output solution file.
 * \return 0 if failed, 1 otherwise.
 *
 *  Set the name of output solution file.
 *
 */
int  MMG5_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solout);
#ifdef SINGUL
/**
 * \param mesh pointer toward the mesh structure.
 * \param sing pointer toward the sing structure.
 * \param singin name for the input singularies file.
 * \return 1.
 *
 * Set the name of input singularities file (only for insertion of
 * singularities mode).
 *
 */
int  MMG5_Set_inputSingulName(MMG5_pMesh mesh,MMG5_pSingul sing, char* singin);
#endif

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
#ifdef SINGUL
/**
 * \param mesh pointer toward the mesh structure.
 * \param sing pointer toward the sing structure.
 * \param np number of singular vertices.
 * \param na number of singular edges.
 * \return 1.
 *
 * Set the number of singular vertices and edges and allocate the
 * associated tables (only for insertion of singularities mode: \a
 * SINGUL preprocessor flag).
 *
 */
int  MMG5_Set_singulSize(MMG5_pMesh mesh,MMG5_pSingul sing, int np, int na);
#endif

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
 * \param pos vertex index.
 * \return 1.
 *
 * Set corner at point \a pos.
 *
 */
int  MMG5_Set_corner(MMG5_pMesh mesh, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param pos vertex index.
 * \return 1.
 *
 * Set point \a pos as required.
 *
 */
int  MMG5_Set_requiredVertex(MMG5_pMesh mesh, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param pos element index.
 * \return 1.
 *
 * Set element \a pos as required.
 *
 */
int  MMG5_Set_requiredTetrahedron(MMG5_pMesh mesh, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param pos triangle index.
 * \return 1.
 *
 * Set triangle \a pos as required.
 *
 */
int  MMG5_Set_requiredTriangle(MMG5_pMesh mesh, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param pos edge index.
 * \return 1.
 *
 * Set ridge at edge \a pos.
 *
 */
int  MMG5_Set_ridge(MMG5_pMesh mesh, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param pos edge index.
 * \return 1.
 *
 * Set edge \a k as required.
 *
 */
int  MMG5_Set_requiredEdge(MMG5_pMesh mesh, int pos);
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
#ifdef SINGUL
/**
 * \param sing pointer toward the sing structure.
 * \param c0 coordinate of the point along the first dimension.
 * \param c1 coordinate of the point along the second dimension.
 * \param c2 coordinate of the point along the third dimension.
 * \param typ unused parameter.
 * \param pos position of the point in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set singular point of coordinates \a c0, \a c1, \a c2 at position
 * \a pos in the singularities structure (only for insertion of
 * singularities mode: SINGUL preprocessor flag).
 *
 */
int  MMG5_Set_singulVertex(MMG5_pSingul sing, double c0,
                           double c1, double c2, int typ,int pos);
/**
 * \param sing pointer toward the sing structure.
 * \param v0 first extremity of the edge.
 * \param v1 second extremity of the edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set singular edge of extremities \a v0,\a v1 and reference \a ref at
 * position \a pos in the singularities structure (only for insertion of
 * singularities mode: \a SINGUL preprocessor flag).
 *
 */
int  MMG5_Set_singulEdge(MMG5_pSingul sing, int v0, int v1, int ref,int pos);
/**
 * \param sing pointer toward the sing structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set corner at singular vertex \a k (only for insertion of
 * singularities mode: \a SINGUL preprocessor flag).
 *
 */
int  MMG5_Set_singulCorner(MMG5_pSingul sing, int pos);
/**
 * \param sing pointer toward the sing structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set required vertex at singular vertex \a k (only for insertion of
 * singularities mode: \a SINGUL preprocessor flag).
 *
 */
int  MMG5_Set_singulRequiredVertex(MMG5_pSingul sing, int pos);
/**
 * \param sing pointer toward the sing structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set ridge at singular edge \a k (only for insertion of
 * singularities mode: \a SINGUL preprocessor flag).
 *
 */
int  MMG5_Set_singulRidge(MMG5_pSingul sing, int pos);
/**
 * \param sing pointer toward the sing structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set required edge at singular edge \a k (only for insertion of
 * singularities mode: \a SINGUL preprocessor flag).
 *
 */
int  MMG5_Set_singulRequiredEdge(MMG5_pSingul sing, int pos);
#endif
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
 * \val value of the parameter.
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

/** input/output functions */
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
int  (*MMG5_saveMesh)(MMG5_pMesh mesh);
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
#ifdef SINGUL
/**
 * \param mesh pointer toward the mesh structure.
 * \param singul pointer toward the singul structure.
 * \return 0 if failed, 1 otherwise
 *
 * Read singul data. Here we suppose that the file contains the
 * singularities (corner, required, ridges....)
 *
 */
int  MMG5_loadSingul(MMG5_pMesh mesh,MMG5_pSingul singul);
#endif

/** deallocations */
#ifdef SINGUL
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param sing pointer toward the sing structure (only for insertion of
 * singularities mode).
 *
 * Deallocations before return.
 *
 */
void MMG5_Free_all(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSingul sing);
#else
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * Deallocations before return.
 *
 */
void MMG5_Free_all(MMG5_pMesh mesh, MMG5_pSol met);
#endif

#ifdef SINGUL
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param sing pointer toward the sing structure (only for insertion of singularities mode).
 *
 * Structure deallocations before return.
 *
 */
void MMG5_Free_structures(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSingul sing);
#else
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * Structure deallocations before return.
 *
 */
void MMG5_Free_structures(MMG5_pMesh mesh, MMG5_pSol met);
#endif

#ifdef SINGUL
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param sing pointer toward the sing structure (only for insertion of
 * singularities mode).
 *
 * File name deallocations before return.
 *
 */
void MMG5_Free_names(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSingul sing);
#else
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * File name deallocations before return.
 *
 */
void MMG5_Free_names(MMG5_pMesh mesh, MMG5_pSol met);
#endif


/** library */
#ifdef SINGUL
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param sing pointer toward the sing structure (only for insertion of
 * singularities mode).
 * \return Return \ref MMG5_SUCCESS if success,
 * \ref MMG5_LOWFAILURE if fail but a conform mesh is saved or
 * \ref MMG5_STRONGFAILURE if fail and we can't save the mesh.
 *
 * Main program for the library.
 *
 */
int  MMG5_mmg3dlib(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSingul singul);
#else
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return Return \ref MMG5_SUCCESS if success,
 * \ref MMG5_LOWFAILURE if fail but a conform mesh is saved or
 * \ref MMG5_STRONGFAILURE if fail and we can't save the mesh.
 *
 * Main program for the library.
 *
 */
int  MMG5_mmg3dlib(MMG5_pMesh mesh, MMG5_pSol sol);
#endif

/** for PAMPA library */
/** Options management */
#ifdef SINGUL
/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param sing pointer toward the sing structure (only for insertion of
 * singularities mode).
 * \return 1.
 * \note Developped for the PaMPA library interface.
 *
 * Store command line arguments.
 *
 */
int  MMG5_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSingul sing);
#else
/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 1.
 * \note Developped for the PaMPA library interface.
 *
 * Store command line arguments.
 *
 */
int  MMG5_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met);
#endif
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 1.
 * \note Developped for the PaMPA library interface.
 *
 * Read local parameters file. This file must have the same name as
 * the mesh with the \a .mmg3d5 extension or must be named \a
 * DEFAULT.mmg3d5.
 *
 */
int  MMG5_parsop(MMG5_pMesh mesh,MMG5_pSol met);
/**
 * \param prog pointer toward the program name.
 * \note Developped for the PaMPA library interface.
 *
 * Print help for mmg3d5 options.
 *
 */
void  MMG5_usage(char *prog);
/**
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 * \return 1.
 * \note Developped for the PaMPA library interface.
 *
 * Store the info structure in the mesh structure.
 *
 */
int  MMG5_stockOptions(MMG5_pMesh mesh, MMG5_Info *info);
/**
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 * \note Developped for the PaMPA library interface.
 *
 * Recover the info structure stored in the mesh structure.
 *
 */
void  MMG5_destockOptions(MMG5_pMesh mesh, MMG5_Info *info);

/** Checks */
#ifdef SINGUL
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param sing pointer toward the sing structure (only for insertion of
 * singularities mode).
 * \param critmin minimum quality for elements.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param eltab table of invalid elements.
 * \note Developped for the PaMPA library interface.
 *
 * Search invalid elements (in term of quality or edge length).
 *
 */
int  MMG5_mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSingul sing,
                     double critmin, double lmin, double lmax, int *eltab);
#else
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param critmin minimum quality for elements.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param eltab table of invalid elements.
 * \note Developped for the PaMPA library interface.
 *
 * Search invalid elements (in term of quality or edge length).
 *
 */
int MMG5_mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol sol,
                    double critmin, double lmin, double lmax, int *eltab);
#endif
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param critmin minimum quality for elements.
 * \param eltab pointer toward the table of invalid elements.
 * \note Developped for the PaMPA library interface.
 *
 * Store elements which have worse quality than \a critmin in \a eltab,
 * \a eltab is allocated and could contain \a mesh->ne elements.
 *
 */
void  MMG5_searchqua(MMG5_pMesh mesh, MMG5_pSol met, double critmin, int *eltab);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param eltab table of invalid elements.
 * \note Developped for the PaMPA library interface.
 *
 * Store in \a eltab elements which have edge lengths shorter than \a lmin
 * or longer than \a lmax, \a eltab is allocated and could contain \a mesh->ne
 * elements.
 *
 */
int  MMG5_searchlen(MMG5_pMesh mesh, MMG5_pSol met, double lmin, double lmax, int *eltab);

/** Utils */
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
 * \note Developped for the PaMPA library interface.
 *
 * Find the indices of the 4 adjacent elements of tetrahedron \a
 * kel. \f$v_i = 0\f$ if the \f$i^{th}\f$ face has no adjacent element
 * (so we are on a boundary face).
 *
 */
int MMG5_Get_adjaTet(MMG5_pMesh mesh,int kel, int* v0, int* v1, int* v2, int* v3);
/**
 * \param ca pointer toward the coordinates of the first edge's extremity.
 * \param cb pointer toward the coordinates of the second edge's extremity.
 * \param ma pointer toward the metric associated to the first edge's extremity.
 * \param mb pointer toward the metric associated to the second edge's extremity.
 * \return edge length.
 * \note Developped for the PaMPA library interface.
 *
 * Compute length of edge \f$[ca,cb]\f$ (with \a ca and \a cb
 * coordinates of edge extremities) according to the size
 * prescription.
 *
 */
double (*MMG5_lenedgCoor)(double *ca,double *cb,double *sa,double *sb);
/**
 * \param mesh pointer toward the mesh structure.
 * \param pack we pack the mesh at function begining if \f$pack=1\f$.
 * \return 0 if failed, 1 otherwise.
 *
 * Create table of adjacency. Set pack variable to 0 for a compact
 * mesh and to 1 for a mesh that need to be packed.
 *
 */
int  (*MMG5_hashTetra)(MMG5_pMesh mesh, int pack);

/** To associate function pointers without calling MMG5_mmg3dlib */
/**
 * \param mesh pointer toward the mesh structure (unused).
 * \note Developped for the PaMPA library interface.
 *
 * Set pointer for MMG5_saveMesh function.
 *
 */
void  MMG5_Set_saveFunc(MMG5_pMesh mesh);
/**
 * \param mesh pointer toward the mesh structure (unused).
 * \param met pointer toward the sol structure (unused).
 * \note Developped for the PaMPA library interface.
 *
 * Set function pointers for caltet, lenedg, defsiz and gradsiz.
 *
 */
void  MMG5_setfunc(MMG5_pMesh mesh,MMG5_pSol met);
/**
 * \param mesh pointer toward the mesh structure (unused).
 * \param met pointer toward the sol structure (unused).
 * \warning Copy of the \ref setfunc function of the \ref mmg3d/shared_func.h
 * file.
 * \note Developped for the PaMPA library interface.
 *
 * Set function pointers for lenedgeCoor, hashTetra and saveMesh.
 *
 */
void  MMG5_pampa_setfunc(MMG5_pMesh mesh,MMG5_pSol met);

#endif
