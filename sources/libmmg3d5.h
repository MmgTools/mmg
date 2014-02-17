#ifndef _MMG3DLIB_H
#define _MMG3DLIB_H

#include "chrono.h"

/* Return Values */
#define MMG5_SUCCESS       0 /**< Return value for success */
#define MMG5_LOWFAILURE    1 /**< Return value if we have problem in the remesh
                                process but we can save a conform mesh */
#define MMG5_STRONGFAILURE 2 /**< Return value if we fail and the mesh is non-conform */

#define SIZE 0.75 /**< Size of the mesh of singularities inside the main mesh */

enum MMG5_type /**< type of solutions */
{
  MMG5_Notype,
  MMG5_Scalar,
  MMG5_Vector,
  MMG5_Tensor
};
enum MMG5_entities /**< entities of MMG5 */
  {
    MMG5_Noentity,
    MMG5_Vertex,
    MMG5_Triangle
  };

enum MMG5_Param /**<  Options for mmg3d2 (integers) */
  {
    MMG5_IPARAM_verbose,           /**<  Tune level of verbosity, [-10..10]         */
    MMG5_IPARAM_mem,               /**<  Set memory size to n Mbytes                */
    MMG5_IPARAM_debug,             /**<  Turn on debug mode                         */
    MMG5_IPARAM_angle,             /**<  Turn on angle detection                    */
    MMG5_IPARAM_iso,               /**<  Level-set meshing                          */
    MMG5_IPARAM_noinsert,          /**<  No point insertion/deletion                */
    MMG5_IPARAM_noswap,            /**<  No edge or face flipping                   */
    MMG5_IPARAM_nomove,            /**<  No point relocation                        */
    MMG5_IPARAM_numberOfLocalParam,/**<  Number of local parameters                 */
    MMG5_IPARAM_renum,             /**<  Turn on point relocation with Scotch       */
    MMG5_IPARAM_sing,              /**<  Turn on the insertion of singularities     */
    MMG5_IPARAM_bucket,            /**<  Specify the size of the bucket per dimension (DELAUNAY)     */
    MMG5_DPARAM_angleDetection,    /**<  Value for angle detection                  */
    MMG5_DPARAM_hmin,              /**<  Minimal mesh size                          */
    MMG5_DPARAM_hmax,              /**<  Maximal mesh size                          */
    MMG5_DPARAM_hausd,             /**<  control global Hausdorff distance          */
    /*                                   (on all the boundary surfaces of the mesh) */
    MMG5_DPARAM_hgrad,             /**<  control gradation                          */
    MMG5_DPARAM_ls,                /**<  Value of level-set (not use for now)       */
    MMG5_PARAM_size,               /**<  Size of table of double parameters         */
  };

typedef struct {
  double   hausd;
  int      ref;
  char     elt;
} MMG5_Par; /** specific parameters */
typedef MMG5_Par * MMG5_pPar;

typedef struct {
  double   c[3]; /**< coordinates of point */
  int      ref; /**< ref of point */
  int      xp; /**< surface point number */
  int      tmp; /**< tmp: numero of points for the saving (we don't count the unused points)*/
  int      flag; /**< flag to know if we have already treated the point */
  char     tag; /**< contains binary flags :
                   if tag=23=16+4+2+1, then the point is MG_REF, MG_GEO,MG_REQ,MG_BDY */
  char     tagdel; /**< tag for delaunay */
} MMG5_Point;
typedef MMG5_Point * MMG5_pPoint;

typedef struct {
  double   n1[3],n2[3]; /**< normals at boundary vertex;
                           n1!=n2 if the vertex belong to a ridge */
  double   t[3]; /** tangeant at vertex */
} MMG5_xPoint;
typedef MMG5_xPoint * MMG5_pxPoint;

typedef struct {
  int      a,b,ref; /**< extremities and ref of edges */
  char     tag; /**< binary flags */
} MMG5_Edge;
typedef MMG5_Edge * MMG5_pEdge;

typedef struct {
  int      v[3],ref; /** vertices and ref of tria */
  int      base;
  int      edg[3]; /**< edg[i] contains the ref of the i^th edge of triangle */
  int      flag;
  char     tag[3]; /**< tag[i] contains the tag associated to th i^th edge of tri */
} MMG5_Tria;
typedef MMG5_Tria * MMG5_pTria;

typedef struct {
  int      v[4],ref; /** vertices and ref of tetra */
  int      base,mark; //CECILE rajout mark pour delaunay
  int      xt;   /**< xt : number of the surfaces xtetra */
  int      flag;
  char     tag;
  double   qual; /**< quality of element */
} MMG5_Tetra;
typedef MMG5_Tetra * MMG5_pTetra;

typedef struct {
  int      ref[4]; /**< ref[i] : ref de la face opp au pt i;*/
  int      edg[6]; /**< edg[i] contains the ref of the i^th edge of tet */
  char     ftag[4]; /**< ftag[i] contains the tag associated to the i^th face of tet */
  char     tag[6]; /**< tag[i] contains the tag associated to the i^th edge of tet */
  char     ori; /**< orientation of tris of the tetra:
                 * i^th bit of ori is set to 0 when the i^th face is bad orientated */
} MMG5_xTetra;
typedef MMG5_xTetra * MMG5_pxTetra;

/** to store geometric edges */
typedef struct {
  int   a,b,ref,nxt;
  char  tag;
} MMG5_hgeom;
typedef struct {
  int         siz,max,nxt;
  MMG5_hgeom  *geom;
} MMG5_HGeom;

typedef struct {
  double        dhd,hmin,hmax,hgrad,hausd,min[3],max[3],delta,ls;
  int           mem,sing,npar,npari;
  int           renum;
  char          imprim,ddebug,badkal,iso,fem;
  unsigned char noinsert, noswap, nomove;
#ifndef PATTERN
  int           bucket;
#endif
  MMG5_pPar     par;
} MMG5_Info;

typedef struct {
  int       ver,dim,type;
  long long memMax; /**< maximum memory available */
  long long memCur; /**< current memory used */
  double    gap; /**< gap for table reallocation */
  int       npi,nti,nai,nei,np,na,nt,ne,npmax,namax,ntmax,nemax,xpmax,xtmax;
  int       base; /**< used with flag to know if an entity has been treated */
  int       mark;//CECILE rajout mark pour delaunay
  int       xp,xt; /**< nb of surfaces points/triangles */
  int       npnil,nenil; /**< nb of first unused point/element */
  int      *adja; /**< tab of tetrahedra adjacency : if adjt[4*i+1+j]=4*k+l then
                     the i^th and k^th tets are adjacent and share their
                     faces j and l (resp.) */
  int      *adjt; /**< tab of triangles adjacency : if adjt[3*i+1+j]=3*k+l then
                     the i^th and k^th triangles are adjacent and share their
                     edges j and l (resp.) */
  char     *namein,*nameout;

  MMG5_pPoint    point;
  MMG5_pxPoint   xpoint;
  MMG5_pTetra    tetra;
  MMG5_pxTetra   xtetra;
  MMG5_pTria     tria;
  MMG5_pEdge     edge;
  MMG5_HGeom     htab;
  MMG5_Info      info;
} MMG5_Mesh;
typedef MMG5_Mesh  * MMG5_pMesh;

typedef struct {
  int       dim,ver,np,npi,npmax,size,type;
  double   *m;
  char     *namein,*nameout;
} MMG5_Sol;
typedef MMG5_Sol * MMG5_pSol;

/* structures only use for insertion of singularities (#ifdef SINGUL) */
typedef struct {
  double         c[3],n[3];
  int            flag,tmp,tet;
  unsigned char  tag;
} MMG5_sPoint;
typedef MMG5_sPoint * MMG5_psPoint; /**< structure to store singular points */

typedef struct {
  char     *namein; /**< name of mesh */
  double   min[3],max[3]; /**< min and max of coordinates for rescaling */
  int      nsi,ns,na; /**< singular vertices and singular edges number */
  MMG5_psPoint  point;
  MMG5_pEdge    edge;
} MMG5_Singul; /**< structure to store the singularities of a mesh */
typedef MMG5_Singul * MMG5_pSingul;
/* end structures for insertion of singularities */

/*----------------------------- functions header -----------------------------*/
/** Initialization functions */
/* init structures */
#ifndef SINGUL
void  MMG5_Init_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol);
void  MMG5_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);
#else
void  MMG5_Init_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol, MMG5_pSingul *sing);
void  MMG5_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSingul sing);
#endif
void  MMG5_Init_parameters(MMG5_pMesh mesh);

/* init file names */
int  MMG5_Set_inputMeshName(MMG5_pMesh mesh, char* meshin);
int  MMG5_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solin);
int  MMG5_Set_outputMeshName(MMG5_pMesh mesh, char* meshout);
int  MMG5_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solout);
#ifdef SINGUL
int  MMG5_Set_inputSingulName(MMG5_pMesh mesh,MMG5_pSingul sing, char* singin);
#endif

/* init structure sizes */
int  MMG5_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, int np, int typSol);
int  MMG5_Set_meshSize(MMG5_pMesh mesh, int np, int ne, int nt, int na);
#ifdef SINGUL
int  MMG5_Set_singulSize(MMG5_pMesh mesh,MMG5_pSingul sing, int np, int na);
#endif

/* init structure datas */
int  MMG5_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                     double c2, int ref,int pos);
int  MMG5_Set_tetrahedra(MMG5_pMesh mesh, int v0, int v1,
                         int v2, int v3, int ref, int pos);
int  MMG5_Set_triangle(MMG5_pMesh mesh, int v0, int v1,
                       int v2, int ref,int pos);
int  MMG5_Set_edges(MMG5_pMesh mesh, int v0, int v1, int ref,int pos);
int  MMG5_Set_corner(MMG5_pMesh mesh, int pos);
int  MMG5_Set_requiredVertex(MMG5_pMesh mesh, int pos);
int  MMG5_Set_requiredTetrahedra(MMG5_pMesh mesh, int pos);
int  MMG5_Set_requiredTriangle(MMG5_pMesh mesh, int pos);
int  MMG5_Set_ridge(MMG5_pMesh mesh, int pos);
int  MMG5_Set_requiredEdge(MMG5_pMesh mesh, int pos);
int  MMG5_Set_scalarSol(MMG5_pSol met, double s,int pos);
#ifndef SINGUL
void MMG5_Set_handGivenMesh(MMG5_pMesh mesh);
#else
int  MMG5_Set_singulVertex(MMG5_pSingul sing, double c0,
                           double c1, double c2, int typ,int pos);
int  MMG5_Set_singulEdge(MMG5_pSingul sing, int v0, int v1, int ref,int pos);
int  MMG5_Set_singulCorner(MMG5_pSingul sing, int pos);
int  MMG5_Set_singulRequiredVertex(MMG5_pSingul sing, int pos);
int  MMG5_Set_singulRidge(MMG5_pSingul sing, int pos);
int  MMG5_Set_singulRequiredEdge(MMG5_pSingul sing, int pos);
void MMG5_Set_handGivenMesh(MMG5_pMesh mesh);
#endif


/* check init */
int MMG5_Chk_meshData(MMG5_pMesh mesh, MMG5_pSol met);

/** functions to set parameters */
int  MMG5_Set_iparameters(MMG5_pMesh mesh,MMG5_pSol sol, int iparam, int val);
int  MMG5_Set_dparameters(MMG5_pMesh mesh,MMG5_pSol sol, int dparam, double val);
int  MMG5_Set_localParameters(MMG5_pMesh mesh, MMG5_pSol sol, int typ, int ref, double val);

/** recover datas */
int  MMG5_Get_meshSize(MMG5_pMesh mesh, int* np, int* ne, int* nt, int* na);
int  MMG5_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, int* np,
                      int* typSol);
int  MMG5_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, int* ref,
                     int* isCorner, int* isRequired);
int  MMG5_Get_tetrahedra(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* v3,
                         int* ref, int* isRequired);
int  MMG5_Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref,
                       int* isRequired);
int  MMG5_Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired);
int  MMG5_Get_scalarSol(MMG5_pSol met, double* s);

/** input/output functions */
int  MMG5_loadMesh(MMG5_pMesh );
int  (*MMG5_saveMesh)(MMG5_pMesh );
int  MMG5_loadMet(MMG5_pMesh,MMG5_pSol );
int  MMG5_saveMet(MMG5_pMesh mesh, MMG5_pSol met);
#ifdef SINGUL
int  MMG5_loadSingul(MMG5_pMesh,MMG5_pSingul singul);
#endif

/** deallocations */
#ifdef SINGUL
void MMG5_Free_all(MMG5_pMesh, MMG5_pSol, MMG5_pSingul);
#else
void MMG5_Free_all(MMG5_pMesh, MMG5_pSol);
#endif

#ifdef SINGUL
void MMG5_Free_structures(MMG5_pMesh, MMG5_pSol, MMG5_pSingul);
#else
void MMG5_Free_structures(MMG5_pMesh, MMG5_pSol);
#endif

#ifdef SINGUL
void MMG5_Free_names(MMG5_pMesh, MMG5_pSol, MMG5_pSingul);
#else
void MMG5_Free_names(MMG5_pMesh, MMG5_pSol);
#endif



/** library */
#ifdef SINGUL
int  MMG5_mmg3dlib(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSingul singul);
#else
int  MMG5_mmg3dlib(MMG5_pMesh mesh, MMG5_pSol sol);
#endif

/** for PAMPA library */
/** Options management */
#ifdef SINGUL
int  MMG5_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSingul sing);
#else
int  MMG5_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met);
#endif
int  MMG5_parsop(MMG5_pMesh mesh,MMG5_pSol met);
void  MMG5_usage(char *prog);
void  MMG5_stockOptions(MMG5_pMesh mesh, MMG5_Info *info);
void  MMG5_destockOptions(MMG5_pMesh mesh, MMG5_Info *info);

/** Checks */
#ifdef SINGUL
int  MMG5_mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSingul sing,
                     double critmin, double lmin, double lmax, int *eltab);
#else
int mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol sol,
               double critmin, double lmin, double lmax, int *eltab);
#endif
void  MMG5_searchqua(MMG5_pMesh mesh, MMG5_pSol met, double critmin, int *eltab);
int  MMG5_searchlen(MMG5_pMesh mesh, MMG5_pSol met, double lmin, double lmax, int *eltab);

/** Utils */
int    MMG5_Get_adjaTet(MMG5_pMesh mesh,int kel, int*, int*, int*, int*);
double (*MMG5_lenedgCoor)(double *ca,double *cb,double *sa,double *sb);

#endif
