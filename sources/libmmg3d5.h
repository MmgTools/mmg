#ifndef _MMG3DLIB_H
#define _MMG3DLIB_H

#include "chrono.h"

/* Return Values */
#define MMG5_SUCCESS       0 /**< Return value for success */
#define MMG5_LOWFAILURE    1 /**< Return value if we have problem in the remesh
                            *   process but we can save a conform mesh */
#define MMG5_STRONGFAILURE 2 /**< Return value if we fail and the mesh is non-conform */

#define SIZE 0.75 /**< Size of the mesh of singularities inside the main mesh */

enum MMG5_optIntCod
  {
    MMG5_imprim,\
    MMG5_mem,\
    MMG5_debug,\
    MMG5_angle,\
    MMG5_iso,\
    MMG5_noinsert,\
    MMG5_noswap,\
    MMG5_nomove,\
    MMG5_renum,\
    MMG5_sing,\
  };
enum MMG5_optDblCod
  {
    MMG5_dhd,
    MMG5_hmin,
    MMG5_hmax,
    MMG5_hausd,
    MMG5_hgrad,
    MMG5_ls,
  };

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
  int           mem,sing;
  int           renum;
  char          imprim,ddebug,badkal,iso,fem;
  unsigned char noinsert, noswap, nomove;
  mytime        ctim[TIMEMAX];
} MMG5_Info;

typedef struct {
  int       ver,dim,type;
  int       npi,nai,nei,np,na,nt,ne,npmax,namax,ntmax,nemax,xpmax,xtmax;
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
} MMG5_Mesh;
typedef MMG5_Mesh  * MMG5_pMesh;

typedef struct {
  int       dim,ver,np,npmax,size,type;
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
  int      ns,na; /**< singular vertices and singular edges number */
  MMG5_psPoint  point;
  MMG5_pEdge    edge;
} MMG5_Singul; /**< structure to store the singularities of a mesh */
typedef MMG5_Singul * MMG5_pSingul;
/* end structures for insertion of singularities */

extern MMG5_Info info;

/** input/output functions */
int  MMG5_loadMesh(MMG5_pMesh );
int  MMG5_saveMesh(MMG5_pMesh );
int  MMG5_loadMet(MMG5_pSol );
int  MMG5_saveMet(MMG5_pMesh mesh,MMG5_pSol met);
#ifdef SINGUL
int  MMG5_loadSingul(MMG5_pSingul singul);
#endif

/** free the pMesh and pSol structures */
#ifdef SINGUL
void MMG5_freeAll(MMG5_pMesh,MMG5_pSol,MMG5_pSingul);
#else
void MMG5_freeAll(MMG5_pMesh,MMG5_pSol);
#endif

/** stock the user options (opt_i and opt_d) in the "info" structure */
void MMG5_stockOption(int opt_i[10],double opt_d[6],MMG5_pMesh mesh);

/** initialize to default values opt_i and opt_d */
void MMG5_mmg3dinit(int opt_i[10],double opt_d[6]);


/** library */
#ifdef SINGUL
int  MMG5_mmg3dlib(int opt_i[10],double opt_d[6],MMG5_pMesh mesh,MMG5_pSol sol,
                   MMG5_pSingul singul);
#else
int  MMG5_mmg3dlib(int opt_i[10],double opt_d[6],MMG5_pMesh mesh,MMG5_pSol sol);
#endif

#endif
