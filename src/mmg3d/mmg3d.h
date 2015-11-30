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

#ifndef _MMG3D_H
#define _MMG3D_H

#include "libmmg3d.h"

#define MG_SMSGN(a,b)  (((double)(a)*(double)(b) > (0.0)) ? (1) : (0))

/** Free allocated pointers of mesh and sol structure and return value val */
#define _MMG5_RETURN_AND_FREE(mesh,met,disp,val)do  \
  {                                                 \
    MMG3D_Free_all(mesh,met,disp);                  \
    return(val);                                    \
  }while(0)

/** Reallocation of point table and sol table and creation
    of point ip with coordinates o and tag tag*/
#define _MMG5_POINT_REALLOC(mesh,sol,ip,wantedGap,law,o,tag ) do        \
  {                                                                     \
    int klink;                                                          \
                                                                        \
    _MMG5_TAB_RECALLOC(mesh,mesh->point,mesh->npmax,wantedGap,MMG5_Point, \
                       "larger point table",law);                       \
                                                                        \
    mesh->npnil = mesh->np+1;                                           \
    for (klink=mesh->npnil; klink<mesh->npmax-1; klink++)               \
      mesh->point[klink].tmp  = klink+1;                                \
                                                                        \
    /* solution */                                                      \
    if ( sol->m ) {                                                     \
      _MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax-sol->npmax))*sizeof(double), \
                    "larger solution",law);                             \
      _MMG5_SAFE_REALLOC(sol->m,sol->size*(mesh->npmax+1),double,"larger solution"); \
    }                                                                   \
    sol->npmax = mesh->npmax;                                           \
                                                                        \
    /* We try again to add the point */                                 \
    ip = _MMG3D_newPt(mesh,o,tag);                                       \
    if ( !ip ) {law;}                                                   \
  }while(0)

/** Reallocation of point table, sol table and bucket table and creation
    of point ip with coordinates o and tag tag*/
#define _MMG5_POINT_AND_BUCKET_REALLOC(mesh,sol,ip,wantedGap,law,o,tag ) do \
  {                                                                     \
    int klink,gap;                                                      \
                                                                        \
    if ( (mesh->memMax-mesh->memCur) <                                  \
         (long long) (wantedGap*mesh->npmax*                            \
                      (sizeof(MMG5_Point)+sol->size*sizeof(int))) ) {   \
      gap = (int)((mesh->memMax-mesh->memCur)/                          \
                  (sizeof(MMG5_Point)+sol->size*sizeof(int)));          \
      if(gap < 1) {                                                     \
        fprintf(stdout,"  ## Error:");                                  \
        fprintf(stdout," unable to allocate %s.\n","larger point/bucket table"); \
        fprintf(stdout,"  ## Check the mesh size or ");                 \
        fprintf(stdout,"increase maximal authorized memory with the -m option.\n"); \
        law;                                                            \
      }                                                                 \
    }                                                                   \
    else                                                                \
      gap = (int)(wantedGap*mesh->npmax);                               \
                                                                        \
    _MMG5_ADD_MEM(mesh,gap*(sizeof(MMG5_Point)+sizeof(int)),            \
                  "point and bucket",law);                              \
    _MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,                      \
                        mesh->npmax+gap+1,MMG5_Point,"larger point table"); \
    _MMG5_SAFE_RECALLOC(bucket->link,mesh->npmax+1,                     \
                        mesh->npmax+gap+1,int,"larger bucket table");   \
    mesh->npmax = mesh->npmax+gap;                                      \
                                                                        \
    mesh->npnil = mesh->np+1;                                           \
    for (klink=mesh->npnil; klink<mesh->npmax-1; klink++)               \
      mesh->point[klink].tmp  = klink+1;                                \
                                                                        \
    /* solution */                                                      \
    if ( sol->m ) {                                                     \
      _MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax-sol->npmax))*sizeof(double), \
                    "larger solution",law);                             \
      _MMG5_SAFE_REALLOC(sol->m,sol->size*(mesh->npmax+1),double,"larger solution"); \
    }                                                                   \
    sol->npmax = mesh->npmax;                                           \
                                                                        \
    /* We try again to add the point */                                 \
    ip = _MMG3D_newPt(mesh,o,tag);                                       \
    if ( !ip ) {law;}                                                   \
  }while(0)

/** Reallocation of tetra table and creation
    of tetra jel */
#define _MMG5_TETRA_REALLOC(mesh,jel,wantedGap,law ) do                 \
  {                                                                     \
    int klink,oldSiz;                                                   \
                                                                        \
    oldSiz = mesh->nemax;                                               \
    _MMG5_TAB_RECALLOC(mesh,mesh->tetra,mesh->nemax,wantedGap,MMG5_Tetra, \
                       "larger tetra table",law);                       \
                                                                        \
    mesh->nenil = mesh->ne+1;                                           \
    for (klink=mesh->nenil; klink<mesh->nemax-1; klink++)               \
      mesh->tetra[klink].v[3]  = klink+1;                               \
                                                                        \
    if ( mesh->adja ) {                                                 \
      /* adja table */                                                  \
      _MMG5_ADD_MEM(mesh,4*(mesh->nemax-oldSiz)*sizeof(int),            \
                    "larger adja table",law);                           \
      _MMG5_SAFE_RECALLOC(mesh->adja,4*mesh->ne+5,4*mesh->nemax+5,int   \
                          ,"larger adja table");                        \
    }                                                                   \
                                                                        \
    /* We try again to add the point */                                 \
    jel = _MMG3D_newElt(mesh);                                           \
    if ( !jel ) {law;}                                                  \
  }while(0)

/* numerical accuracy */
#define _MMG5_ALPHAD    20.7846096908265    //0.04811252243247      /* 12*sqrt(3) */
#define _MMG5_LLONG     2.5//2.0   // 1.414213562373
#define _MMG5_LSHRT     0.3  // 0.707106781186
#define _MMG5_LOPTL     1.3
#define _MMG5_LOPTS     0.6

#define _MMG5_LMAX      10240
#define _MMG5_BADKAL    0.2
#define _MMG5_NULKAL    1.e-30

#define _MMG5_NPMAX  1000000 //200000
#define _MMG5_NAMAX   200000 //40000
#define _MMG5_NTMAX  2000000 //400000
#define _MMG5_NEMAX  6000000 //1200000

#define _MMG5_BOXSIZE 500

#define _MMG5_SHORTMAX     0x7fff

/* Domain refs in iso mode */
#define MG_PLUS    2
#define MG_MINUS   3


/*! \var next vertex of tetra: {1,2,3,0,1,2,3} */
static const unsigned char _MMG5_inxt3[7] = { 1,2,3,0,1,2,3 };
/*! \var previous vertex of tetra: {3,0,1,2,3,0,1} */
static const unsigned char _MMG5_iprv3[7] = { 3,0,1,2,3,0,1 };
/*! \var idir[i] : vertices of face opposite to vertex i */
static const unsigned char _MMG5_idir[4][3] = { {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1} };
/* \var< idirinv[i][j] : num of the jth point in the ith face */
static const          char _MMG5_idirinv[4][4] = {{-1,0,1,2},{0,-1,2,1},{0,1,-1,2},{0,2,1,-1}};
/*! \var iarf[i]: edges of face opposite to vertex i */
static const unsigned char _MMG5_iarf[4][3] = { {5,4,3}, {5,1,2}, {4,2,0}, {3,0,1} };
/*! \var num of the j^th edge in the i^th face */
static const unsigned char _MMG5_iarfinv[4][6] = { {-1,-1,-1,2,1,0}, {-1,1,2,-1,-1,0},{2,-1,1,-1,0,-1},{1,2,-1,0,-1,-1}};
/*! \var vertices of extremities of the edges of the tetra */
static const unsigned char _MMG5_iare[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
/*! \var ifar[i][]: faces sharing the ith edge of the tetra */
static const unsigned char _MMG5_ifar[6][2] = { {2,3}, {1,3}, {1,2}, {0,3}, {0,2}, {0,1} };
/*! \var isar[i][]: vertices of extremities of the edge opposite to the ith edge */
static const unsigned char _MMG5_isar[6][2] = { {2,3}, {3,1}, {1,2}, {0,3}, {2,0}, {0,1} };
/*! \var arpt[i]: edges passing through vertex i */
static const unsigned char _MMG5_arpt[4][3] = { {0,1,2}, {0,4,3}, {1,3,5}, {2,5,4} };


typedef struct {
  int     size;
  int    *head;
  int    *link;
} _MMG5_Bucket;
typedef _MMG5_Bucket * _MMG5_pBucket;

/* bucket */
_MMG5_pBucket _MMG5_newBucket(MMG5_pMesh ,int );
int     _MMG5_addBucket(MMG5_pMesh ,_MMG5_pBucket ,int );
int     _MMG5_delBucket(MMG5_pMesh ,_MMG5_pBucket ,int );
int     _MMG5_buckin_iso(MMG5_pMesh mesh,MMG5_pSol sol,_MMG5_pBucket bucket,int ip);
int     _MMG5_buckin_ani(MMG5_pMesh mesh,MMG5_pSol sol,_MMG5_pBucket bucket,int ip);

/* prototypes */
extern double _MMG5_det3pt1vec(double c0[3],double c1[3],double c2[3],double v[3]);
extern double _MMG5_det4pt(double c0[3],double c1[3],double c2[3],double c3[3]);
extern double _MMG5_orvol(MMG5_pPoint point,int *v);
extern int _MMG5_directsurfball(MMG5_pMesh mesh, int ip, int *list, int ilist, double n[3]);

int  _MMG3D_newPt(MMG5_pMesh mesh,double c[3],char tag);
int  _MMG3D_newElt(MMG5_pMesh mesh);
void _MMG3D_delElt(MMG5_pMesh mesh,int iel);
void _MMG3D_delPt(MMG5_pMesh mesh,int ip);
int  _MMG5_zaldy(MMG5_pMesh mesh);
void _MMG5_freeXTets(MMG5_pMesh mesh);
char _MMG5_chkedg(MMG5_pMesh mesh,MMG5_pTria pt,char ori);
int  _MMG5_chkNumberOfTri(MMG5_pMesh mesh);
void _MMG5_tet2tri(MMG5_pMesh mesh,int k,char ie,MMG5_Tria *ptt);
int    _MMG5_mmg3dBezierCP(MMG5_pMesh mesh,MMG5_Tria *pt,_MMG5_pBezier pb,char ori);
extern int    _MMG5_BezierTgt(double c1[3],double c2[3],double n1[3],double n2[3],double t1[3],double t2[3]);
extern double _MMG5_BezierGeod(double c1[3], double c2[3], double t1[3], double t2[3]);
int  _MMG3D_bezierInt(_MMG5_pBezier pb,double uv[2],double o[3],double no[3],double to[3]);
extern int  _MMG5_BezierReg(MMG5_pMesh mesh,int ip0, int ip1, double s, double v[3], double *o, double *no);
extern int  _MMG5_BezierRef(MMG5_pMesh mesh,int ip0, int ip1, double s, double *o, double *no, double *to);
extern int  _MMG5_BezierEdge(MMG5_pMesh mesh,int ip0, int ip1, double b0[3], double b1[3],char isrid, double v[3]);
extern int  _MMG5_BezierRidge(MMG5_pMesh mesh,int ip0, int ip1, double s, double *o, double *no1, double *no2, double *to);
extern int  _MMG5_BezierNom(MMG5_pMesh mesh,int ip0,int ip1,double s,double *o,double *no,double *to);
extern int  _MMG5_norface(MMG5_pMesh mesh ,int k, int iface, double v[3]);
int  _MMG5_boulernm (MMG5_pMesh mesh, int start, int ip, int *ng, int *nr);
int  _MMG5_boulenm(MMG5_pMesh mesh, int start, int ip, int iface, double n[3],double t[3]);
int  _MMG5_boulevolp(MMG5_pMesh mesh, int start, int ip, int * list);
int  _MMG5_boulesurfvolp(MMG5_pMesh mesh,int start,int ip,int iface,int *listv,
                         int *ilistv,int *lists,int*ilists, int isnm);
int  _MMG5_bouletrid(MMG5_pMesh,int,int,int,int *,int *,int *,int *,int *,int *);
int  _MMG5_startedgsurfball(MMG5_pMesh mesh,int nump,int numq,int *list,int ilist);
int  _MMG5_srcbdy(MMG5_pMesh mesh,int start,int ia);
int  _MMG5_coquil(MMG5_pMesh mesh, int start, int ia, int * list);
int  _MMG5_coquilface(MMG5_pMesh mesh, int start, int ia, int * list, int * it1, int *it2);
void _MMG5_coquilFaceErrorMessage(MMG5_pMesh mesh, int k1, int k2);
char _MMG5_coquilTravel(MMG5_pMesh, int, int, int*, int*, char*, int*);
void _MMG5_openCoquilTravel(MMG5_pMesh, int, int, int*, int*, char*, int*);
extern int  _MMG5_settag(MMG5_pMesh,int,int,int,int);
int  _MMG5_setNmTag(MMG5_pMesh mesh, _MMG5_Hash *hash);
int  _MMG5_chkcol_int(MMG5_pMesh ,MMG5_pSol met,int,char,char,int *,char typchk);
int  _MMG5_chkcol_bdy(MMG5_pMesh,MMG5_pSol met,int,char,char,int *,char typchk);
int  _MMG5_chkmanicoll(MMG5_pMesh,int,int,int,int,int,char,char);
int  _MMG5_chkmani(MMG5_pMesh mesh);
int  _MMG5_colver(MMG5_pMesh,MMG5_pSol,int *,int,char,char);
int  _MMG3D_analys(MMG5_pMesh mesh);
int  _MMG3D_hashTria(MMG5_pMesh mesh, _MMG5_Hash*);
int  _MMG5_hashPop(_MMG5_Hash *hash,int a,int b);
int  _MMG5_hPop(MMG5_HGeom *hash,int a,int b,int *ref,char *tag);
int  _MMG5_hTag(MMG5_HGeom *hash,int a,int b,int ref,char tag);
int  _MMG5_hGet(MMG5_HGeom *hash,int a,int b,int *ref,char *tag);
void _MMG5_hEdge(MMG5_pMesh mesh,int a,int b,int ref,char tag);
int  _MMG5_hNew(MMG5_HGeom *hash,int hsiz,int hmax,int secure);
int  _MMG5_hGeom(MMG5_pMesh mesh);
int  _MMG5_bdryTria(MMG5_pMesh );
int  _MMG5_bdryIso(MMG5_pMesh );
int  _MMG5_bdrySet(MMG5_pMesh );
int  _MMG5_bdryUpdate(MMG5_pMesh );
int  _MMG5_bdryPerm(MMG5_pMesh );
int  _MMG5_chkfemtopo(MMG5_pMesh mesh);
int  _MMG5_cntbdypt(MMG5_pMesh mesh, int nump);
extern double _MMG5_lenedg_ani(MMG5_pMesh ,MMG5_pSol ,int,  MMG5_pTetra);
extern double _MMG5_lenedg_iso(MMG5_pMesh ,MMG5_pSol ,int,  MMG5_pTetra);
long long _MMG5_memSize(void);
void _MMG3D_memOption(MMG5_pMesh mesh);
int  _MMG5_mmg3d1_pattern(MMG5_pMesh ,MMG5_pSol );
int  _MMG5_mmg3d1_delone(MMG5_pMesh ,MMG5_pSol );
int  _MMG5_mmg3d2(MMG5_pMesh ,MMG5_pSol );
int  _MMG5_mmg3dChkmsh(MMG5_pMesh,int,int);
int  _MMG3D_split1_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
void _MMG5_split1(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp);
int  _MMG5_split1b(MMG5_pMesh,MMG5_pSol,int*,int,int,int,char);
int  _MMG5_split2sf_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
void _MMG5_split2sf(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
void _MMG5_split2(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
int  _MMG3D_split3_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
void _MMG5_split3(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
void _MMG5_split3cone(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
void _MMG5_split3op(MMG5_pMesh mesh, MMG5_pSol met, int k, int vx[6],char);
void _MMG5_split4sf(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
void _MMG5_split4op(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
void _MMG5_split5(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
void _MMG5_split6(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
int  _MMG5_split4bar(MMG5_pMesh mesh,MMG5_pSol met,int k,char);
int  _MMG3D_simbulgept(MMG5_pMesh mesh,MMG5_pSol met, int *list, int ilist,int);
void _MMG5_nsort(int ,double *,char *);
extern double _MMG5_orcal(MMG5_pMesh mesh,MMG5_pSol met,int iel);
int    _MMG5_movintpt(MMG5_pMesh ,MMG5_pSol, int *, int , int );
int    _MMG5_movbdyregpt(MMG5_pMesh, MMG5_pSol, int*, int, int*, int);
int    _MMG5_movbdyrefpt(MMG5_pMesh, MMG5_pSol, int*, int, int*, int);
int    _MMG5_movbdynompt(MMG5_pMesh, MMG5_pSol, int*, int, int*, int);
int    _MMG5_movbdyridpt(MMG5_pMesh, MMG5_pSol, int*, int, int*, int);
int  _MMG5_chkswpbdy(MMG5_pMesh, MMG5_pSol,int*, int, int, int,char);
int  _MMG5_swpbdy(MMG5_pMesh,MMG5_pSol,int*,int,int,_MMG5_pBucket,char);
int  _MMG5_swpgen(MMG5_pMesh,MMG5_pSol,int, int, int*,_MMG5_pBucket,char);
int  _MMG5_chkswpgen(MMG5_pMesh,MMG5_pSol,int,int,int*,int*,double,char);
int  _MMG5_srcface(MMG5_pMesh mesh,int n0,int n1,int n2);
int _MMG5_chkptonbdy(MMG5_pMesh,int);
double _MMG5_orcal_poi(double a[3],double b[3],double c[3],double d[3]);
int _MMG5_countelt(MMG5_pMesh mesh,MMG5_pSol sol, double *weightelt, long *npcible);
int MMG3D_opttyp(MMG5_pMesh mesh, MMG5_pSol met,_MMG5_pBucket bucket);
int _MMG5_trydisp(MMG5_pMesh,double *,short);
int _MMG5_dichodisp(MMG5_pMesh,double *);
int _MMG5_lapantilap(MMG5_pMesh,double *);
int _MMG5_ppgdisp(MMG5_pMesh,double *);
int _MMG5_denoisbdy(MMG5_pMesh);
void _MMG3D_inqua(MMG5_pMesh mesh,MMG5_pSol met);
void _MMG3D_outqua(MMG5_pMesh mesh,MMG5_pSol met);
int  _MMG5_badelt(MMG5_pMesh mesh,MMG5_pSol met);
int _MMG3D_prilen(MMG5_pMesh mesh,MMG5_pSol met,char);
int _MMG5_DoSol(MMG5_pMesh mesh,MMG5_pSol met);
void _MMG5_defaultValues(MMG5_pMesh);
int  _MMG5_intridmet(MMG5_pMesh,MMG5_pSol,int,int,double,double*,double*);
int  _MMG5_intregmet(MMG5_pMesh,MMG5_pSol,int,char,double, double*);
int  _MMG5_intvolmet(MMG5_pMesh,MMG5_pSol,int,char,double, double*);
int  _MMG3D_saveAllMesh(MMG5_pMesh mesh);

/* useful functions to debug */
int  _MMG3D_indElt(MMG5_pMesh mesh,int kel);
int  _MMG3D_indPt(MMG5_pMesh mesh,int kp);
void _MMG5_printTetra(MMG5_pMesh mesh,char* fileName);


#ifdef USE_SCOTCH
int _MMG5_mmg3dRenumbering(int vertBoxNbr, MMG5_pMesh mesh, MMG5_pSol sol);
#endif

int    _MMG5_meancur(MMG5_pMesh mesh,int np,double c[3],int ilist,int *list,double h[3]);
double _MMG5_surftri(MMG5_pMesh,int,int);
double _MMG5_timestepMCF(MMG5_pMesh,double);
int    _MMG5_bdyMCF(MMG5_pMesh);
double _MMG5_volint(MMG5_pMesh);

/* Lagrangian mode functions */
double _MMG5_estavglen(MMG5_pMesh);
int   _MMG5_stiffelt(MMG5_pMesh,int,double*,double*);
int  _MMG5_mmg3d3(MMG5_pMesh ,MMG5_pSol, MMG5_pSol );
int  _MMG5_velextLS(MMG5_pMesh ,MMG5_pSol );
int _MMG5_saveDisp(MMG5_pMesh ,MMG5_pSol );

/* Delaunay functions*/
int _MMG5_delone(MMG5_pMesh mesh,MMG5_pSol sol,int ip,int *list,int ilist);
int _MMG5_cavity_iso(MMG5_pMesh mesh,MMG5_pSol sol,int iel,int ip,int *list,int lon);
int _MMG5_cavity_ani(MMG5_pMesh mesh,MMG5_pSol sol,int iel,int ip,int *list,int lon);
int _MMG5_cenrad_iso(MMG5_pMesh mesh,double *ct,double *c,double *rad);
int _MMG5_cenrad_ani(MMG5_pMesh mesh,double *ct,double *m,double *c,double *rad);

/*mmg3d1.c*/
void _MMG5_tet2tri(MMG5_pMesh mesh,int k,char ie,MMG5_Tria *ptt);
int  _MMG3D_dichoto(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  _MMG3D_dichoto1b(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ret,int);
char _MMG5_chkedg(MMG5_pMesh mesh,MMG5_Tria *pt,char ori);
int  _MMG5_anatet(MMG5_pMesh mesh,MMG5_pSol met,char typchk, int patternMode) ;
int  _MMG5_movtet(MMG5_pMesh mesh,MMG5_pSol met,int maxitin);
int  _MMG5_swpmsh(MMG5_pMesh mesh,MMG5_pSol met,_MMG5_pBucket bucket, int);
int  _MMG5_swptet(MMG5_pMesh mesh,MMG5_pSol met,double,_MMG5_pBucket, int);

/* pointers */
/* init structures */
void  _MMG5_Init_parameters(MMG5_pMesh mesh);
/* iso/aniso computations */
extern double _MMG5_caltet_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt);
extern double _MMG5_caltet_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt);
double _MMG5_caltet33_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt);
extern double _MMG5_lenedgCoor_ani(double*, double*, double*, double*);
extern double _MMG5_lenedgCoor_iso(double*, double*, double*, double*);
int    _MMG5_intmet_iso(MMG5_pMesh,MMG5_pSol,int,char,int, double);
int    _MMG5_intmet_ani(MMG5_pMesh,MMG5_pSol,int,char,int, double);
int    _MMG3D_intmet33_ani(MMG5_pMesh,MMG5_pSol,int,char,int, double);
int    _MMG5_interp4bar_ani(MMG5_pMesh,MMG5_pSol,int,int,double *);
int    _MMG5_interp4bar33_ani(MMG5_pMesh,MMG5_pSol,int,int,double *);
int    _MMG5_interp4bar_iso(MMG5_pMesh,MMG5_pSol,int,int,double *);
int    _MMG3D_defsiz_iso(MMG5_pMesh,MMG5_pSol );
int    _MMG3D_defsiz_ani(MMG5_pMesh ,MMG5_pSol );
int    _MMG5_gradsiz_iso(MMG5_pMesh ,MMG5_pSol );
int    _MMG5_gradsiz_ani(MMG5_pMesh ,MMG5_pSol );
extern int    _MMG5_moymet(MMG5_pMesh ,MMG5_pSol ,MMG5_pTetra ,double *);
double _MMG5_lenedgspl_ani(MMG5_pMesh  ,MMG5_pSol , int , MMG5_pTetra );
extern double _MMG5_lenedgspl33_ani(MMG5_pMesh  ,MMG5_pSol , int , MMG5_pTetra );
double _MMG5_lenedgspl_iso(MMG5_pMesh  ,MMG5_pSol , int , MMG5_pTetra );
extern double _MMG5_lenedg33_ani(MMG5_pMesh  ,MMG5_pSol , int , MMG5_pTetra );

double (*_MMG5_lenedg)(MMG5_pMesh ,MMG5_pSol ,int, MMG5_pTetra );
double (*_MMG5_lenedgspl)(MMG5_pMesh ,MMG5_pSol ,int, MMG5_pTetra );
double (*_MMG5_caltet)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt);
double (*_MMG5_caltri)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
int    (*_MMG5_defsiz)(MMG5_pMesh ,MMG5_pSol );
int    (*_MMG5_gradsiz)(MMG5_pMesh ,MMG5_pSol );
int    (*_MMG5_intmet)(MMG5_pMesh,MMG5_pSol,int,char,int, double);
int    (*_MMG5_interp4bar)(MMG5_pMesh,MMG5_pSol,int,int,double *);
int    (*_MMG5_cavity)(MMG5_pMesh ,MMG5_pSol ,int ,int ,int *,int );
int    (*_MMG5_buckin)(MMG5_pMesh ,MMG5_pSol ,_MMG5_pBucket ,int );
int    (*_MMG3D_saveMeshinternal)(MMG5_pMesh mesh);

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Warn user that some tetrahedra of the mesh have been reoriented.
 *
 */
static inline
void _MMG5_warnOrientation(MMG5_pMesh mesh) {
  if ( mesh->xt ) {
    if ( mesh->xt != mesh->ne ) {
      fprintf(stdout,"  ## Warning: %d tetra on %d reoriented.\n",
              mesh->xt,mesh->ne);
      fprintf(stdout,"  Your mesh may be non-conform.\n");
    }
    else {
      fprintf(stdout,"  ## Warning: all tetra reoriented.\n");
    }
  }
  mesh->xt = 0;
}

/**
 * Set common pointer functions between mmgs and mmg3d to the matching mmg3d
 * functions.
 */
static inline
void _MMG3D_Set_commonFunc() {
  _MMG5_bezierCP          = _MMG5_mmg3dBezierCP;
  _MMG5_chkmsh            = _MMG5_mmg3dChkmsh;
#ifdef USE_SCOTCH
  _MMG5_renumbering       = _MMG5_mmg3dRenumbering;
#endif
}

#endif
