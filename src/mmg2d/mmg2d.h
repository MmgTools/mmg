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
#ifndef _MMG2D_H
#define _MMG2D_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <float.h>

#include "libmmg2d.h"
#include "mmgcommon.h"

#ifdef __cplusplus
extern "C" {
#endif

/* constantes */

#define M_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define M_MIN(a,b) (((a) < (b)) ? (a) : (b))

#define M_LAMBDA  0.34
#define M_MU      0.33

#define _MMG2_EPSD   1.e-10 //e-20??
#define _MMG2D_EPSA   1.e-12

#define _MMG2D_PRECI  1.
#define _MMG2D_SIZE   0.75
#define _MMG2D_ALPHA  0.28867513459
#define _MMG2D_ALPHAD 3.464101615137755   /* 6.0 / sqrt(3.0)  */
#define MMG2_LONMAX 1024
#define _MMG2D_BADKAL    0.2
#define _MMG2_NULKAL    1.e-6
#define _MMG2_ANGCORN   -1.e-6
#define _MMG2_SHORTMAX     0x7fff

#define MMG2_LLONG  2.0
#define MMG2_LSHRT  0.3
#define MMG2_LOPTL      1.4
#define MMG2_LOPTS     0.71

#define _MMG2D_NPMAX   50000
#define _MMG2D_NEDMAX  100000
#define _MMG2D_NEMAX   100000

/** Free allocated pointers of mesh and sol structure and return value val */
#define _MMG2D_RETURN_AND_FREE(mesh,met,disp,val)do                 \
  {                                                                 \
    if ( !MMG2D_Free_all(MMG5_ARG_start,                            \
                         MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met, \
                         MMG5_ARG_end) ) {                          \
      return MMG5_LOWFAILURE;                                       \
    }                                                               \
    return(val);                                                    \
  }while(0)

/**
 * \param sigid signal number.
 *
 * Signal handling: specify error messages depending from catched signal.
 *
 */
static inline
void _MMG2_excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
  case SIGABRT:
    fprintf(stdout,"  Abnormal stop\n"); break;
  case SIGFPE:
    fprintf(stdout,"  Floating-point exception\n"); break;
  case SIGILL:
    fprintf(stdout,"  Illegal instruction\n"); break;
  case SIGSEGV:
    fprintf(stdout,"  Segmentation fault\n"); break;
  case SIGTERM:
  case SIGINT:
    fprintf(stdout,"  Program killed\n"); break;
  }
  exit(EXIT_FAILURE);
}

typedef struct squeue {
  int    *stack,cur;
} Queue;
typedef Queue * pQueue;

typedef struct {
  int     size;
  int    *head;
  int    *link;
} Bucket;
typedef Bucket * pBucket;

typedef struct {
  int      min,max,iel,nxt;
} Hedge;

typedef struct {
  int      size,nxtmax,hnxt;
  Hedge    *item;
} HashTable;
typedef HashTable * pHashTable;

static const int MMG2_iare[3][2] = {{1,2},{2,0},{0,1}};
static const int MMG2_iopp[3][2] = {{1,2},{0,2},{0,1}};
static const unsigned int MMG2_idir[5] = {0,1,2,0,1};
static const unsigned int MMG2_inxt[5] = {1,2,0,1,2};


/** Reallocation of point table and sol table and creation
    of point ip with coordinates o and tag tag*/
#define _MMG2D_POINT_REALLOC(mesh,sol,ip,wantedGap,law,o,tag,retval ) do \
  {                                                                     \
    int klink;                                                          \
                                                                        \
    _MMG5_TAB_RECALLOC(mesh,mesh->point,mesh->npmax,wantedGap,MMG5_Point, \
                       "larger point table",law,retval);                \
                                                                        \
    mesh->npnil = mesh->np+1;                                           \
    for (klink=mesh->npnil; klink<mesh->npmax-1; klink++)               \
      mesh->point[klink].tmp  = klink+1;                                \
                                                                        \
    /* solution */                                                      \
    if ( sol->m ) {                                                     \
      _MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax-sol->npmax))*sizeof(double), \
                    "larger solution",law);                             \
      _MMG5_SAFE_REALLOC(sol->m,sol->size*(mesh->npmax+1),              \
                         double,"larger solution",retval);              \
    }                                                                   \
    sol->npmax = mesh->npmax;                                           \
                                                                        \
    /* We try again to add the point */                                 \
    ip = _MMG2D_newPt(mesh,o,tag);                                      \
    if ( !ip ) {law;}                                                   \
  }while(0)

/** Reallocation of tria table and creation
    of tria jel */
#define _MMG2D_TRIA_REALLOC(mesh,jel,wantedGap,law,retval ) do          \
  {                                                                     \
   int klink,oldSiz;                                                    \
                                                                        \
   oldSiz = mesh->ntmax;                                                \
   _MMG5_TAB_RECALLOC(mesh,mesh->tria,mesh->ntmax,wantedGap,MMG5_Tria,  \
                      "larger tria table",law,retval);                  \
                                                                        \
   mesh->nenil = mesh->nt+1;                                            \
   for (klink=mesh->nenil; klink<mesh->ntmax-1; klink++)                \
     mesh->tria[klink].v[2]  = klink+1;                                 \
                                                                        \
   if ( mesh->adja ) {                                                  \
     /* adja table */                                                   \
     _MMG5_ADD_MEM(mesh,3*(mesh->ntmax-oldSiz)*sizeof(int),             \
                   "larger adja table",law);                            \
     _MMG5_SAFE_RECALLOC(mesh->adja,3*oldSiz+5,3*mesh->ntmax+5,int      \
                         ,"larger adja table",retval);                  \
   }                                                                    \
                                                                        \
   /* We try again to add the point */                                  \
   jel = _MMG2D_newElt(mesh);                                           \
   if ( !jel ) {law;}                                                   \
   }while(0)

/** Reallocation of edge table and creation
    of edge jel */
#define _MMG2D_EDGE_REALLOC(mesh,jel,wantedGap,law,retval ) do          \
  {                                                                     \
    int klink;                                                          \
                                                                        \
    _MMG5_TAB_RECALLOC(mesh,mesh->edge,mesh->namax,wantedGap,MMG5_Edge, \
                       "larger edge table",law,retval);                 \
                                                                        \
    mesh->nanil = mesh->na+1;                                           \
    for (klink=mesh->nanil; klink<mesh->namax-1; klink++)               \
      mesh->edge[klink].b  = klink+1;                                   \
                                                                        \
                                                                        \
    /* We try again to add the point */                                 \
    jel = _MMG5_newEdge(mesh);                                          \
    if ( !jel ) {law;}                                                  \
  }while(0)


/* Prototypes */
/*zaldy*/
int _MMG2D_newPt(MMG5_pMesh mesh,double c[2],int16_t tag);
void _MMG2D_delPt(MMG5_pMesh mesh,int ip) ;
int _MMG5_newEdge(MMG5_pMesh mesh);
void _MMG5_delEdge(MMG5_pMesh mesh,int iel);
int _MMG2D_newElt(MMG5_pMesh mesh);
int  _MMG2D_delElt(MMG5_pMesh mesh,int iel);
int _MMG5_getnElt(MMG5_pMesh mesh,int n);
int MMG2D_zaldy(MMG5_pMesh mesh);
long long _MMG5_memSize(void);
int _MMG2D_memOption(MMG5_pMesh mesh);

int MMG2_scaleMesh(MMG5_pMesh ,MMG5_pSol );
int MMG2_unscaleMesh(MMG5_pMesh ,MMG5_pSol );
int MMG2_pack(MMG5_pMesh ,MMG5_pSol );
void MMG2_outqua(MMG5_pMesh ,MMG5_pSol );
//int MMG2_mmg2d0(MMG5_pMesh ,MMG5_pSol );
int MMG2_mmg2d1(MMG5_pMesh ,MMG5_pSol );
//int MMG2_split(MMG5_pMesh ,MMG5_pSol ,int ,int ,int,double );
//int MMG2_splitbdry(MMG5_pMesh ,MMG5_pSol ,int ,int ,int,double*);
//int MMG2_colpoi(MMG5_pMesh ,MMG5_pSol , int ,int ,int ,int ,double );
//int MMG2_colpoibdry(MMG5_pMesh ,MMG5_pSol , int ,int ,int ,int ,double );

int  _MMG2D_Init_mesh_var( va_list argptr );
int  _MMG2D_Free_all_var( va_list argptr );
int  _MMG2D_Free_structures_var( va_list argptr );
int  _MMG2D_Free_names_var( va_list argptr );

void MMG2D_solTruncatureForOptim(MMG5_pMesh mesh, MMG5_pSol met);

int MMG2_mmg2d2(MMG5_pMesh , MMG5_pSol);
int MMG2_mmg2d6(MMG5_pMesh ,MMG5_pSol );
int MMG2_mmg2d9(MMG5_pMesh ,MMG5_pSol ,MMG5_pSol );
//int MMG2_cendel(MMG5_pMesh ,MMG5_pSol ,double ,int );
int _MMG2_swapdelone(MMG5_pMesh ,MMG5_pSol ,int ,char ,double ,int *);
int _MMG5_mmg2dChkmsh(MMG5_pMesh , int, int );
int MMG2_boulep(MMG5_pMesh , int , int , int * );
//int MMG2_markBdry(MMG5_pMesh );
int MMG2_prilen(MMG5_pMesh ,MMG5_pSol );

int MMG2_coorbary(MMG5_pMesh ,MMG5_pTria ,double c[2],double* ,double* ,double* );
int MMG2_isInTriangle(MMG5_pMesh ,int,double c[2]);
int MMG2_cutEdge(MMG5_pMesh ,MMG5_pTria ,MMG5_pPoint ,MMG5_pPoint );
int MMG2_cutEdgeTriangle(MMG5_pMesh ,int ,int ,int );
int MMG2_findTria(MMG5_pMesh ,int );
//int MMG2_findpos(MMG5_pMesh ,MMG5_pTria ,int ,int ,int ,int ,int );
int MMG2_locateEdge(MMG5_pMesh ,int ,int ,int* ,int* ) ;
int MMG2_bdryenforcement(MMG5_pMesh ,MMG5_pSol);
int MMG2_settagtriangles(MMG5_pMesh ,MMG5_pSol );
int MMG2_findtrianglestate(MMG5_pMesh ,int ,int ,int ,int ,int ,int );

pQueue MMG2_kiuini(MMG5_pMesh mesh,int nbel,double declic,int base);
void MMG2_kiufree(pQueue q);
int MMG2_kiudel(pQueue q,int iel);
int MMG2_kiuput(pQueue q,int iel);
int MMG2_kiupop(pQueue q);

int MMG2_hashEdge(pHashTable edgeTable,int iel,int ia, int ib);
//int MMG2_hashel(MMG5_pMesh mesh);
int MMG2_hashNew(HashTable *hash,int hsize,int hmax);
int MMG2_baseBdry(MMG5_pMesh mesh);

int MMG2_invmat(double *m,double *minv);
int simred(double *m1,double *m2,double *m);

//int MMG2_evalgeom(MMG5_pMesh mesh);

int _MMG2_cavity(MMG5_pMesh ,MMG5_pSol ,int ,int *);
int _MMG2_delone(MMG5_pMesh ,MMG5_pSol ,int ,int *,int );
int _MMG2_cenrad_iso(MMG5_pMesh ,double *,double *,double *);

/* Adds Charles */
int _MMG2_getIniRef(MMG5_pMesh ,int );
int _MMG2_isSplit(MMG5_pMesh ,int ,int *,int *);
int MMG2_parsop(MMG5_pMesh ,MMG5_pSol );
int _MMG2_ismaniball(MMG5_pMesh , MMG5_pSol , int , char );
int _MMG2_snapval(MMG5_pMesh ,MMG5_pSol ,double *);
int _MMG2_chkmanimesh(MMG5_pMesh );
int MMG2_hashTria(MMG5_pMesh );
int _MMG2_resetRef(MMG5_pMesh );
int _MMG2_cuttri_ls(MMG5_pMesh ,MMG5_pSol );
int _MMG2_setref_ls(MMG5_pMesh ,MMG5_pSol );
int _MMG2_split1_sim(MMG5_pMesh ,MMG5_pSol ,int ,int vx[3]);
int _MMG2_split2_sim(MMG5_pMesh ,MMG5_pSol ,int ,int vx[3]);
int _MMG2_split3_sim(MMG5_pMesh ,MMG5_pSol ,int ,int vx[3]);
int _MMG2_split1(MMG5_pMesh ,MMG5_pSol ,int ,int vx[3]);
int _MMG2_split2(MMG5_pMesh ,MMG5_pSol ,int ,int vx[3]);
int _MMG2_split3(MMG5_pMesh ,MMG5_pSol ,int ,int vx[3]);
int  MMG2_assignEdge(MMG5_pMesh );
int  MMG2_bdryEdge(MMG5_pMesh );
int _MMG2_setadj(MMG5_pMesh );
int _MMG2_singul(MMG5_pMesh );
int _MMG2_analys(MMG5_pMesh );
int _MMG2_norver(MMG5_pMesh );
int _MMG2_regnor(MMG5_pMesh );
int _MMG2_boulen(MMG5_pMesh , int ,char ,int *,int *,double *);
int MMG2_mmg2d1n(MMG5_pMesh ,MMG5_pSol );
int _MMG2_anatri(MMG5_pMesh ,MMG5_pSol ,char );
int _MMG2_adptri(MMG5_pMesh ,MMG5_pSol );
int _MMG2_defsiz_iso(MMG5_pMesh ,MMG5_pSol );
int _MMG2_defsiz_ani(MMG5_pMesh ,MMG5_pSol );
int _MMG2_defmetbdy_2d(MMG5_pMesh ,MMG5_pSol ,int ,char );
int _MMG2_defaultmet_2d(MMG5_pMesh ,MMG5_pSol ,int ,char );
int _MMG2_grad2met_ani(MMG5_pMesh ,MMG5_pSol ,double *,double *,double );
int _MMG2_gradsiz_iso(MMG5_pMesh ,MMG5_pSol );
int _MMG2_gradsiz_ani(MMG5_pMesh ,MMG5_pSol );
int _MMG2_anaelt(MMG5_pMesh ,MMG5_pSol ,int );
int _MMG2_colelt(MMG5_pMesh ,MMG5_pSol ,int );
int _MMG2_swpmsh(MMG5_pMesh ,MMG5_pSol ,int );
double _MMG2_lencurv_iso(MMG5_pMesh ,MMG5_pSol ,int ,int );
double _MMG2_lencurv_ani(MMG5_pMesh ,MMG5_pSol ,int ,int );
int _MMG2_chkedg(MMG5_pMesh ,int );
int _MMG2_bezierCurv(MMG5_pMesh ,int ,char ,double ,double *,double *);
int _MMG2_dichoto(MMG5_pMesh ,MMG5_pSol ,int ,int *);
double _MMG2_quickcal(MMG5_pMesh , MMG5_pTria );
int _MMG2_chkcol(MMG5_pMesh,MMG5_pSol,int,char,int *,char);
int _MMG2_colver(MMG5_pMesh,int,int*);
int _MMG2_colver3(MMG5_pMesh,int*);
int _MMG2_colver2(MMG5_pMesh,int*);
int _MMG2_boulet(MMG5_pMesh,int,char,int*);
int _MMG2_bouleendp(MMG5_pMesh,int,char,int*,int*);
int _MMG2_savemesh_db(MMG5_pMesh ,char* ,char );
int _MMG2_savemet_db(MMG5_pMesh ,MMG5_pSol ,char* ,char );
int _MMG2_chkswp(MMG5_pMesh , MMG5_pSol ,int ,char ,char );
int _MMG2_swapar(MMG5_pMesh ,int ,char );
int _MMG5_interpmet22(MMG5_pMesh ,double *,double *,double ,double *);
int _MMG2_intmet_iso(MMG5_pMesh ,MMG5_pSol ,int ,char ,int ,double );
int _MMG2_intmet_ani(MMG5_pMesh ,MMG5_pSol ,int ,char ,int ,double );
int _MMG2_adpspl(MMG5_pMesh ,MMG5_pSol );
int _MMG2_adpcol(MMG5_pMesh ,MMG5_pSol );
int _MMG2_movtri(MMG5_pMesh ,MMG5_pSol ,int ,char );
int _MMG2_chkspl(MMG5_pMesh ,MMG5_pSol ,int ,char );
int _MMG2_split1b(MMG5_pMesh ,int ,char ,int );
int _MMG2_movedgpt(MMG5_pMesh ,MMG5_pSol ,int ,int *,char );
int _MMG2_movintpt(MMG5_pMesh ,MMG5_pSol ,int ,int *,char );
int _MMG2_movintpt_ani(MMG5_pMesh ,MMG5_pSol ,int ,int *,char );
int _MMG2_chkmsh(MMG5_pMesh );
int _MMG2_chkor(MMG5_pMesh );
int _MMG2_savenor_db(MMG5_pMesh ,char *,char );
int _MMG2_savedisp_db(MMG5_pMesh mesh,MMG5_pSol ,char *,char );
int _MMG2_velextLS(MMG5_pMesh ,MMG5_pSol );

/* useful functions to debug */
int  _MMG2D_indElt(MMG5_pMesh mesh,int kel);
int  _MMG2D_indPt(MMG5_pMesh mesh,int kp);


/* functions pointers */
double long_ani(double *ca,double *cb,double *ma,double *mb);
double long_iso(double *ca,double *cb,double *ma,double *mb);
double _MMG2_caltri_ani(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pTria );
double _MMG2_caltri_iso(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pTria );
int    optlen_ani(MMG5_pMesh mesh,MMG5_pSol sol,double declic,int base);
int    optlen_iso(MMG5_pMesh mesh,MMG5_pSol sol,double declic,int base);
int    optlen_iso_bar(MMG5_pMesh mesh,MMG5_pSol sol,double declic,int base);
int    interp_ani(double *,double *,double * ,double );
int    interp_iso(double *,double *,double * ,double );
int    lissmet_iso(MMG5_pMesh mesh,MMG5_pSol sol);
int    lissmet_ani(MMG5_pMesh mesh,MMG5_pSol sol);

double (*MMG2D_lencurv)(MMG5_pMesh ,MMG5_pSol ,int ,int );
double (*MMG2D_caltri)(MMG5_pMesh ,MMG5_pSol ,MMG5_pTria );
int    (*MMG2_optlen)(MMG5_pMesh ,MMG5_pSol ,double ,int );
int    (*MMG2D_intmet)(MMG5_pMesh ,MMG5_pSol ,int ,char ,int ,double );
int    (*MMG2D_gradsiz)(MMG5_pMesh ,MMG5_pSol );
int    (*MMG2D_defsiz)(MMG5_pMesh ,MMG5_pSol );

/* init structures */
void  _MMG2_Init_parameters(MMG5_pMesh mesh);

/**
 * Set common pointer functions between mmgs and mmg2d to the matching mmg2d
 * functions.
 */
static inline
void _MMG2D_Set_commonFunc() {
  _MMG5_chkmsh            = _MMG5_mmg2dChkmsh;
  return;
}

#ifdef __cplusplus
}
#endif

#endif
