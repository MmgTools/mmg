#ifndef _MMG3D_H
#define _MMG3D_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "chrono.h"
#include "memory.h"
#include "libmesh5.h"

#define MG_VER   "5.2c"
#define MG_REL   "Jul. 6, 2012"
#define MG_CPY   "Copyright (c) IMB-LJLL, 2004-"
#define MG_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

/* Macros */
#define MG_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MG_MIN(a,b) (((a) < (b)) ? (a) : (b))

#define MG_SMSGN(a,b)  (((double)(a)*(double)(b) > (0.0)) ? (1) : (0))

#define RETURN_AND_FREE(mesh,met,val)do         \
    {                                           \
      freeAll(mesh,met);                        \
      return(val);                              \
    }while(0)

/* numerical accuracy */
#define ALPHAD    20.7846096908265    //0.04811252243247      /* 12*sqrt(3) */
#define LLONG     2.0   // 1.414213562373
#define LSHRT     0.3  // 0.707106781186
#define LOPTL     1.3
#define LOPTS     0.6
#define ANGEDG    0.707106781186548   /*0.573576436351046 */
#define ANGLIM   -0.999999
#define SQR32     0.866025403784
#define ATHIRD    0.333333333333
#define EPSD      1.e-30
#define EPSD2     1.0e-200
#define EPS       1.e-06
#define LMAX      10240
#define BADKAL    0.2
#define NULKAL    1.e-30

#define NPMAX   1000000
#define NAMAX    200000
#define NTMAX   2000000
#define NEMAX   6000000

#ifndef M_PI
#define M_PI            3.14159265358979323846   /* pi   */
#define M_PI_2          1.57079632679489661923   /* pi/2 */
#endif

/* tags */
#define  MG_NOTAG     (0)
#define  MG_REF       (1 << 0)        /* 1 edge reference  */
#define  MG_GEO       (1 << 1)        /* 2 geometric ridge */
#define  MG_REQ       (1 << 2)        /* 4 required entity */
#define  MG_NOM       (1 << 3)        /* 8 non manifold    */
#define  MG_BDY       (1 << 4)        /* 16 boundary entity */
#define  MG_CRN       (1 << 5)        /* 32 corner         */
#define  MG_NUL       (1 << 6)        /* 64 vertex removed */

#define MG_PLUS    2
#define MG_MINUS   3
#define MG_ISO    -1024

#define MG_VOK(ppt)      (ppt && (ppt->tag < MG_NUL))
#define MG_EOK(pt)       (pt && (pt->v[0] > 0))
#define MG_EDG(tag)      ((tag & MG_GEO) || (tag & MG_REF))
#define MG_SIN(tag)      ((tag & MG_CRN) || (tag & MG_REQ))

#define MG_SET(flag,bit) ((flag) |= (1 << (bit)))
#define MG_CLR(flag,bit) ((flag) &= ~(1 << (bit)))
#define MG_GET(flag,bit) ((flag) & (1 << (bit)))

extern unsigned char inxt2[3];
extern unsigned char iprv2[3];
extern unsigned char inxt3[7];
extern unsigned char iprv3[7];
extern unsigned char idir[4][3];
extern          char idirinv[4][4];
extern unsigned char iarf[4][3];
extern unsigned char iarfinv[4][6];
extern unsigned char iare[6][2];
extern unsigned char ifar[6][2];
extern unsigned char isar[6][2];
extern unsigned char arpt[4][3];

typedef struct {
  double   c[3];
  int      ref,xp,tmp,flag;
  char     tag;
} Point;
typedef Point * pPoint;

typedef struct {
  double   n1[3],n2[3],t[3];
} xPoint;
typedef xPoint * pxPoint;

typedef struct {
  int      a,b,ref;
  char     tag;
} Edge;
typedef Edge * pEdge;

typedef struct {
  int      v[3],base,ref,edg[3];
  char     flag,tag[3];
} Tria;
typedef Tria * pTria;

typedef struct {
  int      v[4],ref,base,xt,flag;
  char     tag;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
  int      ref[4],edg[6];
  char     ftag[4],tag[6];
} xTetra;
typedef xTetra * pxTetra;

typedef struct {
  double  b[10][3],n[6][3],t[6][3];
  pPoint  p[3];
} Bezier;
typedef Bezier * pBezier;

typedef struct {
  double    dhd,hmin,hmax,hgrad,hausd,min[3],max[3],delta,ls;
  int       mem;
  char      imprim,ddebug,badkal,iso,fem;
  mytime    ctim[TIMEMAX];
} Info;

/* used to hash edges */
typedef struct {
  int   a,b,s,k,nxt;
} hedge;
typedef struct {
  int     siz,max,nxt;
  hedge  *item;
} Hash;

/* to store geometric edges */
typedef struct {
  int   a,b,ref,nxt;
  char  tag;
} hgeom;
typedef struct {
  int     siz,max,nxt;
  hgeom  *geom;
} HGeom;

typedef struct {
  int       ver,dim,type;
  int       npi,nai,nei,np,na,nt,ne,xp,xt,npmax,namax,ntmax,nemax,xpmax,xtmax,npnil,nenil,base;
  int      *adja,*adjt;
  char     *namein,*nameout;

  pPoint    point;
  pxPoint   xpoint;
  pTetra    tetra;
  pxTetra   xtetra;
  pTria     tria;
  pEdge     edge;
  HGeom     htab;
} Mesh;
typedef Mesh  * pMesh;

typedef struct {
  int       dim,ver,np,npmax,size,type;
  double   *m;
  char     *namein,*nameout;
} Sol;
typedef Sol * pSol;

/* prototypes */
void rotmatrix(double n[3],double r[3][3]);
double det3pt1vec(double c0[3],double c1[3],double c2[3],double v[3]);
double det4pt(double c0[3],double c1[3],double c2[3],double c3[3]);
double orvol(pPoint point,int *v);
int directsurfball(pMesh mesh, int ip, int *list, int ilist, double n[3]);

int  newPt(pMesh mesh,double c[3],char tag);
int  newElt(pMesh mesh);
void delElt(pMesh mesh,int iel);
void delPt(pMesh mesh,int ip);
int  zaldy(pMesh mesh);
char chkedg(pMesh mesh,pTria pt);
void tet2tri(pMesh mesh,int k,char ie,Tria *ptt);
int  bezierCP(pMesh mesh,Tria *pt,pBezier pb);
int  BezierTgt(double c1[3],double c2[3],double n1[3],double n2[3],double t1[3],double t2[3]);
double BezierGeod(double c1[3], double c2[3], double t1[3], double t2[3]);
int  bezierInt(pBezier pb,double uv[2],double o[3],double no[3],double to[3]);
int  BezierReg(pMesh mesh,int ip0, int ip1, double s, double v[3], double *o, double *no);
int  BezierRef(pMesh mesh,int ip0, int ip1, double s, double *o, double *no, double *to);
int  BezierEdge(pMesh mesh,int ip0, int ip1, double b0[3], double b1[3],char isrid, double v[3]);
int  BezierRidge(pMesh mesh,int ip0, int ip1, double s, double *o, double *no1, double *no2, double *to);
int  BezierNom(pMesh mesh,int ip0,int ip1,double s,double *o,double *no,double *to);
int  norface(pMesh mesh ,int k, int iface, double v[3]);
int  boulen(pMesh mesh,int start,int ip,double *nn);
int  bouler(pMesh mesh,int start,int ip,int *list);
int  boulenm(pMesh mesh, int start, int ip, int iface, double n[3],double t[3]);
int  boulec(pMesh mesh,int start,int ip,double *tt);
int  boulevolp(pMesh mesh, int start, int ip, int * list);
int  boulesurfvolp(pMesh mesh,int start,int ip,int iface,int *listv,int *ilistv,int *lists,int*ilists);
int  startedgsurfball(pMesh mesh,int nump,int numq,int *list,int ilist);
int  srcbdy(pMesh mesh,int start,int ia);
int  coquil(pMesh mesh, int start, int ia, int * list);
int  coquilface(pMesh mesh, int start, int ia, int * list, int * it1, int *it2);
int  gettag(pMesh,int,int,int *,int *);
int  settag(pMesh,int,int,int,int);
int  chkcol_int(pMesh ,pSol met,int,char,char,int *,char typchk);
int  chkcol_bdy(pMesh,int,char,char,int *);
int  chkmanicoll(pMesh mesh,int k,int iface,int iedg,int ndepmin,int ndepplus);
int  chkmani(pMesh mesh);
int  colver(pMesh,int *,int,char);
int  loadMesh(pMesh );
int  saveMesh(pMesh );
int  loadMet(pSol );
int  saveMet(pMesh mesh,pSol met);
int  analys(pMesh mesh);
int  hashTetra(pMesh mesh);
int  hashTria(pMesh mesh);
int  hashEdge(Hash *hash,int a,int b,int k);
int  hashGet(Hash *hash,int a,int b);
int  hashPop(Hash *hash,int a,int b);
int  hashNew(Hash *hash,int hsiz,int hmax);
int  hPop(HGeom *hash,int a,int b,int *ref,char *tag);
int  hTag(HGeom *hash,int a,int b,int ref,char tag);
int  hGet(HGeom *hash,int a,int b,int *ref,char *tag);
void hEdge(HGeom *hash,int a,int b,int ref,char tag);
void hNew(HGeom *hash,int hsiz,int hmax);
int  hGeom(pMesh mesh);
int  bdryTria(pMesh );
int  bdryIso(pMesh );
int  bdrySet(pMesh );
int  bdryPerm(pMesh );
int  chkmsh(pMesh,int,int);
int  chkfemtopo(pMesh mesh);
int  cntbdypt(pMesh mesh, int nump);
int  mmg3d1(pMesh ,pSol );
int  mmg3d2(pMesh ,pSol );
int  split1_sim(pMesh mesh,pSol met,int k,int vx[6]);
void split1(pMesh mesh,pSol met,int k,int vx[6]);
void split1b(pMesh mesh,int *list,int ret,int ip);
int  split2sf_sim(pMesh mesh,pSol met,int k,int vx[6]);
void split2sf(pMesh mesh,pSol met,int k,int vx[6]);
void split2(pMesh mesh,pSol met,int k,int vx[6]);
int  split3_sim(pMesh mesh,pSol met,int k,int vx[6]);
void split3(pMesh mesh,pSol met,int k,int vx[6]);
void split3cone(pMesh mesh,pSol met,int k,int vx[6]);
void split3op(pMesh mesh, pSol met, int k, int vx[6]);
void split4sf(pMesh mesh,pSol met,int k,int vx[6]);
void split4op(pMesh mesh,pSol met,int k,int vx[6]);
void split5(pMesh mesh,pSol met,int k,int vx[6]);
void split6(pMesh mesh,pSol met,int k,int vx[6]);
int  split4bar(pMesh mesh,pSol met,int k);
int  simbulgept(pMesh mesh, int *list, int ilist, double o[3]);
int  dichoto1b(pMesh mesh,int *list,int ret,double o[3],double ro[3]);
void nsort(int ,double *,char *);
int  nortri(pMesh mesh,pTria pt,double *n);
double orcal(pMesh mesh,int iel);
int  movintpt(pMesh mesh, int *list, int ilist);
int  movbdyregpt(pMesh mesh, int *listv, int ilistv, int *lists, int ilists);
int  movbdyrefpt(pMesh mesh, int *listv, int ilistv, int *lists, int ilists);
int  movbdynompt(pMesh mesh, int *listv, int ilistv, int *lists, int ilists);
int  movbdyridpt(pMesh mesh, int *listv, int ilistv, int *lists, int ilists);
double caltri(pMesh mesh,pTria ptt);
int  scaleMesh(pMesh mesh,pSol met);
int  unscaleMesh(pMesh mesh,pSol met);
int  chkswpbdy(pMesh mesh,int *list,int ilist,int it1,int it2);
int  swpbdy(pMesh mesh,pSol met,int *list,int ret,int it1);
void swpgen(pMesh mesh,pSol met,int nconf, int ilist, int *list);
int  chkswpgen(pMesh mesh, int start, int ia, int *ilist, int *list);
int  srcface(pMesh mesh,int n0,int n1,int n2);
int  bouleext(pMesh mesh, int start, int ip, int iface, int *listv, int *ilistv, int *lists, int*ilists);
int chkptonbdy(pMesh,int);
int norpts(pMesh,int,int,int,double *);
double orcal_poi(double a[3],double b[3],double c[3],double d[3]);
int trydisp(pMesh,double *,short);
int dichodisp(pMesh,double *);
int lapantilap(pMesh,double *);
int ppgdisp(pMesh,double *);
int denoisbdy(pMesh);
int eigensym(double m[3],double lambda[2],double vp[2][2]);
int sys33sym(double a[6], double b[3], double r[3]);
void outqua(pMesh mesh,pSol met);
int  badelt(pMesh mesh,pSol met);
int prilen(pMesh mesh,pSol met);
int DoSol(pMesh mesh,pSol met);
//int PutMetIn_h(pMesh mesh,pSol met);

int meancur(pMesh mesh,int np,double c[3],int ilist,int *list,double h[3]);
double surftri(pMesh,int,int);
double timestepMCF(pMesh,double);
int bdyMCF(pMesh);
double volint(pMesh);

/* pointers */
double caltet_ani(pMesh mesh,pSol met,int ia,int ib,int ic,int id);
double caltet_iso(pMesh mesh,pSol met,int ia,int ib,int ic,int id);
double lenedg_ani(pMesh ,pSol ,int ,int );
double lenedg_iso(pMesh ,pSol ,int ,int );
int    defsiz_iso(pMesh,pSol );
int    defsiz_ani(pMesh ,pSol );
int    gradsiz_iso(pMesh ,pSol );
int    gradsiz_ani(pMesh ,pSol );

double (*caltet)(pMesh mesh,pSol met,int ia,int ib,int ic,int id);
double (*lenedg)(pMesh ,pSol ,int ,int );
int    (*defsiz)(pMesh ,pSol );
int    (*gradsiz)(pMesh ,pSol );

#endif
