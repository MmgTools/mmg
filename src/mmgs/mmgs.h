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

#ifndef _MMGS_H
#define _MMGS_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <complex.h>

#include "chrono.h"
#include "eigenv.h"
#include "memory.h"

#define MS_VER   "2.0 a"
#define MS_REL   "Sept. 11, 2012"
#define MS_CPY   "Copyright (c) LJLL, 2009-"
#define MS_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

#define MS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MS_MIN(a,b) (((a) < (b)) ? (a) : (b))

/* numerical accuracy */
#define ATHIRD    0.333333333333333
#define A64TH     0.015625
#define A16TH     0.0625
#define A32TH     0.03125
#define ALPHAD    3.464101615137755   /* 6.0 / sqrt(3.0)  */
#define ANGEDG    0.707106781186548   /*0.573576436351046 */
#define SQR32     0.866025403784
#define LOPTL     1.4
#define LOPTS     0.71
#define LLONG     2.0
#define LSHRT     0.3
#define ANGLIM   -0.999999
#define EPSD      1.e-30
#define EPSD2     1.0e-200
#define EPS       1.e-06
#define LMAX      1024
#define BADKAL    2.e-2
#define NULKAL    1.e-4

#define NPMAX     500000
#define NTMAX    1000000
#define NGMAX     500000

#ifndef M_PI
#define M_PI            3.14159265358979323846   /* pi   */
#define M_PI_2          1.57079632679489661923   /* pi/2 */
#endif

/* edge tag */
#define  MS_NOTAG     (0)
#define  MS_REF       (1 << 0)        /* 1 edge reference     */
#define  MS_GEO       (1 << 1)        /* 2 geometric ridge    */
#define  MS_REQ       (1 << 2)        /* 4 required entity    */
#define  MS_NOM       (1 << 3)        /* 8 non manifold      */
/* point tags */
#define  MS_CRN       (1 << 5)        /* 32 corner           */
#define  MS_NUL       (1 << 6)        /* 64 vertex removed   */

#define MS_VOK(ppt)      (ppt && (ppt->tag < MS_NUL))
#define MS_EOK(pt)       (pt && (pt->v[0] > 0))
#define MS_EDG(tag)      ((tag & MS_GEO) || (tag & MS_REF))
#define MS_SIN(tag)      ((tag & MS_CRN) || (tag & MS_REQ) || (tag & MS_NOM))

#define MS_SET(flag,bit) ((flag) |= (1 << (bit)))
#define MS_CLR(flag,bit) ((flag) &= ~(1 << (bit)))
#define MS_GET(flag,bit) ((flag) & (1 << (bit)))

#define MS_Ver       1
#define MS_Tri       2

extern unsigned char inxt[3];
extern unsigned char iprv[3];


typedef struct {
    double  c[3],n[3];
    int     ref,ig,s,tmp;
    unsigned char tag,flag;
} Point;
typedef Point * pPoint;

typedef struct {
    double  b[10][3],n[6][3],t[6][3];
    pPoint  p[3];
} Bezier;
typedef Bezier * pBezier;

typedef struct {
    int   a,b,ref;
    unsigned char  tag;
} Edge;
typedef Edge *  pEdge;

typedef struct {
    int     v[3],edg[3],ref,cc,base;
    char    tag[3],flag;
} Tria;
typedef Tria * pTria;

typedef struct {
    double   n1[3],n2[3];
} Geom;
typedef Geom * pGeom;

/* specific parameters */
typedef struct {
    double   hmin,hmax,hausd;
    int      ref;
    char     elt;
} Par;
typedef Par * pPar;

typedef struct {
    double    dhd,hmin,hmax,hausd,hgrad,min[3],max[3],delta;
    int       ncc,npar,mem;
    char      imprim,ddebug,badkal,nreg,opt,mani;
    mytime    ctim[TIMEMAX];
    pPar      par;
} Info;

typedef struct {
    int       ver,dim,type,base;
    int       npi,nti,np,na,ng,nt,nc1,ngmax,npmax,ntmax,npnil,ntnil;
    int      *adja;
    char     *namein,*nameout;

    pPoint    point;
    pTria     tria;
    pEdge     edge;
    pGeom     geom;
} Mesh;
typedef Mesh  * pMesh;

typedef struct {
    int       dim,ver,np,npmax,size;
    double   *m;
    char     *namein,*nameout;
} Sol;
typedef Sol * pSol;

typedef struct {
    int   a,b,k,nxt;
} hedge;

typedef struct {
    int     siz,max,nxt;
    hedge  *item;
} Hash;


/* prototypes */
int  loadMesh(pMesh );
int  saveMesh(pMesh );
int  loadMet(pSol );
int  saveMet(pMesh ,pSol );
int  zaldy(pMesh mesh);
int  assignEdge(pMesh mesh);
int  scaleMesh(pMesh mesh,pSol met);
int  unscaleMesh(pMesh mesh,pSol met);
int  analys(pMesh mesh);
int  nortri(pMesh ,pTria ,double *n);
int  norpts(pPoint ,pPoint ,pPoint ,double *);
void outqua(pMesh ,pSol );
int  hashTria(pMesh );
int  curvpo(pMesh ,pSol );
int  mmgs1(pMesh ,pSol );
int  boulet(pMesh mesh,int start,int ip,int *list);
int  boulechknm(pMesh mesh,int start,int ip,int *list);
int  boulen(pMesh mesh,int start,int ip,double *nn);
int  boulep(pMesh mesh,int start,int ip,int *list);
int  boulec(pMesh mesh,int k,int i,double *tt);
int  bouler(pMesh mesh,int k,int i,int *list,int *ng,int *nr);
int  bouletrid(pMesh mesh,int start,int ip,int *il1,int *l1,int *il2,int *l2,int *ip0,int *ip1);
int  hashNew(Hash *hash,int hmax);
int  hashGet(Hash *hash,int a,int b);
int  hashEdge(Hash *hash,int a,int b,int k);
int  newPt(pMesh mesh,double c[3],double n[3]);
void delPt(pMesh mesh,int ip);
int  newElt(pMesh mesh);
void delElt(pMesh mesh,int iel);
int  chkedg(pMesh ,int );
int  bezierCP(pMesh ,int ,pBezier );
int  bezierInt(pBezier ,double *,double *,double *,double *);
void bezierEdge(pMesh mesh,int i0,int i1,double b0[3],double b1[3],char isrid,double v[3]);
int  split1(pMesh mesh,pSol met,int k,int i,int *vx);
int  split2(pMesh mesh,pSol met,int k,int *vx);
int  split3(pMesh mesh,pSol met,int k,int *vx);
int  split1b(pMesh mesh,int k,char i,int ip);
int  chkcol(pMesh mesh,pSol met,int k,char i,int *list,char typchk);
int  colver(pMesh mesh,int *list,int ilist);
int  colver3(pMesh mesh,int*list);
int  colver2(pMesh mesh,int *ilist);
int  swapar(pMesh mesh,int k,int i);
int  chkswp(pMesh mesh,pSol met,int k,int i,char typchk);
int  swpedg(pMesh mesh,pSol met,int *list,int ilist,char typchk);
char typelt(pPoint p[3],char *ia);
int  litswp(pMesh mesh,int k,char i,double kal);
int  litcol(pMesh mesh,int k,char i,double kal);
int  intmetsavedir(double *m,double *n,double *mr);
int  eigensym(double m[3], double lambda[2], double vp[2][2]);
int  sys33sym(double a[6],double b[3],double r[3]);
int  intmet33(pMesh mesh,pSol met,int np,int nq,int ip,double s);
int  intextmet(pMesh mesh,pSol met,int np,double me[6]);
int  invmatg(double m[9],double mi[9]);
int  rootDeg2(double complex a[3], double complex r[2]);
int  rootDeg3(double a[4],double complex r[3]);

int  buildridmet(pMesh mesh,pSol met,int np0,double ux,double uy,double uz,double mr[6]);
int  buildridmetnor(pMesh mesh,pSol met,int np0,double nt[3],double mr[6]);
int  buildridmetfic(pMesh mesh,double t[3],double n[3],double dtan,double dv,double m[6]);
int  rotmatrix(double n[3], double r[3][3]);
int  rmtr(double r[3][3],double m[6], double mr[6]);
int  paratmet(double c0[3],double n0[3],double m[6],double c1[3],double n1[3],double mt[6]);
int  intregmet(pMesh mesh,pSol met,int k,char i,double s,double mr[6]);
int  intridmet(pMesh mesh,pSol met,int k,char i,double s,double v[3],double mr[6]);
double surftri_iso(pMesh mesh,pSol met,int iel);
double surftri_ani(pMesh mesh,pSol met,int iel);
int  setref(pMesh,int,int,int);
int  delref(pMesh);
int  chkmet(pMesh,pSol);
int  chknor(pMesh);
void inqua(pMesh mesh,pSol met);

/* function pointers */
double calelt_ani(pMesh mesh,pSol met,int iel);
double calelt_iso(pMesh mesh,pSol met,int iel);
double caleltsig_ani(pMesh mesh,pSol met,int iel);
double caleltsig_iso(pMesh mesh,pSol met,int iel);
int    defsiz_iso(pMesh mesh,pSol met);
int    defsiz_ani(pMesh mesh,pSol met);
int    gradsiz_iso(pMesh mesh,pSol met);
int    gradsiz_ani(pMesh mesh,pSol met);
double lenedg_iso(pMesh mesh,pSol met,int ip1,int ip2,char isedg);
double lenedg_ani(pMesh mesh,pSol met,int ip1,int ip2,char isedg);
void   intmet_iso(pMesh mesh,pSol met,int k,char i,int ip,double s);
void   intmet_ani(pMesh mesh,pSol met,int k,char i,int ip,double s);
int    movridpt_iso(pMesh mesh,pSol met,int *list,int ilist);
int    movintpt_iso(pMesh mesh,pSol met,int *list,int ilist);
int    movridpt_ani(pMesh mesh,pSol met,int *list,int ilist);
int    movintpt_ani(pMesh mesh,pSol met,int *list,int ilist);

double (*calelt)(pMesh mesh,pSol met,int iel);
double (*lenedg)(pMesh ,pSol ,int ,int ,char );
int    (*defsiz)(pMesh mesh,pSol met);
int    (*gradsiz)(pMesh mesh,pSol met);
void   (*intmet)(pMesh mesh,pSol met,int k,char i,int ip,double s);
int    (*movridpt)(pMesh mesh,pSol met,int *list,int ilist);
int    (*movintpt)(pMesh mesh,pSol met,int *list,int ilist);


#endif
