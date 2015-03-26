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

#include <complex.h>

//#include "memory.h"
#include "libmmgs.h"

/* numerical accuracy */
#define A64TH     0.015625
#define A16TH     0.0625
#define A32TH     0.03125
#define ALPHAD    3.464101615137755   /* 6.0 / sqrt(3.0)  */

#define LOPTL     1.4
#define LOPTS     0.71
#define LLONG     2.0
#define LSHRT     0.3

#define LMAX      1024
#define BADKAL    2.e-2
#define NULKAL    1.e-4

#define _MMG5_NPMAX     500000
#define _MMG5_NTMAX    1000000
#define _MMG5_XPMAX     500000


#define MS_SIN(tag)      ((tag & MG_CRN) || (tag & MG_REQ) || (tag & MG_NOM))

#define MS_Ver       1
#define MS_Tri       2

/* prototypes */
int  loadMesh(MMG5_pMesh );
int  saveMesh(MMG5_pMesh );
int  MMG5_loadMet(MMG5_pMesh,MMG5_pSol );
int  saveMet(MMG5_pMesh ,MMG5_pSol );
int  zaldy(MMG5_pMesh mesh);
int  assignEdge(MMG5_pMesh mesh);
int  analys(MMG5_pMesh mesh);
int  nortri(MMG5_pMesh ,MMG5_pTria ,double *n);
int  norpts(MMG5_pPoint ,MMG5_pPoint ,MMG5_pPoint ,double *);
void outqua(MMG5_pMesh ,MMG5_pSol );
int  hashTria(MMG5_pMesh );
int  curvpo(MMG5_pMesh ,MMG5_pSol );
int  mmgs1(MMG5_pMesh ,MMG5_pSol );
int  boulet(MMG5_pMesh mesh,int start,int ip,int *list);
int  boulechknm(MMG5_pMesh mesh,int start,int ip,int *list);
int  boulen(MMG5_pMesh mesh,int start,int ip,double *nn);
int  boulep(MMG5_pMesh mesh,int start,int ip,int *list);
int  boulec(MMG5_pMesh mesh,int k,int i,double *tt);
int  bouler(MMG5_pMesh mesh,int k,int i,int *list,int *xp,int *nr);
int  bouletrid(MMG5_pMesh mesh,int start,int ip,int *il1,int *l1,int *il2,int *l2,int *ip0,int *ip1);
int  newPt(MMG5_pMesh mesh,double c[3],double n[3]);
void delPt(MMG5_pMesh mesh,int ip);
int  newElt(MMG5_pMesh mesh);
void delElt(MMG5_pMesh mesh,int iel);
int  chkedg(MMG5_pMesh ,int );
int  bezierCP(MMG5_pMesh ,int ,_MMG5_pBezier );
int  bezierInt(_MMG5_pBezier ,double *,double *,double *,double *);
void bezierEdge(MMG5_pMesh mesh,int i0,int i1,double b0[3],double b1[3],char isrid,double v[3]);
int  split1(MMG5_pMesh mesh,MMG5_pSol met,int k,int i,int *vx);
int  split2(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  split3(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  split1b(MMG5_pMesh mesh,int k,char i,int ip);
int  chkcol(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int *list,char typchk);
int  colver(MMG5_pMesh mesh,int *list,int ilist);
int  colver3(MMG5_pMesh mesh,int*list);
int  colver2(MMG5_pMesh mesh,int *ilist);
int  swapar(MMG5_pMesh mesh,int k,int i);
int  chkswp(MMG5_pMesh mesh,MMG5_pSol met,int k,int i,char typchk);
int  swpedg(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist,char typchk);
char typelt(MMG5_pPoint p[3],char *ia);
int  litswp(MMG5_pMesh mesh,int k,char i,double kal);
int  litcol(MMG5_pMesh mesh,int k,char i,double kal);
int  intmetsavedir(MMG5_pMesh mesh, double *m,double *n,double *mr);
int  sys33sym(double a[6],double b[3],double r[3]);
int  intmet33(MMG5_pMesh mesh,MMG5_pSol met,int np,int nq,int ip,double s);
int  intextmet(MMG5_pMesh mesh,MMG5_pSol met,int np,double me[6]);
int  invmatg(double m[9],double mi[9]);
int  rootDeg2(double complex a[3], double complex r[2]);
int  rootDeg3(double a[4],double complex r[3]);

int  buildridmetnor(MMG5_pMesh mesh,MMG5_pSol met,int np0,double nt[3],double mr[6]);
int  buildridmetfic(MMG5_pMesh mesh,double t[3],double n[3],double dtan,double dv,double m[6]);
int  rotmatrix(double n[3], double r[3][3]);
int  rmtr(double r[3][3],double m[6], double mr[6]);
int  paratmet(double c0[3],double n0[3],double m[6],double c1[3],double n1[3],double mt[6]);
int  intregmet(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,double s,double mr[6]);
int  intridmet(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,double s,double v[3],double mr[6]);
double surftri_iso(MMG5_pMesh mesh,MMG5_pSol met,int iel);
double surftri_ani(MMG5_pMesh mesh,MMG5_pSol met,int iel);
int  setref(MMG5_pMesh,int,int,int);
int  delref(MMG5_pMesh);
int  chkmet(MMG5_pMesh,MMG5_pSol);
int  chknor(MMG5_pMesh);
void inqua(MMG5_pMesh mesh,MMG5_pSol met);
long long _MMG5_memSize(void);
void _MMG5_memOption(MMG5_pMesh mesh);

/* function pointers */
/* init structures */
void  _MMG5_Init_parameters(MMG5_pMesh mesh);
/* init file names */
int  _MMG5_Set_outputMeshName(MMG5_pMesh mesh, char* meshout);
/* iso/aniso computations */
double calelt_ani(MMG5_pMesh mesh,MMG5_pSol met,int iel);
double calelt_iso(MMG5_pMesh mesh,MMG5_pSol met,int iel);
double caleltsig_ani(MMG5_pMesh mesh,MMG5_pSol met,int iel);
double caleltsig_iso(MMG5_pMesh mesh,MMG5_pSol met,int iel);
int    defsiz_iso(MMG5_pMesh mesh,MMG5_pSol met);
int    defsiz_ani(MMG5_pMesh mesh,MMG5_pSol met);
int    gradsiz_iso(MMG5_pMesh mesh,MMG5_pSol met);
int    gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met);
void   intmet_iso(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
void   intmet_ani(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
int    movridpt_iso(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    movintpt_iso(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    movridpt_ani(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    movintpt_ani(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);

double (*calelt)(MMG5_pMesh mesh,MMG5_pSol met,int iel);
int    (*defsiz)(MMG5_pMesh mesh,MMG5_pSol met);
int    (*gradsiz)(MMG5_pMesh mesh,MMG5_pSol met);
void   (*intmet)(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
int    (*movridpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    (*movintpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);

void _MMG5_Set_APIFunc();

#endif
