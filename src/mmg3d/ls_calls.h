#ifndef __LS_CALLS_H
#define __LS_CALLS_H

enum {Dirichlet=1, Load};
enum {P0=0, P1, P2};
enum {LS_ver=1,LS_edg=2,LS_tri=4,LS_tet=8};


/* data structure */
typedef struct _LSst LSst;

/* prototypes */
LSst *LS_init(int dim, int ver, char typ,char mfree);
int   LS_stop(LSst *lsst);

int   LS_mesh(LSst *lsst,int np,int na,int nt,int ne);
int   LS_iniSol(LSst *lsst,double *u);
int   LS_CreaSol(LSst *lsst);
int   LS_addSol(LSst *lsst,int ip,double *u);
int   LS_addVer(LSst *lsst,int idx,double *c,int ref);
int   LS_addEdg(LSst *lsst,int idx,int *v,int ref);
int   LS_addTri(LSst *lsst,int idx,int *v,int ref);
int   LS_addTet(LSst *lsst,int idx,int *v,int ref);
void  LS_headMesh(LSst *lsst,int *np,int *na,int *nt,int *ne);

void  LS_setPar(LSst *lsst,int imp,int deb,int zip);
int   LS_setBC(LSst *lsst,int typ,int ref,char att,int elt,double *u);
void  LS_setGra(LSst *lsst, double *gr);
int   LS_setLame(LSst *lsst,int ref,double lambda,double mu);
int   LS_elastic(LSst *lsst);

double *LS_getSol(LSst *lsst);



#endif
