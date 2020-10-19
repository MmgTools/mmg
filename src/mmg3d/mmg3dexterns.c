#include "mmgexterns.h"
#include "mmg3d.h"

double (*MMG3D_lenedgCoor)(double *ca,double *cb,double *sa,double *sb)=NULL;
double (*MMG5_lenedg)(MMG5_pMesh ,MMG5_pSol ,int, MMG5_pTetra )=NULL;
double (*MMG5_lenedgspl)(MMG5_pMesh ,MMG5_pSol ,int, MMG5_pTetra )=NULL;
double (*MMG5_caltet)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt)=NULL;
double (*MMG5_caltri)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt)=NULL;
int    (*MMG3D_defsiz)(MMG5_pMesh ,MMG5_pSol )=NULL;
int    (*MMG3D_gradsiz)(MMG5_pMesh ,MMG5_pSol )=NULL;
int    (*MMG3D_gradsizreq)(MMG5_pMesh ,MMG5_pSol )=NULL;
int    (*MMG5_intmet)(MMG5_pMesh,MMG5_pSol,int,int8_t,int, double)=NULL;
int    (*MMG5_interp4bar)(MMG5_pMesh,MMG5_pSol,int,int,double *)=NULL;
int    (*MMG5_movintpt)(MMG5_pMesh ,MMG5_pSol, MMG3D_pPROctree ,int *, int , int )=NULL;
int    (*MMG5_movbdyregpt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int*, int, int*, int, int ,int)=NULL;
int    (*MMG5_movbdyrefpt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int*, int, int*, int ,int)=NULL;
int    (*MMG5_movbdynompt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int*, int, int*, int ,int)=NULL;
int    (*MMG5_movbdyridpt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int*, int, int*, int ,int)=NULL;
int    (*MMG5_cavity)(MMG5_pMesh ,MMG5_pSol ,int ,int ,int *,int ,double)=NULL;
int    (*MMG3D_PROctreein)(MMG5_pMesh ,MMG5_pSol ,MMG3D_pPROctree ,int,double )=NULL;

