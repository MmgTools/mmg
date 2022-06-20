#include "mmgexterns.h"
#include "libmmgs_private.h"

double (*MMG5_calelt)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt)=NULL;
int    (*MMGS_doSol)(MMG5_pMesh mesh,MMG5_pSol met)=NULL;
int    (*MMGS_defsiz)(MMG5_pMesh mesh,MMG5_pSol met)=NULL;
int    (*MMGS_gradsiz)(MMG5_pMesh mesh,MMG5_pSol met)=NULL;
int    (*MMGS_gradsizreq)(MMG5_pMesh mesh,MMG5_pSol met)=NULL;
int    (*intmet)(MMG5_pMesh mesh,MMG5_pSol met,int k,int8_t i,int ip,double s)=NULL;
int    (*movintpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist)=NULL;
int    (*movridpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist)=NULL;
