#include "mmgexterns.h"
#include "mmgs.h"

extern int    (*movintpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
extern int    (*MMGS_defsiz)(MMG5_pMesh mesh,MMG5_pSol met);
extern double (*MMG5_calelt)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
extern int    (*MMGS_gradsiz)(MMG5_pMesh mesh,MMG5_pSol met);
extern int    (*MMGS_gradsizreq)(MMG5_pMesh mesh,MMG5_pSol met);
extern int    (*intmet)(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
extern int    (*movridpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
