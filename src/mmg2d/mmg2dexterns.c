#include "mmgexterns.h"
#include "mmg2d.h"

int    (*MMG2D_defsiz)(MMG5_pMesh ,MMG5_pSol )=NULL;
int    (*MMG2D_intmet)(MMG5_pMesh ,MMG5_pSol ,int ,int8_t ,int ,double )=NULL;
double (*MMG2D_lencurv)(MMG5_pMesh ,MMG5_pSol ,int ,int )=NULL;
int    (*MMG2D_gradsizreq)(MMG5_pMesh ,MMG5_pSol )=NULL;
double (*MMG2D_caltri)(MMG5_pMesh ,MMG5_pSol ,MMG5_pTria )=NULL;
int    (*MMG2D_gradsiz)(MMG5_pMesh ,MMG5_pSol )=NULL;

