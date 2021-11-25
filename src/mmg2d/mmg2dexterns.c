#include "mmgexterns.h"
#include "mmg2d.h"

MMG_int    (*MMG2D_defsiz)(MMG5_pMesh ,MMG5_pSol )=NULL;
MMG_int    (*MMG2D_intmet)(MMG5_pMesh ,MMG5_pSol ,MMG_int ,int8_t ,MMG_int ,double )=NULL;
double (*MMG2D_lencurv)(MMG5_pMesh ,MMG5_pSol ,MMG_int ,MMG_int )=NULL;
MMG_int    (*MMG2D_gradsizreq)(MMG5_pMesh ,MMG5_pSol )=NULL;
double (*MMG2D_caltri)(MMG5_pMesh ,MMG5_pSol ,MMG5_pTria )=NULL;
MMG_int    (*MMG2D_gradsiz)(MMG5_pMesh ,MMG5_pSol )=NULL;

