#include "mmgexterns.h"
#include "mmg2d.h"

extern int    (*MMG2D_defsiz)(MMG5_pMesh ,MMG5_pSol );
extern int    (*MMG2D_intmet)(MMG5_pMesh ,MMG5_pSol ,MMG_int ,int8_t ,MMG_int ,double );
extern double (*MMG2D_lencurv)(MMG5_pMesh ,MMG5_pSol ,MMG_int ,MMG_int );
extern int    (*MMG2D_gradsizreq)(MMG5_pMesh ,MMG5_pSol );
extern double (*MMG2D_caltri)(MMG5_pMesh ,MMG5_pSol ,MMG5_pTria );
extern int    (*MMG2D_gradsiz)(MMG5_pMesh ,MMG5_pSol );

