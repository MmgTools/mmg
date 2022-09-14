#include "libmmgtypes.h"
#include "mmgs_export.h"

#define MMG_EXTERN
#define MMG_ASSIGN_NULL =NULL

#include  "mmgsexterns.h"

LIBMMGS_EXPORT int (*MMGS_doSol)(MMG5_pMesh mesh,MMG5_pSol met)=NULL;
