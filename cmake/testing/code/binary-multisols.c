/**
 * Test binary multi-solution round-trips through the public API.
 */

#include <stdio.h>
#include <stdlib.h>

#if defined(MMG_MULTISOL_2D)
#include "mmg/mmg2d/libmmg2d.h"
#define API(name) MMG2D_##name
#elif defined(MMG_MULTISOL_3D)
#include "mmg/mmg3d/libmmg3d.h"
#define API(name) MMG3D_##name
#elif defined(MMG_MULTISOL_S)
#include "mmg/mmgs/libmmgs.h"
#define API(name) MMGS_##name
#else
#error "An MMG_MULTISOL_* test target must be selected"
#endif

static int set_mesh(MMG5_pMesh mesh) {
#if defined(MMG_MULTISOL_2D)
  if ( API(Set_meshSize)(mesh,4,2,0,0) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,0.0,0.0,0,1) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,1.0,0.0,0,2) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,0.0,1.0,0,3) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,1.0,1.0,0,4) != 1 )  return 0;
  if ( API(Set_triangle)(mesh,1,2,3,0,1) != 1 )  return 0;
  if ( API(Set_triangle)(mesh,2,4,3,0,2) != 1 )  return 0;
#elif defined(MMG_MULTISOL_3D)
  if ( API(Set_meshSize)(mesh,4,1,0,0,0,0) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,0.0,0.0,0.0,0,1) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,1.0,0.0,0.0,0,2) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,0.0,1.0,0.0,0,3) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,0.0,0.0,1.0,0,4) != 1 )  return 0;
  if ( API(Set_tetrahedron)(mesh,1,2,3,4,0,1) != 1 )  return 0;
#else
  if ( API(Set_meshSize)(mesh,4,2,0) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,0.0,0.0,0.0,0,1) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,1.0,0.0,0.0,0,2) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,0.0,1.0,0.0,0,3) != 1 )  return 0;
  if ( API(Set_vertex)(mesh,0.0,0.0,1.0,0,4) != 1 )  return 0;
  if ( API(Set_triangle)(mesh,1,2,3,0,1) != 1 )  return 0;
  if ( API(Set_triangle)(mesh,1,2,4,0,2) != 1 )  return 0;
#endif
  return 1;
}

static int same_values(const double *actual,const double *expected) {
  int i;

  for ( i=0; i<4; ++i ) {
    if ( actual[i] != expected[i] )  return 0;
  }
  return 1;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh;
  MMG5_pSol  met,written,loaded;
  MMG5_int   np;
  double     first[4] = {1.0,2.0,3.0,4.0};
  double     second[4] = {5.0,6.0,7.0,8.0};
  double     actual[4];
  int        nsols,status,types[2],loaded_types[MMG5_NSOLS_MAX];

  if ( argc != 2 ) {
    fprintf(stderr,"Usage: %s output.solb\n",argv[0]);
    return EXIT_FAILURE;
  }

  mesh = NULL;
  met = written = loaded = NULL;
  status = EXIT_FAILURE;
  types[0] = types[1] = MMG5_Scalar;

  if ( API(Init_mesh)(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,
                      MMG5_ARG_ppMet,&met,MMG5_ARG_end) != 1 )  goto cleanup;
  if ( !set_mesh(mesh) )  goto cleanup;
  if ( API(Set_solsAtVerticesSize)(mesh,&written,2,4,types) != 1 )
    goto cleanup;
  if ( API(Set_ithSols_inSolsAtVertices)(written,1,first) != 1 )
    goto cleanup;
  if ( API(Set_ithSols_inSolsAtVertices)(written,2,second) != 1 )
    goto cleanup;
  if ( API(saveAllSols)(mesh,&written,argv[1]) != 1 )  goto cleanup;
  if ( API(loadAllSols)(mesh,&loaded,argv[1]) != 1 )  goto cleanup;
  if ( API(Get_solsAtVerticesSize)(mesh,&loaded,&nsols,&np,loaded_types) != 1 )
    goto cleanup;
  if ( nsols != 2 || np != 4 || loaded_types[0] != MMG5_Scalar ||
       loaded_types[1] != MMG5_Scalar )  goto cleanup;
  if ( API(Get_ithSols_inSolsAtVertices)(loaded,1,actual) != 1 )
    goto cleanup;
  if ( !same_values(actual,first) )  goto cleanup;
  if ( API(Get_ithSols_inSolsAtVertices)(loaded,2,actual) != 1 )
    goto cleanup;
  if ( !same_values(actual,second) )  goto cleanup;

  status = EXIT_SUCCESS;

cleanup:
  if ( loaded )  API(Free_allSols)(mesh,&loaded);
  if ( written )  API(Free_allSols)(mesh,&written);
  API(Free_all)(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,
                MMG5_ARG_ppMet,&met,MMG5_ARG_end);
  return status;
}
