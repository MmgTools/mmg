/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file mmgs_obj_io.c
 *  \brief Test native OBJ input/output through the public MMGS API.
 */

#include "mmg/mmgs/libmmgs.h"

#include <stdio.h>

static int writeInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"# polygon triangulation and OBJ index syntax\n");
  fprintf(out,"v 0 0 0\n");
  fprintf(out,"v 1 0 0\n");
  fprintf(out,"v 1 1 0\n");
  fprintf(out,"v 0 1 0\n");
  fprintf(out,"v 0 0 1\n");
  fprintf(out,"g skin\n");
  fprintf(out,"f 1/1 2/2 3/3 4/4\n");
  fprintf(out,"g mmg_ref_42\n");
  fprintf(out,"f -5//1 -2//1 -1//1\n");
  return !fclose(out);
}

static int checkMesh(MMG5_pMesh mesh) {
  MMG5_int np,nt,na,v0,v1,v2,ref;
  int      required,nref1 = 0,nref42 = 0;

  if ( !MMGS_Get_meshSize(mesh,&np,&nt,&na) || np != 5 || nt != 3 || na ) {
    return 0;
  }
  while ( nt-- ) {
    if ( !MMGS_Get_triangle(mesh,&v0,&v1,&v2,&ref,&required) ) return 0;
    if ( ref == 1 ) ++nref1;
    else if ( ref == 42 ) ++nref42;
    else return 0;
  }
  return nref1 == 2 && nref42 == 1;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh = NULL;
  int        ier = 1;

  if ( argc != 3 || !writeInput(argv[1]) ) return 1;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);

  if ( MMGS_loadGenericMesh(mesh,NULL,NULL,argv[1]) != 1 ) goto cleanup;
  if ( !checkMesh(mesh) ) goto cleanup;
  if ( MMGS_saveGenericMesh(mesh,NULL,argv[2]) != 1 ) goto cleanup;

  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadObjMesh(mesh,argv[2]) != 1 || !checkMesh(mesh) ) goto cleanup;
  ier = 0;

cleanup:
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return ier;
}
