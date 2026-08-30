/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file mmgs_ply_io.c
 *  \brief Test native PLY input/output through the public MMGS API.
 */

#include "mmg/mmgs/libmmgs.h"

#include <stdio.h>

static MMG5_int wideReference(void) {
  return sizeof(MMG5_int) > 4 ? (MMG5_int)INT64_C(5000000000) : 42;
}

static int writeInput(const char *filename) {
  FILE *out = fopen(filename,"wb");

  if ( !out ) return 0;
  /* Deliberately use CRLF and nonstandard property ordering. */
  fprintf(out,"ply\r\nformat ascii 1.0\r\ncomment property order test\r\n");
  fprintf(out,"element vertex 5\r\nproperty uchar red\r\n");
  fprintf(out,"property double z\r\nproperty double x\r\n");
  fprintf(out,"property double y\r\n");
  fprintf(out,"element face 2\r\nproperty %s ref\r\n",
          sizeof(MMG5_int) > 4 ? "int64" : "int");
  fprintf(out,"property list uchar int vertex_indices\r\nend_header\r\n");
  fprintf(out,"255 0 0 0\r\n255 0 1 0\r\n255 0 1 1\r\n");
  fprintf(out,"255 0 0 1\r\n255 1 0 0\r\n");
  fprintf(out,"7 4 0 1 2 3\r\n%" MMG5_PRId " 3 0 3 4\r\n",
          wideReference());
  return !fclose(out);
}

static int writeRejectedInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"ply\nformat binary_big_endian 1.0\n");
  fprintf(out,"element vertex 1\nproperty float x\nproperty float y\n");
  fprintf(out,"property float z\nelement face 1\n");
  fprintf(out,"property list uchar int vertex_indices\nend_header\n");
  return !fclose(out);
}

static int readsLargeFloatCoordinates(const char *filename) {
  MMG5_pMesh mesh = NULL;
  FILE       *out = fopen(filename,"w");
  int        accepted;

  if ( !out ) return 0;
  fprintf(out,"ply\nformat ascii 1.0\nelement vertex 3\n");
  fprintf(out,"property double x\nproperty double y\nproperty double z\n");
  fprintf(out,"element face 1\nproperty list uchar int vertex_indices\n");
  fprintf(out,"end_header\n1e30 0 0\n1e30 1 0\n1e30 0 1\n3 0 1 2\n");
  if ( fclose(out) ) return 0;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  accepted = MMGS_loadPlyMesh(mesh,filename) == 1;
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return accepted;
}

static int checkMesh(MMG5_pMesh mesh) {
  MMG5_int np,nt,na,v0,v1,v2,ref;
  int      required,nref7 = 0,nrefWide = 0;

  if ( !MMGS_Get_meshSize(mesh,&np,&nt,&na) || np != 5 || nt != 3 || na ) {
    return 0;
  }
  while ( nt-- ) {
    if ( !MMGS_Get_triangle(mesh,&v0,&v1,&v2,&ref,&required) ) return 0;
    if ( ref == 7 ) ++nref7;
    else if ( ref == wideReference() ) ++nrefWide;
    else return 0;
  }
  return nref7 == 2 && nrefWide == 1;
}

static int rejectsBigEndian(const char *filename) {
  MMG5_pMesh mesh = NULL;
  int        rejected;

  if ( !writeRejectedInput(filename) ) return 0;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  rejected = MMGS_loadPlyMesh(mesh,filename) < 1;
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return rejected;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh = NULL;
  int        ier = 1;

  if ( argc != 5 || !writeInput(argv[1]) || !rejectsBigEndian(argv[4]) ||
       !readsLargeFloatCoordinates(argv[4]) ) {
    return 1;
  }
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadGenericMesh(mesh,NULL,NULL,argv[1]) != 1 || !checkMesh(mesh) ) {
    goto cleanup;
  }
  if ( MMGS_saveGenericMesh(mesh,NULL,argv[2]) != 1 ||
       MMGS_savePlyMesh(mesh,argv[3]) != 1 ) goto cleanup;

  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadPlyMesh(mesh,argv[2]) != 1 || !checkMesh(mesh) ) goto cleanup;
  ier = 0;

cleanup:
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return ier;
}
