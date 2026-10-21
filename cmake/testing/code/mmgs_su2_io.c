/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file mmgs_su2_io.c
 *  \brief Test native SU2 surface input/output through the public MMGS API.
 */

#include "mmg/mmgs/libmmgs.h"

#include <stdio.h>

static int writeBoundaryInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"NDIME= 3\nNELEM= 1\n10 0 1 2 4 0\n");
  fprintf(out,"NPOIN= 5\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n0 0 1\n");
  fprintf(out,"NMARK= 2\n");
  fprintf(out,"MARKER_TAG= wall\nMARKER_ELEMS= 1\n9 0 1 2 3\n");
  /* Reserve ref 1 before assigning the fallback for the earlier wall marker. */
  fprintf(out,"MARKER_TAG= mmg_ref_1\nMARKER_ELEMS= 1\n5 0 3 4\n");
  return !fclose(out);
}

static int writeDomainInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  /* Exercise the NPOIN-first ordering used by third-party writers. */
  fprintf(out,"NDIME= 3\nNPOIN= 5\n");
  fprintf(out,"0 0 0\n1 0 0\n1 1 0\n0 1 0\n0 0 1\n");
  fprintf(out,"NELEM= 2\n5 0 3 4\n9 0 1 2 3\nNMARK= 0\n");
  return !fclose(out);
}

static int writeRejectedInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"NDIME= 3\nNELEM= 1\n10 0 1 2 3\n");
  fprintf(out,"NPOIN= 4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\nNMARK= 0\n");
  return !fclose(out);
}

static int checkMesh(MMG5_pMesh mesh,int domain) {
  MMG5_int np,nt,na,a,b,c,ref;
  int      required,nref0=0,nref1=0,nref2=0;

  if ( !MMGS_Get_meshSize(mesh,&np,&nt,&na) || np != 5 || nt != 3 || na ) {
    return 0;
  }
  while ( nt-- ) {
    if ( !MMGS_Get_triangle(mesh,&a,&b,&c,&ref,&required) ) return 0;
    if ( ref == 0 ) ++nref0;
    else if ( ref == 1 ) ++nref1;
    else if ( ref == 2 ) ++nref2;
    else return 0;
  }
  return domain ? nref0 == 3 : (nref1 == 1 && nref2 == 2);
}

static int rejectsVolumeWithoutBoundary(const char *filename) {
  MMG5_pMesh mesh = NULL;
  int        rejected;

  if ( !writeRejectedInput(filename) ) return 0;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  rejected = MMGS_loadSu2Mesh(mesh,filename) < 1;
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return rejected;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh = NULL;
  int        ier = 1;

  if ( argc != 6 || !writeBoundaryInput(argv[1]) || !writeDomainInput(argv[4]) ||
       !rejectsVolumeWithoutBoundary(argv[5]) ) return 1;

  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadGenericMesh(mesh,NULL,NULL,argv[1]) != 1 ||
       !checkMesh(mesh,0) || MMGS_saveGenericMesh(mesh,NULL,argv[2]) != 1 ||
       MMGS_saveSu2Mesh(mesh,argv[3]) != 1 ) goto cleanup;

  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadSu2Mesh(mesh,argv[2]) != 1 || !checkMesh(mesh,0) ) goto cleanup;

  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadSu2Mesh(mesh,argv[3]) != 1 || !checkMesh(mesh,0) ) goto cleanup;

  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadSu2Mesh(mesh,argv[4]) != 1 || !checkMesh(mesh,1) ) goto cleanup;
  ier = 0;

cleanup:
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return ier;
}
