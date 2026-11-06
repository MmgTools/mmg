/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file mmg2d_su2_io.c
 *  \brief Test native SU2 input/output through the public MMG2D API.
 */

#include "mmg/mmg2d/libmmg2d.h"

#include <stdio.h>

static int writeInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"%% Native triangle and quadrilateral with optional IDs\n");
  /* meshio writes the point section before the element section. */
  fprintf(out,"NDIME= 2\nNPOIN= 7\n");
  fprintf(out,"0 0 101\n1 0 102\n0 1 103\n");
  fprintf(out,"2 0 104\n3 0 105\n3 1 106\n2 1 107\n");
  fprintf(out,"NELEM= 2\n");
  fprintf(out,"5 0 1 2 17\n");
  fprintf(out,"9 3 4 5 6 23\n");
  fprintf(out,"NMARK= 2\n");
  fprintf(out,"MARKER_TAG= wall\nMARKER_ELEMS= 2\n");
  fprintf(out,"3 0 1\n3 1 2\n");
  /* This explicit value must be reserved before `wall` gets its fallback. */
  fprintf(out,"MARKER_TAG= mmg_ref_1\nMARKER_ELEMS= 2\n");
  fprintf(out,"3 3 4\n3 4 5\n");
  return !fclose(out);
}

static int writeRejectedInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"NDIME= 2\nNELEM= 1\n10 0 1 2 3\n");
  fprintf(out,"NPOIN= 4\n0 0\n1 0\n0 1\n1 1\nNMARK= 0\n");
  return !fclose(out);
}

static int checkMesh(MMG5_pMesh mesh) {
  MMG5_int np,nt,nquad,na,a,b,ref;
  int      isRidge,isRequired,nref1=0,nref2=0;

  if ( !MMG2D_Get_meshSize(mesh,&np,&nt,&nquad,&na) ||
       np != 7 || nt != 1 || nquad != 1 || na != 4 ) return 0;
  while ( na-- ) {
    if ( !MMG2D_Get_edge(mesh,&a,&b,&ref,&isRidge,&isRequired) ) return 0;
    if ( ref == 1 ) ++nref1;
    else if ( ref == 2 ) ++nref2;
    else return 0;
  }
  return nref1 == 2 && nref2 == 2;
}

static int rejectsUnsupportedVolume(const char *filename) {
  MMG5_pMesh mesh = NULL;
  int        rejected;

  if ( !writeRejectedInput(filename) ) return 0;
  MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  rejected = MMG2D_loadSu2Mesh(mesh,filename) < 1;
  MMG2D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return rejected;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh = NULL;
  int        ier = 1;

  if ( argc != 5 || !writeInput(argv[1]) ||
       !rejectsUnsupportedVolume(argv[4]) ) return 1;

  MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMG2D_loadGenericMesh(mesh,NULL,NULL,argv[1]) != 1 ||
       !checkMesh(mesh) ) goto cleanup;
  if ( MMG2D_saveGenericMesh(mesh,NULL,argv[2]) != 1 ||
       MMG2D_saveSu2Mesh(mesh,argv[3]) != 1 ) goto cleanup;

  MMG2D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMG2D_loadSu2Mesh(mesh,argv[2]) != 1 || !checkMesh(mesh) ) {
    goto cleanup;
  }
  ier = 0;

cleanup:
  MMG2D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return ier;
}
