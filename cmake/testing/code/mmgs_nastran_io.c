/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file mmgs_nastran_io.c
 *  \brief Test native Nastran input/output through the public MMGS API.
 */

#include "mmg/mmgs/libmmgs.h"

#include <stdio.h>

static int writeInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"$MMG_REF,1,-7\nBEGIN BULK\n");
  fprintf(out,"GRID,10,,0.0,0.0,0.0\n");
  fprintf(out,"%-8s%8d%8s%8s%8s%8s\n","GRID",20,"","1.0","0.0","0.0");
  fprintf(out,"%-8s%8d%8s%8s%8s%8s\n","GRID",30,"","0.0","1.0","0.0");
  fprintf(out,"%-8s%16d%16s%16s%16s\n","GRID*",40,"","2.0","0.0");
  fprintf(out,"%-8s%16s\n","*","0.0");
  /* Exercise Nastran's exponent-without-E real syntax. */
  fprintf(out,"%-8s%8d%8s%8s%8s%8s\n","GRID",50,"","2.0+0","1.0","0.0");
  fprintf(out,"CTRIAR,1,1,10,20,30\n");
  fprintf(out,"%-8s%8d%8d%8d%8d%8d%8d\n",
          "CQUADR",2,42,20,40,50,30);
  fprintf(out,"ENDDATA\n");
  return !fclose(out);
}

static int writeRejectedInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"BEGIN BULK\nGRID,1,9,0.0,0.0,0.0\nENDDATA\n");
  return !fclose(out);
}

static int writeHighOrderInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"BEGIN BULK\nGRID,1,,0,0,0\nGRID,2,,1,0,0\n");
  fprintf(out,"GRID,3,,0,1,0\nGRID,4,,.5,0,0\n");
  fprintf(out,"GRID,5,,.5,.5,0\nGRID,6,,0,.5,0\n");
  fprintf(out,"CTRIA6,1,1,1,2,3,4,5,6\nENDDATA\n");
  return !fclose(out);
}

static int checkMesh(MMG5_pMesh mesh) {
  MMG5_int np,nt,na,v0,v1,v2,ref;
  int      required,nrefMinus7=0,nref42=0;

  if ( !MMGS_Get_meshSize(mesh,&np,&nt,&na) || np != 5 || nt != 3 || na ) {
    return 0;
  }
  while ( nt-- ) {
    if ( !MMGS_Get_triangle(mesh,&v0,&v1,&v2,&ref,&required) ) return 0;
    if ( ref == -7 ) ++nrefMinus7;
    else if ( ref == 42 ) ++nref42;
    else return 0;
  }
  return nrefMinus7 == 1 && nref42 == 2;
}

static int rejectsNonBasicCoordinates(const char *filename) {
  MMG5_pMesh mesh = NULL;
  int        rejected,highOrderRejected;

  if ( !writeRejectedInput(filename) ) return 0;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  rejected = MMGS_loadNastranMesh(mesh,filename) < 1;
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( !rejected || !writeHighOrderInput(filename) ) return 0;
  mesh = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  highOrderRejected = MMGS_loadNastranMesh(mesh,filename) < 1;
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return highOrderRejected;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh = NULL;
  int        ier = 1;

  if ( argc != 5 || !writeInput(argv[1]) ||
       !rejectsNonBasicCoordinates(argv[4]) ) return 1;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadGenericMesh(mesh,NULL,NULL,argv[1]) != 1 || !checkMesh(mesh) ||
       MMGS_saveGenericMesh(mesh,NULL,argv[2]) != 1 ||
       MMGS_saveNastranMesh(mesh,argv[3]) != 1 ) goto cleanup;
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadNastranMesh(mesh,argv[2]) != 1 || !checkMesh(mesh) ) {
    goto cleanup;
  }
  ier = 0;

cleanup:
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return ier;
}
