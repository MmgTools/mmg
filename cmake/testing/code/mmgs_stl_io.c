/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

#include "mmg/mmgs/libmmgs.h"

#include <stdio.h>

static int writeFacet(FILE *out,const double coordinates[9]) {
  int i;

  fprintf(out,"facet normal 0 0 0\nouter loop\n");
  for ( i=0; i<3; ++i ) {
    fprintf(out,"vertex %.17g %.17g %.17g\n",coordinates[3*i],
            coordinates[3*i+1],coordinates[3*i+2]);
  }
  fprintf(out,"endloop\nendfacet\n");
  return !ferror(out);
}

static int writeInput(const char *filename) {
  const double facets[4][9] = {
    {0,0,0, 0,1,0, 1,0,0},
    {1e-15,0,0, 1,0,0, 0,0,1},
    {0,0,0, 0,0,1, 0,1,0},
    {1,0,0, 0,1,0, 0,0,1}
  };
  FILE *out = fopen(filename,"w");
  int  i;

  if ( !out ) return 0;
  fprintf(out,"solid tetrahedron\n");
  for ( i=0; i<4; ++i ) {
    if ( !writeFacet(out,facets[i]) ) { fclose(out); return 0; }
  }
  fprintf(out,"endsolid tetrahedron\n");
  return !fclose(out);
}

static int writeEmptyBinary(const char *filename) {
  unsigned char header[84] = {0};
  FILE          *out = fopen(filename,"wb");

  if ( !out ) return 0;
  return fwrite(header,sizeof(header),1,out) == 1 && !fclose(out);
}

static int rejectsEmptyBinary(const char *filename) {
  MMG5_pMesh mesh = NULL;
  int        rejected;

  if ( !writeEmptyBinary(filename) ) return 0;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  rejected = MMGS_loadStlMesh(mesh,filename) < 1;
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return rejected;
}

static int writeTranslatedInput(const char *filename) {
  const double origin = 1.e15;
  const double facets[2][9] = {
    {origin,origin,origin, origin+4.,origin,origin,
     origin+4.,origin+4.,origin},
    {origin,origin,origin, origin+4.,origin+4.,origin,
     origin,origin+4.,origin}
  };
  FILE *out = fopen(filename,"w");
  int  i;

  if ( !out ) return 0;
  fprintf(out,"solid translated_square\n");
  for ( i=0; i<2; ++i ) {
    if ( !writeFacet(out,facets[i]) ) { fclose(out); return 0; }
  }
  fprintf(out,"endsolid translated_square\n");
  return !fclose(out);
}

static int readsTranslatedInput(const char *filename) {
  MMG5_pMesh mesh = NULL;
  MMG5_int   np,nt,na;
  int        valid = 0;

  if ( !writeTranslatedInput(filename) ) return 0;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadStlMesh(mesh,filename) == 1 &&
       MMGS_Get_meshSize(mesh,&np,&nt,&na) && np == 4 && nt == 2 && !na ) {
    valid = 1;
  }
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return valid;
}

static int readsTinyInput(const char *filename) {
  const double coordinates[9] = {
    0.,0.,0., 1.e-15,0.,0., 0.,1.e-15,0.
  };
  MMG5_pMesh mesh = NULL;
  MMG5_int   np,nt,na;
  FILE       *out = fopen(filename,"w");
  int        valid = 0;

  if ( !out ) return 0;
  fprintf(out,"solid tiny_triangle\n");
  if ( !writeFacet(out,coordinates) ) { fclose(out); return 0; }
  fprintf(out,"endsolid tiny_triangle\n");
  if ( fclose(out) ) return 0;

  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadStlMesh(mesh,filename) == 1 &&
       MMGS_Get_meshSize(mesh,&np,&nt,&na) && np == 3 && nt == 1 && !na ) {
    valid = 1;
  }
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return valid;
}

static int checkMesh(MMG5_pMesh mesh) {
  MMG5_int np,nt,na;

  return MMGS_Get_meshSize(mesh,&np,&nt,&na) && np == 4 && nt == 4 && !na;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh = NULL;
  int        ier = 1;

  if ( argc != 6 || !writeInput(argv[1]) ||
       !rejectsEmptyBinary(argv[4]) || !readsTinyInput(argv[4]) ||
       !readsTranslatedInput(argv[5]) ) {
    return 1;
  }
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadGenericMesh(mesh,NULL,NULL,argv[1]) != 1 || !checkMesh(mesh) ) {
    goto cleanup;
  }
  if ( MMGS_saveGenericMesh(mesh,NULL,argv[2]) != 1 ||
       MMGS_saveStlMesh(mesh,argv[3]) != 1 ) goto cleanup;

  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadStlMesh(mesh,argv[2]) != 1 || !checkMesh(mesh) ) goto cleanup;
  ier = 0;

cleanup:
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return ier;
}
