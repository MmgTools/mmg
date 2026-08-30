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

static int checkMesh(MMG5_pMesh mesh) {
  MMG5_int np,nt,na;

  return MMGS_Get_meshSize(mesh,&np,&nt,&na) && np == 4 && nt == 4 && !na;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh = NULL;
  int        ier = 1;

  if ( argc != 4 || !writeInput(argv[1]) ) return 1;
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
