/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file mmg3d_su2_io.c
 *  \brief Test native SU2 input/output through the public MMG3D API.
 */

#include "mmg/mmg3d/libmmg3d.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int writeNativeInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  /* meshio writes the point section before the element section. */
  fprintf(out,"NDIME= 3\nNPOIN= 10\n");
  fprintf(out,"0 0 0 100\n1 0 0 101\n0 1 0 102\n0 0 1 103\n");
  fprintf(out,"2 0 0 104\n3 0 0 105\n2 1 0 106\n");
  fprintf(out,"2 0 1 107\n3 0 1 108\n2 1 1 109\n");
  fprintf(out,"NELEM= 2\n");
  fprintf(out,"10 0 1 2 3 81\n");
  /* SU2 orders each triangular prism face opposite to Mmg. */
  fprintf(out,"13 4 6 5 7 9 8 82\n");
  fprintf(out,"NMARK= 2\n");
  fprintf(out,"MARKER_TAG= wall\nMARKER_ELEMS= 1\n5 0 2 1\n");
  fprintf(out,"MARKER_TAG= mmg_ref_42\nMARKER_ELEMS= 1\n9 4 5 8 7\n");
  return !fclose(out);
}

static int writeMixedInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"NDIME= 3\nNELEM= 4\n");
  /* Two adjacent hexahedra deliberately use opposite shared-face ordering. */
  fprintf(out,"12 0 1 2 3 4 5 6 7 0\n");
  fprintf(out,"12 1 8 9 2 5 10 11 6 1\n");
  fprintf(out,"14 12 13 14 15 16 2\n");
  fprintf(out,"13 17 18 19 20 21 22 3\n");
  fprintf(out,"NPOIN= 23\n");
  fprintf(out,"0 0 0\n1 0 0\n1 1 0\n0 1 0\n");
  fprintf(out,"0 0 1\n1 0 1\n1 1 1\n0 1 1\n");
  fprintf(out,"2 0 0\n2 1 0\n2 0 1\n2 1 1\n");
  fprintf(out,"0 3 0\n1 3 0\n1 4 0\n0 4 0\n0.5 3.5 1\n");
  fprintf(out,"0 6 0\n1 6 0\n0 7 0\n0 6 1\n1 6 1\n0 7 1\n");
  fprintf(out,"NMARK= 1\n");
  fprintf(out,"MARKER_TAG= mmg_ref_7\nMARKER_ELEMS= 1\n");
  fprintf(out,"9 0 3 7 4\n");
  return !fclose(out);
}

static int writeRejectedInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"NDIME= 3\nNELEM= 1\n11 0 1 2 3 4 5 6 7\n");
  fprintf(out,"NPOIN= 8\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n");
  fprintf(out,"0 0 1\n1 0 1\n1 1 1\n0 1 1\nNMARK= 0\n");
  return !fclose(out);
}

static int checkNativeMesh(MMG5_pMesh mesh) {
  MMG5_int np,ne,nprism,nt,nquad,na,v[4],ref;
  int      required;

  if ( !MMG3D_Get_meshSize(mesh,&np,&ne,&nprism,&nt,&nquad,&na) ||
       np != 10 || ne != 1 || nprism != 1 || nt != 1 || nquad != 1 || na ) {
    return 0;
  }
  if ( !MMG3D_Get_triangle(mesh,&v[0],&v[1],&v[2],&ref,&required) ||
       ref != 1 ) return 0;
  return MMG3D_Get_quadrilateral(mesh,&v[0],&v[1],&v[2],&v[3],&ref,
                                 &required) && ref == 42;
}

static void sort3(MMG5_int value[3]) {
  MMG5_int tmp;

  if ( value[0] > value[1] ) {
    tmp = value[0]; value[0] = value[1]; value[1] = tmp;
  }
  if ( value[1] > value[2] ) {
    tmp = value[1]; value[1] = value[2]; value[2] = tmp;
  }
  if ( value[0] > value[1] ) {
    tmp = value[0]; value[0] = value[1]; value[1] = tmp;
  }
}

static int countFace(const MMG5_int *tetra,MMG5_int ne,
                     MMG5_int a,MMG5_int b,MMG5_int c) {
  static const int localFace[4][3] = {
    {0,1,2},{0,1,3},{0,2,3},{1,2,3}
  };
  MMG5_int target[3] = {a,b,c};
  MMG5_int face[3],i;
  int      count=0,j,k;

  sort3(target);
  for ( i=0; i<ne; ++i ) {
    for ( j=0; j<4; ++j ) {
      for ( k=0; k<3; ++k ) face[k] = tetra[4*i+localFace[j][k]];
      sort3(face);
      if ( face[0] == target[0] && face[1] == target[1] &&
           face[2] == target[2] ) ++count;
    }
  }
  return count;
}

static double tetraVolume(const double *points,const MMG5_int *tetra) {
  const double *a = &points[3*(tetra[0]-1)];
  const double *b = &points[3*(tetra[1]-1)];
  const double *c = &points[3*(tetra[2]-1)];
  const double *d = &points[3*(tetra[3]-1)];
  double ab[3],ac[3],ad[3],det;
  int    i;

  for ( i=0; i<3; ++i ) {
    ab[i] = b[i]-a[i];
    ac[i] = c[i]-a[i];
    ad[i] = d[i]-a[i];
  }
  det = ab[0]*(ac[1]*ad[2]-ac[2]*ad[1])
        -ab[1]*(ac[0]*ad[2]-ac[2]*ad[0])
        +ab[2]*(ac[0]*ad[1]-ac[1]*ad[0]);
  return fabs(det)/6.;
}

static int checkConvertedMesh(MMG5_pMesh mesh) {
  MMG5_int *tetra = NULL;
  double   *points = NULL;
  MMG5_int np,ne,nprism,nt,nquad,na,i;
  MMG5_int quad[4],ref;
  double   volume=0.;
  int      required,valid=0;

  if ( !MMG3D_Get_meshSize(mesh,&np,&ne,&nprism,&nt,&nquad,&na) ||
       np != 27 || ne != 38 || nprism || nt || nquad != 1 || na ||
       !MMG3D_Get_quadrilateral(mesh,&quad[0],&quad[1],&quad[2],&quad[3],
                                &ref,&required) || ref != 7 ) return 0;
  tetra = (MMG5_int *)malloc((size_t)4*ne*sizeof(MMG5_int));
  points = (double *)malloc((size_t)3*np*sizeof(double));
  if ( !tetra || !points ||
       !MMG3D_Get_tetrahedra(mesh,tetra,NULL,NULL) ||
       !MMG3D_Get_vertices(mesh,points,NULL,NULL,NULL) ) goto cleanup;

  for ( i=0; i<ne; ++i ) volume += tetraVolume(points,&tetra[4*i]);
  valid = fabs(volume-17./6.) < 1.e-12 &&
          countFace(tetra,ne,2,3,7) == 2 &&
          countFace(tetra,ne,2,6,7) == 2;

cleanup:
  free(tetra);
  free(points);
  return valid;
}

static int rejectsUnsupportedVolume(const char *filename) {
  MMG5_pMesh mesh = NULL;
  int        rejected;

  if ( !writeRejectedInput(filename) ) return 0;
  MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  rejected = MMG3D_loadSu2Mesh(mesh,filename) < 1;
  MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return rejected;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh = NULL;
  int        ier = 1;

  if ( argc != 7 || !writeNativeInput(argv[1]) || !writeMixedInput(argv[4]) ||
       !rejectsUnsupportedVolume(argv[6]) ) return 1;

  MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMG3D_loadGenericMesh(mesh,NULL,NULL,argv[1]) != 1 ||
       !checkNativeMesh(mesh) ||
       MMG3D_saveGenericMesh(mesh,NULL,argv[2]) != 1 ||
       MMG3D_saveSu2Mesh(mesh,argv[3]) != 1 ) goto cleanup;

  MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMG3D_loadSu2Mesh(mesh,argv[2]) != 1 || !checkNativeMesh(mesh) ) {
    goto cleanup;
  }

  MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMG3D_loadSu2Mesh(mesh,argv[4]) != 1 || !checkConvertedMesh(mesh) ||
       MMG3D_saveSu2Mesh(mesh,argv[5]) != 1 ) goto cleanup;

  MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMG3D_loadGenericMesh(mesh,NULL,NULL,argv[5]) != 1 ||
       !checkConvertedMesh(mesh) ) goto cleanup;
  ier = 0;

cleanup:
  MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return ier;
}
