/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file mmg3d_nastran_io.c
 *  \brief Test native Nastran input/output through the public MMG3D API.
 */

#include "mmg/mmg3d/libmmg3d.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int writeNativeInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"BEGIN BULK\n");
  fprintf(out,"GRID,10,,0,0,0\nGRID,20,,1,0,0\n");
  fprintf(out,"GRID,30,,0,1,0\nGRID,40,,0,0,1\n");
  fprintf(out,"GRID,50,,2,0,0\nGRID,60,,3,0,0\nGRID,70,,2,1,0\n");
  fprintf(out,"GRID,80,,2,0,1\nGRID,90,,3,0,1\nGRID,100,,2,1,1\n");
  /* A blank PID is common in mesh-only exchange decks and maps to ref 0. */
  fprintf(out,"CTETRA,1,,10,20,30,40\n");
  fprintf(out,"%-8s%8d%8d%8d%8d%8d%8d%8d%8d\n",
          "CPENTA",2,22,50,60,70,80,90,100);
  fprintf(out,"CTRIA3,3,31,10,30,20\n");
  fprintf(out,"CQUAD4,4,32,50,60,90,80\nENDDATA\n");
  return !fclose(out);
}

static void writeHex(FILE *out,int eid,int pid,const int vertex[8],
                     const char *continuation) {
  fprintf(out,"%-8s%8d%8d%8d%8d%8d%8d%8d%8d%-8s\n","CHEXA",eid,pid,
          vertex[0],vertex[1],vertex[2],vertex[3],vertex[4],vertex[5],
          continuation);
  fprintf(out,"%-8s%8d%8d\n",continuation,vertex[6],vertex[7]);
}

static int writeMixedInput(const char *filename) {
  static const int hex0[8] = {1,2,3,4,5,6,7,8};
  static const int hex1[8] = {2,9,10,3,6,11,12,7};
  FILE *out = fopen(filename,"w");
  int  i;
  const double point[23][3] = {
    {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},
    {1,1,1},{0,1,1},{2,0,0},{2,1,0},{2,0,1},{2,1,1},
    {0,3,0},{1,3,0},{1,4,0},{0,4,0},{.5,3.5,1},
    {0,6,0},{1,6,0},{0,7,0},{0,6,1},{1,6,1},{0,7,1}
  };

  if ( !out ) return 0;
  fprintf(out,"$MMG_REF,1,101\n$MMG_REF,2,102\n");
  fprintf(out,"$MMG_REF,3,103\n$MMG_REF,4,104\nBEGIN BULK\n");
  for ( i=0; i<23; ++i ) {
    fprintf(out,"GRID,%d,,%.17g,%.17g,%.17g\n",i+1,point[i][0],
            point[i][1],point[i][2]);
  }
  /* Shared face order is opposite in these adjacent fixed-field CHEXA cards. */
  writeHex(out,1,1,hex0,"+H1");
  writeHex(out,2,2,hex1,"+H2");
  fprintf(out,"CPYRAM,3,3,13,14,15,16,17\n");
  fprintf(out,"CPENTA,4,4,18,19,20,21,22,23\n");
  fprintf(out,"CQUAD4,5,7,1,4,8,5\nENDDATA\n");
  return !fclose(out);
}

static int writeRejectedInput(const char *filename) {
  FILE *out = fopen(filename,"w");
  int  i;

  if ( !out ) return 0;
  fprintf(out,"BEGIN BULK\n");
  for ( i=1; i<=10; ++i ) fprintf(out,"GRID,%d,,%d,0,0\n",i,i);
  fprintf(out,"CTETRA,1,1,1,2,3,4,5,6,7,8,9,10\nENDDATA\n");
  return !fclose(out);
}

static int checkNativeMesh(MMG5_pMesh mesh) {
  MMG5_int np,ne,nprism,nt,nquad,na,v[6],ref;
  int      required;

  if ( !MMG3D_Get_meshSize(mesh,&np,&ne,&nprism,&nt,&nquad,&na) ||
       np != 10 || ne != 1 || nprism != 1 || nt != 1 || nquad != 1 || na ) {
    return 0;
  }
  if ( !MMG3D_Get_tetrahedron(mesh,&v[0],&v[1],&v[2],&v[3],&ref,
                              &required) || ref ) return 0;
  if ( !MMG3D_Get_prism(mesh,&v[0],&v[1],&v[2],&v[3],&v[4],&v[5],&ref,
                        &required) || ref != 22 ) return 0;
  if ( !MMG3D_Get_triangle(mesh,&v[0],&v[1],&v[2],&ref,&required) ||
       ref != 31 ) return 0;
  return MMG3D_Get_quadrilateral(mesh,&v[0],&v[1],&v[2],&v[3],&ref,
                                 &required) && ref == 32;
}

static void sort3(MMG5_int value[3]) {
  MMG5_int tmp;

  if ( value[0] > value[1] ) { tmp=value[0]; value[0]=value[1]; value[1]=tmp; }
  if ( value[1] > value[2] ) { tmp=value[1]; value[1]=value[2]; value[2]=tmp; }
  if ( value[0] > value[1] ) { tmp=value[0]; value[0]=value[1]; value[1]=tmp; }
}

static int countFace(const MMG5_int *tetra,MMG5_int ne,
                     MMG5_int a,MMG5_int b,MMG5_int c) {
  static const int localFace[4][3] = {
    {0,1,2},{0,1,3},{0,2,3},{1,2,3}
  };
  MMG5_int target[3] = {a,b,c},face[3],i;
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
  const double *a=&points[3*(tetra[0]-1)],*b=&points[3*(tetra[1]-1)];
  const double *c=&points[3*(tetra[2]-1)],*d=&points[3*(tetra[3]-1)];
  double ab[3],ac[3],ad[3],det;
  int    i;

  for ( i=0; i<3; ++i ) { ab[i]=b[i]-a[i]; ac[i]=c[i]-a[i]; ad[i]=d[i]-a[i]; }
  det = ab[0]*(ac[1]*ad[2]-ac[2]*ad[1])
        -ab[1]*(ac[0]*ad[2]-ac[2]*ad[0])
        +ab[2]*(ac[0]*ad[1]-ac[1]*ad[0]);
  return fabs(det)/6.;
}

static int checkConvertedMesh(MMG5_pMesh mesh) {
  MMG5_int *tetra=NULL,*refs=NULL,np,ne,nprism,nt,nquad,na,i;
  double   *points=NULL,volume=0.;
  int      n101=0,n102=0,n103=0,n104=0,valid=0;

  if ( !MMG3D_Get_meshSize(mesh,&np,&ne,&nprism,&nt,&nquad,&na) ||
       np != 27 || ne != 38 || nprism || nt || nquad != 1 || na ) return 0;
  tetra = (MMG5_int *)malloc((size_t)4*ne*sizeof(MMG5_int));
  refs = (MMG5_int *)malloc((size_t)ne*sizeof(MMG5_int));
  points = (double *)malloc((size_t)3*np*sizeof(double));
  if ( !tetra || !refs || !points ||
       !MMG3D_Get_tetrahedra(mesh,tetra,refs,NULL) ||
       !MMG3D_Get_vertices(mesh,points,NULL,NULL,NULL) ) goto cleanup;
  for ( i=0; i<ne; ++i ) {
    volume += tetraVolume(points,&tetra[4*i]);
    if ( refs[i] == 101 ) ++n101;
    else if ( refs[i] == 102 ) ++n102;
    else if ( refs[i] == 103 ) ++n103;
    else if ( refs[i] == 104 ) ++n104;
    else goto cleanup;
  }
  valid = fabs(volume-17./6.) < 1.e-12 && n101 == 12 && n102 == 12 &&
          n103 == 6 && n104 == 8 && countFace(tetra,ne,2,3,7) == 2 &&
          countFace(tetra,ne,2,6,7) == 2;

cleanup:
  free(tetra); free(refs); free(points);
  return valid;
}

static int rejectsQuadraticSolid(const char *filename) {
  MMG5_pMesh mesh = NULL;
  int        rejected;

  if ( !writeRejectedInput(filename) ) return 0;
  MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  rejected = MMG3D_loadNastranMesh(mesh,filename) < 1;
  MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return rejected;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh = NULL;
  int        ier = 1;

  if ( argc != 7 || !writeNativeInput(argv[1]) || !writeMixedInput(argv[4]) ||
       !rejectsQuadraticSolid(argv[6]) ) return 1;
  MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMG3D_loadGenericMesh(mesh,NULL,NULL,argv[1]) != 1 ||
       !checkNativeMesh(mesh) || MMG3D_saveGenericMesh(mesh,NULL,argv[2]) != 1 ||
       MMG3D_saveNastranMesh(mesh,argv[3]) != 1 ) goto cleanup;
  MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMG3D_loadNastranMesh(mesh,argv[2]) != 1 || !checkNativeMesh(mesh) ) {
    goto cleanup;
  }
  MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMG3D_loadNastranMesh(mesh,argv[4]) != 1 || !checkConvertedMesh(mesh) ||
       MMG3D_saveNastranMesh(mesh,argv[5]) != 1 ) goto cleanup;
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
