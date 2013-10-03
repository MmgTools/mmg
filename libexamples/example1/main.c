/*Authors Cécile Dobrzynski

Example for using mmg3dlib

*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#include "libmmg3d5.h"

int main(int argc,char *argv[]) {
  pMesh      mmgMesh;
  pSol       mmgSol;
  int        opt_i[10],k,ier;
  double     opt_d[6];

  fprintf(stdout,"  -- TEST MMG3DLIB \n");

  mmgMesh = (pMesh)calloc(1,sizeof(Mesh));
  assert(mmgMesh);

  /* allocation */
  mmgMesh->ver    = 2;
  mmgMesh->dim    = 3;

  mmgMesh->np     = 12;
  mmgMesh->nt     = 0;
  mmgMesh->ne     = 12;


  mmgMesh->npmax  = 500000;
  mmgMesh->ntmax  = 1000000;
  mmgMesh->nemax  = 3000000;

  mmgMesh->point = (pPoint)calloc(mmgMesh->npmax+1,sizeof(Point));
  assert(mmgMesh->point);
  mmgMesh->tetra = (pTetra)calloc(mmgMesh->nemax+1,sizeof(Tetra));
  assert(mmgMesh->tetra);
  mmgMesh->tria  = (pTria)calloc(mmgMesh->ntmax+1,sizeof(Tria));
  assert(mmgMesh->tria);

  /*coordinates vertices*/
  mmgMesh->point[1].c[0]  = 0.;  mmgMesh->point[1].c[1]  = 0.; mmgMesh->point[1].c[2]  = 0.; mmgMesh->point[1].ref  = 0;
  mmgMesh->point[2].c[0]  = 0.5; mmgMesh->point[2].c[1]  = 0;  mmgMesh->point[2].c[2]  = 0;  mmgMesh->point[2].ref  = 0;
  mmgMesh->point[3].c[0]  = 0.5; mmgMesh->point[3].c[1]  = 0;  mmgMesh->point[3].c[2]  = 1;  mmgMesh->point[3].ref  = 0;
  mmgMesh->point[4].c[0]  = 0;   mmgMesh->point[4].c[1]  = 0;  mmgMesh->point[4].c[2]  = 1;  mmgMesh->point[4].ref  = 0;
  mmgMesh->point[5].c[0]  = 0;   mmgMesh->point[5].c[1]  = 1;  mmgMesh->point[5].c[2]  = 0;  mmgMesh->point[5].ref  = 0;
  mmgMesh->point[6].c[0]  = 0.5; mmgMesh->point[6].c[1]  = 1;  mmgMesh->point[6].c[2]  = 0;  mmgMesh->point[6].ref  = 0;
  mmgMesh->point[7].c[0]  = 0.5; mmgMesh->point[7].c[1]  = 1;  mmgMesh->point[7].c[2]  = 1;  mmgMesh->point[7].ref  = 0;
  mmgMesh->point[8].c[0]  = 0;   mmgMesh->point[8].c[1]  = 1;  mmgMesh->point[8].c[2]  = 1;  mmgMesh->point[8].ref  = 0;
  mmgMesh->point[9].c[0]  = 1;   mmgMesh->point[9].c[1]  = 0;  mmgMesh->point[9].c[2]  = 0;  mmgMesh->point[9].ref  = 0;
  mmgMesh->point[10].c[0] = 1;   mmgMesh->point[10].c[1] = 1;  mmgMesh->point[10].c[2] = 0;  mmgMesh->point[10].ref = 0;
  mmgMesh->point[11].c[0] = 1;   mmgMesh->point[11].c[1] = 0;  mmgMesh->point[11].c[2] = 1;  mmgMesh->point[11].ref = 0;
  mmgMesh->point[12].c[0] = 1;   mmgMesh->point[12].c[1] = 1;  mmgMesh->point[12].c[2] = 1;  mmgMesh->point[12].ref = 0;

  /* tetra*/
  /* warning: here we suppose that tetras are positively oriented */
  mmgMesh->tetra[1].v[0]  = 1;  mmgMesh->tetra[1].v[1]  = 4;  mmgMesh->tetra[1].v[2]  = 2;  mmgMesh->tetra[1].v[3]  = 8;  mmgMesh->tetra[1].ref  = 1;
  mmgMesh->tetra[2].v[0]  = 8;  mmgMesh->tetra[2].v[1]  = 3;  mmgMesh->tetra[2].v[2]  = 2;  mmgMesh->tetra[2].v[3]  = 7;  mmgMesh->tetra[2].ref  = 1;
  mmgMesh->tetra[3].v[0]  = 5;  mmgMesh->tetra[3].v[1]  = 2;  mmgMesh->tetra[3].v[2]  = 6;  mmgMesh->tetra[3].v[3]  = 8;  mmgMesh->tetra[3].ref  = 1;
  mmgMesh->tetra[4].v[0]  = 5;  mmgMesh->tetra[4].v[1]  = 8;  mmgMesh->tetra[4].v[2]  = 1;  mmgMesh->tetra[4].v[3]  = 2;  mmgMesh->tetra[4].ref  = 1;
  mmgMesh->tetra[5].v[0]  = 7;  mmgMesh->tetra[5].v[1]  = 2;  mmgMesh->tetra[5].v[2]  = 8;  mmgMesh->tetra[5].v[3]  = 6;  mmgMesh->tetra[5].ref  = 1;
  mmgMesh->tetra[6].v[0]  = 2;  mmgMesh->tetra[6].v[1]  = 4;  mmgMesh->tetra[6].v[2]  = 3;  mmgMesh->tetra[6].v[3]  = 8;  mmgMesh->tetra[6].ref  = 1;
  mmgMesh->tetra[7].v[0]  = 9;  mmgMesh->tetra[7].v[1]  = 2;  mmgMesh->tetra[7].v[2]  = 3;  mmgMesh->tetra[7].v[3]  = 7;  mmgMesh->tetra[7].ref  = 2;
  mmgMesh->tetra[8].v[0]  = 7;  mmgMesh->tetra[8].v[1]  = 11; mmgMesh->tetra[8].v[2]  = 9;  mmgMesh->tetra[8].v[3]  = 12; mmgMesh->tetra[8].ref  = 2;
  mmgMesh->tetra[9].v[0]  = 6;  mmgMesh->tetra[9].v[1]  = 9;  mmgMesh->tetra[9].v[2]  = 10; mmgMesh->tetra[9].v[3]  = 7;  mmgMesh->tetra[9].ref  = 2;
  mmgMesh->tetra[10].v[0] = 6;  mmgMesh->tetra[10].v[1] = 7;  mmgMesh->tetra[10].v[2] = 2;  mmgMesh->tetra[10].v[3] = 9;  mmgMesh->tetra[10].ref = 2;
  mmgMesh->tetra[11].v[0] = 12;  mmgMesh->tetra[11].v[1] = 9; mmgMesh->tetra[11].v[2] = 7;  mmgMesh->tetra[11].v[3] = 10; mmgMesh->tetra[11].ref = 2;
  mmgMesh->tetra[12].v[0] = 9;  mmgMesh->tetra[12].v[1] = 3;  mmgMesh->tetra[12].v[2] = 11; mmgMesh->tetra[12].v[3] = 7;  mmgMesh->tetra[12].ref = 2;


  /*metric*/
  mmgSol           = (pSol)calloc(1,sizeof(Sol));
  assert(mmgSol);
  mmgSol->size = 1;

  /*scalaire size*/
  mmgSol->np = mmgMesh->np;
  mmgSol->npmax = mmgMesh->npmax;
  mmgSol->m    = (double*)calloc(mmgSol->npmax+1,mmgSol->size*sizeof(double));
  assert(mmgSol->m);
  for(k=1 ; k<=mmgMesh->np ; k++) {
    mmgSol->m[k] = 0.5;
  }

  MMG5_mmg3dinit(opt_i, opt_d);

  opt_i[imprim]=5;
  opt_d[dhd] =50;
  mmgMesh->nameout = "result0.mesh";
  mmgSol->nameout = "result0.sol";


  ier = MMG5_mmg3dlib(opt_i,opt_d,mmgMesh,mmgSol);
  if ( ier == MG_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MG_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /*save result*/
  MMG5_saveMesh(mmgMesh);

  /*save metric*/
  MMG5_saveMet(mmgMesh,mmgSol);

  mmgMesh->nameout = "result2.mesh";
  mmgSol->nameout = "result2.sol";

  ier = MMG5_mmg3dlib(opt_i,opt_d,mmgMesh,mmgSol);
  if ( ier == MG_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MG_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /*save result*/
  MMG5_saveMesh(mmgMesh);

  /*save metric*/
  MMG5_saveMet(mmgMesh,mmgSol);

  opt_i[analysis] = 1;
  opt_d[dhd] = 40;

  mmgMesh->nameout = "result1.mesh";
  mmgSol->nameout = "result1.sol";

  ier = MMG5_mmg3dlib(opt_i,opt_d,mmgMesh,mmgSol);
  if ( ier == MG_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  }
  else if ( ier == MG_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /*save result*/
  MMG5_saveMesh(mmgMesh);

  /*save metric*/
  MMG5_saveMet(mmgMesh,mmgSol);


  /* free mem */
  free(mmgMesh->point);
  mmgMesh->point = NULL;
  free(mmgMesh->tria);
  mmgMesh->tria = NULL;
  free(mmgMesh->tetra);
  mmgMesh->tetra = NULL;
  free(mmgMesh);
  if ( mmgSol->np ) {
    free(mmgSol->m);
    mmgSol->m = NULL;
  }
  free(mmgSol);


  return(ier);
}
