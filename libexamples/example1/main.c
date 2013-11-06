/** Authors Cécile Dobrzynski */
/* Example for using mmg3dlib */

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
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  int             opt_i[10],k,ier;
  double          opt_d[6];

  fprintf(stdout,"  -- TEST MMG3DLIB \n");

  mmgMesh = (MMG5_pMesh)calloc(1,sizeof(MMG5_Mesh));
  assert(mmgMesh);

  /* allocation */
  mmgMesh->np     = 12;
  mmgMesh->nt     = 0;
  mmgMesh->ne     = 12;


  mmgMesh->npmax  = 500000;
  mmgMesh->ntmax  = 1000000;
  mmgMesh->nemax  = 3000000;

  mmgMesh->point = (MMG5_pPoint)calloc(mmgMesh->npmax+1,sizeof(MMG5_Point));
  assert(mmgMesh->point);
  mmgMesh->tetra = (MMG5_pTetra)calloc(mmgMesh->nemax+1,sizeof(MMG5_Tetra));
  assert(mmgMesh->tetra);
  mmgMesh->tria  = (MMG5_pTria)calloc(mmgMesh->ntmax+1,sizeof(MMG5_Tria));
  assert(mmgMesh->tria);

  /*coordinates vertices*/
  mmgMesh->point[1].c[0]  = 0.;  mmgMesh->point[1].c[1]  = 0;
  mmgMesh->point[1].c[2]  = 0;
  mmgMesh->point[2].c[0]  = 0.5; mmgMesh->point[2].c[1]  = 0;
  mmgMesh->point[2].c[2]  = 0;
  mmgMesh->point[3].c[0]  = 0.5; mmgMesh->point[3].c[1]  = 0;
  mmgMesh->point[3].c[2]  = 1;
  mmgMesh->point[4].c[0]  = 0;   mmgMesh->point[4].c[1]  = 0;
  mmgMesh->point[4].c[2]  = 1;
  mmgMesh->point[5].c[0]  = 0;   mmgMesh->point[5].c[1]  = 1;
  mmgMesh->point[5].c[2]  = 0;
  mmgMesh->point[6].c[0]  = 0.5; mmgMesh->point[6].c[1]  = 1;
  mmgMesh->point[6].c[2]  = 0;
  mmgMesh->point[7].c[0]  = 0.5; mmgMesh->point[7].c[1]  = 1;
  mmgMesh->point[7].c[2]  = 1;
  mmgMesh->point[8].c[0]  = 0;   mmgMesh->point[8].c[1]  = 1;
  mmgMesh->point[8].c[2]  = 1;
  mmgMesh->point[9].c[0]  = 1;   mmgMesh->point[9].c[1]  = 0;
  mmgMesh->point[9].c[2]  = 0;
  mmgMesh->point[10].c[0] = 1;   mmgMesh->point[10].c[1] = 1;
  mmgMesh->point[10].c[2] = 0;
  mmgMesh->point[11].c[0] = 1;   mmgMesh->point[11].c[1] = 0;
  mmgMesh->point[11].c[2] = 1;
  mmgMesh->point[12].c[0] = 1;   mmgMesh->point[12].c[1] = 1;
  mmgMesh->point[12].c[2] = 1;

  mmgMesh->point[1].ref  = 0;
  mmgMesh->point[2].ref  = 0;
  mmgMesh->point[3].ref  = 0;
  mmgMesh->point[4].ref  = 0;
  mmgMesh->point[5].ref  = 0;
  mmgMesh->point[6].ref  = 0;
  mmgMesh->point[7].ref  = 0;
  mmgMesh->point[8].ref  = 0;
  mmgMesh->point[9].ref  = 0;
  mmgMesh->point[10].ref = 0;
  mmgMesh->point[11].ref = 0;
  mmgMesh->point[12].ref = 0;

  /* tria */
  /* not mandatory */
  mmgMesh->nt = 20;
  mmgMesh->tria[ 1].v[0] =  1; mmgMesh->tria[ 1].v[1] = 4;
  mmgMesh->tria[ 1].v[2] =  8; mmgMesh->tria[ 1].ref  = 2;
  mmgMesh->tria[ 2].v[0] =  1; mmgMesh->tria[ 2].v[1] = 2;
  mmgMesh->tria[ 2].v[2] =  4; mmgMesh->tria[ 2].ref  = 2;
  mmgMesh->tria[ 3].v[0] =  8; mmgMesh->tria[ 3].v[1] = 3;
  mmgMesh->tria[ 3].v[2] =  7; mmgMesh->tria[ 3].ref  = 0;
  mmgMesh->tria[ 4].v[0] =  5; mmgMesh->tria[ 4].v[1] = 8;
  mmgMesh->tria[ 4].v[2] =  6; mmgMesh->tria[ 4].ref  = 0;
  mmgMesh->tria[ 5].v[0] =  5; mmgMesh->tria[ 5].v[1] = 6;
  mmgMesh->tria[ 5].v[2] =  2; mmgMesh->tria[ 5].ref  = 0;
  mmgMesh->tria[ 6].v[0] =  5; mmgMesh->tria[ 6].v[1] = 2;
  mmgMesh->tria[ 6].v[2] =  1; mmgMesh->tria[ 6].ref  = 1;
  mmgMesh->tria[ 7].v[0] =  5; mmgMesh->tria[ 7].v[1] = 1;
  mmgMesh->tria[ 7].v[2] =  8; mmgMesh->tria[ 7].ref  = 0;
  mmgMesh->tria[ 8].v[0] =  7; mmgMesh->tria[ 8].v[1] = 6;
  mmgMesh->tria[ 8].v[2] =  8; mmgMesh->tria[ 8].ref  = 0;
  mmgMesh->tria[ 9].v[0] =  4; mmgMesh->tria[ 9].v[1] = 3;
  mmgMesh->tria[ 9].v[2] =  8; mmgMesh->tria[ 9].ref  = 0;
  mmgMesh->tria[10].v[0] =  2; mmgMesh->tria[10].v[1] = 3;
  mmgMesh->tria[10].v[2] =  4; mmgMesh->tria[10].ref  = 0;
  mmgMesh->tria[11].v[0] =  9; mmgMesh->tria[11].v[1] = 3;
  mmgMesh->tria[11].v[2] =  2; mmgMesh->tria[11].ref  = 0;
  mmgMesh->tria[12].v[0] = 11; mmgMesh->tria[12].v[1] = 9;
  mmgMesh->tria[12].v[2] = 12; mmgMesh->tria[12].ref  = 0;
  mmgMesh->tria[13].v[0] =  7; mmgMesh->tria[13].v[1] = 11;
  mmgMesh->tria[13].v[2] = 12; mmgMesh->tria[13].ref  = 0;
  mmgMesh->tria[14].v[0] =  6; mmgMesh->tria[14].v[1] = 7;
  mmgMesh->tria[14].v[2] = 10; mmgMesh->tria[14].ref  = 0;
  mmgMesh->tria[15].v[0] =  6; mmgMesh->tria[15].v[1] = 10;
  mmgMesh->tria[15].v[2] =  9; mmgMesh->tria[15].ref  = 0;
  mmgMesh->tria[16].v[0] =  6; mmgMesh->tria[16].v[1] = 9;
  mmgMesh->tria[16].v[2] =  2; mmgMesh->tria[16].ref  = 0;
  mmgMesh->tria[17].v[0] = 12; mmgMesh->tria[17].v[1] = 10;
  mmgMesh->tria[17].v[2] =  7; mmgMesh->tria[17].ref  = 0;
  mmgMesh->tria[18].v[0] = 12; mmgMesh->tria[18].v[1] = 9;
  mmgMesh->tria[18].v[2] = 10; mmgMesh->tria[18].ref  = 0;
  mmgMesh->tria[19].v[0] =  3; mmgMesh->tria[19].v[1] = 11;
  mmgMesh->tria[19].v[2] =  7; mmgMesh->tria[19].ref  = 0;
  mmgMesh->tria[20].v[0] =  9; mmgMesh->tria[20].v[1] = 11;
  mmgMesh->tria[20].v[2] =  3; mmgMesh->tria[20].ref  = 0;

  /* tetra*/
  /* warning: here we suppose that tetras are positively oriented */
  mmgMesh->tetra[1].v[0]  = 1;  mmgMesh->tetra[1].v[1]  = 4;
  mmgMesh->tetra[1].v[2]  = 2;  mmgMesh->tetra[1].v[3]  = 8;
  mmgMesh->tetra[2].v[0]  = 8;  mmgMesh->tetra[2].v[1]  = 3;
  mmgMesh->tetra[2].v[2]  = 2;  mmgMesh->tetra[2].v[3]  = 7;
  mmgMesh->tetra[3].v[0]  = 5;  mmgMesh->tetra[3].v[1]  = 2;
  mmgMesh->tetra[3].v[2]  = 6;  mmgMesh->tetra[3].v[3]  = 8;
  mmgMesh->tetra[4].v[0]  = 5;  mmgMesh->tetra[4].v[1]  = 8;
  mmgMesh->tetra[4].v[2]  = 1;  mmgMesh->tetra[4].v[3]  = 2;
  mmgMesh->tetra[5].v[0]  = 7;  mmgMesh->tetra[5].v[1]  = 2;
  mmgMesh->tetra[5].v[2]  = 8;  mmgMesh->tetra[5].v[3]  = 6;
  mmgMesh->tetra[6].v[0]  = 2;  mmgMesh->tetra[6].v[1]  = 4;
  mmgMesh->tetra[6].v[2]  = 3;  mmgMesh->tetra[6].v[3]  = 8;
  mmgMesh->tetra[7].v[0]  = 9;  mmgMesh->tetra[7].v[1]  = 2;
  mmgMesh->tetra[7].v[2]  = 3;  mmgMesh->tetra[7].v[3]  = 7;
  mmgMesh->tetra[8].v[0]  = 7;  mmgMesh->tetra[8].v[1]  = 11;
  mmgMesh->tetra[8].v[2]  = 9;  mmgMesh->tetra[8].v[3]  = 12;
  mmgMesh->tetra[9].v[0]  = 6;  mmgMesh->tetra[9].v[1]  = 9;
  mmgMesh->tetra[9].v[2]  = 10; mmgMesh->tetra[9].v[3]  = 7;
  mmgMesh->tetra[10].v[0] = 6;  mmgMesh->tetra[10].v[1] = 7;
  mmgMesh->tetra[10].v[2] = 2;  mmgMesh->tetra[10].v[3] = 9;
  mmgMesh->tetra[11].v[0] = 12; mmgMesh->tetra[11].v[1] = 9;
  mmgMesh->tetra[11].v[2] = 7;  mmgMesh->tetra[11].v[3] = 10;
  mmgMesh->tetra[12].v[0] = 9;  mmgMesh->tetra[12].v[1] = 3;
  mmgMesh->tetra[12].v[2] = 11; mmgMesh->tetra[12].v[3] = 7;

  mmgMesh->tetra[1].ref  = 0;
  mmgMesh->tetra[2].ref  = 0;
  mmgMesh->tetra[3].ref  = 0;
  mmgMesh->tetra[4].ref  = 0;
  mmgMesh->tetra[5].ref  = 0;
  mmgMesh->tetra[6].ref  = 0;
  mmgMesh->tetra[7].ref  = 0;
  mmgMesh->tetra[8].ref  = 0;
  mmgMesh->tetra[9].ref  = 0;
  mmgMesh->tetra[10].ref = 0;
  mmgMesh->tetra[11].ref = 0;
  mmgMesh->tetra[12].ref = 0;

  /*metric*/
  mmgSol           = (MMG5_pSol)calloc(1,sizeof(MMG5_Sol));
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

  /** first wave of refinment with a detection of angles (between normals */
  /*  at 2 adjacent surfaces) smallest than 90 and a maximal size of 0.2 */
  opt_i[MMG5_imprim]   = 5;
  opt_d[MMG5_dhd]      = 90;
  opt_d[MMG5_hmax]     = 0.2;

  mmgMesh->nameout = "result0.mesh";
  mmgSol->nameout = "result0.sol";

  ier = MMG5_mmg3dlib(opt_i,opt_d,mmgMesh,mmgSol);
  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");


  /** Second wave of refinment with a smallest maximal size */
  opt_d[MMG5_hmax]  = 0.1;

  ier = MMG5_mmg3dlib(opt_i,opt_d,mmgMesh,mmgSol);
  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /*save result*/
  MMG5_saveMesh(mmgMesh);

  /*save metric*/
  MMG5_saveMet(mmgMesh,mmgSol);

  /* free mem */
  MMG5_freeAll(mmgMesh,mmgSol);
  free(mmgSol);
  free(mmgMesh);

  return(ier);
}
