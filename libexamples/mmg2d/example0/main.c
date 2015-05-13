/*Authors Cécile Dobrzynski

  Example for using mmg2dlib

*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#include "mmg2d.h"
void tryinsert(int num,int t1,int t2,int t3,int t4) {

  //printf("nouveau point %d %d %d %d %d\n",num,t1,t2,t3,t4);
  return;
}

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  double          optdbl[2];
  int             opt[7],k;

  fprintf(stdout,"  -- TEST MMG2DLIB \n");

  mmgMesh = (MMG5_pMesh)calloc(1,sizeof(MMG5_Mesh));
  assert(mmgMesh);

  /* allocation */
  mmgMesh->np      = 4; /*nb of nodes*/
  mmgMesh->nt      = 2; /*nb of triangles*/
  mmgMesh->na      = 4; /*nb of edges*/


  mmgMesh->npmax   = 500000; /*nb maximal of nodes*/
  mmgMesh->ntmax   = 1000000;/*nb maximal of triangles*/
  mmgMesh->namax   = 3000000;/*nb maximal of edges*/

  mmgMesh->point = (MMG5_pPoint)calloc(mmgMesh->npmax+1,sizeof(MMG5_Point));
  assert(mmgMesh->point);
  mmgMesh->tria  = (MMG5_pTria)calloc(mmgMesh->ntmax+1,sizeof(MMG5_Tria));
  assert(mmgMesh->tria);
  mmgMesh->edge  = (MMG5_pEdge)calloc(mmgMesh->namax+1,sizeof(MMG5_Edge));
  assert(mmgMesh->edge);
  mmgMesh->adja = (int*)calloc(3*mmgMesh->ntmax+3,sizeof(int));
  assert(mmgMesh->adja);

  /*coordinates vertices*/
  mmgMesh->point[1].c[0]  = 0.;  mmgMesh->point[1].c[1]  = 0.; mmgMesh->point[1].ref  = 0;
  mmgMesh->point[2].c[0]  = 1.;  mmgMesh->point[2].c[1]  = 0.;  mmgMesh->point[2].ref  = 0;
  mmgMesh->point[3].c[0]  = 1.;  mmgMesh->point[3].c[1]  = 1.;  mmgMesh->point[3].ref  = 0;
  mmgMesh->point[4].c[0]  = 0.;  mmgMesh->point[4].c[1]  = 1.;  mmgMesh->point[4].ref  = 0;

  /*triangles*/
  mmgMesh->tria[1].v[0]  = 1;  mmgMesh->tria[1].v[1]  = 2;  mmgMesh->tria[1].v[2]  = 4; mmgMesh->tria[1].ref  = 1;
  mmgMesh->tria[2].v[0]  = 2;  mmgMesh->tria[2].v[1]  = 3;  mmgMesh->tria[2].v[2]  = 4; mmgMesh->tria[2].ref  = 1;

  /*edges*/
  mmgMesh->edge[1].a  = 1;  mmgMesh->edge[1].b  = 2;  mmgMesh->edge[1].ref  = 1;
  mmgMesh->edge[2].a  = 2;  mmgMesh->edge[2].b  = 3;  mmgMesh->edge[2].ref  = 2;
  mmgMesh->edge[3].a  = 3;  mmgMesh->edge[3].b  = 4;  mmgMesh->edge[3].ref  = 3;
  mmgMesh->edge[4].a  = 4;  mmgMesh->edge[4].b  = 1;  mmgMesh->edge[4].ref  = 4;

  /*save init mesh*/
  MMG2_saveMesh(mmgMesh,"init.mesh");

  /*metric*/
  mmgSol           = (MMG5_pSol)calloc(1,sizeof(MMG5_Sol));
  assert(mmgSol);
  mmgSol->size = 1; /*1 : isotropic size*/

  /*scalaire size*/
  mmgSol->np = mmgMesh->np;
  mmgSol->m    = (double*)calloc(mmgMesh->npmax+1,mmgSol->size*sizeof(double));
  assert(mmgSol->m);
  for(k=1 ; k<=mmgMesh->np ; k++) {
    mmgSol->m[k] = 0.01;
  }

  /*save init size*/
  //MMG2_saveSol(mmgMesh,mmgSol,"init");

  opt[0]=1; //adaptation
  opt[1]=0; //no debug
  opt[2]=1;//noswap OFF
  opt[3]=0;//noinsert OFF
  opt[4]=1;//nomove OFF
  opt[5]=5; //imprim
  opt[6]=0;  //no ridge OFF
  optdbl[0] = -1; //gradation
  optdbl[1] = 30; //angle for ridge detection

  if(MMG2_mmg2dlib(opt,optdbl,mmgMesh,mmgSol,(void*)tryinsert)) {
    fprintf(stdout,"BAD ENDING OF MMG2DLIB\n");
  }

  /*save result*/
  MMG2_saveMesh(mmgMesh,"result.mesh");

  /*save metric*/
  MMG2_saveSol(mmgMesh,mmgSol,"result");

  /* free mem */
  free(mmgMesh->point);
  free(mmgMesh->tria);
  free(mmgMesh->edge);
  free(mmgMesh->adja);

  free(mmgMesh);
  free(mmgSol);


  return(0);
}
