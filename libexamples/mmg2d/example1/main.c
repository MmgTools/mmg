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

/** Include the mmg2d library hader file */
// if the header file is in the "include" directory
#include "libmmg2d.h"
// if the header file is in "include/mmg/mmg2d"
// #include "mmg/mmg2d/libmmg2d.h"

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  char            *pwd,*filename;

  int             k,np,nt,na,ier;

  fprintf(stdout,"  -- TEST MMG2DLIB \n");

  /* Name and path of the mesh file */
  pwd = getenv("PWD");
  filename = (char *) calloc(strlen(pwd) + 47, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  sprintf(filename, "%s%s%s", pwd, "/../libexamples/mmg2d/example1/", "dom");

  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh: mesh=&mmgMesh, sol=&mmgSol, input mesh name, input sol name,
     output mesh name */
  mmgMesh = NULL;
  mmgSol  = NULL;
  MMG2D_Init_mesh(&mmgMesh,&mmgSol);

 /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG2D_Set* functions */

  /** with MMG2D_loadMesh function */
  /** a) (not mandatory): give the mesh name
     (by default, the "mesh.mesh" file is oppened)*/
  if ( MMG2D_Set_inputMeshName(mmgMesh,filename) != 1 )
    exit(EXIT_FAILURE);
  /** b) function calling */
  if ( MMG2D_loadMesh(mmgMesh,filename) != 1 )  exit(EXIT_FAILURE);


  /*save init mesh*/
  if ( MMG2D_saveMesh(mmgMesh,"init.mesh") != 1 )  exit(EXIT_FAILURE);

  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG2D_loadMet function that will read a .sol(b)
      file formatted or manually set your sol using the MMG2D_Set* functions */

  /** Manually set of the sol */
  /** a) Get np the number of vertex */
  if ( MMG2D_Get_meshSize(mmgMesh,&np,&nt,&na) != 1 )
    exit(EXIT_FAILURE);

  /** b) give info for the sol structure: sol applied on vertex entities,
      number of vertices=np, the sol is scalar*/
  if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,np,MMG5_Scalar) != 1 )
    exit(EXIT_FAILURE);

  /** c) give solutions values and positions */
  for(k=1 ; k<=np ; k++) {
    if ( MMG2D_Set_scalarSol(mmgSol,0.01,k) != 1 ) exit(EXIT_FAILURE);
  }

  /** 4) (not mandatory): check if the number of given entities match with mesh size */
  if ( MMG2D_Chk_meshData(mmgMesh,mmgSol) != 1 ) exit(EXIT_FAILURE);


  /*save init size*/
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"init") != 1 )  exit(EXIT_FAILURE);

  ier = MMG2D_mmg2dlib(mmgMesh,mmgSol,NULL);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /*save result*/
  if ( MMG2D_saveMesh(mmgMesh,"result.mesh") != 1 )  exit(EXIT_FAILURE);

  /*save metric*/
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"result") != 1 )  exit(EXIT_FAILURE);

  /** 5) Free the MMG3D5 structures */
  MMG2D_Free_all(mmgMesh,mmgSol);


  return(0);
}
