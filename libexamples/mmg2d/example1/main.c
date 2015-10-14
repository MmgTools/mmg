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
  MMG2_Init_mesh(&mmgMesh,&mmgSol);

 /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG5_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG5_Set* functions */

  /** with MMG5_loadMesh function */
  /** a) (not mandatory): give the mesh name
     (by default, the "mesh.mesh" file is oppened)*/
  if ( !MMG5_Set_inputMeshName(mmgMesh,filename) )
    exit(EXIT_FAILURE);
  /** b) function calling */
  if ( !MMG2_loadMesh(mmgMesh,filename) )  exit(EXIT_FAILURE);


  /*save init mesh*/
  MMG2_saveMesh(mmgMesh,"init.mesh");

  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG5_loadMet function that will read a .sol(b)
      file formatted or manually set your sol using the MMG5_Set* functions */
  
  /** Manually set of the sol */
  /** a) Get np the number of vertex */
  if (!MMG5_Get_meshSize(mmgMesh,&np,&nt,&na) )
    exit(EXIT_FAILURE);
  
  /** b) give info for the sol structure: sol applied on vertex entities,
      number of vertices=np, the sol is scalar*/
  if ( !MMG5_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,np,MMG5_Scalar) )
    exit(EXIT_FAILURE);

  /** c) give solutions values and positions */
  for(k=1 ; k<=np ; k++) {
    if ( !MMG5_Set_scalarSol(mmgSol,0.01,k) ) exit(EXIT_FAILURE);
  }

  /** 4) (not mandatory): check if the number of given entities match with mesh size */
  if ( !MMG5_Chk_meshData(mmgMesh,mmgSol) ) exit(EXIT_FAILURE);
  

  /*save init size*/
  MMG2_saveSol(mmgMesh,mmgSol,"init");

  ier = MMG2_mmg2dlib(mmgMesh,mmgSol,NULL);
  
  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /*save result*/
  MMG2_saveMesh(mmgMesh,"result.mesh");

  /*save metric*/
  MMG2_saveSol(mmgMesh,mmgSol,"result");

  /** 5) Free the MMG3D5 structures */
  MMG5_Free_all(mmgMesh,mmgSol);


  return(0);
}
