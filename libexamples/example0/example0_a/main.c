/** Authors Cécile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly */
/** \include Example for using mmg3dlib (basic use) */

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
  int             ier;
  char            *pwd,*filename;

  fprintf(stdout,"  -- TEST MMG3DLIB \n");

  /* Name and path of the mesh file */
  pwd = getenv("PWD");
  filename = (char *) calloc(strlen(pwd) + 42, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  sprintf(filename, "%s%s%s", pwd, "/../libexamples/example0/example0_a/", "cube");

  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh: mesh=&mmgMesh, sol=&mmgSol, input mesh name, input sol name,
     output mesh name */
  mmgMesh = NULL;
  mmgSol  = NULL;
  MMG5_Init_mesh(&mmgMesh,&mmgSol);

  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG5_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG5_Set* functions */

  /** with MMG5_loadMesh function */
  /** a) (not mandatory): give the mesh name
     (by default, the "mesh.mesh" file is oppened)*/
  if ( !MMG5_Set_inputMeshName(mmgMesh,filename) )
    exit(EXIT_FAILURE);
  /** b) function calling */
  if ( !MMG5_loadMesh(mmgMesh) )  exit(EXIT_FAILURE);

  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG5_loadMet function that will read a .sol(b)
      file formatted or manually set your sol using the MMG5_Set* functions */

  /** With MMG5_loadMet function */
  /** a) (not mandatory): give the sol name
     (by default, the "mesh.sol" file is oppened)*/
  if ( !MMG5_Set_inputSolName(mmgMesh,mmgSol,filename) )
    exit(EXIT_FAILURE);
  /** b) function calling */
  if ( !MMG5_loadMet(mmgMesh,mmgSol) )
    exit(EXIT_FAILURE);

  /** 4) (not mandatory): check if the number of given entities match with mesh size */
  if ( !MMG5_Chk_meshData(mmgMesh,mmgSol) ) exit(EXIT_FAILURE);

  /** ------------------------------ STEP  II -------------------------- */
  /** library call */
  ier = MMG5_mmg3dlib(mmgMesh,mmgSol);
  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");


  /** ------------------------------ STEP III -------------------------- */
  /** get results */
  /** Two solutions: just use the MMG5_saveMesh/MMG5_saveMet functions
      that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
      using the MMG5_getMesh/MMG5_getSol functions */

  /** 1) Automatically save the mesh */
  /** a)  (not mandatory): give the ouptut mesh name using MMG5_Set_outputMeshName
     (by default, the mesh is saved in the "mesh.o.mesh" file */
  // MMG5_Set_outputMeshName(mmgMesh,"output.mesh");
  /** b) function calling */
  MMG5_saveMesh(mmgMesh);

  /** 2) Automatically save the solution */
  /** a)  (not mandatory): give the ouptut sol name using MMG5_Set_outputSolName
     (by default, the mesh is saved in the "mesh.o.sol" file */
  // MMG5_Set_outputSolName(mmgSol,"output.sol");
  /** b) function calling */
  MMG5_saveMet(mmgMesh,mmgSol);

  /** 3) Free the MMG3D5 structures */
  MMG5_Free_all(mmgMesh,mmgSol);

  free(filename);
  filename = NULL;

  return(ier);
}
