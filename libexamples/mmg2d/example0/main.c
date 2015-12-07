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
  int             ier,k;

  fprintf(stdout,"  -- TEST MMG2DLIB \n");

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

  /** Manually set of the mesh */
  /** a) give the size of the mesh: 4 vertices, 2 triangles, 4 edges */
  /* allocation */
  if ( !MMG5_Set_meshSize(mmgMesh,4,2,5) )  exit(EXIT_FAILURE);

/** b) give the vertices: for each vertex, give the coordinates, the reference
      and the position in mesh of the vertex */
  if ( !MMG5_Set_vertex(mmgMesh,0  ,0  ,0  ,  1) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,1  ,0  ,0  ,  2) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,1  ,1  ,0  ,  3) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0  ,1  ,0  ,  4) )  exit(EXIT_FAILURE);

 /** c) give the triangles: for each triangle,
      give the vertices index, the reference and the position of the triangle */
  if ( !MMG5_Set_triangle(mmgMesh,  1,  2,  4, 1, 1) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  2,  3,  4, 1, 2) )  exit(EXIT_FAILURE);
 
  /** d) give the edges (not mandatory): for each edge,
      give the vertices index, the reference and the position of the edge */
  if ( !MMG5_Set_edge(mmgMesh,  1,  2, 1, 1) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_edge(mmgMesh,  2,  3, 2, 2) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_edge(mmgMesh,  3,  4, 3, 3) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_edge(mmgMesh,  4,  1, 4, 4) )  exit(EXIT_FAILURE);
 

  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG5_loadSol function that will read a .sol(b)
      file formatted or manually set your sol using the MMG5_Set* functions */

  /** Manually set of the sol */
  /** a) give info for the sol structure: sol applied on vertex entities,
      number of vertices=4, the sol is scalar*/
  if ( !MMG5_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,4,MMG5_Scalar) )
    exit(EXIT_FAILURE);

  /** b) give solutions values and positions */
  for(k=1 ; k<=4 ; k++) {
    if ( !MMG5_Set_scalarSol(mmgSol,0.01,k) ) exit(EXIT_FAILURE);
  }

 /** 4) (not mandatory): check if the number of given entities match with mesh size */
  if ( !MMG5_Chk_meshData(mmgMesh,mmgSol) ) exit(EXIT_FAILURE);

  /*save init mesh*/
  MMG2_saveMesh(mmgMesh,"init.mesh");
  MMG2_saveSol(mmgMesh,mmgSol,"init.sol");
 
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

  /** 3) Free the MMG2D5 structures */
  MMG5_Free_all(mmgMesh,mmgSol);



  return(0);
}
