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
  MMG2D_Init_mesh(&mmgMesh,&mmgSol);

 /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG2D_Set* functions */

  /** Manually set of the mesh */
  /** a) give the size of the mesh: 4 vertices, 2 triangles, 4 edges */
  /* allocation */
  if ( MMG2D_Set_meshSize(mmgMesh,4,2,5) != 1 )  exit(EXIT_FAILURE);

/** b) give the vertices: for each vertex, give the coordinates, the reference
      and the position in mesh of the vertex */
  if ( MMG2D_Set_vertex(mmgMesh,0  ,0  ,0  ,  1) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_vertex(mmgMesh,1  ,0  ,0  ,  2) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_vertex(mmgMesh,1  ,1  ,0  ,  3) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_vertex(mmgMesh,0  ,1  ,0  ,  4) != 1 )  exit(EXIT_FAILURE);

 /** c) give the triangles: for each triangle,
      give the vertices index, the reference and the position of the triangle */
  if ( MMG2D_Set_triangle(mmgMesh,  1,  2,  4, 1, 1) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_triangle(mmgMesh,  2,  3,  4, 1, 2) != 1 )  exit(EXIT_FAILURE);

  /** d) give the edges (not mandatory): for each edge,
      give the vertices index, the reference and the position of the edge */
  if ( MMG2D_Set_edge(mmgMesh,  1,  2, 1, 1) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_edge(mmgMesh,  2,  3, 2, 2) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_edge(mmgMesh,  3,  4, 3, 3) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_edge(mmgMesh,  4,  1, 4, 4) != 1 )  exit(EXIT_FAILURE);

  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG2D_loadSol function that will read a .sol(b)
      file formatted or manually set your sol using the MMG2D_Set* functions */

  /** Manually set of the sol */
  /** a) give info for the sol structure: sol applied on vertex entities,
      number of vertices=4, the sol is scalar*/
  if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,4,MMG5_Scalar) != 1 )
    exit(EXIT_FAILURE);

  /** b) give solutions values and positions */
  for(k=1 ; k<=4 ; k++) {
    if ( MMG2D_Set_scalarSol(mmgSol,0.01,k) != 1 ) exit(EXIT_FAILURE);
  }

 /** 4) (not mandatory): check if the number of given entities match with mesh size */
  if ( MMG2D_Chk_meshData(mmgMesh,mmgSol) != 1 ) exit(EXIT_FAILURE);

  /*save init mesh*/
  if ( MMG2D_saveMesh(mmgMesh,"init.mesh") != 1 )
    exit(EXIT_FAILURE);
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"init.sol") != 1 )
    exit(EXIT_FAILURE);

  ier = MMG2D_mmg2dlib(mmgMesh,mmgSol,NULL);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /*save result*/
  if ( MMG2D_saveMesh(mmgMesh,"result.mesh") != 1 )
    exit(EXIT_FAILURE);

  /*save metric*/
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"result") != 1 )
    exit(EXIT_FAILURE);

  /** 3) Free the MMG2D5 structures */
  MMG2D_Free_all(mmgMesh,mmgSol);



  return(0);
}
