/** Authors Cécile Dobrzynski */
/** \include Example for using mmg3dlib */

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
  int             *opt_i,k,ier;
  double          *opt_d;

  fprintf(stdout,"  -- TEST MMG3DLIB \n");

  /** Step 1: Initialisation of mesh and sol structures */
  /* args of InitMesh: mesh=&mmgMesh, sol=&mmgSol, input mesh name, input sol name,
   output mesh name */
  MMG5_Init_mesh(&mmgMesh,&mmgSol);

  /** Step 2: Build mesh in MMG5 format */
  /** Two solutions: just use the MMG5_loadMesh function that will read a .mesh(b)
     file formatted or manually set your mesh using the MMG5_Set* functions */

  /** Manually set of the mesh */
  /** a) give the size of the mesh: 12 vertices, 12 tetra, 20 triangles, 0 edges */
  if ( !MMG5_Set_meshSize(mmgMesh,12,12,20,0) )  exit(EXIT_FAILURE);

  /** b) give the vertices: for each vertex, give the coordinates and the reference of vertex */
  if ( !MMG5_Set_vertex(mmgMesh,0  ,0  ,0  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0.5,0  ,0  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0.5,0  ,1  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0  ,0  ,1  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0  ,1  ,0  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0.5,1  ,0  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0.5,1  ,1  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0  ,1  ,1  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,1  ,0  ,0  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,1  ,1  ,0  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,1  ,0  ,1  ,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,1  ,1  ,1  ,0) )  exit(EXIT_FAILURE);

  /** c) give the tetrahedras: for each tetrahedra,
      give the vertices index and the reference of the tetra */
  /* warning: here we suppose that tetras are positively oriented */
  if ( !MMG5_Set_tetrahedra(mmgMesh,  1,  4,  2,  8,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh,  8,  3,  2,  7,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh,  5,  2,  6,  8,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh,  5,  8,  1,  2,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh,  7,  2,  8,  6,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh,  2,  4,  3,  8,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh,  9,  2,  3,  7,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh,  7, 11,  9, 12,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh,  6,  9, 10,  7,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh,  6,  7,  2,  9,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh, 12,  9,  7, 10,0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedra(mmgMesh,  9,  3, 11,  7,0) )  exit(EXIT_FAILURE);

  /** d) give the triangles (not mandatory): for each triangle,
      give the vertices index and the reference of the triangle */
  if ( !MMG5_Set_triangle(mmgMesh,  1,  4,  8, 2) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  1,  2,  4, 2) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  8,  3,  7, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  5,  8,  6, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  5,  6,  2, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  5,  2,  1, 1) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  5,  1,  8, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  7,  6,  8, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  4,  3,  8, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  2,  3,  4, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  9,  3,  2, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh, 11,  9, 12, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  7, 11, 12, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  6,  7, 10, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  6, 10,  9, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  6,  9,  2, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh, 12, 10,  7, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh, 12,  9, 10, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  3, 11,  7, 0) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  9, 11,  3, 0) )  exit(EXIT_FAILURE);


  /** Step 3: Build sol in MMG5 format */
  /** Two solutions: just use the MMG5_loadMet function that will read a .sol(b)
      file formatted or manually set your sol using the MMG5_Set* functions */

  /** Manually set of the sol */
  /** a) give info for the sol structure: sol applied on vertex entities,
       number of vertices=12, the sol is scalar*/
  if ( !MMG5_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,12,MMG5_Scalar) )
    exit(EXIT_FAILURE);

  /** b) give solutions values */
  for(k=1 ; k<=12 ; k++) {
    if ( !MMG5_Set_scalarSol(mmgSol,0.5) ) exit(EXIT_FAILURE);
  }

  /** Step 4 (not mandatory): check if the number of given entities match with mesh size */
  if ( !MMG5_Chk_meshData(mmgMesh,mmgSol) ) exit(EXIT_FAILURE);

  /** Step 5 (not mandatory): set your global parameters using the MMG5_Set_iparameters and
   MMG5_Set_dparameters function (resp. for integer parameters and double param)*/

  /* first wave of refinment with a detection of angles (between normals
   *  at 2 adjacent surfaces) smallest than 90 and a maximal size of 0.2 */
  if ( !MMG5_Set_iparameters(mmgMesh,MMG5_IPARAM_verbose, 5) )
    exit(EXIT_FAILURE);
  if ( !MMG5_Set_dparameters(mmgMesh,MMG5_DPARAM_angleDetection,90) )
    exit(EXIT_FAILURE);
  if ( !MMG5_Set_dparameters(mmgMesh,MMG5_DPARAM_hmax,0.2) )
    exit(EXIT_FAILURE);

  /** Step 6 (not mandatory): set your local parameters */
  /* use 2 local hausdorff numbers on ref 2 (hausd = 0.001) and 0 (hausd = 0.005) */
  if ( !MMG5_Set_iparameters(mmgMesh,MMG5_IPARAM_numberOfLocalParam,2) )
    exit(EXIT_FAILURE);

  /** for each local parameter: give the type and reference of the element on which you
      will apply a particular hausdorff number and the hausdorff number to apply */
  if ( !MMG5_Set_localParameters(mmgMesh,MMG5_Triangle,2,0.001) )
    exit(EXIT_FAILURE);
  if ( !MMG5_Set_localParameters(mmgMesh,MMG5_Triangle,0,0.005) )
    exit(EXIT_FAILURE);

  if ( !MMG5_Set_outputMeshName(mmgMesh,"result0.mesh") )
    exit(EXIT_FAILURE);
  if ( !MMG5_Set_outputSolName(mmgMesh,mmgSol,"result0.sol") )
    exit(EXIT_FAILURE);

  /** Step 7: library call */
  ier = MMG5_mmg3dlib(mmgMesh,mmgSol);
  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");


  /* Second wave of refinment with a smallest maximal size and with a different
     second local parameter */
   if ( !MMG5_Set_dparameters(mmgMesh,MMG5_DPARAM_hmax,0.1) )
     exit(EXIT_FAILURE);
  if ( !MMG5_Set_localParameters(mmgMesh,MMG5_Triangle,0,0.01) )
    exit(EXIT_FAILURE);

  ier = MMG5_mmg3dlib(mmgMesh,mmgSol);
  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /** Step 8: get results */
  /** Two solutions: just use the MMG5_saveMesh/MMG5_saveMet functions
      that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
      using the MMG5_getMesh/MMG5_getSol functions */

  /* Automatically save the mesh */
  MMG5_saveMesh(mmgMesh);

  /* Automatically save the solution */
  MMG5_saveMet(mmgMesh,mmgSol);

  /* Step 9: free the MMG3D5 structures */
  MMG5_Free_all(mmgMesh,mmgSol);

  return(ier);
}
