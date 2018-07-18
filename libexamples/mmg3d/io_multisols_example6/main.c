/**
 * Example of input output for the mmg3d library for multiple solutions at mesh
 * vertices
 *
 * \author Algiane Froehly (InriaSoft)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/** Include the mmg3d library hader file */
// if the header file is in the "include" directory
// #include "libmmg3d.h"
// if the header file is in "include/mmg/mmg3d"
#include "mmg/mmg3d/libmmg3d.h"

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol,tmpSol;
  int             i;

  /* To manually recover the mesh */
  int             nsol,np,typEntity[MMG5_NSOL_MAX],typSol[MMG5_NSOL_MAX];
  double          *sols;

  /* Filenames */
  char            *filename, *fileout;

  fprintf(stdout,"  -- TEST MMG3DLIB \n");

  if ( argc != 3 ) {
    printf(" Usage: %s filein fileout\n",argv[0]);
    return(1);
  }

  /* Name and path of the mesh file */
  filename = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filename,argv[1]);

  fileout = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[2]);

  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh:
   * MMG5_ARG_start: we start to give the args of a variadic func
   * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
   * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
   * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
   * &mmgSol: pointer toward your MMG5_pSol (that store your metric) */

  mmgMesh = NULL;
  mmgSol  = NULL;
  tmpSol  = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /** 2) Build initial mesh and solutions in MMG5 format */
  /** Two solutions: just use the MMG3D_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG3D_Set* functions */

  /** Automatic loading of the mesh and multiple solutions */
  if ( MMG3D_loadMesh(mmgMesh,filename) != 1 )  exit(EXIT_FAILURE);

  if ( MMG3D_loadAllSols(mmgMesh,&mmgSol,filename) != 1 )
    exit(EXIT_FAILURE);

  /** ------------------------------ STEP II --------------------------- */

  /** 3) Transfer the solutions in a new solutions array */
  /** a) Get the solutions sizes */
  if ( MMG3D_Get_allSolsSizes(mmgMesh,&mmgSol,&nsol,typEntity,&np,typSol) != 1 )
    exit(EXIT_FAILURE);

  /** b) Manually set the size of the new solution: give info for the sol
      structure: number of solutions, type of entities on which applied the
      solutions, number of vertices, type of the solution */
  if ( MMG3D_Set_allSolsSizes(mmgMesh,&tmpSol,nsol,typEntity,np,typSol) != 1 )
    exit(EXIT_FAILURE);

  /** c) Get each solution and set it in the new structure */

  /** b) give solutions values and positions */
  /* Get the entire field of a given solution */
  for ( i=0; i<nsol; ++i ) {
    assert ( typEntity[i] == MMG5_Vertex );

    /* Get the ith solution array */
    if ( typSol[i] == MMG5_Scalar )
      sols = (double*) calloc(np, sizeof(double));
    else if ( typSol[i] == MMG5_Vector )
      sols = (double*) calloc(np*3, sizeof(double));
    else if ( typSol[i] == MMG5_Tensor ) {
      sols = (double*) calloc(np*6, sizeof(double));
    }

    if ( MMG3D_Get_ithSols_inAllSols(mmgSol,i,sols) !=1 ) exit(EXIT_FAILURE);

    /* Set the ith solution in the new structure */
    if ( MMG3D_Set_ithSols_inAllSols(tmpSol,i,sols) !=1 ) exit(EXIT_FAILURE);

    free(sols); sols = NULL;
  }


  /** ------------------------------ STEP III -------------------------- */
  /** Save the new data */
  /** Use the MMG3D_saveMesh/MMG3D_saveAllSols functions */
  /* save the mesh */
  if ( MMG3D_saveMesh(mmgMesh,fileout) != 1 )
    exit(EXIT_FAILURE);

  /*s ave the solutions array */
  if ( MMG3D_saveAllSols(mmgMesh,&tmpSol,fileout) != 1 )
    exit(EXIT_FAILURE);

  /** 3) Free the MMG3D structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppSols,&tmpSol,
                 MMG5_ARG_ppSols,&mmgSol,
                 MMG5_ARG_end);

  free(filename);
  filename = NULL;

  free(fileout);
  fileout = NULL;

  return(0);
}
