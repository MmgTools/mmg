/**
 * Example of input output for the mmg2d library for multiple solutions at mesh
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

/** Include the mmg2d library hader file */
// if the header file is in the "include" directory
// #include "libmmg2d.h"
// if the header file is in "include/mmg/mmg2d"
#include "mmg/mmg2d/libmmg2d.h"

#define MAX0(a,b)     (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX0(a,b)) > (MAX0(c,d))) ? (MAX0(a,b)) : (MAX0(c,d)))

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       *mmgSol,*tmpSol;
  int             ier,k,i;

  /* To manually recover the mesh */
  int             nsol,np,typEntity[MMG5_NSOL_MAX],typSol[MMG5_NSOL_MAX];
  double          *sols;
  double          scalarSol,vectorialSol[2];
  char            *fileout,*solout;

  fprintf(stdout,"  -- TEST MMG2DLIB \n");

  if ( argc != 2 ) {
    printf(" Usage: %s fileout\n",argv[0]);
    return(1);
  }

  /* Name and path of the mesh file */
  fileout = (char *) calloc(strlen(argv[1]) + 6, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[1]);
  strcat(fileout,".mesh");

  solout = (char *) calloc(strlen(argv[1]) + 5, sizeof(char));
  if ( solout == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(solout,argv[1]);
  strcat(solout,".sol");

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
  MMG2D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /** 2) Build initial mesh and solutions in MMG5 format */
  /** Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG2D_Set* functions */

  /** Automaric loading of the mesh and multiple solutions */
  if ( MMG2D_loadMesh(mmgMesh,filename) != 1 )  exit(EXIT_FAILURE);

  if ( MMG2D_loadAllSols(mmgMesh,mmgSol,filename) != 1 )
    exit(EXIT_FAILURE);

  /** ------------------------------ STEP II --------------------------- */

  /** 3) Transfer the solutions in a new solutions array */
  /** a) Get the solutions sizes */
  if ( MMG2D_Get_allSolsSizes(mmgMesh,mmgSol,&nsol,typEntity,&np,typSol) != 1 )
    exit(EXIT_FAILURE);

  /** b) Manually set the size of the new solution: give info for the sol
      structure: number of solutions, type of entities on which applied the
      solutions, number of vertices, type of the solution */
  if ( MMG2D_Set_allSolsSizes(mmgMesh,tmpSol,nsol,typEntity,np,typSol) != 1 )
    exit(EXIT_FAILURE);

  /** c) Get each solution and set it in the new structure */

  /** b) give solutions values and positions */
  /* Get the entire field of a given solution */
  for ( i=0; i<*nsol; ++i ) {
    assert ( typEntity[i] == MMG5_Vertex );

    /* Get the ith solution array */
    if ( typSol[i] = MMG5_Scalar )
      sols = (double*) calloc(np, sizeof(double));
    else if ( typSol[i] = MMG5_Vector )
      sols = (double*) calloc(np*2, sizeof(double));
    else if ( typSol[i] = MMG5_Tensor )
      sols = (double*) calloc(np*3, sizeof(double));

    if ( MMG2D_Get_ithSols_inAllSols(mmgSol,i,sols) !=1 ) exit(EXIT_FAILURE);

    /* Set the ith solution in the new structure */
    if ( MMG2D_Set_ithSols_inAllSols(tmpSol,i,sols) !=1 ) exit(EXIT_FAILURE);
  }


  /** ------------------------------ STEP III -------------------------- */
  /** Save the new data */
  /** Use the MMG2D_saveMesh/MMG2D_saveAllSols functions */
  /* save the mesh */
  if ( MMG2D_saveMesh(mmgMesh,fileout) != 1 )
    exit(EXIT_FAILURE);

  /*s ave the solutions array */
  if ( MMG2D_saveAllSols(mmgMesh,tmpSol,fileout) != 1 )
    exit(EXIT_FAILURE);

  /** 3) Free the MMG2D structures */
  MMG2D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  free(fileout);
  fileout = NULL;

  free(solout);
  solout = NULL;

  return(0);
}
