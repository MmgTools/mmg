/**
 * Test matrix algebra in 3D.
 *
 * \author Luca Cirrottola (Inria)
 * \version 1
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

/** Include the mmg3d library header file */
// if the header file is in the "include" directory
// #include "libmmg3d.h"
// if the header file is in "include/mmg/mmg3d"
#include "mmg/mmg3d/libmmg3d.h"
#include "mmgcommon.h"

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;

  fprintf(stdout,"  -- TEST MATRIX ALGEBRA IN 3D \n");

  if ( argc != 1 ) {
    printf(" Usage: %s\n",argv[0]);
    return(1);
  }


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
  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);


  /** ------------------------------ STEP  II -------------------------- */

  if( !MMG5_test_invmat33() )
    return(EXIT_FAILURE);

  if( !MMG5_test_eigenvmatsym3d() )
    return(EXIT_FAILURE);

  if( !MMG5_test_eigenvmatnonsym3d() )
    return(EXIT_FAILURE);

  /** ------------------------------ STEP III -------------------------- */

  /** 3) Free the MMG3D structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  return(0);
}
