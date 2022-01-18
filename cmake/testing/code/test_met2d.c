/**
 * Test matrix algebra in 2D.
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

/** Include the mmg2d library header file */
// if the header file is in the "include" directory
// #include "libmmg2d.h"
// if the header file is in "include/mmg/mmg2d"
#include "mmg/mmg2d/libmmg2d.h"
#include "mmgcommon.h"

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;

  fprintf(stdout,"  -- TEST MATRIX ALGEBRA IN 2D \n");

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
  MMG2D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /* Set inoffensive hmin and hmax for metric intersection */
  if( !MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmin,1.e-6) )
    return EXIT_FAILURE;
  if( !MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmax,1.e+6) )
    return EXIT_FAILURE;


  /** ------------------------------ STEP  II -------------------------- */

  if( !MMG5_test_invmat22() )
    return(EXIT_FAILURE);

  if( !MMG5_test_eigenvmatsym2d() )
    return EXIT_FAILURE;

  if( !MMG5_test_eigenvmatnonsym2d() )
    return EXIT_FAILURE;

  if( !MMG5_test_simred2d() )
    return EXIT_FAILURE;

  if( !MMG5_test_updatemet2d_ani() )
    return EXIT_FAILURE;

  if( !MMG5_test_intersecmet22(mmgMesh) )
    return EXIT_FAILURE;

  /** ------------------------------ STEP III -------------------------- */

  /** 3) Free the MMG2D structures */
  MMG2D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  return(0);
}
