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

  /* Set inoffensive hmin and hmax for metric intersection */
  if( !MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hmin,1.e-6) )
    return EXIT_FAILURE;
  if( !MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hmax,1.e+6) )
    return EXIT_FAILURE;


  /** ------------------------------ STEP  II -------------------------- */

  /* matrix inversion test */
  if( !MMG5_test_invmat33() )
    return(EXIT_FAILURE);

  /* symmetric matrix eigendecomposition test */
  double m_sym[2][6] = {{2.,0.,0.,3.,4.,9.},
                        {1.,0.,0.,50.5,49.5,50.5}}; /* Test matrices */
  double lambda_sym[2][3] = {{1.,2.,11.},
                             {1.,1.,100.}}; /* Exact eigenvalues */
  double vp_sym[2][3][3] = {{{0.,-2./sqrt(5.),1./sqrt(5.)},
                             {1.,0.,0.},
                             {0.,1./sqrt(5.),2./sqrt(5.)}},
                            {{1.,0.,0.},
                             {0.,1./sqrt(2.),-1./sqrt(2.)},
                             {0.,1./sqrt(2.),1./sqrt(2.)}}}; /* Exact eigenvectors */
  for( int8_t i = 0; i < 2; i++ )
    if( !MMG5_test_eigenvmatsym3d(mmgMesh,m_sym[i],lambda_sym[i],vp_sym[i]) )
      return(EXIT_FAILURE);

  /* non-symmetric matrix eigendecomposition test */
  double m_nonsym[2][9] = {{500.5,-499.5,499.5,
                            -49.5,  50.5, 49.5,
                            450., -450., 550.},
                           {50.5,-49.5,49.5,
                             0.,   1.,  0.,
                            49.5,-49.5,50.5}}; /* Test matrices */
  double lambda_nonsym[2][3] = {{1.,100.,1000.},
                                {1.,  1., 100.}}; /* Exact eigenvalues */
  double vp_nonsym[3][3] = {{1./sqrt(2.),1./sqrt(2.),0.},
                            {0.,         1./sqrt(2.),1./sqrt(2.)},
                            {1./sqrt(2.),         0.,1./sqrt(2.)}}; /* Exact right eigenvectors */
  double ivp_nonsym[3][3] = {{ 1./sqrt(2.),-1./sqrt(2.), 1./sqrt(2.)},
                             { 1./sqrt(2.), 1./sqrt(2.),-1./sqrt(2.)},
                             {-1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.)}}; /* Exact right eigenvectors inverse */
  for( int8_t i = 0; i < 2; i++ )
    if( !MMG5_test_eigenvmatnonsym3d(mmgMesh,m_nonsym[i],lambda_nonsym[i],vp_nonsym,ivp_nonsym) )
      return(EXIT_FAILURE);

  /* symmetric matrix multiplication test */
  if( !MMG5_test_mn() )
    return(EXIT_FAILURE);

  /* matrix linear transformation test */
  if( !MMG5_test_rmtr() )
    return(EXIT_FAILURE);

  /* rotation matrix test */
  if( !MMG5_test_rotmatrix() )
    return(EXIT_FAILURE);

  /* simultaneous reduction test */
  if( !MMG5_test_simred3d() )
    return(EXIT_FAILURE);

  /* matrix inverse transformation test */
  if( !MMG5_test_updatemet3d_ani() )
    return EXIT_FAILURE;

  /* metrics intersection test */
  if( !MMG5_test_intersecmet33(mmgMesh) )
    return EXIT_FAILURE;

  /** ------------------------------ STEP III -------------------------- */

  /** 3) Free the MMG3D structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  return(0);
}
