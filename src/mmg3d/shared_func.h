/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file mmg3d/shared_func.h
 * \brief Common functions between MMG3D5 library and executable.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

/* global variables */
unsigned char _MMG5_inxt2[3] = {1,2,0};
unsigned char _MMG5_iprv2[3] = {2,0,1};
unsigned char _MMG5_idir[4][3] = { {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1} };
char _MMG5_idirinv[4][4] = {{-1,0,1,2},{0,-1,2,1},{0,1,-1,2},{0,2,1,-1}};
unsigned char _MMG5_iarf[4][3] = { {5,4,3}, {5,1,2}, {4,2,0}, {3,0,1} };
unsigned char _MMG5_iarfinv[4][6] = { {-1,-1,-1,2,1,0}, {-1,1,2,-1,-1,0},{2,-1,1,-1,0,-1},{1,2,-1,0,-1,-1}};
unsigned char _MMG5_inxt3[7] = { 1,2,3,0,1,2,3 };
unsigned char _MMG5_iprv3[7] = { 3,0,1,2,3,0,1 };
unsigned char _MMG5_iare[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
unsigned char _MMG5_ifar[6][2] = { {2,3}, {1,3}, {1,2}, {0,3}, {0,2}, {0,1} };
unsigned char _MMG5_isar[6][2] = { {2,3}, {3,1}, {1,2}, {0,3}, {2,0}, {0,1} };
unsigned char _MMG5_arpt[4][3] = { {0,1,2}, {0,4,3}, {1,3,5}, {2,5,4} };

/* functions shared by executable and library versions of MMG3D5 */

#ifdef USE_SCOTCH
/** Warn user that we overflow asked memory during scotch call */
static inline
void _MMG5_warnScotch(MMG5_pMesh mesh) {
  if ( mesh->info.imprim > 4 || mesh->info.ddebug ) {
    if ( mesh->info.mem >= 0 ) {
      fprintf(stdout,"  ## Warning: we will overflow the memory asked with \"-m\"");
      fprintf(stdout," option during Scotch call.\n" );
    }
  }
}
#endif

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Warn user that some tetrahedra of the mesh have been reoriented.
 *
 */
static inline
void _MMG5_warnOrientation(MMG5_pMesh mesh) {
  if ( mesh->xt ) {
    if ( mesh->xt != mesh->ne ) {
      fprintf(stdout,"  ## Warning: %d tetra on %d reoriented.\n",
              mesh->xt,mesh->ne);
      fprintf(stdout,"  Your mesh may be non-conform.\n");
    }
    else {
      fprintf(stdout,"  ## Warning: all tetra reoriented.\n");
    }
  }
  mesh->xt = 0;
}

/**
 * \param sigid signal number.
 *
 * Signal handling: specify error messages depending from catched signal.
 *
 */
static inline
void _MMG5_excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
  case SIGABRT:
    fprintf(stdout,"  *** potential lack of memory.\n");  break;
  case SIGFPE:
    fprintf(stdout,"  Floating-point exception\n"); break;
  case SIGILL:
    fprintf(stdout,"  Illegal instruction\n"); break;
  case SIGSEGV:
    fprintf(stdout,"  Segmentation fault\n");  break;
  case SIGTERM:
  case SIGINT:
    fprintf(stdout,"  Program killed\n");  break;
  }
  exit(EXIT_FAILURE);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * Set function pointers.
 *
 */
void _MMG5_setfunc(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( met->size == 1 || ( met->size == 3 && mesh->info.lag >= 0 ) ) {
    _MMG5_caltet     = _MMG5_caltet_iso;
    _MMG5_caltri     = _MMG5_caltri_iso;
    _MMG5_lenedg     = _MMG5_lenedg_iso;
    _MMG5_intmetvol  = _MMG5_interp_iso;
    _MMG5_interp4bar = _MMG5_interp4bar_iso;
    _MMG5_defsiz     = _MMG5_defsiz_iso;
    _MMG5_gradsiz    = _MMG5_gradsiz_iso;
  }
  else if ( met->size == 6 ) {
    _MMG5_caltet     = _MMG5_caltet_ani;
    _MMG5_caltri     = _MMG5_caltri_ani;
    _MMG5_lenedg     = _MMG5_lenedg_ani;
    _MMG5_intmetvol  = _MMG5_intmetvol_ani;
    _MMG5_interp4bar = _MMG5_interp4bar_ani;
    _MMG5_defsiz     = _MMG5_defsiz_ani;
    _MMG5_gradsiz    = _MMG5_gradsiz_ani;
  }
}

/**
 * Set common pointer functions between mmgs and mmg3d to the matching mmg3d
 * functions.
 */
void _MMG5_Set_commonFunc() {
  MMG5_Init_parameters    = _MMG5_Init_parameters;
  _MMG5_bezierCP          = _MMG5_mmg3dBezierCP;
  _MMG5_chkmsh            = _MMG5_mmg3dChkmsh;
#ifdef USE_SCOTCH
  _MMG5_renumbering       = _MMG5_mmg3dRenumbering;
#endif
}
