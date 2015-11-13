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
 * \file mmg2d/libmmg2df.c
 * \brief Private Fortran API functions for MMG2D library.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \note Please, refer to the \ref mmg2d/libmmg2d.h file for functions
 * documentation.
 *
 * Define the private Fortran API functions for MMG2D library
 * (incompatible functions with the main binary): adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */

#include "libmmg2d.h"

/**
 * \def FORTRAN_NAME(nu,nl,pl,pc)
 * \brief Adds function definitions.
 * \param nu function name in upper case.
 * \param nl function name in lower case.
 * \param pl type of arguments.
 * \param pc name of arguments.
 * \note Macro coming from Scotch library.
 *
 * Adds function definitions with upcase, underscore and double
 * underscore to match any fortran compiler.
 *
 */
#define FORTRAN_NAME(nu,nl,pl,pc)               \
  void nu pl;                                   \
  void nl pl                                    \
  { nu pc; }                                    \
  void nl##_ pl                                 \
  { nu pc; }                                    \
  void nl##__ pl                                \
  { nu pc; }                                    \
  void nu pl

/**
 * See \ref MMG5_Free_all function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG5_FREE_ALL,mmg5_free_all,(MMG5_pMesh *mesh,MMG5_pSol *met
               ),(mesh,met
                 )){

  MMG5_Free_all(*mesh,*met);

  return;
}

/**
 * See \ref MMG2_saveMesh function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2_SAVEMESH,mmg2_savemesh,(MMG5_pMesh *mesh,char *meshin,int *strlen, int* retval),
             (mesh,meshin,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMG2_saveMesh(*mesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}
/**
 * See \ref MMG2_mmg2dlib function in \ref mmg2d/libmmg2d.h file.
 */
//#warning todo : add callbackfunction argument in Fortran
FORTRAN_NAME(MMG2_MMG2DLIB,mmg2_mmg2dlib,(MMG5_pMesh *mesh,MMG5_pSol *met
                                          ,int* retval),(mesh,met
                                                         ,retval)){

  *retval = MMG2_mmg2dlib(*mesh,*met,NULL);

  return;
}
