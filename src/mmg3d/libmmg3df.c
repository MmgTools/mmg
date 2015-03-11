/* =============================================================================
**  This file is part of the Mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
**
**  Mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with Mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the Mmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file mmg3d/libmmg3d.c
 * \brief Private Fortran API functions for MMG3D library.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \note Please, refer to the \ref mmg3d/libmmg3d5.h file for functions
 * documentation.
 *
 * Define the private Fortran API functions for MMG3D library
 * (incompatible functions with the main binary): adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */

#include "libmmg3d5.h"

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
    void nu pl;                                 \
    void nl pl                                  \
    { nu pc; }                                  \
    void nl##_ pl                               \
    { nu pc; }                                  \
    void nl##__ pl                              \
    { nu pc; }                                  \
    void nu pl

/**
 * See \ref MMG5_Free_all function in \ref mmg3d/libmmg3d5.h file.
 */
FORTRAN_NAME(MMG5_FREE_ALL,mmg5_free_all,(MMG5_pMesh *mesh,MMG5_pSol *met
#ifdef SINGUL
                                          ,MMG5_pSingul *singul
#endif
                 ),(mesh,met
#ifdef SINGUL
                    ,singul
#endif
                     )){

#ifdef SINGUL
    MMG5_Free_all(*mesh,*met,*singul);
#else
    MMG5_Free_all(*mesh,*met);
#endif
    return;
}

/**
 * See \ref MMG5_saveMesh function in \ref mmg3d/libmmg3d5.h file.
 */
FORTRAN_NAME(MMG5_SAVEMESH,mmg5_savemesh,(MMG5_pMesh *mesh, int* retval),
             (mesh,retval)){
    *retval = MMG5_saveMesh(*mesh);
    return;
}

/**
 * See \ref MMG5_mmg3dlib function in \ref mmg3d/libmmg3d5.h file.
 */
FORTRAN_NAME(MMG5_MMG3DLIB,mmg5_mmg3dlib,(MMG5_pMesh *mesh,MMG5_pSol *met
#ifdef SINGUL
                                          ,MMG5_pSingul *sing
#endif
                                          ,int* retval),(mesh,met
#ifdef SINGUL
                                                         ,sing
#endif
                                                         ,retval)){

#ifdef SINGUL
    *retval = MMG5_mmg3dlib(*mesh,*met,*sing);
#else
    *retval = MMG5_mmg3dlib(*mesh,*met);
#endif
    return;
}
