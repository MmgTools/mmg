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
 * \file mmg2d/libmmg2d_toolsf.c
 * \brief Fortran API functions for MMG2D library.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \note Please, refer to the \ref mmgs/libmmgs.h file for functions
 * documentation.
 *
 * Define the private Fortran API functions for MMG2D library
 * (incompatible functions with the main binary): adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */

#include "mmg2d.h"

/**
 * See \ref MMG2D_setfunc function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_SETFUNC,mmg2d_setfunc,
             (MMG5_pMesh *mesh,MMG5_pSol *met),
             (mesh,met)) {
  MMG2D_setfunc(*mesh,*met);
  return;
}


/**
 * See \ref MMG2D_Get_adjaTri function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_GET_ADJATRI,mmg2d_get_adjatri,
             (MMG5_pMesh *mesh,int* kel, int* listri, int* retval),
             (mesh,kel,listri,retval)) {
  *retval =  MMG2D_Get_adjaTri(*mesh,*kel,listri);
  return;
}

/**
 * See \ref MMG2D_Get_adjaVertices function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_GET_ADJAVERTICES,mmg2d_get_adjavertices,
             (MMG5_pMesh *mesh,int* ip, int* lispoi, int* retval),
             (mesh,ip,lispoi,retval)) {
  *retval =  MMG2D_Get_adjaVertices(*mesh, *ip,lispoi);
  return;
}

/**
 * See \ref MMG2D_Get_adjaVerticesFast function in \ref mmg2d/libmmg2d.h file.
 */
FORTRAN_NAME(MMG2D_GET_ADJAVERTICESFAST,mmg2d_get_adjaverticesfast,
             (MMG5_pMesh *mesh,int* ip, int *start, int* lispoi, int* retval),
             (mesh,ip,start,lispoi,retval)) {
  *retval =  MMG2D_Get_adjaVerticesFast(*mesh,*ip, *start,lispoi);
  return;
}
