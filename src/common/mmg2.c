/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
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
 * \file mmg/mmg2.c
 * \author Algiane Froehly (Bx INP/Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * Common functions for ls discretization.
 *
 */

#include "mmgcommon.h"

/**
 * \param mesh   pointer toward the mesh structure.
 * \param ref    initial reference.
 * \param refint internal reference after ls discretization.
 * \param refint internal reference after ls discretization.
 * \return 1 if entity can be splitted, 0 if cannot be splitted.
 *
 * Identify whether an entity with reference ref should be split, and the
 * labels of the resulting entities.
 *
 */
int MMG5_isSplit(MMG5_pMesh mesh,int ref,int *refint,int *refext) {
  MMG5_pMat    pm;
  int          k;

  /* Check in the info->mat table if reference ref is supplied by the user */
  for (k=0; k<mesh->info.nmat; k++) {
    pm = &mesh->info.mat[k];
    if ( pm->ref == ref ) {
      if ( !pm->dospl ) {
        return 0;
      }
      else {
        *refint = pm->rin;
        *refext = pm->rex;
        return 1;
      }
    }
  }

  /* Default case: split with references MG_MINUS, MG_PLUS */
  *refint = MG_MINUS;
  *refext = MG_PLUS;
  return 1;

}

/**
 * \param mesh pointer toward the mesh
 * \param ref  final reference for which we are searching the initial one
 * \return initial reference associated to \a ref if founded, \a ref if not founded.
 *
 * Retrieve the initial domain reference associated to the (split) reference ref.
 *
 */
int MMG5_getIniRef(MMG5_pMesh mesh,int ref) {
  MMG5_pMat     pm;
  int           k;

  for (k=0; k<mesh->info.nmat; k++) {
    pm = &mesh->info.mat[k];
    if ( pm->ref == ref && !pm->dospl ) return pm->ref;
    if ( ref == pm->rin || ref == pm->rex ) return pm->ref;
  }

  if ( k==0 ) {
    /* No materials */
    return 0;
  }

  return ref;
}
