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
 * \file mmg3d/anisosiz.c
 * \brief Fonctions for anisotropic size map computation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmg3d.h"


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric stucture.
 * \return 0 if fail, 1 otherwise.
 *
 * Define size at points by intersecting the surfacic metric and the
 * physical metric.
 *
 */
int defsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  printf("defsiz_ani\n");
  exit(0);

  return(1);
}
/* Enforces mesh gradation by truncating metric field */
int gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  printf("gradsize_ani not implemented\n");
  return(1);
}
