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
 * \file mmgs/mettools.c
 * \brief Algebraic tools for metric handling.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

extern char ddb;


/* Build metric tensor at a fictitious ridge point, whose normal and tangent are provided */
inline int buildridmetfic(MMG5_pMesh mesh,double t[3],double n[3],double dtan,double dv,double m[6]) {
  double u[3],r[3][3];

  u[0] = n[1]*t[2] - n[2]*t[1];
  u[1] = n[2]*t[0] - n[0]*t[2];
  u[2] = n[0]*t[1] - n[1]*t[0];

  /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(p0->m[0],dv,0)*/
  r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n[0];
  r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n[1];
  r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n[2];

  m[0] = dtan*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1];
  m[1] = dtan*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1];
  m[2] = dtan*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1];
  m[3] = dtan*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1];
  m[4] = dtan*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1];
  m[5] = dtan*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1];

  return(1);
}
