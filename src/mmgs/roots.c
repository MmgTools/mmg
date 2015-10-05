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
 * \file mmgs/roots.c
 * \brief Functions to compute roots of degree 2 and 3 polynomial.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmgs.h"
#define EPSRO     1.e-3

/* Returns the 2 complex roots of a degree 2 polynomial with real coefficients a[2]T^2 + ... + a[0]
   By convention, the real roots are stored first (same thing for multiple roots) :
   return value = number of roots, counted with multiplicity */
int rootDeg2(double complex a[3], double complex r[2]){
  double complex Delta,delta,r1,r2;

  if( cabs(a[2])<_MMG5_EPSD ) {
    if( cabs(a[1])<_MMG5_EPSD ) {
      r[0] = r[1] = 0.0;
      return(0);
    }

    else{
      r[0] = r[1] = -a[0]/a[1];
      return(1);
    }
  }

  Delta = a[1]*a[1] - 4.0*a[2]*a[0];
  /* ONE square root of Delta */
  delta = cpow(Delta,0.5);
  r1 = 0.5/a[2]*(-a[1] - delta);
  r2 = 0.5/a[2]*(-a[1] + delta);

  if( fabs(cimag(r1)) < EPSRO ){
    r[0] = creal(r1);
    r[1] = ( fabs(cimag(r2)) < EPSRO ) ? creal(r2) : r2;
  }
  else {
    r[0] = ( fabs(cimag(r2)) < EPSRO ) ? creal(r2) : r2;
    r[1] = r1;
  }

  return(2);
}
