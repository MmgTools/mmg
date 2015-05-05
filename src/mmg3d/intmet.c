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
 * \file mmg3d/intmet.c
 * \brief Metric interpolations.
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
 * \param ma pointer on a metric
 * \param mb pointer on a metric
 * \param mp pointer on the computed interpolated metric
 * \param t interpolation parameter (comprise between 0 and 1)
 *
 *
 * Linear interpolation of anisotropic sizemap along an internal edge
 *
 */
int _MMG5_intmetvol_ani(double *ma,double *mb,double *mp,double t) {
  double        dma[6],dmb[6],mai[6],mbi[6],mi[6];
  int           i;

  for (i=0; i<6; i++) {
    dma[i] = ma[i];
    dmb[i] = mb[i];
  }
  if ( !_MMG5_invmat(dma,mai) || !_MMG5_invmat(dmb,mbi) ) {
    fprintf(stderr,"  ## INTERP INVALID METRIC.\n");
    return(0);
  }
  for (i=0; i<6; i++)
    mi[i] = (1.0-t)*mai[i] + t*mbi[i];

  if ( !_MMG5_invmat(mi,mai) ) {
    fprintf(stderr,"  ## INTERP INVALID METRIC.\n");
    return(0);
  }

  for (i=0; i<6; i++)  mp[i] = mai[i];

  return 1;
}
