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
 * \file common/mmg_quality.c
 * \brief Functions to compute elements quality and edge lengths.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the meric structure.
 * \param ptt pointer toward the triangle structure.
 * \return The computed quality.
 *
 * Compute the quality of the surface triangle \a ptt with respect to
 * an anisotropic metric.
 *
 */
inline double _MMG5_caltri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) {
  double        rap,anisurf,l[3];
  int           ia,ib,ic;

  ia = ptt->v[0];
  ib = ptt->v[1];
  ic = ptt->v[2];

  anisurf = _MMG5_surftri_ani(mesh,met,ptt);

  l[0] = _MMG5_lenedg_ani(mesh,met,ib,ic,( ptt->tag[0] & MG_GEO ));
  l[1] = _MMG5_lenedg_ani(mesh,met,ia,ic,( ptt->tag[1] & MG_GEO ));
  l[2] = _MMG5_lenedg_ani(mesh,met,ia,ib,( ptt->tag[2] & MG_GEO ));

  rap = l[0]*l[0] + l[1]*l[1] + l[2]*l[2];

  if ( rap < _MMG5_EPSD ) return(0.0);

  return (anisurf / rap);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the meric structure.
 * \param ptt pointer toward the triangle structure.
 * \return The computed quality.
 *
 * Compute the quality of the surface triangle \a ptt with respect to
 * an isotropic metric.
 *
 */
inline double _MMG5_caltri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) {
  double   *a,*b,*c,cal,abx,aby,abz,acx,acy,acz,bcx,bcy,bcz,rap;

  a = &mesh->point[ptt->v[0]].c[0];
  b = &mesh->point[ptt->v[1]].c[0];
  c = &mesh->point[ptt->v[2]].c[0];

  /* area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  cal  = (aby*acz - abz*acy) * (aby*acz - abz*acy);
  cal += (abz*acx - abx*acz) * (abz*acx - abx*acz);
  cal += (abx*acy - aby*acx) * (abx*acy - aby*acx);

  if ( cal < _MMG5_EPSD2 )  return(0.0);

  /* qual = 2.*surf / length */
  rap  = abx*abx + aby*aby + abz*abz;
  rap += acx*acx + acy*acy + acz*acz;
  rap += bcx*bcx + bcy*bcy + bcz*bcz;

  if ( rap < _MMG5_EPSD2 )  return(0.0);

  return(sqrt(cal) / rap);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param ned edges number.
 * \param avlen pointer toward the average edges lengths.
 * \param amin index of first extremity of the smallest edge.
 * \param bmin index of second extremity of the smallest edge.
 * \param lmin smallest edge length.
 * \param amax index of first extremity of the largest edge.
 * \param bmax index of second extremity of the largest edge.
 * \param lmax largest edge length.
 * \param bd pointer toward the table of the quality span.
 * \param hl pointer toward the table that store the number of edges for each
 * span of quality
 *
 * Display histogram of edge length.
 *
 */
void _MMG5_displayHisto(MMG5_pMesh mesh, int ned, double *avlen,
                        int amin, int bmin, double lmin,
                        int amax, int bmax, double lmax, double *bd, int *hl )
{
  double dned;
  int    k;

  dned     = (double)ned;
  (*avlen) = (*avlen) / dned;

  fprintf(stdout,"\n  -- RESULTING EDGE LENGTHS  %d\n",ned);
  fprintf(stdout,"     AVERAGE LENGTH         %12.4f\n",(*avlen));
  fprintf(stdout,"     SMALLEST EDGE LENGTH   %12.4f   %6d %6d\n",
          lmin,amin,bmin);
  fprintf(stdout,"     LARGEST  EDGE LENGTH   %12.4f   %6d %6d \n",
          lmax,amax,bmax);
  if ( abs(mesh->info.imprim) < 3 ) return;

  /* if ( hl[3]+hl[4]+hl[5] ) */
  /*   fprintf(stdout,"   %6.2f < L <%5.2f  %8d   %5.2f %%  \n", */
  /*           bd[3],bd[6],hl[3]+hl[4]+hl[5],100.*(hl[3]+hl[4]+hl[5])/(double)ned); */
  if ( hl[2]+hl[3]+hl[4] )
    fprintf(stdout,"   %6.2f < L <%5.2f  %8d   %5.2f %%  \n",
            bd[2],bd[5],hl[2]+hl[3]+hl[4],100.*(hl[2]+hl[3]+hl[4])/(double)ned);

  if ( abs(mesh->info.imprim) < 4 ) return;

  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"\n     HISTOGRAMM:\n");
    if ( hl[0] )
      fprintf(stdout,"     0.00 < L < 0.30  %8d   %5.2f %%  \n",
              hl[0],100.*(hl[0]/(float)ned));
    if ( lmax > 0.2 ) {
      for (k=2; k<9; k++) {
        if ( hl[k-1] > 0 )
          fprintf(stdout,"   %6.2f < L <%5.2f  %8d   %5.2f %%  \n",
                  bd[k-1],bd[k],hl[k-1],100.*(hl[k-1]/(float)ned));
      }
      if ( hl[8] )
        fprintf(stdout,"     5.   < L         %8d   %5.2f %%  \n",
                hl[8],100.*(hl[8]/(float)ned));
    }
  }
}
