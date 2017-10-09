/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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
 * \file mmg2d/isosiz_2d.c
 * \brief Interpolation of metrics
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "mmg2d.h"

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return 0 if fail, 1 otherwise
 *
 * New version for the definition of a size map; takes into account the
 * curvature of the external and internal curves present in the mesh
 *
 */
int _MMG2_defsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria       pt;
  MMG5_pPoint      p1,p2;
  double           t1[2],t2[2],b1[2],b2[2],gpp1[2],gpp2[2],pv,M1,M2;
  double           ps1,ps2,ux,uy,ll,li,lm,hmax,hausd,hmin;
  int              k,ip1,ip2;
  unsigned char    i,i1,i2;


  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Defining isotropic map\n");

  if ( mesh->info.hmax < 0.0 )  mesh->info.hmax = 0.5 * mesh->info.delta;

  hmax = mesh->info.hmax;
  hausd = mesh->info.hausd;
  hmin = mesh->info.hmin;

  /* Allocate the structure */
  if ( !met->np ) {

    /* Allocate and store the header informations for each solution */
    if ( !MMG2D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,1) ) {
      return 0;
    }

    /* Initialize metric with a constant size in the case met->np = 0 (meaning that no metric was supplied) */
    for (k=1; k<=mesh->np; k++)
      met->m[k] = hmax;
  }

  /* Only the boundary edges impose a minimal size feature */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      if ( !MG_EDG(pt->tag[i]) ) continue;
      i1 = _MMG5_inxt2[i];
      i2 = _MMG5_iprv2[i];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];

      p1 = &mesh->point[ip1];
      p2 = &mesh->point[ip2];

      ux = p2->c[0] - p1->c[0];
      uy = p2->c[1] - p1->c[1];
      ll = ux*ux + uy*uy;
      if ( ll < _MMG5_EPSD ) continue;
      li = 1.0 / sqrt(ll);

      /* Recovery of the two tangent vectors associated to points p1,p2; they need not be oriented in the same fashion */
      if ( MG_SIN(p1->tag) || (p1->tag & MG_NOM) ) {
        t1[0] = li*ux;
        t1[1] = li*uy;
      }
      else {
        t1[0] = -p1->n[1];
        t1[1] = p1->n[0];
      }

      if ( MG_SIN(p2->tag) || (p2->tag & MG_NOM) ) {
        li = 1.0 / sqrt(ll);
        t2[0] = li*ux;
        t2[1] = li*uy;
      }
      else {
        t2[0] = -p2->n[1];
        t2[1] = p2->n[0];
      }

      /* Calculation of the two Bezier coefficients along the curve */
      ps1   = ux*t1[0] + uy*t1[1];
      b1[0] = p1->c[0] + _MMG5_ATHIRD*ps1*t1[0];
      b1[1] = p1->c[1] + _MMG5_ATHIRD*ps1*t1[1];

      ps2   = ux*t2[0]+uy*t2[1];
      b2[0] = p2->c[0] - _MMG5_ATHIRD*ps2*t2[0];
      b2[1] = p2->c[1] - _MMG5_ATHIRD*ps2*t2[1];

      ps1 *= ps1;
      ps2 *= ps2;

      if ( ps1 < _MMG5_EPSD || ps2 < _MMG5_EPSD ) continue;

      /* \gamma^{\prime\prime}(0); \gamma^\prime(0) = ps*t1 by construction */
      gpp1[0] = 6.0*(p1->c[0] - 2.0*b1[0] + b2[0]);
      gpp1[1] = 6.0*(p1->c[1] - 2.0*b1[1] + b2[1]);

      /* Vector product gpp1 ^ t1 */
      pv = gpp1[0]*t1[1] - gpp1[1]*t1[0];
      M1 = fabs(pv)/ps1;

      /* \gamma^{\prime\prime}(1); \gamma^\prime(1) = -ps*t2 by construction */
      gpp2[0] = 6.0*(p2->c[0] - 2.0*b2[0] + b1[0]);
      gpp2[1] = 6.0*(p2->c[1] - 2.0*b2[1] + b1[1]);

      /* Vector product gpp2 ^ t2 */
      pv = gpp2[0]*t2[1] - gpp2[1]*t2[0];
      M2 = fabs(pv)/ps2;

      M1 = MG_MAX(M1,M2);
      if ( M1 < _MMG5_EPSD)
        lm = hmax;
      else {
        lm = 8.0*hausd / M1;
        lm = sqrt(lm);
      }
      met->m[ip1] = MG_MAX(hmin,MG_MIN(met->m[ip1],lm));
      met->m[ip2] = MG_MAX(hmin,MG_MIN(met->m[ip2],lm));
    }
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return 0 if fail, 1 otherwise
 *
 * Isotropic mesh gradation routine
 *
 */
int _MMG2_gradsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            hgrad,ll,h1,h2,hn;
  int               k,it,ip1,ip2,maxit,nup,nu;
  unsigned char     i,i1,i2;


  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Grading mesh\n");

  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  hgrad = log(mesh->info.hgrad);
  it = nup = 0;
  maxit = 100;

  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<3; i++) {
        i1  = _MMG5_inxt2[i];
        i2  = _MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        p1 = &mesh->point[ip1];
        p2 = &mesh->point[ip2];
        if ( p1->flag < mesh->base-1 && p2->flag < mesh->base-1 )  continue;

        ll = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]);
        ll = sqrt(ll);

        h1 = met->m[ip1];
        h2 = met->m[ip2];
        if ( h1 < h2 ) {
          if ( h1 < _MMG5_EPSD )  continue;
          hn  = h1 + hgrad*ll;
          if ( h2 > hn ) {
            met->m[ip2] = hn;
            p2->flag    = mesh->base;
            nu++;
          }
        }
        else {
          if ( h2 < _MMG5_EPSD )  continue;
          hn = h2 + hgrad*ll;
          if ( h1 > hn ) {
            met->m[ip1] = hn;
            p1->flag    = mesh->base;
            nu++;
          }
        }
      }
    }
    nup += nu;
  }
  while ( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"     gradation: %7d updated, %d iter.\n",nup,it);

  return(1);
}
