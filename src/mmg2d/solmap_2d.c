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
 * \file mmg2d/solmap_2d.c
 * \brief  Compute isotropic size map according to the mean of the length of the edges
 * passing through a point.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg2d.h"

/**
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward the sol structure
 * \return 1 if success
 *
 * Compute isotropic size map according to the mean of the length of the edges
 * passing through a point.
 *
 */
int MMG2D_doSol_iso(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria      ptt,pt;
  MMG5_pPoint     p1,p2;
  double          ux,uy,dd;
  int             i,k,ib,ipa,ipb;
  int             MMG_inxtt[5] = {0,1,2,0,1};

  /* Memory alloc */
  if ( sol->size!=1 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,sol->size);
    return 0;
  }

  if ( !MMG2D_Set_solSize(mesh,sol,MMG5_Vertex,mesh->np,sol->size) )
    return 0;

  /* tagdel will be used to count the number of edges passing through each
   * point */
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    p1->tagdel = 0;
  }

  /* Travel the triangles edges and add the edge contribution to edges
   * extermities */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !ptt->v[0] )  continue;

    for (i=0; i<3; i++) {
      ib  = MMG_inxtt[i+1];
      ipa = ptt->v[i];
      ipb = ptt->v[ib];
      p1  = &mesh->point[ipa];
      p2  = &mesh->point[ipb];

      ux  = p1->c[0] - p2->c[0];
      uy  = p1->c[1] - p2->c[1];
      dd  = sqrt(ux*ux + uy*uy);

      sol->m[ipa] += dd;
      p1->tagdel++;
      sol->m[ipb] += dd;
      p2->tagdel++;
    }
  }

  /* if hmax is not specified, compute it from the metric */
  if ( mesh->info.hmax < 0. ) {
    dd = 0.;
    for (k=1; k<=mesh->np; k++) {
      p1 = &mesh->point[k];
      if ( !p1->tagdel ) continue;
      dd = MG_MAX(dd,sol->m[k]);
    }
    assert ( dd > 0. );
    mesh->info.hmax = 10.*dd;
  }

  /* vertex size */
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !p1->tagdel )  {
      sol->m[k] = mesh->info.hmax;
      continue;
    }

    sol->m[k] = sol->m[k] / (double)p1->tagdel;
    p1->tagdel = 0;
  }

  /* compute quality */
  if ( MMG2D_caltri ) {
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      pt->qual = MMG2D_caltri_iso(mesh,sol,pt);
    }
  }

  if ( mesh->info.imprim < -4 )
    fprintf(stdout,"   HMAX %f\n",mesh->info.hmax);
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward the sol structure
 * \return 1 if success
 *
 * Compute unit mesh anisotropic size map using statistical concept of
 * length distribution tensors (formula 5 of \cite COUPEZ20112391).
 *
 */
int MMG2D_doSol_ani(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria      ptt,pt;
  MMG5_pPoint     p1,p2;
  double          ux,uy,dd,tensordot[3];
  int             i,k,ib,iadr,ipa,ipb;
  int             MMG_inxtt[5] = {0,1,2,0,1};

  /* Memory alloc */
  if ( sol->size!=3 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,sol->size);
    return 0;
  }

  if ( !MMG2D_Set_solSize(mesh,sol,MMG5_Vertex,mesh->np,sol->size) )
    return 0;

  /* tagdel will be used to count the number of edges passing through each
   * point */
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    p1->tagdel = 0;
  }

  /* Travel the triangles edges and add the edge contribution to edges
   * extermities */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !ptt->v[0] )  continue;

    for (i=0; i<3; i++) {
      ib  = MMG_inxtt[i+1];
      ipa = ptt->v[i];
      ipb = ptt->v[ib];
      p1  = &mesh->point[ipa];
      p2  = &mesh->point[ipb];

      ux  = p1->c[0] - p2->c[0];
      uy  = p1->c[1] - p2->c[1];

      tensordot[0] = ux*ux;
      tensordot[1] = ux*uy;
      tensordot[2] = uy*uy;

      iadr = 3*ipa;
      sol->m[iadr]   += tensordot[0];
      sol->m[iadr+1] += tensordot[1];
      sol->m[iadr+2] += tensordot[2];
      p1->tagdel++;

      iadr = 3*ipb;
      sol->m[iadr]   += tensordot[0];
      sol->m[iadr+1] += tensordot[1];
      sol->m[iadr+2] += tensordot[2];
      p2->tagdel++;
    }
  }

  /* Compute metric tensor and hmax if not specified */
  double hmax = FLT_MAX;
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !p1->tagdel ) {
      continue;
    }

    iadr = 3*k;
    /* Metric = nedges/dim * inv (sum(tensor_dot(edges,edges))).
     * sum(tensor_dot) is stored in sol->m so reuse tensordot to
     * compute M.  */
    dd = 1./(sol->m[iadr]*sol->m[iadr+2] - sol->m[iadr+1]*sol->m[iadr+1]);
    dd *= (double)p1->tagdel*0.5;

    tensordot[0] = sol->m[iadr+2];
    tensordot[1] = -sol->m[iadr+1];
    tensordot[2] = sol->m[iadr];

    sol->m[iadr]   = dd*tensordot[0];
    sol->m[iadr+1] = dd*tensordot[1];
    sol->m[iadr+2] = dd*tensordot[2];

    /* Check metric */
    double lambda[2],vp[2][2];
    MMG5_eigensym(sol->m+iadr,lambda,vp);

    assert (lambda[0] > 0. && lambda[1] > 0. && "Negative eigenvalue");

    /* Normally the case where the point belongs to only 2 colinear points is
    impossible */
    assert (isfinite(lambda[0]) && isfinite(lambda[1]) && "wrong eigenvalue");

    hmax = MG_MIN(hmax,lambda[0]);
    hmax = MG_MIN(hmax,lambda[1]);

  }
  if ( mesh->info.hmax < 0.) {
    assert ( hmax > 0. && hmax < FLT_MAX && "Wrong hmax value" );
    mesh->info.hmax = 10./sqrt(hmax);
  }

  /* vertex size: impose hmax size */
  hmax = 1./(mesh->info.hmax*mesh->info.hmax);
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !p1->tagdel )  {
      iadr = 3*k;
      sol->m[iadr]   = hmax;
      sol->m[iadr+2] = sol->m[iadr];
      continue;
    }
    p1->tagdel = 0;
  }

  /* compute quality */
  if ( MMG2D_caltri ) {
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      pt->qual = MMG2D_caltri_ani(mesh,sol,pt);
    }
  }

  if ( mesh->info.imprim < -4 )
    fprintf(stdout,"   HMAX %f\n",mesh->info.hmax);
  return 1;
}
