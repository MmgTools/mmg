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
 * \file mmg3d/anisomovpt.c
 * \brief Functions to move a point in the mesh.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmg3d.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param list pointer toward the volumic ball of the point.
 * \param ilist size of the volumic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 0.9 of the old minimum element quality.
 *
 * \return 0 if fail, 1 if success.
 *
 * Move internal point whose volumic is passed.
 *
 */
int _MMG5_movintpt_ani(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist,
                       int improve) {
  printf("TO IMPLEMENT: movintpt_ani\n");
  return(0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param listv pointer toward the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer toward the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \return 0 if fail, 1 if success.
 *
 * Move boundary regular point, whose volumic and surfacic balls are passed.
 *
 */
int _MMG5_movbdyregpt_ani(MMG5_pMesh mesh, MMG5_pSol met,int *listv,
                          int ilistv,int *lists,int ilists) {
  MMG5_pTetra       pt,pt0;
  MMG5_pxTetra      pxt;
  MMG5_pPoint       p0,p1,p2,ppt0;
  MMG5_Tria         tt;
  MMG5_pxPoint      pxp;
  _MMG5_Bezier      b;
  double            *n,r[3][3],lispoi[3*_MMG5_LMAX+1],ux,uy,uz,det2d;
  double            detloc,oppt[2],step,lambda[3];
  double            ll,m[2],uv[2],o[3],no[3],to[3];
  double            calold,calnew,caltmp,callist[ilistv];
  int               k,kel,iel,l,n0,na,nb,ntempa,ntempb,ntempc,nut,nxp;
  unsigned char     i0,iface,i;

  step = 0.1;
  nut    = 0;
  oppt[0] = 0.0;
  oppt[1] = 0.0;
  if ( ilists < 2 )      return(0);

  k      = listv[0] / 4;
  i0 = listv[0] % 4;
  pt = &mesh->tetra[k];
  n0 = pt->v[i0];
  p0 = &mesh->point[n0];
  assert( p0->xp && !MG_EDG(p0->tag) );

  n = &(mesh->xpoint[p0->xp].n1[0]);

  /** Step 1 : rotation matrix that sends normal n to the third coordinate vector of R^3 */
  if ( !_MMG5_rotmatrix(n,r) ) return(0);

  /** Step 2 : rotation of the oriented surfacic ball with r : lispoi[k] is the common edge
      between faces lists[k-1] and lists[k] */
  k                     = lists[0] / 4;
  iface = lists[0] % 4;
  pt            = &mesh->tetra[k];
  na = nb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != n0 ) {
      if ( !na )
        na = pt->v[_MMG5_idir[iface][i]];
      else
        nb = pt->v[_MMG5_idir[iface][i]];
    }
  }

  for (l=1; l<ilists; l++) {
    k                   = lists[l] / 4;
    iface = lists[l] % 4;
    pt          = &mesh->tetra[k];
    ntempa = ntempb = 0;
    for (i=0; i<3; i++) {
      if ( pt->v[_MMG5_idir[iface][i]] != n0 ) {
        if ( !ntempa )
          ntempa = pt->v[_MMG5_idir[iface][i]];
        else
          ntempb = pt->v[_MMG5_idir[iface][i]];
      }
    }
    if ( ntempa == na )
      p1 = &mesh->point[na];
    else if ( ntempa == nb )
      p1 = &mesh->point[nb];
    else if ( ntempb == na )
      p1 = &mesh->point[na];
    else {
      assert(ntempb == nb);
      p1 = &mesh->point[nb];
    }
    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*l+1] =      r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*l+2] =      r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*l+3] =      r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

    na = ntempa;
    nb = ntempb;
  }

  /* Finish with point 0 */;
  k      = lists[0] / 4;
  iface  = lists[0] % 4;
  pt     = &mesh->tetra[k];
  ntempa = ntempb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != n0 ) {
      if ( !ntempa )
        ntempa = pt->v[_MMG5_idir[iface][i]];
      else
        ntempb = pt->v[_MMG5_idir[iface][i]];
    }
  }
  if ( ntempa == na )
    p1 = &mesh->point[na];
  else if ( ntempa == nb )
    p1 = &mesh->point[nb];
  else if ( ntempb == na )
    p1 = &mesh->point[na];
  else {
    assert(ntempb == nb);
    p1 = &mesh->point[nb];
  }

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  lispoi[1] =    r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
  lispoi[2] =    r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
  lispoi[3] =    r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

  /* list goes modulo ilist */
  lispoi[3*ilists+1] = lispoi[1];
  lispoi[3*ilists+2] = lispoi[2];
  lispoi[3*ilists+3] = lispoi[3];

  /* At this point, lispoi contains the oriented surface ball of point p0, that has been rotated
     through r, with the convention that triangle l has edges lispoi[l]; lispoi[l+1] */

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( det2d < 0.0 )  return(0);
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  if ( det2d < 0.0 )    return(0);

  /** Step 3 : Compute optimal position to make current triangle equilateral, and average of
      these positions*/
  for (k=0; k<ilists; k++) {
    m[0] = 0.5*(lispoi[3*(k+1)+1] + lispoi[3*k+1]);
    m[1] = 0.5*(lispoi[3*(k+1)+2] + lispoi[3*k+2]);
    ux = lispoi[3*(k+1)+1] - lispoi[3*k+1];
    uy = lispoi[3*(k+1)+2] - lispoi[3*k+2];
    ll = ux*ux + uy*uy;
    if ( ll < _MMG5_EPSD )    continue;
    nut++;
    oppt[0] += (m[0]-_MMG5_SQR32*uy);
    oppt[1] += (m[1]+_MMG5_SQR32*ux);
  }
  oppt[0] *= (1.0 / nut);
  oppt[1] *= (1.0 / nut);

  /** Step 4 : locate new point in the ball, and compute its barycentric coordinates */
  det2d = lispoi[1]*oppt[1] - lispoi[2]*oppt[0];
  kel = 0;
  if ( det2d >= 0.0 ) {
    for (k=0; k<ilists; k++) {
      detloc = oppt[0]*lispoi[3*(k+1)+2] - oppt[1]*lispoi[3*(k+1)+1];
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == ilists ) return(0);
  }
  else {
    for (k=ilists-1; k>=0; k--) {
      detloc = lispoi[3*k+1]*oppt[1] - lispoi[3*k+2]*oppt[0];
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == -1 ) return(0);
  }

  /* Sizing of time step : make sure point does not go out corresponding triangle. */
  det2d = -oppt[1]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]) + \
    oppt[0]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]);
  if ( fabs(det2d) < _MMG5_EPSD ) return(0);

  det2d = 1.0 / det2d;
  step *= det2d;

  det2d = lispoi[3*(kel)+1]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]) - \
    lispoi[3*(kel)+2 ]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]);
  step *= det2d;
  step  = fabs(step);
  oppt[0] *= step;
  oppt[1] *= step;

  /* Computation of the barycentric coordinates of the new point in the corresponding triangle. */
  det2d = lispoi[3*kel+1]*lispoi[3*(kel+1)+2] - lispoi[3*kel+2]*lispoi[3*(kel+1)+1];
  if ( det2d < _MMG5_EPSD )    return(0);
  det2d = 1.0 / det2d;
  lambda[1] = lispoi[3*(kel+1)+2]*oppt[0] - lispoi[3*(kel+1)+1]*oppt[1];
  lambda[2] = -lispoi[3*(kel)+2]*oppt[0] + lispoi[3*(kel)+1]*oppt[1];
  lambda[1]*= (det2d);
  lambda[2]*= (det2d);
  lambda[0] = 1.0 - lambda[1] - lambda[2];

  /** Step 5 : come back to original problem, and compute patch in triangle iel */
  iel    = lists[kel] / 4;
  iface  = lists[kel] % 4;
  pt     = &mesh->tetra[iel];
  pxt    = &mesh->xtetra[pt->xt];

  _MMG5_tet2tri(mesh,iel,iface,&tt);

  if(!_MMG5_bezierCP(mesh,&tt,&b,MG_GET(pxt->ori,iface))){
    fprintf(stdout,"%s:%d: Error: function _MMG5_bezierCP return 0\n",
            __FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }

  /* Now, for Bezier interpolation, one should identify which of i,i1,i2 is 0,1,2
     recall uv[0] = barycentric coord associated to pt->v[1], uv[1] associated to pt->v[2] and
     1 - uv[0] - uv[1] is associated to pt->v[0]. For this, use the fact that kel, kel + 1 is
     positively oriented with respect to n */
  na = nb = 0;
  for( i=0 ; i<4 ; i++ ){
    if ( (pt->v[i] != n0) && (pt->v[i] != pt->v[iface]) ) {
      if ( !na )
        na = pt->v[i];
      else
        nb = pt->v[i];
    }
  }
  p1 = &mesh->point[na];
  p2 = &mesh->point[nb];
  detloc = _MMG5_det3pt1vec(p0->c,p1->c,p2->c,n);

  /* ntempa = point to which is associated 1 -uv[0] - uv[1], ntempb = uv[0], ntempc = uv[1] */
  ntempa = pt->v[_MMG5_idir[iface][0]];
  ntempb = pt->v[_MMG5_idir[iface][1]];
  ntempc = pt->v[_MMG5_idir[iface][2]];

  /* na = lispoi[kel] -> lambda[1], nb = lispoi[kel+1] -> lambda[2] */
  if ( detloc > 0.0 ) {
    if ( ntempb == na )
      uv[0] = lambda[1];
    else if (ntempb == nb )
      uv[0] = lambda[2];
    else {
      assert(ntempb == n0);
      uv[0] = lambda[0];
    }
    if ( ntempc == na )
      uv[1] = lambda[1];
    else if (ntempc == nb )
      uv[1] = lambda[2];
    else {
      assert(ntempc == n0);
      uv[1] = lambda[0];
    }
  }
  /* nb = lispoi[kel] -> lambda[1], na = lispoi[kel+1] -> lambda[2] */
  else {
    if ( ntempb == na )
      uv[0] = lambda[2];
    else if (ntempb == nb )
      uv[0] = lambda[1];
    else {
      assert(ntempb == n0);
      uv[0] = lambda[0];
    }
    if ( ntempc == na )
      uv[1] = lambda[2];
    else if (ntempc == nb )
      uv[1] = lambda[1];
    else {
      assert(ntempc == n0);
      uv[1] = lambda[0];
    }
  }
  if(!_MMG3D_bezierInt(&b,uv,o,no,to)){
    fprintf(stdout,"%s:%d: Error: function _MMG3D_bezierInt return 0\n",
            __FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }

  /* Test : make sure that geometric approximation has not been degraded too much */
  ppt0 = &mesh->point[0];
  ppt0->c[0] = o[0];
  ppt0->c[1] = o[1];
  ppt0->c[2] = o[2];

  ppt0->tag      = p0->tag;
  ppt0->ref      = p0->ref;


  nxp = mesh->xp + 1;
  if ( nxp > mesh->xpmax ) {
    _MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                       "larger xpoint table",
                       return(0));
    n = &(mesh->xpoint[p0->xp].n1[0]);
  }
  ppt0->xp = nxp;
  pxp = &mesh->xpoint[nxp];
  memcpy(pxp,&(mesh->xpoint[p0->xp]),sizeof(MMG5_xPoint));
  pxp->n1[0] = no[0];
  pxp->n1[1] = no[1];
  pxp->n1[2] = no[2];

  memcpy(&(met->m[0]),&(met->m[met->size*n0]),met->size*sizeof(double));

  /* For each surfacic triangle, build a virtual displaced triangle for check purposes */
  calold = calnew = DBL_MAX;
  for (l=0; l<ilists; l++) {
    k                   = lists[l] / 4;
    iface = lists[l] % 4;
    pt          = &mesh->tetra[k];
    _MMG5_tet2tri(mesh,k,iface,&tt);
    calold = MG_MIN(calold,_MMG5_caltri(mesh,met,&tt));
    for( i=0 ; i<3 ; i++ )
      if ( tt.v[i] == n0 )      break;
    assert(i<3);
    tt.v[i] = 0;
    caltmp = _MMG5_caltri(mesh,met,&tt);
    if ( caltmp < _MMG5_EPSD )        return(0.0);
    calnew = MG_MIN(calnew,caltmp);
  }
  if ( calold < _MMG5_NULKAL && calnew <= calold )    return(0);
  else if (calnew < _MMG5_NULKAL) return(0);
  else if ( calnew < 0.3*calold )        return(0);
  memset(pxp,0,sizeof(MMG5_xPoint));

  /* Test : check whether all volumes remain positive with new position of the point */
  calold = calnew = DBL_MAX;
  for (l=0; l<ilistv; l++) {
    k    = listv[l] / 4;
    i0 = listv[l] % 4;
    pt = &mesh->tetra[k];
    pt0 = &mesh->tetra[0];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[i0] = 0;
    calold = MG_MIN(calold, pt->qual);
    callist[l]=_MMG5_orcal(mesh,met,0);
    if ( callist[l] < _MMG5_EPSD )        return(0);
    calnew = MG_MIN(calnew,callist[l]);
  }
  if ( calold < _MMG5_NULKAL && calnew <= calold )    return(0);
  else if (calnew < _MMG5_NULKAL) return(0);
  else if ( calnew < 0.3*calold )        return(0);

  /* When all tests have been carried out, update coordinates and normals */
  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  n[0] = no[0];
  n[1] = no[1];
  n[2] = no[2];

//#warning update metric?

  for(l=0; l<ilistv; l++){
    (&mesh->tetra[listv[l]/4])->qual= callist[l];
  }
  return(1);
}
