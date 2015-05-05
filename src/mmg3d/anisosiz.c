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
 * \param met pointer toward the metric structure.
 * \param kel index of the tetra in which we work.
 * \param iface face of the tetra on which we work.
 * \param ip index of the point on which we want to compute the metric
 * (in tetra \a kel).
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a SINGULARITY of the geometry, associated to
 * the geometric approx of the surface. metric \f$=\alpha*Id\f$, \f$\alpha =\f$
 * size.
 *
 */
static int _MMG5_defmetsin(MMG5_pMesh mesh,MMG5_pSol met,int kel, int iface, int ip) {
  MMG5_pTetra        pt;
  MMG5_pxTetra       pxt;
  MMG5_pPoint        p0;
  MMG5_pPar          par;
  double             *m,n[3],isqhmin,isqhmax,b0[3],b1[3],ps1,tau[3];
  double             ntau2,gammasec[3];
  double             c[3],kappa,maxkappa,alpha, hausd;
  int                lists[_MMG5_LMAX+2],listv[_MMG5_LMAX+2],ilist,ilists,ilistv;
  int                k,iel,idp,ifac;
  unsigned char      i,i0,i1,i2;

  pt  = &mesh->tetra[kel];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  ilist = _MMG5_boulesurfvolp(mesh,kel,ip,iface,
                              listv,&ilistv,lists,&ilists,(p0->tag & MG_NOM));

  if ( ilist!=1 ) {
    printf("Error; unable to compute the ball af the point %d.\n", idp);
    printf("Exit program.\n");
    exit(EXIT_FAILURE);
  }

  isqhmin  = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax  = 1.0 / (mesh->info.hmax*mesh->info.hmax);
  maxkappa = 0.0;
  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    ifac  = lists[k] % 4;
    pt    = &mesh->tetra[iel];
    assert(pt->xt && (mesh->xtetra[pt->xt].ftag[ifac] & MG_BDY) );
    pxt = &mesh->xtetra[pt->xt];

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    i0  = _MMG5_idir[ifac][i];
    i1  = _MMG5_idir[ifac][_MMG5_inxt2[i]];
    i   = _MMG5_iprv2[i];
    i2  = _MMG5_idir[ifac][i];

    /* Computation of the two control points associated to edge p0p1 with
     * p1=mesh->point[pt->v[i1]]: p0 is singular */
    _MMG5_norpts(mesh,pt->v[i0],pt->v[i1],pt->v[i2],n);

    _MMG5_bezierEdge(mesh,idp,pt->v[i1],b0,b1,
                     MG_EDG(pxt->tag[_MMG5_iarf[ifac][i]]),n);

    /* tangent vector */
    tau[0] = 3.0*(b0[0] - p0->c[0]);
    tau[1] = 3.0*(b0[1] - p0->c[1]);
    tau[2] = 3.0*(b0[2] - p0->c[2]);
    ntau2  = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];

    /* 2nd order derivative */
    gammasec[0] = 6.0*p0->c[0] - 12.0*b0[0] + 6.0*b1[0];
    gammasec[1] = 6.0*p0->c[1] - 12.0*b0[1] + 6.0*b1[1];
    gammasec[2] = 6.0*p0->c[2] - 12.0*b0[2] + 6.0*b1[2];
    if ( ntau2 < _MMG5_EPSD )  continue;
    ntau2 = 1.0 / ntau2;

    /* derivative via the normal parametrization */
    ps1  = gammasec[0]*tau[0] + gammasec[1]*tau[1] + gammasec[2]*tau[2];
    c[0] = gammasec[0] - ps1*tau[0]*ntau2;
    c[1] = gammasec[1] - ps1*tau[1]*ntau2;
    c[2] = gammasec[2] - ps1*tau[2]*ntau2;

    kappa = ntau2 * sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

    /* local hausdorff for triangle */
    hausd = mesh->info.hausd;
    for (i=0; i<mesh->info.npar; i++) {
      par = &mesh->info.par[i];
      if ( (par->elt == MMG5_Triangle) && (pxt->ref[ifac] == par->ref ) )
        hausd = par->hausd;
    }
    maxkappa = MG_MAX(kappa,maxkappa/hausd);
  }

  alpha = 1.0 / 8.0 * maxkappa;
  alpha = MG_MIN(alpha,isqhmin);
  alpha = MG_MAX(alpha,isqhmax);

  m = &met->m[6*idp];
  memset(m,0,6*sizeof(double));
  m[0] = m[3] = m[5] = alpha;

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param kel index of the tetra in which we work.
 * \param iface face of the tetra on which we work.
 * \param ip index of the point on which we want to compute the metric
 * (in tetra \a kel).
 * \return 1 if success, 0 otherwise.
 *
 * Compute metric tensor associated to a ridge point : convention is a bit weird
 * here :
 * \a p->m[0] is the specific size in direction \a t,
 * \a p->m[1] is the specific size in direction \f$ u_1 = n_1^t\f$
 * \a p->m[2] is the specific size in direction \f$ u_2 = n_2^t\f$,
 * and at each time, metric tensor has to be recomputed, depending on the side.
 *
 */
static int _MMG5_defmetrid(MMG5_pMesh mesh,MMG5_pSol met,int kel,
                           int iface, int ip)
{
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_Tria      ptt;
  MMG5_pPoint    p0,p1,p2;
  _MMG5_Bezier   b;
  int            k,iel,idp,ilist1,ilist2,ilist,*list;
  int            list1[_MMG5_LMAX+2],list2[_MMG5_LMAX+2],iprid[2],ier;
  double         *m,isqhmin,isqhmax,*n1,*n2,*n,*t;
  double         trot[2],u[2],ux,uy,uz,det,bcu[3];
  double         r[3][3],lispoi[3*_MMG5_LMAX+1];
  double         detg,detd;
  int            i,i0,i1,i2,ifac;

  pt  = &mesh->tetra[kel];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  isqhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  n1 = &mesh->xpoint[p0->xp].n1[0];
  n2 = &mesh->xpoint[p0->xp].n2[0];
  t  = p0->n;

  m = &met->m[6*idp];
  memset(m,0,6*sizeof(double));
  m[0] = isqhmax;
  m[1] = isqhmax;
  m[2] = isqhmax;

  // Call bouletrid that construct the surfacic ball
  ier = _MMG5_bouletrid(mesh,kel,iface,ip,&ilist1,list1,&ilist2,list2,
                        &iprid[0],&iprid[1] );
  if ( !ier ) {
    printf("Error: unable to compute the two balls at the ridge point %d.\n", idp);
    return(0);
  }

  // Check the ball orientation
  // If _MMG5_directsurfball return 1 it is useless to call this function,
  // thus it is valid here to call it inside the assert.
  assert(_MMG5_directsurfball(mesh, idp,list1,ilist1,n1) == 1);
  assert(_MMG5_directsurfball(mesh, idp,list2,ilist2,n2) == 1);

  /* Specific size in direction of t */
  m[0] = _MMG5_ridSizeInTangentDir(mesh,p0,idp,iprid,isqhmin,isqhmax);

  /* Characteristic sizes in directions u1 and u2 */
  for (i=0; i<2; i++) {
    if ( i==0 ) {
      n = n1;
      ilist = ilist1;
      list  = &list1[0];
    }
    else {
      n = n2;
      ilist = ilist2;
      list  = &(list2[0]);
    }
    _MMG5_rotmatrix(n,r);

    /* Apply rotation to the half-ball under consideration */
    i1 = 0;
    for (k=0; k<ilist; k++) {
      iel  = list[k] / 4;
      ifac = list[k] % 4;
      pt = &mesh->tetra[iel];
      for ( i0=0; i0!=3; ++i0 ) {
        if ( pt->v[_MMG5_idir[ifac][i0]]==idp ) break;
      }
      assert(i0!=3);
      i1   = _MMG5_inxt2[i0];
      p1 = &mesh->point[pt->v[_MMG5_idir[ifac][i1]]];

      ux = p1->c[0] - p0->c[0];
      uy = p1->c[1] - p0->c[1];
      uz = p1->c[2] - p0->c[2];

      lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
      lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
      lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
    }

    /* last point : the half-ball is open : ilist tria, and ilist +1 points ;
       lists are enumerated in direct order */
    i2 = _MMG5_inxt2[i1];
    p2 = &mesh->point[pt->v[_MMG5_idir[ifac][i2]]];

    ux = p2->c[0] - p0->c[0];
    uy = p2->c[1] - p0->c[1];
    uz = p2->c[2] - p0->c[2];

    lispoi[3*ilist+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*ilist+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*ilist+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

    /* At this point, lispoi contains all the points of the half-ball of p0, rotated
       so that t_{p_0}S = [z = 0] */

    /* Rotated tangent vector (trot[2] = 0), and third direction */
    trot[0] = r[0][0]*t[0] + r[0][1]*t[1] + r[0][2]*t[2];
    trot[1] = r[1][0]*t[0] + r[1][1]*t[1] + r[1][2]*t[2];

    u[0] = -trot[1];
    u[1] =  trot[0];

    /* Travel half-ball at p0 and stop at first triangle containing u */
    for (k=0; k<ilist; k++) {
      detg = lispoi[3*k+1]*u[1] - lispoi[3*k+2]*u[0];
      detd = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
      if ( detg > 0.0 && detd > 0.0 )  break;
    }

    /* If triangle not found, try with -u */
    if ( k == ilist ) {
      u[0] *= -1.0;
      u[1] *= -1.0;

      for (k=0; k<ilist; k++) {
        detg = lispoi[3*k+1]*u[1] - lispoi[3*k+2]*u[0];
        detd = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
        if ( detg > 0.0 && detd > 0.0 )  break;
      }
    }
    if ( k == ilist )  continue;

    iel  = list[k] / 4;
    ifac = list[k] % 4;
    pt = &mesh->tetra[iel];
    for ( i0=0; i0!=3; ++i0 ) {
      if ( pt->v[_MMG5_idir[ifac][i0]]==idp ) break;
    }
    assert(i0!=3);
    i1  = _MMG5_inxt2[i0];
    i2  = _MMG5_iprv2[i0];

    _MMG5_tet2tri(mesh,iel,ifac,&ptt);
    assert(pt->xt);
    pxt = &mesh->xtetra[pt->xt];
    if ( !_MMG5_bezierCP(mesh,&ptt,&b,MG_GET(pxt->ori,i)) )  continue;

    /* Barycentric coordinates of vector u in tria iel */
    detg = lispoi[3*k+1]*u[1] - lispoi[3*k+2]*u[0];
    detd = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
    det = detg + detd;
    if ( det < _MMG5_EPSD )  continue;

    det = 1.0 / det;
    bcu[0] = 0.0;
    bcu[1] = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
    bcu[1] *= det;
    assert(bcu[1] <= 1.0);
    bcu[2] = 1.0 - bcu[1];

    /* Computation of tangent vector and second derivative of curve t \mapsto
       b(tbcu) (not in rotated frame) */
    m[i+1] = MG_MAX(m[i+1],
                    _MMG5_ridSizeInNormalDir(mesh,i0,bcu,&b,isqhmin,isqhmax));
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param kel index of the triangle in which we work.
 * \param iface face of the tetra on which we work.
 * \param ip index of the point on which we want to compute the metric
 * (in tetra \a kel).
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a REF vertex of the mesh, associated to the
 * geometric approx of the surface.
 *
 */
static int _MMG5_defmetref(MMG5_pMesh mesh,MMG5_pSol met,int kel, int iface, int ip) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_Tria     ptt;
  MMG5_pPoint   p0,p1;
  MMG5_pxPoint  px0;
  _MMG5_Bezier  b;
  MMG5_pPar     par;
  int           lists[_MMG5_LMAX+2],listv[_MMG5_LMAX+2],ilists,ilistv,ilist;
  int           k,iel,ipref[2],idp,ifac;
  double        *m,isqhmin,isqhmax,*n,r[3][3],lispoi[3*_MMG5_LMAX+1];
  double        ux,uy,uz,det2d,c[3];
  double        tAA[6],tAb[3], hausd, hausdloc;
  unsigned char i1,i2,itri1,itri2,i;

  ipref[0] = ipref[1] = 0;
  pt  = &mesh->tetra[kel];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  ilist = _MMG5_boulesurfvolp(mesh,kel,ip,iface,listv,&ilistv,lists,&ilists,0);

  if ( ilist!=1 ) {
    printf("Error; unable to compute the ball af the point %d.\n", idp);
    printf("Exit program.\n");
    exit(EXIT_FAILURE);
  }

  isqhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  /* Computation of the rotation matrix T_p0 S -> [z = 0] */
  assert( p0->xp && !MG_SIN(p0->tag) && MG_EDG(p0->tag) && !(MG_NOM & p0->tag) );
  px0 = &mesh->xpoint[p0->xp];

  n  = &px0->n1[0];
  assert ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] > _MMG5_EPSD2 );

  _MMG5_rotmatrix(n,r);
  m = &met->m[6*idp];

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    ifac  = lists[k] % 4;
    pt    = &mesh->tetra[iel];
    assert(pt->xt && (mesh->xtetra[pt->xt].ftag[ifac] & MG_BDY) );
    pxt   = &mesh->xtetra[pt->xt];

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0    = _MMG5_idir[ifac][i];
    itri1 = _MMG5_inxt2[i];
    i1    = _MMG5_idir[ifac][itri1];
    itri2 = _MMG5_iprv2[i];
    i2    = _MMG5_idir[ifac][itri2];
    p1    = &mesh->point[pt->v[i1]];

    /* Store the two ending points of ref curves */
    if ( MG_REF & pxt->tag[_MMG5_iarf[ifac][itri1]] ) {
      if ( !ipref[0] ) {
        ipref[0] = pt->v[i2];
      }
      else if ( !ipref[1] && (pt->v[i2] != ipref[0]) ) {
        ipref[1] = pt->v[i2];
      }
      else if ( (pt->v[i2] != ipref[0]) && (pt->v[i2] != ipref[1]) ) {
        printf("Problem (func defmetref) : three adjacent ref at a non singular point\n");
        exit(EXIT_FAILURE);
      }
    }

    if ( MG_REF & pxt->tag[_MMG5_iarf[ifac][itri2]] ) {
      if ( !ipref[0] ) {
        ipref[0] = pt->v[i1];
      }
      else if ( !ipref[1] && (pt->v[i1] != ipref[0]) ) {
        ipref[1] = pt->v[i1];
      }
      else if ( (pt->v[i1] != ipref[0]) && (pt->v[i1] != ipref[1]) ) {
        printf("Problem (func defmetref) : three adjacent ref at a non singular point\n");
        exit(EXIT_FAILURE);
      }
    }

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
  }

  /* list goes modulo ilist */
  lispoi[3*ilists+1] =  lispoi[1];
  lispoi[3*ilists+2] =  lispoi[2];
  lispoi[3*ilists+3] =  lispoi[3];

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    assert(det2d);
    if ( det2d <= 0.0 ) {
      printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
      return(0);
    }
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  assert(det2d);
  if ( det2d <= 0.0 ) {
    printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
    return(0);
  }
  assert(ipref[0] && ipref[1]);

  /* At this point, lispoi contains all the points of the ball of p0, rotated so
     that t_{p_0}S = [z = 0], ipref1 and ipref2 are the indices of other ref
     points. */

  /* Second step : reconstitution of the curvature tensor at p0 in the tangent
     plane, with a quadric fitting approach */
  memset(tAA,0.0,6*sizeof(double));
  memset(tAb,0.0,3*sizeof(double));

  hausd = -1.;
  for (k=0; k<ilists; k++) {
    /* Approximation of the curvature in the normal section associated to tau :
       by assumption, p1 is either regular, either on a ridge (or a
       singularity), but p0p1 is not ridge*/
    iel  = lists[k] / 4;
    ifac = lists[k] % 4;
    pt  = &mesh->tetra[iel];
    assert(pt->xt);
    pxt = &mesh->xtetra[pt->xt];

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0  = _MMG5_idir[ifac][i];
    // i1  = _MMG5_idir[ifac][_MMG5_inxt2[i]];

    _MMG5_tet2tri(mesh,iel,ifac,&ptt);

    _MMG5_bezierCP(mesh,&ptt,&b,MG_GET(pxt->ori,ifac));

    /* 1. Fill matrice tAA and second member tAb with A=(\sum X_{P_i}^2 \sum
     * Y_{P_i}^2 \sum X_{P_i}Y_{P_i}) and b=\sum Z_{P_i} with P_i the physical
     * points at edge [i0;i1] extremities and middle.
     * 2. Compute the physical coor \a c of the curve edge's
     * mid-point.
     */
    _MMG5_fillDefmetregSys(k,p0,i,b,r,c,lispoi,tAA,tAb);

    /* local hausdorff */
    hausdloc = -1.;
    for (i=0; i<mesh->info.npar; i++) {
      par = &mesh->info.par[i];
      if ( (par->elt == MMG5_Triangle) && (pxt->ref[ifac] == par->ref ) )
        hausdloc = par->hausd;
    }
    if ( hausdloc > 0. ) {
      if ( hausd > 0. )
        hausd = MG_MIN(hausd,hausdloc);
      else
        hausd = hausdloc;
    }
  }
  if ( hausd <= 0. ) hausd = mesh->info.hausd;

  /* Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) */
  return(_MMG5_solveDefmetrefSys( mesh, p0, ipref, r, c, tAA, tAb, m,
                                  isqhmin, isqhmax, hausd));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param kel index of the triangle in which we work.
 * \param ip index of the point on which we want to compute the metric
 * in (tetra \a kel).
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a REGULAR vertex of the mesh, associated to
 * the geometric approx of the surface.
 *
 */
static int _MMG5_defmetreg(MMG5_pMesh mesh,MMG5_pSol met,int kel,int iface, int ip) {
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_Tria      ptt;
  MMG5_pPoint    p0,p1;
  MMG5_pxPoint   px0;
  _MMG5_Bezier   b;
  MMG5_pPar      par;
  int            lists[_MMG5_LMAX+2],listv[_MMG5_LMAX+2],ilists,ilistv,ilist;
  int            k,iel,idp,ifac;
  double         *n,*m,r[3][3],ux,uy,uz,lispoi[3*_MMG5_LMAX+1];
  double         det2d,c[3],isqhmin,isqhmax;
  double         tAA[6],tAb[3],hausd, hausdloc;
  unsigned char  i1,i;

  pt  = &mesh->tetra[kel];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  ilist = _MMG5_boulesurfvolp(mesh,kel,ip,iface,listv,&ilistv,lists,&ilists,0);

  if ( ilist!=1 ) {
    printf("Error; unable to compute the ball af the point %d.\n", idp);
    printf("Exit program.\n");
    exit(EXIT_FAILURE);
  }

  isqhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  /* Computation of the rotation matrix T_p0 S -> [z = 0] */
  assert( p0->xp && !MG_SIN(p0->tag) && !MG_EDG(p0->tag) && !(MG_NOM & p0->tag) );
  px0 = &mesh->xpoint[p0->xp];

  n  = &px0->n1[0];
  assert ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] > _MMG5_EPSD2 );

  _MMG5_rotmatrix(n,r);
  m = &met->m[6*idp];

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    ifac  = lists[k] % 4;
    pt = &mesh->tetra[iel];
    assert(pt->xt && (mesh->xtetra[pt->xt].ftag[ifac] & MG_BDY) );

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0 = _MMG5_idir[ifac][i];
    i1 = _MMG5_idir[ifac][_MMG5_inxt2[i]];
    p1 = &mesh->point[pt->v[i1]];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
  }

  /* list goes modulo ilist */
  lispoi[3*ilists+1] =  lispoi[1];
  lispoi[3*ilists+2] =  lispoi[2];
  lispoi[3*ilists+3] =  lispoi[3];

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    assert(det2d);
    if ( det2d <= 0.0 ) {
      printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
      return(0);
    }
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  assert(det2d);
  if ( det2d <= 0.0 ) {
    printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
    return(0);
  }

  /* At this point, lispoi contains all the points of the ball of p0, rotated
     so that t_{p_0}S = [z = 0] */

  /* Second step : reconstitution of the curvature tensor at p0 in the tangent
     plane, with a quadric fitting approach */
  memset(tAA,0.0,6*sizeof(double));
  memset(tAb,0.0,3*sizeof(double));

  hausd = -1.;
  for (k=0; k<ilists; k++) {
    /* Approximation of the curvature in the normal section associated to tau :
       by assumption, p1 is either regular, either on a ridge (or a
       singularity), but p0p1 is not ridge*/
    iel  = lists[k] / 4;
    ifac = lists[k] % 4;
    pt  = &mesh->tetra[iel];
    assert(pt->xt);
    pxt = &mesh->xtetra[pt->xt];

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0  = _MMG5_idir[ifac][i];
    // i1  = _MMG5_idir[ifac][_MMG5_inxt2[i]];

    _MMG5_tet2tri(mesh,iel,ifac,&ptt);

    _MMG5_bezierCP(mesh,&ptt,&b,MG_GET(pxt->ori,ifac));

    /* 1. Fill matrice tAA and second member tAb with A=(\sum X_{P_i}^2 \sum
     * Y_{P_i}^2 \sum X_{P_i}Y_{P_i}) and b=\sum Z_{P_i} with P_i the physical
     * points at edge [i0;i1] extremities and middle.
     * 2. Compute the physical coor \a c of the curve edge's
     * mid-point.
     */
    _MMG5_fillDefmetregSys(k,p0,i,b,r,c,lispoi,tAA,tAb);

    /* local hausdorff */
    hausdloc = -1.;
    for (i=0; i<mesh->info.npar; i++) {
      par = &mesh->info.par[i];
      if ( (par->elt == MMG5_Triangle) && (pxt->ref[ifac] == par->ref ) )
        hausdloc = par->hausd;
    }
    if ( hausdloc > 0. ) {
      if ( hausd > 0. )
        hausd = MG_MIN(hausd,hausdloc);
      else
        hausd = hausdloc;
    }
  }
  if ( hausd <= 0. ) hausd = mesh->info.hausd;

  /* Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) */
  return(_MMG5_solveDefmetregSys( mesh,r, c, tAA, tAb, m, isqhmin, isqhmax,
                                  hausd));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric stucture.
 * \return 0 if fail, 1 otherwise.
 *
 * Define size at points by intersecting the surfacic metric and the
 * physical metric.
 *
 */
int _MMG5_defsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   ppt;
  double        mm[6];
  int           k,l,iploc;
  char          i,ismet;

 if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Defining anisotropic map\n");

  if ( mesh->info.hmax < 0.0 ) {
    //  mesh->info.hmax = 0.5 * mesh->info.delta;
    fprintf(stdout,"%s:%d:Error: negative hmax value.\n",__FILE__,__LINE__);
    return(0);
  }

  ismet = (met->m > 0);

  if ( !met->m ) {
    met->np    = mesh->np;
    met->npmax = mesh->npmax;
    met->size  = 6;
    met->dim   = 3;
    _MMG5_ADD_MEM(mesh,(6*(met->npmax+1))*sizeof(double),"solution",return(0));
    _MMG5_SAFE_MALLOC(met->m,(6*(mesh->npmax+1)),double);

  }

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->flag = 0;
  }

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    // Warning: why are we skipped the tetra with negative refs ?
    if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
    else if ( !pt->xt )  continue;

    pxt = &mesh->xtetra[pt->xt];
    for (l=0; l<4; l++) {
      if ( !(pxt->ftag[l] & MG_BDY) ) continue;
      // In multidomain case, acces the face through a tetra for which it is
      // well oriented.
      if ( !(MG_GET(pxt->ori,l)) ) continue;

      for (i=0; i<3; i++) {
        iploc = _MMG5_idir[l][i];
        ppt   = &mesh->point[pt->v[iploc]];

        if ( ppt->flag || !MG_VOK(ppt) )  continue;
        if ( ismet )  memcpy(mm,&met->m[6*(pt->v[iploc])],6*sizeof(double));

        if ( MG_SIN(ppt->tag) || (ppt->tag & MG_NOM) ) {
          if ( !_MMG5_defmetsin(mesh,met,k,l,iploc) )  continue;
        }
        else if ( ppt->tag & MG_GEO ) {
          if ( !_MMG5_defmetrid(mesh,met,k,l,iploc))  continue;
        }
        else if ( ppt->tag & MG_REF ) {
          if ( !_MMG5_defmetref(mesh,met,k,l,iploc) )  continue;
        }
        else {
          if ( !_MMG5_defmetreg(mesh,met,k,l,iploc) )  continue;
        }
/* A FAIRE */
        // if ( ismet )  intextmet(mesh,met,pt->v[iploc],mm);
        ppt->flag = 1;
      }
    }
  }

  /* search for unintialized metric */
  _MMG5_defUninitSize(mesh,met,ismet);

  return(1);
}

/* Enforces mesh gradation by truncating metric field */
int _MMG5_gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  printf("gradsize_ani not implemented\n");
  return(1);
}
