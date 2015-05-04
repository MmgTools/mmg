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
 * \file mmgs/anisosiz.c
 * \brief Fonctions for anisotropic size map computation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param it index of the triangle in which we work.
 * \param ip index of the point on which we want to compute the metric in \a it.
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a SINGULARITY of the geometry, associated to
 * the geometric approx of the surface. metric \f$=\alpha*Id\f$, \f$\alpha =\f$
 * size.
 *
 */
static int _MMG5_defmetsin(MMG5_pMesh mesh,MMG5_pSol met,int it,int ip) {
  MMG5_pTria         pt;
  MMG5_pPoint        p0;
  double             *m,n[3],isqhmin,isqhmax,b0[3],b1[3],ps1,tau[3];
  double             ntau2,gammasec[3];
  double             c[3],kappa,maxkappa,alpha;
  int                ilist,list[_MMG5_LMAX+2],k,iel,idp;
  unsigned char      i0,i1,i2;

  pt  = &mesh->tria[it];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  ilist = boulet(mesh,it,ip,list);
  if ( !ilist ) {
    printf("Error: unable to compute the ball af the point %d.\n", idp);
    return(0);
  }

  isqhmin  = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax  = 1.0 / (mesh->info.hmax*mesh->info.hmax);
  maxkappa = 0.0;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    i2  = _MMG5_iprv2[i0];
    pt  = &mesh->tria[iel];

    /* Computation of the two control points associated to edge p0p1 with
     * p1=mesh->point[pt->v[i1]]: p0 is singular */
    _MMG5_nortri(mesh,pt,n);
    if ( MG_EDG(pt->tag[i2]) )
      _MMG5_bezierEdge(mesh,idp,pt->v[i1],b0,b1,1,n);
    else
      _MMG5_bezierEdge(mesh,idp,pt->v[i1],b0,b1,0,n);

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
    maxkappa = MG_MAX(kappa,maxkappa);
  }
  alpha = 1.0 / 8.0 * maxkappa / mesh->info.hausd;
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
 * \param it index of the triangle in which we work.
 * \param ip index of the point on which we want to compute the metric in \a it.
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
static int _MMG5_defmetrid(MMG5_pMesh mesh,MMG5_pSol met,int it,int ip) {
  MMG5_pTria     pt;
  MMG5_pPoint    p0,p1,p2;
  _MMG5_Bezier   b;
  int            k,iel,idp,ilist1,ilist2,ilist,*list,list1[_MMG5_LMAX+2];
  int            list2[_MMG5_LMAX+2],iprid[2],ier;
  double         *m,isqhmin,isqhmax,*n1,*n2,*n,*t,trot[2],u[2];
  double         r[3][3],lispoi[3*_MMG5_LMAX+1],ux,uy,uz,det,bcu[3];
  double         detg,detd;
  unsigned char  i,i0,i1,i2;

  pt  = &mesh->tria[it];
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

  ier = bouletrid(mesh,it,ip,&ilist1,list1,&ilist2,list2,&iprid[0],&iprid[1]);
  if ( !ier ) {
    printf("Error: unable to compute the two balls at the ridge point %d.\n", idp);
    return(0);
  }

  /* Specific size in direction of t */
  m[0] = MG_MAX(m[0],_MMG5_ridSizeInTangentDir(mesh,p0,idp,iprid,isqhmin,isqhmax));

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
      iel = list[k] / 3;
      i0  = list[k] % 3;
      i1  = _MMG5_inxt2[i0];
      pt = &mesh->tria[iel];
      p1 = &mesh->point[pt->v[i1]];

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
    p2 = &mesh->point[pt->v[i2]];

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

    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    i2  = _MMG5_iprv2[i0];
    pt = &mesh->tria[iel];
    if ( !_MMG5_bezierCP(mesh,pt,&b,1) )  continue;

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
 * \param it index of the triangle in which we work.
 * \param ip index of the point on which we want to compute the metric in \a it.
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a REF vertex of the mesh, associated to the
 * geometric approx of the surface.
 *
 */
static int _MMG5_defmetref(MMG5_pMesh mesh,MMG5_pSol met,int it,int ip) {
  MMG5_pTria         pt;
  MMG5_pPoint        p0,p1;
  _MMG5_Bezier       b;
  int                ilist,list[_MMG5_LMAX+2],k,iel,ipref[2],idp;
  double             *m,isqhmin,isqhmax,*n,r[3][3],lispoi[3*_MMG5_LMAX+1];
  double             ux,uy,uz,det2d,intm[3],c[3];
  double             tAA[6],tAb[3];
  unsigned char      i0,i1,i2;

  ipref[0] = ipref[1] = 0;
  pt  = &mesh->tria[it];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];
  ilist = boulet(mesh,it,ip,list);
  if ( !ilist ) {
    printf("Error: unable to compute the ball af the point %d.\n", idp);
    return(0);
  }

  isqhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  /* Computation of the rotation matrix T_p0 S -> [z = 0] */
  n  = &mesh->xpoint[p0->xp].n1[0];
  assert ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] > _MMG5_EPSD2 );

  _MMG5_rotmatrix(n,r);
  m = &met->m[6*idp];

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    i2  = _MMG5_iprv2[i0];
    pt = &mesh->tria[iel];
    p1 = &mesh->point[pt->v[i1]];

    /* Store the two ending points of ref curves */
    if ( MG_REF & pt->tag[i1] ) {
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

    if ( MG_REF & pt->tag[i2] ) {
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
  lispoi[3*ilist+1] =  lispoi[1];
  lispoi[3*ilist+2] =  lispoi[2];
  lispoi[3*ilist+3] =  lispoi[3];

  /* Check all projections over tangent plane. */
  for (k=0; k<ilist-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( det2d < 0.0 ) {
      printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
      return(0);
    }
  }
  det2d = lispoi[3*(ilist-1)+1]*lispoi[3*0+2] - lispoi[3*(ilist-1)+2]*lispoi[3*0+1];
  if ( det2d < 0.0 ) {
    printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
    return(0);
  }
  assert(ipref[0] && ipref[1]);

  /* At this point, lispoi contains all the points of the ball of p0, rotated
     so that t_{p_0}S = [z = 0], ipref1 and ipref2 are the indices of other ref points. */

  /* Second step : reconstitution of the curvature tensor at p0 in the tangent plane,
     with a quadric fitting approach */
  memset(intm,0.0,3*sizeof(double));
  memset(tAA,0.0,6*sizeof(double));
  memset(tAb,0.0,3*sizeof(double));

  for (k=0; k<ilist; k++) {
    /* Approximation of the curvature in the normal section associated to tau : by assumption,
       p1 is either regular, either on a ridge (or a singularity), but p0p1 is not ridge*/
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    pt = &mesh->tria[iel];
    _MMG5_bezierCP(mesh,pt,&b,1);


    /* 1. Fill matrice tAA and second member tAb with A=(\sum X_{P_i}^2 \sum
     * Y_{P_i}^2 \sum X_{P_i}Y_{P_i}) and b=\sum Z_{P_i} with P_i the physical
     * points at edge [i0;i1] extremities and middle.
     * 2. Compute the physical coor \a c of the curve edge's
     * mid-point.
     */
    _MMG5_fillDefmetregSys(k,p0,i0,b,r,c,lispoi,tAA,tAb);
  }
  return( _MMG5_solveDefmetrefSys(mesh,p0,ipref,r,c,tAA,tAb,m,
                                  isqhmin,isqhmax,mesh->info.hausd) );
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param it index of the triangle in which we work.
 * \param ip index of the point on which we want to compute the metric in \a it.
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a REGULAR vertex of the mesh, associated to
 * the geometric approx of the surface.
 *
 */
static int _MMG5_defmetreg(MMG5_pMesh mesh,MMG5_pSol met,int it,int ip) {
  MMG5_pTria          pt;
  MMG5_pPoint         p0,p1;
  _MMG5_Bezier        b;
  int                 ilist,list[_MMG5_LMAX+2],k,iel,idp;
  double              *n,*m,r[3][3],ux,uy,uz,lispoi[3*_MMG5_LMAX+1];
  double              det2d,c[3],isqhmin,isqhmax;
  double              tAA[6],tAb[3];
  unsigned char       i0,i1;

  pt  = &mesh->tria[it];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  ilist = boulet(mesh,it,ip,list);
  if ( !ilist ) {
    printf("Error: unable to compute the ball af the point %d.\n", idp);
    return(0);
  }

  isqhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  /* Computation of the rotation matrix T_p0 S -> [z = 0] */
  n  = &p0->n[0];
  assert ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] > _MMG5_EPSD2 );

  if ( !_MMG5_rotmatrix(n,r) ) return(0);
  m = &met->m[6*idp];

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    pt = &mesh->tria[iel];
    p1 = &mesh->point[pt->v[i1]];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
  }

  /* list goes modulo ilist */
  lispoi[3*ilist+1] =  lispoi[1];
  lispoi[3*ilist+2] =  lispoi[2];
  lispoi[3*ilist+3] =  lispoi[3];

  /* Check all projections over tangent plane. */
  for (k=0; k<ilist-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( det2d <= 0.0 ) {
      printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
      return(0);
    }
  }
  det2d = lispoi[3*(ilist-1)+1]*lispoi[3*0+2] - lispoi[3*(ilist-1)+2]*lispoi[3*0+1];
  if ( det2d <= 0.0 ) {
    printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
    return(0);
  }

  /* At this point, lispoi contains all the points of the ball of p0, rotated
     so that t_{p_0}S = [z = 0] */

  /* Second step : reconstitution of the curvature tensor at p0 in the tangent plane,
     with a quadric fitting approach */
  memset(tAA,0.0,6*sizeof(double));
  memset(tAb,0.0,3*sizeof(double));

  for (k=0; k<ilist; k++) {
    /* Approximation of the curvature in the normal section associated to tau : by assumption,
       p1 is either regular, either on a ridge (or a singularity), but p0p1 is not ridge*/
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    pt = &mesh->tria[iel];
    _MMG5_bezierCP(mesh,pt,&b,1);

    /* 1. Fill matrice tAA and second member tAb with A=(\sum X_{P_i}^2 \sum
     * Y_{P_i}^2 \sum X_{P_i}Y_{P_i}) and b=\sum Z_{P_i} with P_i the physical
     * points at edge [i0;i1] extremities and middle.
     * 2. Compute the physical coor \a c of the curve edge's
     * mid-point.
     */
    _MMG5_fillDefmetregSys(k,p0,i0,b,r,c,lispoi,tAA,tAb);
  }

  /* 2. Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) */
  return(_MMG5_solveDefmetregSys( mesh,r, c, tAA, tAb, m, isqhmin, isqhmax,
                                  mesh->info.hausd));
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
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  double        mm[6];
  int           k;
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
    _MMG5_SAFE_CALLOC(met->m,6*(mesh->npmax+1),double);
  }

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->flag = 0;
  }

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( ppt->flag || !MG_VOK(ppt) )  continue;
      if ( ismet )  memcpy(mm,&met->m[6*(pt->v[i])],6*sizeof(double));

      if ( MS_SIN(ppt->tag) ) {
        if ( !_MMG5_defmetsin(mesh,met,k,i) )  continue;
      }
      else if ( ppt->tag & MG_GEO ) {
        if ( !_MMG5_defmetrid(mesh,met,k,i))  continue;
      }
      else if ( ppt->tag & MG_REF ) {
        if ( !_MMG5_defmetref(mesh,met,k,i) )  continue;
      }
      else if ( ppt->tag )  continue;
      else {
        if ( !_MMG5_defmetreg(mesh,met,k,i) )  continue;
      }
/* A FAIRE */
      // if ( ismet )  intextmet(mesh,met,pt->v[i],mm);
      ppt->flag = 1;
    }
  }

  /* search for unintialized metric */
  _MMG5_defUninitSize(mesh,met,ismet);

  return(1);
}

/* Enforces gradation of metric in one extremity of edge i in tria k with respect to the other,
   along the direction of the associated support curve
   Return -1 if no gradation is needed, else index of graded point */
static int grad2met(MMG5_pMesh mesh, MMG5_pSol met, int iel, int i){
  MMG5_pTria    pt;
  MMG5_pPoint   p1,p2;
  double   *mm1,*mm2,*nn1,*nn2,ps1,ps2,ux,uy,uz,m1[6],m2[6],n1[3],n2[3],nt[3];
  double   r1[3][3],r2[3][3],t1[3],t2[3],c[3],mtan1[3],mtan2[3],mr[6]/*,l1,l2*/,l,dd;
  double   lambda[2],vp[2][2],alpha,beta,mu;
  int      np1,np2;
  char     i1,i2,ichg;

  pt = &mesh->tria[iel];

  i1 = _MMG5_inxt2[i];
  i2 = _MMG5_iprv2[i];
  np1 = pt->v[i1];
  np2 = pt->v[i2];

  p1 = &mesh->point[np1];
  p2 = &mesh->point[np2];

  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  uz = p2->c[2] - p1->c[2];

  mm1 = &met->m[6*np1];
  mm2 = &met->m[6*np2];

  if( !_MMG5_nortri(mesh,pt,nt) )
    return(-1);

  /* Recover normal and metric associated to p1 */
  if( MS_SIN(p1->tag) ){
    memcpy(n1,nt,3*sizeof(double));
    memcpy(m1,mm1,6*sizeof(double));
  }
  else if( MG_GEO & p1->tag ){
    nn1 = &mesh->xpoint[p1->xp].n1[0];
    nn2 = &mesh->xpoint[p1->xp].n2[0];
    ps1 = nt[0]*nn1[0] + nt[1]*nn1[1] + nt[2]*nn1[2];
    ps2 = nt[0]*nn2[0] + nt[1]*nn2[1] + nt[2]*nn2[2];

    if( fabs(ps1) < fabs(ps2))
      memcpy(n1,nn2,3*sizeof(double));
    else
      memcpy(n1,nn1,3*sizeof(double));

    if( !_MMG5_buildridmet(mesh,met,np1,ux,uy,uz,m1) )
      return(-1);
  }
  else if( MG_REF & p1->tag ){
    memcpy(n1,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
    memcpy(m1,mm1,6*sizeof(double));
  }
  else{
    memcpy(n1,p1->n,3*sizeof(double));
    memcpy(m1,mm1,6*sizeof(double));
  }

  /* Recover normal and metric associated to p2 */
  if ( MS_SIN(p2->tag) ) {
    memcpy(n2,nt,3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }
  else if ( MG_GEO & p2->tag ) {
    nn1 = &mesh->xpoint[p2->xp].n1[0];
    nn2 = &mesh->xpoint[p2->xp].n2[0];
    ps1 = nt[0]*nn1[0] + nt[1]*nn1[1] + nt[2]*nn1[2];
    ps2 = nt[0]*nn2[0] + nt[1]*nn2[1] + nt[2]*nn2[2];

    if( fabs(ps1) < fabs(ps2))
      memcpy(n2,nn2,3*sizeof(double));
    else
      memcpy(n2,nn1,3*sizeof(double));

    if( !_MMG5_buildridmet(mesh,met,np2,ux,uy,uz,m2) )
      return(-1);
  }
  else if( MG_REF & p2->tag ){
    memcpy(n2,&(mesh->xpoint[p2->xp].n1[0]),3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }
  else{
    memcpy(n2,p2->n,3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }

  /* Rotation matrices mapping n1/n2 to e_3 */
  _MMG5_rotmatrix(n1,r1);
  _MMG5_rotmatrix(n2,r2);

  /* Geodesic length of support curve to edge i */
  ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
  t1[0] = ux - ps1*n1[0];
  t1[1] = uy - ps1*n1[1];
  t1[2] = uz - ps1*n1[2];

  ps2 = - (ux*n2[0] + uy*n2[1] + uz*n2[2]);
  t2[0] = -ux - ps1*n2[0];
  t2[1] = -uy - ps1*n2[1];
  t2[2] = -uz - ps1*n2[2];

  // Commentated because the next line overwrite it... to check!
  /* l1 = m1[0]*t1[0]*t1[0] + m1[3]*t1[1]*t1[1] + m1[5]*t1[2]*t1[2] \ */
  /*   + 2.0 * ( m1[1]*t1[0]*t1[1] + m1[2]*t1[0]*t1[2] + m1[4]*t1[1]*t1[2] ) ; */
  /* l2 = m2[0]*t2[0]*t2[0] + m2[3]*t2[1]*t2[1] + m2[5]*t2[2]*t2[2] \ */
  /*   + 2.0 * ( m2[1]*t2[0]*t2[1] + m2[2]*t2[0]*t2[2] + m2[4]*t2[1]*t2[2] ) ; */
  /* l = 0.5* ( sqrt(l1) + sqrt(l2) ) ; */

  l = sqrt(ux*ux+uy*uy+uz*uz);

  /* Characteristic sizes in direction of support curve */
  _MMG5_rmtr(r1,m1,mr);
  mtan1[0] = mr[0];
  mtan1[1] = mr[1];
  mtan1[2] = mr[3];
  c[0] = r1[0][0]*ux + r1[0][1]*uy + r1[0][2]*uz;
  c[1] = r1[1][0]*ux + r1[1][1]*uy + r1[1][2]*uz;
  c[2] = r1[2][0]*ux + r1[2][1]*uy + r1[2][2]*uz;
  memcpy(t1,c,3*sizeof(double));
  dd = t1[0]*t1[0] + t1[1]*t1[1];
  if(dd < _MMG5_EPSD2)
    return(-1);

  dd = 1.0/sqrt(dd);
  t1[0] *= dd;
  t1[1] *= dd;
  ps1 = mtan1[0]*t1[0]*t1[0] + 2.0*mtan1[1]*t1[0]*t1[1] + mtan1[2]*t1[1]*t1[1];
  ps1 = sqrt(ps1);

  _MMG5_rmtr(r2,m2,mr);
  mtan2[0] = mr[0];
  mtan2[1] = mr[1];
  mtan2[2] = mr[3];
  c[0] = - ( r2[0][0]*ux + r2[0][1]*uy + r2[0][2]*uz );
  c[1] = - ( r2[1][0]*ux + r2[1][1]*uy + r2[1][2]*uz );
  c[2] = - ( r2[2][0]*ux + r2[2][1]*uy + r2[2][2]*uz );
  memcpy(t2,c,3*sizeof(double));

  dd = t2[0]*t2[0] + t2[1]*t2[1];
  if(dd < _MMG5_EPSD2)
    return(-1);

  dd = 1.0/sqrt(dd);
  t2[0] *= dd;
  t2[1] *= dd;
  ps2 = mtan2[0]*t2[0]*t2[0] + 2.0*mtan2[1]*t2[0]*t2[1] + mtan2[2]*t2[1]*t2[1];
  ps2 = sqrt(ps2);

  /* Metric in p1 has to be changed */
  if( ps2 > ps1 ){
    alpha = ps2 /(1.0+mesh->info.hgrad*l*ps2);
    if( ps1 >= alpha -_MMG5_EPS )
      return(-1);

    _MMG5_eigensym(mtan1,lambda,vp);
    c[0] = t1[0]*vp[0][0] + t1[1]*vp[0][1];
    c[1] = t1[0]*vp[1][0] + t1[1]*vp[1][1];

    if( fabs(c[0]) > fabs(c[1]) ){
      ichg  = 0;
      beta = (alpha*alpha - ps1*ps1)/(c[0]*c[0]);
      mu = lambda[0] + beta ;
      mtan1[0] = mu*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
      mtan1[1] = mu*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
      mtan1[2] = mu*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];
    }
    else{
      ichg = 1;
      beta = (alpha*alpha - ps1*ps1)/(c[1]*c[1]);
      mu = lambda[1] + beta;
      mtan1[0] = lambda[0]*vp[0][0]*vp[0][0] + mu*vp[1][0]*vp[1][0];
      mtan1[1] = lambda[0]*vp[0][0]*vp[0][1] + mu*vp[1][0]*vp[1][1];
      mtan1[2] = lambda[0]*vp[0][1]*vp[0][1] + mu*vp[1][1]*vp[1][1];
    }

    /* Metric update */
    if( MS_SIN(p1->tag) ){
      mm1[0] += 0.5*beta;
      mm1[3] += 0.5*beta;
      mm1[5] += 0.5*beta;
    }
    else if( p1->tag & MG_GEO ){
      c[0] = fabs(mm1[0]-lambda[ichg]);
      c[1] = fabs(mm1[1]-lambda[ichg]);
      c[2] = fabs(mm1[2]-lambda[ichg]);
      if( c[0] < c[1] ){
        if( c[0] < c[2] ){
          mm1[0] += beta;
        }
        else{
          mm1[2] += beta;
        }
      }
      else{
        if( c[1] < c[2] ){
          mm1[1] += beta;
        }
        else{
          mm1[2] += beta;
        }
      }
    }
    else{
      /* Reuse t1 and t2 */
      t1[0] = mtan1[0]*r1[0][0] + mtan1[1]*r1[1][0];  t1[1] = mtan1[0]*r1[0][1] + mtan1[1]*r1[1][1];  t1[2] = mtan1[0]*r1[0][2] + mtan1[1]*r1[1][2];
      t2[0] = mtan1[1]*r1[0][0] + mtan1[2]*r1[1][0];  t2[1] = mtan1[1]*r1[0][1] + mtan1[2]*r1[1][1];  t2[2] = mtan1[1]*r1[0][2] + mtan1[2]*r1[1][2];

      m1[0] = r1[0][0]*t1[0] + r1[1][0]*t2[0];
      m1[1] = r1[0][0]*t1[1] + r1[1][0]*t2[1];
      m1[2] = r1[0][0]*t1[2] + r1[1][0]*t2[2];
      m1[3] = r1[0][1]*t1[1] + r1[1][1]*t2[1];
      m1[4] = r1[0][1]*t1[2] + r1[1][1]*t2[2];
      m1[5] = r1[0][2]*t1[2] + r1[1][2]*t2[2];

      memcpy(mm1,m1,6*sizeof(double));
    }
    return(i1);
  }
  /* Metric in p2 has to be changed */
  else{
    alpha = ps1 /(1.0+mesh->info.hgrad*l*ps1);
    if( ps2 >= alpha - _MMG5_EPS)
      return(-1);

    _MMG5_eigensym(mtan2,lambda,vp);
    c[0] = t2[0]*vp[0][0] + t2[1]*vp[0][1];
    c[1] = t2[0]*vp[1][0] + t2[1]*vp[1][1];

    if( fabs(c[0]) > fabs(c[1]) ){
      ichg = 0;
      beta = (alpha*alpha - ps2*ps2)/(c[0]*c[0]);
      mu = lambda[0] + beta;
      mtan2[0] = mu*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
      mtan2[1] = mu*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
      mtan2[2] = mu*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];
    }
    else{
      ichg = 1;
      beta = (alpha*alpha - ps2*ps2)/(c[1]*c[1]);
      mu = lambda[1] + beta;
      mtan2[0] = lambda[0]*vp[0][0]*vp[0][0] + mu*vp[1][0]*vp[1][0];
      mtan2[1] = lambda[0]*vp[0][0]*vp[0][1] + mu*vp[1][0]*vp[1][1];
      mtan2[2] = lambda[0]*vp[0][1]*vp[0][1] + mu*vp[1][1]*vp[1][1];
    }

    /* Metric update */
    if( MS_SIN(p2->tag) ){
      mm2[0] += 0.5*beta;
      mm2[3] += 0.5*beta;
      mm2[5] += 0.5*beta;
    }
    else if( p2->tag & MG_GEO ){
      c[0] = fabs(mm2[0]-lambda[ichg]);
      c[1] = fabs(mm2[1]-lambda[ichg]);
      c[2] = fabs(mm2[2]-lambda[ichg]);
      if( c[0] < c[1] ){
        if( c[0] < c[2] ){
          mm2[0] += beta;
        }
        else{
          mm2[2] += beta;
        }
      }
      else{
        if( c[1] < c[2] ){
          mm2[1] += beta;
        }
        else{
          mm2[2] += beta;
        }
      }
    }
    else{
      /* Reuse t1 and t2 */
      t1[0] = mtan2[0]*r2[0][0] + mtan2[1]*r2[1][0];  t1[1] = mtan2[0]*r2[0][1] + mtan2[1]*r2[1][1];  t1[2] = mtan2[0]*r2[0][2] + mtan2[1]*r2[1][2];
      t2[0] = mtan2[1]*r2[0][0] + mtan2[2]*r2[1][0];  t2[1] = mtan2[1]*r2[0][1] + mtan2[2]*r2[1][1];  t2[2] = mtan2[1]*r2[0][2] + mtan2[2]*r2[1][2];

      m2[0] = r2[0][0]*t1[0] + r2[1][0]*t2[0];
      m2[1] = r2[0][0]*t1[1] + r2[1][0]*t2[1];
      m2[2] = r2[0][0]*t1[2] + r2[1][0]*t2[2];
      m2[3] = r2[0][1]*t1[1] + r2[1][1]*t2[1];
      m2[4] = r2[0][1]*t1[2] + r2[1][1]*t2[2];
      m2[5] = r2[0][2]*t1[2] + r2[1][2]*t2[2];

      memcpy(mm2,m2,6*sizeof(double));
    }

    return(i2);
  }
}

/* Enforces mesh gradation by truncating metric field */
int gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria   pt;
  MMG5_pPoint  p1,p2;
  double  *m,mv;
  int     k,it,nup,nu,maxit;
  char    i,ier,i1,i2;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Anisotropic mesh gradation\n");

  mesh->base = 0;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  /* First step : make ridges iso */
  for (k=1; k<= mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !MG_VOK(p1) ) continue;
    if ( MS_SIN(p1->tag) ) continue;
    if ( !(p1->tag & MG_GEO) ) continue;

    m = &met->m[6*k];
    mv = MG_MAX(m[0],MG_MAX(m[1],m[2]));
    m[0] = mv;
    m[1] = mv;
    m[2] = mv;
  }

  /* Second step : standard gradation procedure */
  it = nup = 0;
  maxit = 100;
  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<3; i++) {
        i1 = _MMG5_inxt2[i];
        i2 = _MMG5_iprv2[i];
        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt->v[i2]];

        if ( p1->flag < mesh->base-1 && p2->flag < mesh->base-1 )  continue;
        ier = grad2met(mesh,met,k,i);
        if ( ier == i1 ) {
          p1->flag = mesh->base;
          nu++;
        }
        else if ( ier == i2 ) {
          p2->flag = mesh->base;
          nu++;
        }
      }
    }
    nup += nu;
  }
  while( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 )  fprintf(stdout,"     gradation: %7d updated, %d iter.\n",nup,it);
  return(1);
}

