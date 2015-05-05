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




/* Compute the intersected (2 x 2) metric from metrics m and n : take simultaneous reduction,
   and proceed to truncation in sizes */
static int intersecmet22(MMG5_pMesh mesh, double *m,double *n,double *mr) {
  double  det,imn[4],dd,sqDelta,trimn,lambda[2],vp0[2],vp1[2],dm[2],dn[2],vnorm,d0,d1,ip[4];
  double  isqhmin,isqhmax;

  isqhmin  = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax  = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  /* Compute imn = M^{-1}N */
  det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < _MMG5_EPS*_MMG5_EPS ) {
    printf("  ## Function intersecmet : null metric det : %E \n",det);
    return(0);
  }
  det = 1.0 / det;

  imn[0] = det * ( m[2]*n[0] - m[1]*n[1]);
  imn[1] = det * ( m[2]*n[1] - m[1]*n[2]);
  imn[2] = det * (-m[1]*n[0] + m[0]*n[1]);
  imn[3] = det * (-m[1]*n[1] + m[0]*n[2]);
  dd = imn[0] - imn[3];
  sqDelta = sqrt(fabs(dd*dd + 4.0*imn[1]*imn[2]));
  trimn = imn[0] + imn[3];

  lambda[0] = 0.5 * (trimn - sqDelta);
  if ( lambda[0] < 0.0 ) {
    printf(" ## Eigenvalues : %f \n",lambda[0]);
    return(0);
  }

  /* First case : matrices m and n are homothetic : n = lambda0*m */
  if ( sqDelta < _MMG5_EPS ) {
    /* Diagonalize m and truncate eigenvalues : trimn, det, etc... are reused */
    if (fabs(m[1]) < _MMG5_EPS) {
      dm[0]   = m[0];
      dm[1]   = m[2];
      vp0[0] = 1;
      vp0[1] = 0;
      vp1[0] = 0;
      vp1[1] = 1;
    } else {
      dd    = m[0] - m[2];
      trimn = m[0] + m[2];
      det   = m[0]*m[2] - m[1]*m[1];
      
      sqDelta = sqrt(fabs(dd*dd +4*0*m[1]*m[1]));
      dm[0]   = 0.5 * (trimn + sqDelta);
      dm[1]   = 0.5 * (trimn - sqDelta);
      
      vp0[0] = m[1];
      vp0[1] = (dm[0]-m[0]);
      vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);
      if ( vnorm < _MMG5_EPS ) {
        vp0[0] = (dm[0] - m[2]);
        vp0[1] = m[1];
        vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);
        
        if ( vnorm < _MMG5_EPS ) return(0);
      }
      
      vnorm   = 1.0 / vnorm;
      vp0[0] *= vnorm;
      vp0[1] *= vnorm;
      
      vp1[0] = m[1];
      vp1[1] = (dm[1]-m[0]);
      vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);
      
      if ( vnorm < _MMG5_EPS ) {
        vp1[0] = (dm[1] - m[2]);
        vp1[1] = m[1];
        vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);

        if ( vnorm < _MMG5_EPS ) return(0);
      }

      vnorm   = 1.0 / vnorm;
      vp1[0] *= vnorm;
      vp1[1] *= vnorm;
    }
    /* Eigenvalues of the resulting matrix*/
    dn[0] = MG_MAX(dm[0],lambda[0]*dm[0]);
    dn[0] = MG_MIN(isqhmin,MG_MAX(isqhmax,dn[0]));
    dn[1] = MG_MAX(dm[1],lambda[0]*dm[1]);
    dn[1] = MG_MIN(isqhmin,MG_MAX(isqhmax,dn[1]));

    /* Intersected metric = P diag(d0,d1){^t}P, P = (vp0, vp1) stored in columns */
    mr[0] = dn[0]*vp0[0]*vp0[0] + dn[1]*vp1[0]*vp1[0];
    mr[1] = dn[0]*vp0[0]*vp0[1] + dn[1]*vp1[0]*vp1[1];
    mr[2] = dn[0]*vp0[1]*vp0[1] + dn[1]*vp1[1]*vp1[1];

    return(1);
  }

  /* Second case : both eigenvalues of imn are distinct ; theory says qf associated to m and n
     are diagonalizable in basis (vp0, vp1) - the coreduction basis */
  else {
    lambda[1] = 0.5 * (trimn + sqDelta);
    assert(lambda[1] >= 0.0);

    vp0[0] = imn[1];
    vp0[1] = (lambda[0] - imn[0]);
    vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);

    if ( vnorm < _MMG5_EPS ) {
      vp0[0] = (lambda[0] - imn[3]);
      vp0[1] = imn[2];
      vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);
    }

    vnorm   = 1.0 / vnorm;
    vp0[0] *= vnorm;
    vp0[1] *= vnorm;

    vp1[0] = imn[1];
    vp1[1] = (lambda[1] - imn[0]);
    vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);

    if ( vnorm < _MMG5_EPS ) {
      vp1[0] = (lambda[1] - imn[3]);
      vp1[1] = imn[2];
      vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);
    }

    vnorm   = 1.0 / vnorm;
    vp1[0] *= vnorm;
    vp1[1] *= vnorm;

    /* Compute diagonal values in simultaneous reduction basis */
    dm[0] = m[0]*vp0[0]*vp0[0] + 2.0*m[1]*vp0[0]*vp0[1] + m[2]*vp0[1]*vp0[1];
    dm[1] = m[0]*vp1[0]*vp1[0] + 2.0*m[1]*vp1[0]*vp1[1] + m[2]*vp1[1]*vp1[1];
    dn[0] = n[0]*vp0[0]*vp0[0] + 2.0*n[1]*vp0[0]*vp0[1] + n[2]*vp0[1]*vp0[1];
    dn[1] = n[0]*vp1[0]*vp1[0] + 2.0*n[1]*vp1[0]*vp1[1] + n[2]*vp1[1]*vp1[1];

    /* Diagonal values of the intersected metric */
    d0 = MG_MAX(dm[0],dn[0]);
    d0 = MG_MIN(isqhmin,MG_MAX(d0,isqhmax));

    d1 = MG_MAX(dm[1],dn[1]);
    d1 = MG_MIN(isqhmin,MG_MAX(d1,isqhmax));

    /* Intersected metric = tP^-1 diag(d0,d1)P^-1, P = (vp0, vp1) stored in columns */
    det = vp0[0]*vp1[1] - vp0[1]*vp1[0];
    if ( fabs(det) < _MMG5_EPS )  return(0);
    det = 1.0 / det;

    ip[0] =  vp1[1]*det;
    ip[1] = -vp1[0]*det;
    ip[2] = -vp0[1]*det;
    ip[3] =  vp0[0]*det;

    mr[0] = d0*ip[0]*ip[0] + d1*ip[2]*ip[2];
    mr[1] = d0*ip[0]*ip[1] + d1*ip[2]*ip[3];
    mr[2] = d0*ip[1]*ip[1] + d1*ip[3]*ip[3];
  }
  return(1);
}

/* Intersect metric held in np (supported in tangent plane of np) with 3*3 metric in me */
int intextmet(MMG5_pMesh mesh,MMG5_pSol met,int np,double me[6]) {
  MMG5_pPoint         p0;
  MMG5_pxPoint          go;
  double         hu,isqhmin,isqhmax,dd,alpha1,alpha2,alpha3;
  double        *m,*n,*n1,*n2,*t,r[3][3],mrot[6],mr[3],mtan[3],metan[3],u[3],a[4];
  double complex ro[3];
  char           i;

  isqhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  p0 = &mesh->point[np];
  m  = &met->m[6*np];

  /* Case of a singular point : take smallest size prescribed by met, or me in every direction */
  if ( MS_SIN(p0->tag) ) {
    /* Characteristic polynomial of me */
    a[3] = -1.0;
    a[2] = me[0]+me[3]+me[5];
    a[1] = -(me[0]*me[3]+me[0]*me[5]+me[3]*me[5]) + (me[1]*me[1]+me[2]*me[2]+me[4]*me[4]);
    a[0] = me[0]*(me[3]*me[5]-me[4]*me[4]) -me[1]*(me[1]*me[5]-me[2]*me[4]) \
      + me[2]*(me[1]*me[4]-me[2]*me[3]);

    rootDeg3(a,ro);
    hu = m[0];
    for(i=0; i<3; i++) {
      if( cimag(ro[i]) != 0.0 )
        break;
      else
        hu = MG_MAX(hu,creal(ro[i]));
    }
    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[0] = hu;
    m[3] = hu;
    m[5] = hu;
  }
  /* Case of a ridge point : take sizes in 3 directions t,n1,u */
  else if ( p0->tag & MG_GEO ) {
    /* Size prescribed by metric me in direction t */
    t = &p0->n[0];
    hu = me[0]*t[0]*t[0] + me[3]*t[1]*t[1] + me[5]*t[2]*t[2] \
      + 2.0*me[1]*t[0]*t[1] + 2.0*me[2]*t[0]*t[2] + 2.0*me[4]*t[1]*t[2];

    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[0] = MG_MAX(m[0],hu);

    /* Size prescribed by metric me in direction u1 = n1 ^ t */
    assert ( p0->xp );
    go = &mesh->xpoint[p0->xp];
    n1 = &go->n1[0];
    n2 = &go->n2[0];

    u[0] = n1[1]*t[2] - n1[2]*t[1];
    u[1] = n1[2]*t[0] - n1[0]*t[2];
    u[2] = n1[0]*t[1] - n1[1]*t[0];
    dd = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    if ( dd < _MMG5_EPSD ) return(0);
    dd = 1.0 / sqrt(dd);

    u[0] *= dd;
    u[1] *= dd;
    u[2] *= dd;

    hu = me[0]*u[0]*u[0] + me[3]*u[1]*u[1] + me[5]*u[2]*u[2] \
      + 2.0*me[1]*u[0]*u[1] + 2.0*me[2]*u[0]*u[2] + 2.0*me[4]*u[1]*u[2];

    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[1] = MG_MAX(m[1],hu);

    /* Size prescribed by metric me in direction u1 = n1 ^ t */
    u[0] = n2[1]*t[2] - n2[2]*t[1];
    u[1] = n2[2]*t[0] - n2[0]*t[2];
    u[2] = n2[0]*t[1] - n2[1]*t[0];
    dd = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    if ( dd < _MMG5_EPSD ) return(0);
    dd = 1.0 / sqrt(dd);

    u[0] *= dd;
    u[1] *= dd;
    u[2] *= dd;

    hu =     me[0]*u[0]*u[0] +     me[3]*u[1]*u[1] +     me[5]*u[2]*u[2] \
      + 2.0*me[1]*u[0]*u[1] + 2.0*me[2]*u[0]*u[2] + 2.0*me[4]*u[1]*u[2];

    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[2] = MG_MAX(m[2],hu);
  }
  /* Case of a ref, or regular point : intersect metrics in tangent plane */
  else {
    n = &p0->n[0];
    _MMG5_rotmatrix(n,r);

    /* Expression of both metrics in tangent plane */
    _MMG5_rmtr(r,m,mrot);
    mtan[0] = mrot[0];
    mtan[1] = mrot[1];
    mtan[2] = mrot[3];


    _MMG5_rmtr(r,me,mrot);
    metan[0] = mrot[0];
    metan[1] = mrot[1];
    metan[2] = mrot[3];
    /* printf("metsurf %e %e %e\n",mtan[0],mtan[1],mtan[2]); */
    /* printf("metphys %e %e %e\n",metan[0],metan[1],metan[2]); */

    /* Intersection of metrics in the tangent plane */
    if ( !intersecmet22(mesh,mtan,metan,mr) ) {
      fprintf(stdout,"WARNING IMPOSSIBLE INTERSECTION : SURFACIC METRIC SKIPPED \n");
      m[0] = me[0];
      m[1] = me[1];
      m[2] = me[2];
      m[3] = me[3];
      m[4] = me[4];
      m[5] = me[5];

      return(0);
    }

    /* Back to the canonical basis of \mathbb{R}^3 : me = ^tR*mr*R : mtan and metan are reused */
    mtan[0]  = mr[0]*r[0][0] + mr[1]*r[1][0] + r[2][0]*mrot[2];
    mtan[1]  = mr[0]*r[0][1] + mr[1]*r[1][1] + r[2][1]*mrot[2];
    mtan[2]  = mr[0]*r[0][2] + mr[1]*r[1][2] + r[2][2]*mrot[2];
    metan[0] = mr[1]*r[0][0] + mr[2]*r[1][0] + r[2][0]*mrot[4];
    metan[1] = mr[1]*r[0][1] + mr[2]*r[1][1] + r[2][1]*mrot[4];
    metan[2] = mr[1]*r[0][2] + mr[2]*r[1][2] + r[2][2]*mrot[4];
    alpha1 = r[0][0]*mrot[2] + r[1][0]*mrot[4] + r[2][0]*mrot[5];
    alpha2 = r[0][1]*mrot[2] + r[1][1]*mrot[4] + r[2][1]*mrot[5];
    alpha3 = r[0][2]*mrot[2] + r[1][2]*mrot[4] + r[2][2]*mrot[5];

    m[0] = r[0][0] * mtan[0] + r[1][0] * metan[0] + r[2][0]*alpha1;
    m[1] = r[0][0] * mtan[1] + r[1][0] * metan[1] + r[2][0]*alpha2;
    m[2] = r[0][0] * mtan[2] + r[1][0] * metan[2] + r[2][0]*alpha3;
    m[3] = r[0][1] * mtan[1] + r[1][1] * metan[1] + r[2][1]*alpha2;
    m[4] = r[0][1] * mtan[2] + r[1][1] * metan[2] + r[2][1]*alpha3;
    m[5] = r[0][2] * mtan[2] + r[1][2] * metan[2] + r[2][2]*alpha3;

  }

  return(1);
}

