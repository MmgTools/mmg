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
 * \file mmg2d/anisosiz_2d.c
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

/* Impose default metric (isotropic, with size hmax) at vertex i in triangle k */
int _MMG2_defaultmet_2d(MMG5_pMesh mesh,MMG5_pSol met,int k,char i) {
  MMG5_pTria       pt;
  double           *m,isqhmax;
  int              ip;
  
  isqhmax = mesh->info.hmax;
  isqhmax = 1.0 / (isqhmax*isqhmax);
  pt = &mesh->tria[k];
  ip = pt->v[i];
  m = &met->m[3*ip];
  
  m[0] = isqhmax;
  m[1] = 0.0;
  m[2] = isqhmax;
  
  return(1);
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param k index of the tria in which we work
 * \param i index of the point on which we want to compute the metric
 *
 * \return 1 if success, 0 if fail
 *
 * Calculate anisotropic metric tensor at (boundary) vertex i in triangle k on
 * account of geometric approximation of the corresponding curve
 *
 */
int _MMG2_defmetbdy_2d(MMG5_pMesh mesh,MMG5_pSol met,int k,char i) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0,p1,p2;
  double          hausd,sqhmin,sqhmax,ux,uy,ll,li,ps1,ps2,lm,ltmp,pv,M1,M2,t1[2],t2[2],b1[2],b2[2],*n,*m;
  double          gpp1[2],gpp2[2];
  int             ilist,iel,ip,ip1,ip2,it[2],l,list[MMG2_LONMAX+2];
  char            i0,i1,i2,j;
  static char     mmgWarn0=0,mmgWarn1=0;
  
  sqhmin   = mesh->info.hmin*mesh->info.hmin;
  sqhmax   = mesh->info.hmax*mesh->info.hmax;
  hausd  = mesh->info.hausd;
  
  pt = &mesh->tria[k];
  ip = pt->v[i];
  p0 = &mesh->point[ip];
  m = &met->m[3*ip];
  
  ip1 = ip2 = 0;
  ilist = _MMG2_boulet(mesh,k,i,list);
  
  /* Recover the two boundary edges meeting at ip */
  for (l=0; l<ilist; l++) {
    iel = list[l] / 3;
    pt = &mesh->tria[iel];
    
    i0 = list[l] % 3;
    i1 = _MMG5_inxt2[i0];
    i2 = _MMG5_iprv2[i0];
    
    if ( MG_EDG(pt->tag[i1]) ) {
      if ( ip1 == 0 ) {
        ip1 = pt->v[i2];
        it[0] = 3*iel+i1;
      }
      else if ( ip1 != pt->v[i2] ) {
        if ( ip2 == 0 ) {
          ip2 = pt->v[i2];
          it[1] = 3*iel+i1;
        }
        else if ( ip2 != pt->v[i2] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Warning: %s: at least 1 point at the"
                    " intersection of 3 edges. abort.\n",__func__);
          }
          return 0;
        }
      }
    }
    
    if ( MG_EDG(pt->tag[i2]) ) {
      if ( ip1 == 0 ) {
        ip1 = pt->v[i1];
        it[0] = 3*iel+i2;
      }
      else if ( ip1 != pt->v[i1] ) {
        if ( ip2 == 0 ) {
          ip2 = pt->v[i1];
          it[1] = 3*iel+i2;
        }
        else if ( ip2 != pt->v[i1] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Warning: %s: at least 1 point at the"
                    " intersection of 3 edges. abort.\n",__func__);
          }
          return 0;
        }
      }
    }
  }
  
  /* Check that there are exactly two boundary points connected at p0 */
  if ( ip1 == 0 || ip2 == 0 ) {
    if ( !mmgWarn1 ) {
      mmgWarn1 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 point that is not"
              "at the intersection of 2 edges. abort.\n",__func__);
    }
    return 0;
  }
  
  lm = sqhmax;
  
  /* Curvature of both boundary edges meeting at ip */
  for (j=0; j<2; j++) {
    iel = it[j] / 3;
    pt = &mesh->tria[iel];
    i0 = it[j] % 3;
    i1 = _MMG5_inxt2[i0];
    i2 = _MMG5_iprv2[i0];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];
    
    p1 = &mesh->point[ip1];
    p2 = &mesh->point[ip2];
    
    ux = p2->c[0] - p1->c[0];
    uy = p2->c[1] - p1->c[1];
    ll = ux*ux + uy*uy;
    if ( ll < _MMG5_EPSD ) continue;
    li = 1.0 / sqrt(ll);
    
    /* Tangent vector at p1 */
    if ( MG_SIN(p1->tag) || p1->tag & MG_NOM ) {
      t1[0] = li*ux;
      t1[1] = li*uy;
    }
    else {
      t1[0] = p1->n[1];
      t1[1] = -p1->n[0];
    }
    
    /* Tangent vector at p2 */
    if ( MG_SIN(p2->tag) || p2->tag & MG_NOM ) {
      t2[0] = li*ux;
      t2[1] = li*uy;
    }
    else {
      t2[0] = p2->n[1];
      t2[1] = -p2->n[0];
    }
    
    /* Calculation of the two Bezier coefficients along the boundary curve */
    ps1   = ux*t1[0] + uy*t1[1];
    b1[0] = p1->c[0] + _MMG5_ATHIRD*ps1*t1[0];
    b1[1] = p1->c[1] + _MMG5_ATHIRD*ps1*t1[1];
    
    ps2   = ux*t2[0] + uy*t2[1];
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
    if ( M1 < _MMG5_EPSD) continue;
    else {
      ltmp = 8.0*hausd / M1;
      lm = MG_MAX(sqhmin,MG_MIN(ltmp,lm));
    }
  }
  
  /* Expression of the metric tensor diag(lm,hmax) (in the (tau, n) basis) in the canonical basis */
  n = &p0->n[0];
  sqhmax = 1.0 / sqhmax;
  lm = 1.0 / lm;
  
  m[0] = lm*n[1]*n[1] + sqhmax*n[0]*n[0];
  m[1] = n[0]*n[1]*(sqhmax-lm);
  m[2] = lm*n[0]*n[0] + sqhmax*n[1]*n[1];
  
  return(1);
}

/* Definition of an anisotropic metric tensor field based on the geometry of the domain;
 this tensor field is intersected by a user-defined tensor field */
int _MMG2_defsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt;
  double         mm[3],mr[3];
  int            k,ip;
  char           ismet,isdef,i;
  
  
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Defining isotropic map\n");
  
  if ( mesh->info.hmax < 0.0 )  mesh->info.hmax = 0.5 * mesh->info.delta;
  
  /* Allocate the structure */
  if ( met->m )
    ismet = 1;
  else {
    ismet = 0;
    met->npmax = mesh->npmax;
    met->np    = mesh->np;
    met->dim = 2;
    met->size  = 3;
    
    _MMG5_ADD_MEM(mesh,3*(met->npmax+1)*sizeof(double),"solution",return(0));
    _MMG5_SAFE_MALLOC(met->m,3*(mesh->npmax+1),double,0);
  }
  
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;
  
  /* Travel all the points (via triangles) in the mesh and set metric tensor */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 ) continue;
    
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      ppt = &mesh->point[ip];
      if ( !MG_VOK(ppt) || ppt->flag ) continue;
      if ( ismet )
        memcpy(mm,&met->m[3*ip],3*sizeof(double));
      
      isdef = 0;
      /* Calculation of a metric tensor depending on the anisotropic features of the mesh */
      /* At a singular point, an isotropic metric with size hmax is defined */
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM ) {
        if ( _MMG2_defaultmet_2d(mesh,met,k,i) ) isdef = 1;
      }
      else if ( MG_EDG(ppt->tag) ) {
        if ( _MMG2_defmetbdy_2d(mesh,met,k,i) ) isdef = 1;
      }
      
      /* If ppt is an interior point, or if it is a boundary point and the special definition of
       a metric tensor has failed, define a default isotropic metric at ppt */
      if ( !isdef ) _MMG2_defaultmet_2d(mesh,met,k,i);
      
      /* If a metric is supplied by the user, intersect it with the geometric one */
      if ( ismet && _MMG5_intersecmet22(mesh,&met->m[3*ip],mm,mr) )
        memcpy(&met->m[3*ip],mr,3*sizeof(double));
      
      ppt->flag = 1;
    }
  }
  
  return(1);
}

/* Perform simultaneous reduction of matrices m1 and m2, and truncate characteristic sizes so that
 the difference between two corresponding sizes is less than difsiz. Return ier, where (ier & 1) if metric m
 is altered, and (ier & 2) if metric n is altered */
int _MMG2_grad2met_ani(MMG5_pMesh mesh,MMG5_pSol met,double *m,double *n,double difsiz) {
  double       det,dd,sqDelta,trimn,vnorm,hm,hn,lambda[2],dm[2],dn[2],imn[4],vp0[2],vp1[2],ip[4];
  char         ier;
  static char  mmgWarn0=0,mmgWarn1=0;

  ier = 0;
  
  /* Simultaneous reduction of m1 and m2 */
  /* Compute imn = M^{-1}N */
  det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < _MMG5_EPS*_MMG5_EPS ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 null metric det : %E \n",
              __func__,det);
    }
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
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 metric with a "
              "negative eigenvalue: %f \n",__func__,lambda[0]);
    }
    return(0);
  }
  
  /* First case : matrices m and n are homothetic: n = lambda0*m */
  if ( sqDelta < _MMG5_EPS ) {
    /* Subase where m is diaonal */
    if ( fabs(m[1]) < _MMG5_EPS ) {
      dm[0]   = m[0];
      dm[1]   = m[2];
      vp0[0] = 1;
      vp0[1] = 0;
      vp1[0] = 0;
      vp1[1] = 1;
    }
    /* Subcase where m is not diagonal; dd,trimn,... are reused */
    else {
      dd    = m[0] - m[2];
      trimn = m[0] + m[2];
      
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
    
    /* Eigenvalues of metric n */
    dn[0] = lambda[0]*dm[0];
    dn[1] = lambda[0]*dm[1];
    
  }
  /* Second case: both eigenvalues of imn are distinct ; theory says qf associated to m and n
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
  }
  
  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the first direction */
  assert ( dm[0] > _MMG5_EPS && dn[0] > _MMG5_EPS );
  hm = 1.0 / sqrt(dm[0]);
  hn = 1.0 / sqrt(dn[0]);
  
  if ( hn > hm+difsiz ) {
    hn = hm+difsiz;
    dn[0] = 1.0 / (hn*hn);
    ier = ier | 2;
  }
  else if ( hm > hn+difsiz ) {
    hm = hn+difsiz;
    dm[0] = 1.0 / (hm*hm);
    ier = ier | 1;
  }
  
  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the second direction */
  assert ( dm[1] > _MMG5_EPS && dn[1] > _MMG5_EPS );
  hm = 1.0 / sqrt(dm[1]);
  hn = 1.0 / sqrt(dn[1]);
  
  if ( hn > hm+difsiz ) {
    hn = hm+difsiz;
    dn[1] = 1.0 / (hn*hn);
    ier = ier | 2;
  }
  else if ( hm > hn+difsiz ) {
    hm = hn+difsiz;
    dm[1] = 1.0 / (hm*hm);
    ier = ier | 1;
  }
  
  /* Update of the metrics = tP^-1 diag(d0,d1)P^-1, P = (vp0, vp1) stored in columns */
  det = vp0[0]*vp1[1] - vp0[1]*vp1[0];
  if ( fabs(det) < _MMG5_EPS )  return(0);
  det = 1.0 / det;
  
  ip[0] =  vp1[1]*det;
  ip[1] = -vp1[0]*det;
  ip[2] = -vp0[1]*det;
  ip[3] =  vp0[0]*det;
  
  if ( ier | 1 ) {
    m[0] = dm[0]*ip[0]*ip[0] + dm[1]*ip[2]*ip[2];
    m[1] = dm[0]*ip[0]*ip[1] + dm[1]*ip[2]*ip[3];
    m[2] = dm[0]*ip[1]*ip[1] + dm[1]*ip[3]*ip[3];
  }
  if ( ier | 2 ) {
    n[0] = dn[0]*ip[0]*ip[0] + dn[1]*ip[2]*ip[2];
    n[1] = dn[0]*ip[0]*ip[1] + dn[1]*ip[2]*ip[3];
    n[2] = dn[0]*ip[1]*ip[1] + dn[1]*ip[3]*ip[3];
  }
  
  return(ier);
  
}

/* Anisotropic mesh gradation routine */
int _MMG2_gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            hgrad,ll,*m1,*m2,difsiz;
  int               k,it,ip1,ip2,maxit,nup,nu;
  char              i,i1,i2,ier;
  
  
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
        
        /* Maximum allowed difference between the prescribed sizes in p1 and p2 */
        difsiz = ll*hgrad;
        
        m1 = &met->m[3*ip1];
        m2 = &met->m[3*ip2];
        
        /* bit 0 of ier = 0 if metric m1 is untouched, 1 otherwise ;  bit 1 of ier = 0 if metric m2 is untouched, 1 otherwise*/
        ier = _MMG2_grad2met_ani(mesh,met,m1,m2,difsiz);
        
        if ( ier & 1 ) {
          p1->flag = mesh->base;
          nu++;
        }
        if ( ier & 2 ) {
          p2->flag = mesh->base;
          nu++;
        }
      }
    }
    nup += nu;
  }
  while ( ++it < maxit && nu > 0 );
  
  if ( abs(mesh->info.imprim) > 4 )  fprintf(stdout,"     gradation: %7d updated, %d iter.\n",nup,it);
  return(1);
  
}
