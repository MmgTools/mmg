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
 * \file mmg3d/elaslag.c
 * \brief Local solverfor linear elasticity for Lagrangian meshing.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmg3d.h"
#define  _MMG5_LAMBDA    142.85
#define  _MMG5_MU        35.71

extern char  ddb;

/** Calculate an estimate of the average (isotropic) length of edges in the mesh */
double _MMG5_estavglen(MMG5_pMesh mesh) {
  MMG5_pTetra    pt;
  MMG5_pPoint    p1,p2;
  int       k,na;
  double    len,lent,dna;
  char      i,i1,i2;
  
  na = 0;
  lent = 0.0;
  
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    for (i=0; i<6; i++) {
      i1 = _MMG5_iare[i][0];
      i2 = _MMG5_iare[i][1];
      
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];
      
      len = (p2->c[0] - p1->c[0])*(p2->c[0] - p1->c[0])
      + (p2->c[1] - p1->c[1])*(p2->c[1] - p1->c[1])
      + (p2->c[2] - p1->c[2])*(p2->c[2] - p1->c[2]);
      
      lent += sqrt(len);
      na++;
    }
  }
  
  dna = (double)na;
  dna = 1.0 / dna;
  lent *= dna;
  
  return(lent);
}

/** Stiffness of the material associated to element k */
int _MMG5_stiffelt(MMG5_pMesh mesh,int k,double* lambda,double* mu) {
  MMG5_pTetra   pt;
  double        cal;
  
  pt = &mesh->tetra[k];
  cal = MG_MAX(_MMG5_EPS,_MMG5_ALPHAD*pt->qual);
  cal = 1.0/cal;
  
  *lambda = _MMG5_LAMBDA;
  *mu = _MMG5_MU;
  
  //*lambda = _MMG5_LAMBDA*cal;
  //*mu = _MMG5_MU*cal;
  
  return(1);
}

/** Invert 3x3 non-symmetric matrix */
int _MMG5_invmatg(double m[9],double mi[9]) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;
  
  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<9; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return(0);
  
  /* compute sub-dets */
  aa = m[4]*m[8] - m[5]*m[7];
  bb = m[5]*m[6] - m[3]*m[8];
  cc = m[3]*m[7] - m[4]*m[6];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < _MMG5_EPSD )  return(0);
  det = 1.0 / det;
  
  mi[0] = aa*det;
  mi[3] = bb*det;
  mi[6] = cc*det;
  mi[1] = (m[2]*m[7] - m[1]*m[8])*det;
  mi[4] = (m[0]*m[8] - m[2]*m[6])*det;
  mi[7] = (m[1]*m[6] - m[0]*m[7])*det;
  mi[2] = (m[1]*m[5] - m[2]*m[4])*det;
  mi[5] = (m[2]*m[3] - m[0]*m[5])*det;
  mi[8] = (m[0]*m[4] - m[1]*m[3])*det;
  
  return(1);
}

/** Perform one iteration of the explicit elastodynamic solver */
int _MMG5_onestep(MMG5_pMesh mesh,double* un,double* unm1,double* tmp,double* mass,double dt) {
  MMG5_pTetra   pt;
  MMG5_pPoint   ppt;
  double        Dp[3][4],m[9],im[9],*a,*b,*c,*d,lambda,mu,div,vol,aux,Aeueloc,dt2,imass;
  double        Aeu[9][12],eu[9][12],Aeue[12];
  int           k,ndl,ndla;
  char          i,j,s,ia,ja,ndlloc,ndlaloc;
  
char jj;
double lm[3][4];
  
  dt2 = dt*dt;
  
  memset(tmp,0,3*(mesh->np+1)*sizeof(double));
  
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    
    a = &mesh->point[pt->v[0]].c[0];
    b = &mesh->point[pt->v[1]].c[0];
    c = &mesh->point[pt->v[2]].c[0];
    d = &mesh->point[pt->v[3]].c[0];
    
    /* Material parameters associated to tetra k */
    if ( !_MMG5_stiffelt(mesh,k,&lambda,&mu) ) return(0);
    
    /* Volume of k */
    vol = 0.5*_MMG5_ATHIRD*_MMG5_det4pt(a,b,c,d);
    
    /* Calculation of base functions of tetra k */
    Dp[0][0] = 1.0 ; Dp[0][1] = 0.0 ; Dp[0][2] = 0.0 ; Dp[0][3] = -1.0;
    Dp[1][0] = 0.0 ; Dp[1][1] = 1.0 ; Dp[1][2] = 0.0 ; Dp[1][3] = -1.0;
    Dp[2][0] = 0.0 ; Dp[2][1] = 0.0 ; Dp[2][2] = 1.0 ; Dp[2][3] = -1.0;
    
    m[0] = a[0] - d[0] ; m[1] = a[1] - d[1] ; m[2] = a[2] - d[2];
    m[3] = b[0] - d[0] ; m[4] = b[1] - d[1] ; m[5] = b[2] - d[2];
    m[6] = c[0] - d[0] ; m[7] = c[1] - d[1] ; m[8] = c[2] - d[2];
    
    if ( !_MMG5_invmatg(m,im) )  return(0);
    
    /* eu[.][j] = strain of the jth basis function in tetra k 
       Numbering of basis functions: j= 0,1,2,3 = component x, j=4,5,6,7 = y, j=8,9,10,11 = z */
    memset(eu,0,9*12*sizeof(double));
    
    /* Basis functions*/
    memset(lm,0,3*4*sizeof(double));
    for (j=0; j<4; j++) {
      for(i=0; i<3; i++) {
        for(s=0; s<3; s++){
          lm[i][j] += im[3*i+s]*Dp[s][j];
        }
      }
    }
    
    /*printf("TEST\n");
    printf("1 0 0 : ---> %f %f %f \n",lm[0][0]*(a[0]-d[0])+ lm[1][0]*(a[1]-d[1]) + lm[2][0]*(a[2]-d[2]),lm[0][0]*(b[0]-d[0])+ lm[1][0]*(b[1]-d[1]) + lm[2][0]*(b[2]-d[2]), lm[0][0]*(c[0]-d[0])+ lm[1][0]*(c[1]-d[1]) + lm[2][0]*(c[2]-d[2]) );
    printf("0 1 0 : ---> %f %f %f \n",lm[0][1]*(a[0]-d[0])+ lm[1][1]*(a[1]-d[1]) + lm[2][1]*(a[2]-d[2]),lm[0][1]*(b[0]-d[0])+ lm[1][1]*(b[1]-d[1]) + lm[2][1]*(b[2]-d[2]), lm[0][1]*(c[0]-d[0])+ lm[1][1]*(c[1]-d[1]) + lm[2][1]*(c[2]-d[2]) );
    printf("0 0 1 : ---> %f %f %f \n",lm[0][2]*(a[0]-d[0])+ lm[1][2]*(a[1]-d[1]) + lm[2][2]*(a[2]-d[2]),lm[0][2]*(b[0]-d[0])+ lm[1][2]*(b[1]-d[1]) + lm[2][2]*(b[2]-d[2]), lm[0][2]*(c[0]-d[0])+ lm[1][2]*(c[1]-d[1]) + lm[2][2]*(c[2]-d[2]) );
    printf("-1 -1 -1 : ---> %f %f %f \n",lm[0][3]*(a[0]-d[0])+ lm[1][3]*(a[1]-d[1]) + lm[2][3]*(a[2]-d[2]),lm[0][3]*(b[0]-d[0])+ lm[1][3]*(b[1]-d[1]) + lm[2][3]*(b[2]-d[2]), lm[0][3]*(c[0]-d[0])+ lm[1][3]*(c[1]-d[1]) + lm[2][3]*(c[2]-d[2]) );
    */
    
    /*for (j=0; j<4; j++) {
      for(s=0; s<3; s++) {
        eu[0][j]   += im[0+s]*Dp[s][j];
        eu[4][j+4] += im[3+s]*Dp[s][j];
        eu[8][j+8] += im[6+s]*Dp[s][j];
      }
      
      eu[1][j+4] = eu[3][j+4] = eu[2][j+8] = eu[6][j+8] = 0.5*eu[0][j];
      eu[1][j] = eu[3][j] = eu[5][j+8] = eu[7][j+8] = 0.5*eu[4][j+4];
      eu[2][j] = eu[6][j] = eu[5][j+4] = eu[7][j+4] = 0.5*eu[8][j+8];
    }*/
    
    for(j=0; j<4; j++) {
      eu[0][j] = lm[0][j];
      eu[1][j] = eu[3][j] = 0.5*lm[1][j];
      eu[2][j] = eu[6][j] = 0.5*lm[2][j];
      
      eu[4][j+4] = lm[1][j];
      eu[3][j+4] = eu[1][j+4] = 0.5*lm[0][j];
      eu[5][j+4] = eu[7][j+4] = 0.5*lm[2][j];
      
      eu[8][j+8] = lm[2][j];
      eu[7][j+8] = eu[5][j+8] = 0.5*lm[1][j];
      eu[6][j+8] = eu[2][j+8] = 0.5*lm[0][j];
    }
    
    /* Aeu[.][j] = components of Aeu(\phi_j) */
    for (j=0; j<12; j++ ) {
      div = eu[0][j] + eu[4][j] + eu[8][j];
      div *= lambda;
      for (s=0; s<9; s++)
        Aeu[s][j] = 2.0*mu*eu[s][j];
      
      Aeu[0][j] += div;
      Aeu[4][j] += div;
      Aeu[8][j] += div;
    }
    
    /* Aeue[l] = \int_K { Ae(u_n):e(\phi_l)}; 
       ia,ja = indices at which the value is added 
       i,j  = summation indices */
    for (ia=0; ia<4; ia++) {
      for (ja=0; ja<3; ja++) {
        ndla = 3*pt->v[ia]+ja;
        ndlaloc = ia+4*ja;
        
        aux = 0.0;
        for (i=0; i<4; i++) {
          for (j=0; j<3; j++) {
            ndl = 3*pt->v[i]+j;
            ndlloc = i+4*j;
            
            Aeueloc = 0.0;
            for (s=0; s<9; s++)
              Aeueloc += Aeu[s][ndlloc]*eu[s][ndlaloc];
            
            Aeueloc *= un[ndl];
            aux += Aeueloc;
          }
        }
        tmp[ndla] += vol*aux;
      }
    }
  }
  
  /* Update: un <- u_{n+1}, unm1 <- un */
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) continue;
    if ( ( ppt->tag & MG_BDY ) ) continue;
    if ( mass[k] < _MMG5_EPSD )
      return(0);
    
    imass = 1.0 / mass[k];
    
    for (j=0; j<3; j++) {
      ndl = 3*k+j;
      aux = unm1[ndl];
      unm1[ndl] = un[ndl];
      un[ndl] = 2.0*un[ndl] - dt2*imass*tmp[ndl] - aux;
    }
  }
  
  return(1);
}

/** Initialize vectors un and unm1 from the input displacement 
    and put mark 1 on imposed Dirichlet boundary conditions */
int _MMG5_inidispelas(MMG5_pMesh mesh,MMG5_pSol disp,double* un,double* unm1) {
  MMG5_pTetra     pt;
  MMG5_pPoint     ppt;
  MMG5_pxTetra    pxt;
  double          *m;
  int             k,ip;
  char            i,j;

  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;
  
  /* For now, it is maybe safer to access to references of boundary points through the faces -> TO BE UPDATED */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];
    
    for (i=0; i<4; i++) {
      if ( (pxt->ftag[i] & MG_BDY) && (pxt->ref[i] == _MMG5_DISPREF) ) {
        for(j=0; j<3; j++) {
          ip = pt->v[_MMG5_idir[i][j]];
          ppt = &mesh->point[ip];
          ppt->flag = 1;
          m = &disp->m[3*ip+0];
          memcpy(&un[3*ip+0],m,3*sizeof(double));
          memcpy(&unm1[3*ip+0],m,3*sizeof(double));
        }
      }
    }
  }
  
  return(1);
}

/** (Lumped) mass matrix; the coefficients are the same for every 3 dl associated
 to a given vertex, and the common value is stored only once */
int _MMG5_lumpmass(MMG5_pMesh mesh,double* mass) {
  MMG5_pTetra     pt;
  MMG5_pPoint     p0,p1;
  int             k,ip0,ip1;
  double          *a,*b,*c,*d,vols10,vols20;
  char            i,j;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    
    a = &mesh->point[pt->v[0]].c[0];
    b = &mesh->point[pt->v[1]].c[0];
    c = &mesh->point[pt->v[2]].c[0];
    d = &mesh->point[pt->v[3]].c[0];
    
    /* Volume / 10 and / 20 */
    vols10 = 0.05*_MMG5_ATHIRD*_MMG5_det4pt(a,b,c,d);
    vols20 = 0.5*vols10;
    
    for (i=0; i<4; i++) {
      ip0 = pt->v[i];
      p0 = &mesh->point[ip0];
      if ( p0->tag & MG_BDY ) continue;
      
      for (j=0; j<4; j++) {
        if ( i == j ) {
          mass[ip0] += vols10;
          continue;
        }
        ip1 = pt->v[j];
        p1 = &mesh->point[ip1];
        if ( p1->tag & MG_BDY ) continue;
        mass[ip0] += vols20;
      }
    }
  }
  return(1);
}

/** Perform the explicit resolution of the elastodynamic extension equation */
int _MMG5_elaslag(MMG5_pMesh mesh,MMG5_pSol disp) {
  MMG5_pTetra  pt;
  double       *un,*unm1,*tmp,*mass,dt,avlen,cp;
  int          k,it,maxit;
  char         i,ier,data[128];
  
  it    = 1;
  maxit = 10000;
  ier   = 1;
  
  /* Initialization of displacement vectors */
  mass   = (double*)calloc(mesh->np+1,sizeof(double));
  tmp    = (double*)calloc(3*(mesh->np+1),sizeof(double));
  un     = (double*)calloc(3*(mesh->np+1),sizeof(double));
  unm1   = (double*)calloc(3*(mesh->np+1),sizeof(double));
  
  /* Initialize displacement */
  if ( !_MMG5_inidispelas(mesh,disp,un,unm1) ) {
    printf("** Error function inidispelas: unable to initialize displacement.\n");
    return(0);
  }
  
  /* Creation of (lumped) mass matrix; only the diagonal coefficients are stored */
  if ( !_MMG5_lumpmass(mesh,mass) ) {
    printf("** Error function lumpmass: unable to create mass matrix.\n");
    return(0);
  }
  
  /* Time step dimensioning */
  avlen = _MMG5_estavglen(mesh);
  cp = 2.0*_MMG5_MU + _MMG5_LAMBDA;  // velocity of a pressure wave
  cp = sqrt(cp);
  dt = 0.01*avlen / cp;
  
  /* Main loop */
  while ( it <= maxit ) {
    if ( !_MMG5_onestep(mesh,un,unm1,tmp,mass,dt) ) {
      ier = -1;
      break;
    }
    it++;
  }
  
  free(mass);
  free(tmp);
  free(unm1);
  
  if ( !ier ) {
    free (un);
    return(0);
  }
  
  memset(&disp->m[0],0,3*(disp->npmax+1)*sizeof(double));
  memcpy(&disp->m[0],un,3*(mesh->np+1)*sizeof(double));
  
  /* Save displacement */
  strcpy(data,"test.sol");
  disp->nameout = &data[0];
  
  if ( !MMG5_saveMet(mesh,disp) ) {
    printf("** Unable to save displacement.");
    return(0);
  }
  
  free(un);
  exit(0);
  
  
  return(1);
}

