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
 * \file common/mmg_tools.c
 * \brief Various tools for the mmg applications.
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
 * \param m pointer toward the first metric to intersect.
 * \param n pointer toward the second metric to intersect.
 * \param mr pointer toward the computed intersected metric.
 * \return 1.
 *
 * Compute the intersected (2 x 2) metric between metrics \a m and \a n,
 * PRESERVING the directions of \a m. Result is stored in \a mr.
 *
 */
int _MMG5_intmetsavedir(MMG5_pMesh mesh, double *m,double *n,double *mr) {
  int    i;
  double lambda[2],vp[2][2],siz,isqhmin;

  isqhmin = 1.0 / (mesh->info.hmin * mesh->info.hmin);
  _MMG5_eigensym(m,lambda,vp);

  for (i=0; i<2; i++) {
    siz = n[0]*vp[i][0]*vp[i][0] + 2.0*n[1]*vp[i][0]*vp[i][1]
      + n[2]*vp[i][1]*vp[i][1];
    lambda[i] = MG_MAX(lambda[i],siz);
    lambda[i] = MG_MIN(lambda[i],isqhmin);
  }
  mr[0] = lambda[0]*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
  mr[1] = lambda[0]*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
  mr[2] = lambda[0]*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param np0 index of edge's extremity.
 * \param ux distance \f$[p0;p1]\f$ along x axis.
 * \param uy distance \f$[p0;p1]\f$ along y axis.
 * \param uz distance \f$[p0;p1]\f$ along z axis.
 * \param mr computed metric tensor.
 * \return 1.
 *
 * Build metric tensor at ridge point p0, when computations with respect to p1
 * are to be held.
 *
 */
int _MMG5_buildridmet(MMG5_pMesh mesh,MMG5_pSol met,int np0,
                      double ux,double uy,double uz,double mr[6]) {
  MMG5_pPoint p0;
  MMG5_pxPoint  go;
  double ps1,ps2,*n1,*n2,*t,*m,dv,u[3],r[3][3];

  p0 = &mesh->point[np0];
  if ( !(MG_GEO & p0->tag) )  return(0);
  m = &met->m[6*np0];
  go = &mesh->xpoint[p0->xp];
  t = &p0->n[0];

  /* Decide between the two possible configurations */
  n1 = &go->n1[0];
  n2 = &go->n2[0];

  ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
  ps2 = ux*n2[0] + uy*n2[1] + uz*n2[2];

  if ( fabs(ps2)<fabs(ps1) ) {
    n1 = &go->n2[0];
    dv = m[2];
  }
  else{
    dv = m[1];
  }

  u[0] = n1[1]*t[2] - n1[2]*t[1];
  u[1] = n1[2]*t[0] - n1[0]*t[2];
  u[2] = n1[0]*t[1] - n1[1]*t[0];

  /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(m[0],dv,0)*/
  r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n1[0];
  r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n1[1];
  r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n1[2];

  mr[0] = m[0]*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1];
  mr[1] = m[0]*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1];
  mr[2] = m[0]*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1];
  mr[3] = m[0]*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1];
  mr[4] = m[0]*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1];
  mr[5] = m[0]*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1];
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param np0 index of edge's extremity.
 * \param nt.
 * \param mr.
 *
 * \return 1.
 *
 * Build metric tensor at ridge point \a p0, when the 'good' normal direction is given by \a nt
 *
 */
/*  */
int _MMG5_buildridmetnor(MMG5_pMesh mesh,MMG5_pSol met,int np0,double nt[3],double mr[6]) {
  MMG5_pPoint p0;
  MMG5_pxPoint  go;
  double ps1,ps2,*n1,*n2,*t,*m,dv,u[3],r[3][3];

  p0 = &mesh->point[np0];
  if ( !(MG_GEO & p0->tag) )  return(0);
  m = &met->m[6*np0];
  t = &p0->n[0];
  go = &mesh->xpoint[p0->xp];

  /* Decide between the two possible configurations */
  n1 = &go->n1[0];
  n2 = &go->n2[0];

  ps1 = nt[0]*n1[0] + nt[1]*n1[1] + nt[2]*n1[2];
  ps2 = nt[0]*n2[0] + nt[1]*n2[1] + nt[2]*n2[2];

  if ( fabs(ps2) > fabs(ps1) ) {
    n1 = &go->n2[0];
    dv = m[2];
  }
  else{
    dv = m[1];
  }

  u[0] = n1[1]*t[2] - n1[2]*t[1];
  u[1] = n1[2]*t[0] - n1[0]*t[2];
  u[2] = n1[0]*t[1] - n1[1]*t[0];

  /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(m[0],dv,0)*/
  r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n1[0];
  r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n1[1];
  r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n1[2];

  mr[0] = m[0]*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1];
  mr[1] = m[0]*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1];
  mr[2] = m[0]*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1];
  mr[3] = m[0]*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1];
  mr[4] = m[0]*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1];
  mr[5] = m[0]*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1];

  return(1);
}

/**
 * \param mesh pointer toward the mesh stucture.
 * \param ip1 first point of face.
 * \param ip2 second point of face.
 * \param ip3 third point of face.
 * \param n pointer to store the computed normal.
 * \return 1
 *
 * Compute face normal given three points on the surface.
 *
 */
inline int _MMG5_norpts(MMG5_pMesh mesh,int ip1,int ip2, int ip3,double *n) {
  MMG5_pPoint   p1,p2,p3;
  double   dd,abx,aby,abz,acx,acy,acz,det;

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  p3 = &mesh->point[ip3];

  /* area */
  abx = p2->c[0] - p1->c[0];
  aby = p2->c[1] - p1->c[1];
  abz = p2->c[2] - p1->c[2];

  acx = p3->c[0] - p1->c[0];
  acy = p3->c[1] - p1->c[1];
  acz = p3->c[2] - p1->c[2];

  n[0] = aby*acz - abz*acy;
  n[1] = abz*acx - abx*acz;
  n[2] = abx*acy - aby*acx;
  det  = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];

  if ( det < _MMG5_EPSD2 )  return(0);

  dd = 1.0 / sqrt(det);
  n[0] *= dd;
  n[1] *= dd;
  n[2] *= dd;

  return(1);
}

/**
 * \param mesh pointer toward the mesh stucture.
 * \param pt pointer toward the triangle structure.
 * \param n pointer to store the computed normal.
 * \return 1
 *
 * Compute triangle normal.
 *
 */
inline int _MMG5_nortri(MMG5_pMesh mesh,MMG5_pTria pt,double *n) {

  return(_MMG5_norpts(mesh,pt->v[0],pt->v[1],pt->v[2],n));

}

/* Compute product R*M*tR when M is symmetric */
inline int _MMG5_rmtr(double r[3][3],double m[6], double mr[6]){
  double n[3][3];

  n[0][0] = m[0]*r[0][0] + m[1]*r[0][1] + m[2]*r[0][2];
  n[1][0] = m[1]*r[0][0] + m[3]*r[0][1] + m[4]*r[0][2];
  n[2][0] = m[2]*r[0][0] + m[4]*r[0][1] + m[5]*r[0][2];

  n[0][1] = m[0]*r[1][0] + m[1]*r[1][1] + m[2]*r[1][2];
  n[1][1] = m[1]*r[1][0] + m[3]*r[1][1] + m[4]*r[1][2];
  n[2][1] = m[2]*r[1][0] + m[4]*r[1][1] + m[5]*r[1][2];

  n[0][2] = m[0]*r[2][0] + m[1]*r[2][1] + m[2]*r[2][2];
  n[1][2] = m[1]*r[2][0] + m[3]*r[2][1] + m[4]*r[2][2];
  n[2][2] = m[2]*r[2][0] + m[4]*r[2][1] + m[5]*r[2][2];

  mr[0] = r[0][0]*n[0][0] + r[0][1]*n[1][0] + r[0][2]*n[2][0];
  mr[1] = r[0][0]*n[0][1] + r[0][1]*n[1][1] + r[0][2]*n[2][1];
  mr[2] = r[0][0]*n[0][2] + r[0][1]*n[1][2] + r[0][2]*n[2][2];
  mr[3] = r[1][0]*n[0][1] + r[1][1]*n[1][1] + r[1][2]*n[2][1];
  mr[4] = r[1][0]*n[0][2] + r[1][1]*n[1][2] + r[1][2]*n[2][2];
  mr[5] = r[2][0]*n[0][2] + r[2][1]*n[1][2] + r[2][2]*n[2][2];

  return(1);
}

/**
 * \param n pointer toward the vector that we want to send on the third vector
 * of canonical basis.
 * \param r computed rotation matrix.
 *
 * Compute rotation matrix that sends vector \a n to the third vector of
 * canonical basis.
 *
 */
inline int _MMG5_rotmatrix(double n[3],double r[3][3]) {
  double aa,bb,ab,ll,l,cosalpha,sinalpha;

  aa = n[0]*n[0];
  bb = n[1]*n[1];
  ab = n[0]*n[1];
  ll = aa+bb;
  cosalpha = n[2];
  sinalpha = sqrt(1.0- MG_MIN(1.0,cosalpha*cosalpha));

  /* No rotation needed in this case */
  if ( ll < _MMG5_EPS ) {
    if ( n[2] > 0.0 ) {
      r[0][0] = 1.0 ; r[0][1] = 0.0 ; r[0][2] = 0.0;
      r[1][0] = 0.0 ; r[1][1] = 1.0 ; r[1][2] = 0.0;
      r[2][0] = 0.0 ; r[2][1] = 0.0 ; r[2][2] = 1.0;
    }
    else {
      r[0][0] = -1.0 ; r[0][1] = 0.0 ; r[0][2] = 0.0;
      r[1][0] = 0.0 ; r[1][1] = 1.0 ; r[1][2] = 0.0;
      r[2][0] = 0.0 ; r[2][1] = 0.0 ; r[2][2] = -1.0;
    }
  }
  else {
    l = sqrt(ll);

    r[0][0] = (aa*cosalpha + bb)/ll;
    r[0][1] = ab*(cosalpha-1)/ll;
    r[0][2] = -n[0]*sinalpha/l;
    r[1][0] = r[0][1];
    r[1][1] = (bb*cosalpha + aa)/ll;
    r[1][2] = -n[1]*sinalpha/l;
    r[2][0] = n[0]*sinalpha/l;
    r[2][1] = n[1]*sinalpha/l;
    r[2][2] = cosalpha;
  }
  return(1);
}

/**
 * \param m pointer toward a 3x3 matrix
 * \param mi pointer toward the computed 3x3 matrix.
 *
 * Invert \a m (3x3 matrix) and store the result on \a mi
 *
 */
int _MMG5_invmat(double *m,double *mi) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;

  /* check diagonal matrices */
  vmax = fabs(m[1]);
  maxx = fabs(m[2]);
  if( maxx > vmax ) vmax = maxx;
  maxx = fabs(m[4]);
  if( maxx > vmax ) vmax = maxx;
  if ( vmax < _MMG5_EPS ) {
    mi[0]  = 1./m[0];
    mi[3]  = 1./m[3];
    mi[5]  = 1./m[5];
    mi[1] = mi[2] = mi[4] = 0.0;
    return(1);
  }

  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<6; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return(0);
  /* compute sub-dets */
  aa  = m[3]*m[5] - m[4]*m[4];
  bb  = m[4]*m[2] - m[1]*m[5];
  cc  = m[1]*m[4] - m[2]*m[3];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < _MMG5_EPS3 )  return(0);
  det = 1.0 / det;

  mi[0] = aa*det;
  mi[1] = bb*det;
  mi[2] = cc*det;
  mi[3] = (m[0]*m[5] - m[2]*m[2])*det;
  mi[4] = (m[1]*m[2] - m[0]*m[4])*det;
  mi[5] = (m[0]*m[3] - m[1]*m[1])*det;

  return(1);
}
/**
 * \param a matrix to invert.
 * \param b last member.
 * \param r vector of unknowns.
 * \return 0 if fail, 1 otherwise.
 *
 * Solve \f$ 3\times 3\f$ symmetric system \f$ A . r = b \f$.
 *
 */
inline int _MMG5_sys33sym(double a[6], double b[3], double r[3]){
  double ia[6],as[6],det,m;
  int    i;

  /* Multiply matrix by a constant coefficient for stability purpose (because of the scaling) */
  m = fabs(a[0]);
  for(i=1;i<6;i++){
    if(fabs(a[i])<m){
      m = fabs(a[i]);
    }
  }

  if(m < _MMG5_EPSD)
    return(0);

  m = 1.0/m;

  for(i=0;i<6;i++){
    as[i] = a[i]*m;
  }

  det = as[0]*(as[3]*as[5]-as[4]*as[4]) - as[1]*(as[1]*as[5]-as[2]*as[4]) \
    + as[2]*(as[1]*as[4]-as[2]*as[3]);

  if(fabs(det) < _MMG5_EPSD)
    return(0);

  det = 1.0/det;

  ia[0] = (as[3]*as[5]-as[4]*as[4]);
  ia[1] = - (as[1]*as[5]-as[2]*as[4]);
  ia[2] = (as[1]*as[4]-as[2]*as[3]);
  ia[3] = (as[0]*as[5]-as[2]*as[2]);
  ia[4] = -(as[0]*as[4]-as[2]*as[1]);
  ia[5] = (as[0]*as[3]-as[1]*as[1]);

  r[0] = ia[0]*b[0] + ia[1]*b[1] + ia[2]*b[2];
  r[1] = ia[1]*b[0] + ia[3]*b[1] + ia[4]*b[2];
  r[2] = ia[2]*b[0] + ia[4]*b[1] + ia[5]*b[2];

  r[0]*=(det*m);
  r[1]*=(det*m);
  r[2]*=(det*m);

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename pointer toward the file name.
 *
 * Debug function (not use in clean code): write mesh->tria structure in file.
 *
 */
void _MMG5_printTria(MMG5_pMesh mesh,char* fileName) {
  MMG5_pTria ptt;
  int   k;
  FILE  *inm;

  inm = fopen(fileName,"w");

  fprintf(inm,"----------> %d TRIANGLES <----------\n",mesh->nt);
  for(k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    fprintf(inm,"num %d -> %d %d %d\n",k,ptt->v[0],ptt->v[1],
            ptt->v[2]);
    fprintf(inm,"ref   -> %d\n",ptt->ref);
    fprintf(inm,"tag   -> %d %d %d\n",ptt->tag[0],ptt->tag[1],ptt->tag[2]);
    fprintf(inm,"edg   -> %d %d %d\n",ptt->edg[0],ptt->edg[1],ptt->edg[2]);
    fprintf(inm,"\n");
  }
  fprintf(inm,"---------> END TRIANGLES <--------\n");
  fclose(inm);
}

/**
 * \param c0 table of the coordinates of the starting point.
 * \param n0 normal at the starting point.
 * \param m metric to be transported.
 * \param c1 table of the coordinates of the ending point.
 * \param n1 normal at the ending point.
 * \param mt computed metric.
 * \return 0 if fail, 1 otherwise.
 *
 * Parallel transport of a metric tensor field, attached to point \a c0, with
 * normal \a n0, to point \a c1, with normal \a n1.
 *
 */

int _MMG5_paratmet(double c0[3],double n0[3],double m[6],double c1[3],double n1[3],double mt[6]) {
  double  r[3][3],mrot[6],mtan[3],lambda[2],vp[2][2],u[3],ps,ll;

  /* Take the induced metric tensor in the tangent plane by change of basis : R * M * {^t}R*/
  if ( !_MMG5_rotmatrix(n0,r) )  return(0);
  _MMG5_rmtr(r,m,mrot);
  mtan[0] = mrot[0];
  mtan[1] = mrot[1];
  mtan[2] = mrot[3];

  /* Take eigenvectors of metric tensor in tangent plane */
  _MMG5_eigensym(mtan,lambda,vp);

  /* Eigenvector in canonical basis = {t}R*vp[0] */
  u[0] = r[0][0]*vp[0][0] + r[1][0]*vp[0][1];
  u[1] = r[0][1]*vp[0][0] + r[1][1]*vp[0][1];
  u[2] = r[0][2]*vp[0][0] + r[1][2]*vp[0][1];

  /* Projection in the tangent plane of c1 */
  ps = u[0]*n1[0] + u[1]*n1[1] + u[2]*n1[2];
  u[0] -= ps*n1[0];
  u[1] -= ps*n1[1];
  u[2] -= ps*n1[2];
  ll = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
  if ( ll < _MMG5_EPSD )  return(0);
  ll = 1.0 / sqrt(ll);
  u[0] *= ll;
  u[1] *= ll;
  u[2] *= ll;

  /* And the transported metric is diag(lambda[0], lambda[1], 0) in basis (u,n1^u,n1) */
  r[0][0] = u[0];
  r[1][0] = u[1];
  r[2][0] = u[2];

  r[0][1] = n1[1]*u[2] - n1[2]*u[1];
  r[1][1] = n1[2]*u[0] - n1[0]*u[2];
  r[2][1] = n1[0]*u[1] - n1[1]*u[0];

  ll = r[0][1]*r[0][1] + r[1][1]*r[1][1] + r[2][1]*r[2][1];
  if ( ll < _MMG5_EPSD )  return(0);
  ll = 1.0 / sqrt(ll);
  r[0][1] *= ll;
  r[1][1] *= ll;
  r[2][1] *= ll;

  r[0][2] = n1[0];
  r[1][2] = n1[1];
  r[2][2] = n1[2];

  /*mt = R * diag(lambda[0], lambda[1], 0)*{^t}R */
  mt[0] = lambda[0]*r[0][0]*r[0][0] + lambda[1]*r[0][1]*r[0][1];
  mt[1] = lambda[0]*r[0][0]*r[1][0] + lambda[1]*r[0][1]*r[1][1];
  mt[2] = lambda[0]*r[0][0]*r[2][0] + lambda[1]*r[0][1]*r[2][1];
  mt[3] = lambda[0]*r[1][0]*r[1][0] + lambda[1]*r[1][1]*r[1][1];
  mt[4] = lambda[0]*r[2][0]*r[1][0] + lambda[1]*r[2][1]*r[1][1];
  mt[5] = lambda[0]*r[2][0]*r[2][0] + lambda[1]*r[2][1]*r[2][1];

  return(1);
}

/**
 * \return the available memory size of the computer.
 *
 * Compute the available memory size of the computer.
 *
 */
long long _MMG5_memSize (void) {
  long long mem;

#if (defined(__APPLE__) && defined(__MACH__))
  size_t size;

  size = sizeof(mem);
  if ( sysctlbyname("hw.memsize",&mem,&size,NULL,0) == -1)
    return(0);

#elif defined(__unix__) || defined(__unix) || defined(unix)
  mem = ((long long)sysconf(_SC_PHYS_PAGES))*
    ((long long)sysconf(_SC_PAGE_SIZE));
#else
  printf("  ## WARNING: UNKNOWN SYSTEM, RECOVER OF MAXIMAL MEMORY NOT AVAILABLE.\n");
  return(0);
#endif

  return(mem);
}
