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
 * \file mmgs/gentool.c
 * \brief Generic algebraic and algorithmic tools.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgs.h"

/* Compute product R*M*tR when M is symmetric */
inline int rmtr(double r[3][3],double m[6], double mr[6]){
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

/* Compute the intersected (2 x 2) metric between metrics m and n, PRESERVING the directions
   of m. Result is stored in mr*/
int intmetsavedir(MMG5_pMesh mesh, double *m,double *n,double *mr) {
  int    i;
  double lambda[2],vp[2][2],siz,isqhmin;

  isqhmin = 1.0 / (mesh->info.hmin * mesh->info.hmin);
  _MMG5_eigensym(m,lambda,vp);

  for (i=0; i<2; i++) {
    siz = n[0]*vp[i][0]*vp[i][0] + 2.0*n[1]*vp[i][0]*vp[i][1] + n[2]*vp[i][1]*vp[i][1];
    lambda[i] = MG_MAX(lambda[i],siz);
    lambda[i] = MG_MIN(lambda[i],isqhmin);
  }
  mr[0] = lambda[0]*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
  mr[1] = lambda[0]*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
  mr[2] = lambda[0]*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];

  return(1);
}

/* Delete all triangle references in mesh */
int delref(MMG5_pMesh mesh) {
  MMG5_pTria    pt;
  int      k;

  for(k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    pt->ref = 0;
  }

  return(1);
}

/* Start from triangle start, and pile up triangles by adjacency, till a GEO or REF curve is met ;
   pass all references of travelled faces to ref ; putreq = 1 if boundary edges met must
   be set to MG_REQ, 0 otherwise. */
int setref(MMG5_pMesh mesh,int start,int ref,int putreq) {
  MMG5_pTria      pt,pt1;
  int        *list,*adja,cur,base,k,iel,jel,ilist;
  char       j,voy;

  ilist = cur = 0;
  _MMG5_SAFE_CALLOC(list,mesh->nt+1,int);
  base = ++mesh->base;

  /* Pile up triangles from start, till a GEO boundary is met */
  pt = &mesh->tria[start];
  list[ilist] = start;
  ilist++;
  assert( ilist <= mesh->nt );
  pt->flag = base;

  do {
    iel = list[cur];
    pt = &mesh->tria[iel];
    adja = &mesh->adja[3*(iel-1)+1];

    for(j=0; j<3; j++) {
      if( MG_EDG(pt->tag[j]) ) {
        if( putreq ) {
          pt->tag[j] |= MG_REQ;
          jel = adja[j] / 3;
          voy = adja[j] % 3;
          if( !jel ) continue;
          pt1 = &mesh->tria[jel];
          pt1->tag[voy] |= MG_REQ;
        }
        continue;
      }
      jel = adja[j] / 3;
      assert(jel);
      pt1 = &mesh->tria[jel];
      if ( pt1->flag == base )  continue;

      list[ilist] = jel;
      ilist++;
      assert( ilist <= mesh->nt );
      pt1->flag = base;
    }
    cur++;
  }
  while( cur < ilist );

  /* Set all references of triangles of list to ref */
  for (k=0; k<ilist; k++) {
    iel = list[k];
    pt  = &mesh->tria[iel];
    pt->ref = ref;
    printf("Le tria %d passe a %d \n",k,ref);
  }

  return(1);
}

/* invert 3x3 non-symmetric matrix */
int invmatg(double m[9],double mi[9]) {
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

/** find the element number in packed numerotation */
int _MMG5_indElt(MMG5_pMesh mesh, int kel) {
  MMG5_pTria pt;
  int    ne, k;

  ne = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( MG_EOK(pt) ) {
      ne++;
      if ( k == kel )  return(ne);
    }
  }
  return(0);
}

/** find the point number in packed numerotation */
int _MMG5_indPt(MMG5_pMesh mesh, int kp) {
  MMG5_pPoint ppt;
  int    np, k;

  np = 0;
  for (k=1; k<=mesh->nt; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      np++;
      if ( k == kp )  return(np);
    }
  }
  return(0);
}
