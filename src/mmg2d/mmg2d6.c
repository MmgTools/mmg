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
 * \file mmg2d/mmg2d6.c
 * \brief Isosurface discretization.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "mmg2d.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param ip0 First vertex of the triangle
 * \param ip1 Second vertex of the triangle
 * \param ip2 Third vertex of the triangle
 * \return area of the triangle
 *
 * Calculate the area of a triangle given by its vertices
 *
 **/
inline double MMG2D_voltri(MMG5_pMesh mesh,int ip0,int ip1,int ip2) {
  MMG5_pPoint    p0,p1,p2;
  double         vol;

  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];

  vol = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);
  vol = 0.5*fabs(vol);

  return vol;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward the ls function
 * \param k index of the triangle
 * \return volfrac
 *
 * Calculate the area of the positive (if pm == 1) or negative (if pm == -1) subdomain
 * inside triangle k defined by the ls function in sol
 *
 **/
double MMG2D_vfrac(MMG5_pMesh mesh,MMG5_pSol sol,int k,int pm) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt[3];
  double        v[3],vfp,vfm,lam,area,eps,o1[2],o2[2];
  int           ip[3],nplus,nminus,nzero;
  int8_t        i,i0,i1,i2,imin1,iplus1,iz;

  eps = MMG5_EPS*MMG5_EPS;
  pt = &mesh->tria[k];

  ip[0] = pt->v[0];
  ip[1] = pt->v[1];
  ip[2] = pt->v[2];

  ppt[0] = &mesh->point[ip[0]];
  ppt[1] = &mesh->point[ip[1]];
  ppt[2] = &mesh->point[ip[2]];

  v[0] = sol->m[ip[0]];
  v[1] = sol->m[ip[1]];
  v[2] = sol->m[ip[2]];

  /* Identify number of zero, positive and negative vertices, and corresponding indices */
  nplus = nminus = nzero = 0;
  imin1 = iplus1 = iz = -1;

  for (i=0; i<3; i++) {
    if ( fabs(v[i]) < eps ) {
      nzero++;
      if ( iz < 0 ) iz = i;
    }
    else if ( v[i] >= eps ) {
      nplus++;
      if ( iplus1 < 0 ) iplus1 = i;
    }
    else {
      nminus++;
      if ( imin1 < 0 ) imin1 = i;
    }
  }

  /* Degenerate case */
  if ( nzero == 3 ) return 0.0;

  /* Whole triangle is positive */
  if ( nminus == 0 ) {
    vfp = (ppt[1]->c[0]-ppt[0]->c[0])*(ppt[2]->c[1]-ppt[0]->c[1]) - (ppt[1]->c[1]-ppt[0]->c[1])*(ppt[2]->c[0]-ppt[0]->c[0]);
    vfp = 0.5*fabs(vfp);
    if ( pm == 1 ) return vfp;
    else           return 0.0;
  }

  /* Whole triangle is negative */
  if ( nplus == 0 ) {
    vfm = (ppt[1]->c[0]-ppt[0]->c[0])*(ppt[2]->c[1]-ppt[0]->c[1]) - (ppt[1]->c[1]-ppt[0]->c[1])*(ppt[2]->c[0]-ppt[0]->c[0]);
    vfm = 0.5*fabs(vfm);
    if ( pm == -1 ) return vfm;
    else            return 0.0;
  }

  /* Exactly one vertex is negative */
  if ( nminus == 1 ) {
    i0 = imin1;
    i1 = MMG5_inxt2[i0];
    i2 = MMG5_iprv2[i0];

    lam = v[i0] / (v[i0]-v[i1]);
    o1[0] = ppt[i0]->c[0] + lam*(ppt[i1]->c[0]-ppt[i0]->c[0]);
    o1[1] = ppt[i0]->c[1] + lam*(ppt[i1]->c[1]-ppt[i0]->c[1]);

    lam = v[i0] / (v[i0]-v[i2]);
    o2[0] = ppt[i0]->c[0] + lam*(ppt[i2]->c[0]-ppt[i0]->c[0]);
    o2[1] = ppt[i0]->c[1] + lam*(ppt[i2]->c[1]-ppt[i0]->c[1]);

    vfm = (o1[0]-ppt[i0]->c[0])*(o2[1]-ppt[i0]->c[1]) - (o1[1]-ppt[i0]->c[1])*(o2[0]-ppt[i0]->c[0]);
    vfm = 0.5*fabs(vfm);

    if ( pm == -1 ) return vfm;
    else {
      area = (ppt[1]->c[0]-ppt[0]->c[0])*(ppt[2]->c[1]-ppt[0]->c[1]) - (ppt[1]->c[1]-ppt[0]->c[1])*(ppt[2]->c[0]-ppt[0]->c[0]);
      area = 0.5*fabs(area);
      vfp = area-vfm;
      return vfp;
    }
  }

  /* Exactly one vertex is positive */
  if ( nplus == 1 ) {
    i0 = iplus1;
    i1 = MMG5_inxt2[i0];
    i2 = MMG5_iprv2[i0];

    lam = v[i0] / (v[i0]-v[i1]);
    o1[0] = ppt[i0]->c[0] + lam*(ppt[i1]->c[0]-ppt[i0]->c[0]);
    o1[1] = ppt[i0]->c[1] + lam*(ppt[i1]->c[1]-ppt[i0]->c[1]);

    lam = v[i0] / (v[i0]-v[i2]);
    o2[0] = ppt[i0]->c[0] + lam*(ppt[i2]->c[0]-ppt[i0]->c[0]);
    o2[1] = ppt[i0]->c[1] + lam*(ppt[i2]->c[1]-ppt[i0]->c[1]);

    vfp = (o1[0]-ppt[i0]->c[0])*(o2[1]-ppt[i0]->c[1]) - (o1[1]-ppt[i0]->c[1])*(o2[0]-ppt[i0]->c[0]);
    vfp = 0.5*fabs(vfp);

    if ( pm == 1 ) return vfp;
    else {
      area = (ppt[1]->c[0]-ppt[0]->c[0])*(ppt[2]->c[1]-ppt[0]->c[1]) - (ppt[1]->c[1]-ppt[0]->c[1])*(ppt[2]->c[0]-ppt[0]->c[0]);
      area = 0.5*fabs(area);
      vfm = area-vfp;
      return vfm;
    }
  }

  /* Should not pass here */
  return(0.0);
}

/**
 * \param mesh pointer toward the mesh
 *
 * Reset mesh->info.isoref vertex and edge references to 0.
 *
 */
int MMG2D_resetRef(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0;
  int             k,ref;
  int8_t          i;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;

    for (i=0; i<3; i++) {
      p0 = &mesh->point[pt->v[i]];
      if ( pt->edg[i] == mesh->info.isoref ) pt->edg[i] = 0;
      if ( p0->ref == mesh->info.isoref ) p0->ref = 0;
    }
  }

  /* Reset the triangle references to their initial distribution */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;
    if( !MMG5_getStartRef(mesh,pt->ref,&ref) ) return 0;
    pt->ref = ref;
  }

  return 1;
}

/* Check whether snapping the value of vertex i of k to 0 exactly leads to a non manifold situation
 assumption: the triangle k has vertex i with value 0 and the other two with changing values */
int MMG2D_ismaniball(MMG5_pMesh mesh, MMG5_pSol sol, int start, int8_t istart) {
  MMG5_pTria       pt;
  double           v1, v2;
  int              *adja,k,ip1,ip2,end1,refstart;
  int8_t           i,i1,smsgn;
  static int8_t    mmgWarn=0;

  k = start;
  refstart = mesh->tria[k].ref;
  i = MMG5_inxt2[istart];

  /* First loop: stop if an external boundary, or a change in signs (or a 0) is met
     recall that MG_SMGSGN(a,b) = 1 provided a*b >0 */
  do{
    adja = &mesh->adja[3*(k-1)+1];
    k = adja[i] / 3;
    i1 = adja[i] % 3;
    i = MMG5_iprv2[i1];

    if ( k==0 ) break;

    pt = &mesh->tria[k];

    ip1 = pt->v[i1];
    ip2 = pt->v[i];

    v1 = sol->m[ip1];
    v2 = sol->m[ip2];

    if ( (fabs(v1) < MMG5_EPS) && (fabs(v2) < MMG5_EPS) ) {
      /* Do not authorize a snap that leads to a triangle with only 0 vertices */
      return 0;
    }

    /* Authorize change of references only provided the boundary reference is mesh->info.isoref */
    if ( pt->ref != refstart && pt->edg[i1] != mesh->info.isoref ) {
      smsgn = 0;
      k = 0;
    } else
      smsgn = (fabs(v1) < MMG5_EPS) || ( (fabs(v2) > MMG5_EPS) && MG_SMSGN(v1,v2) ) ? 1 : 0;
    // smsgn =  MG_SMSGN(v1,v2) ? 1 : 0;
  }
  while ( smsgn && (k != start) );

  if ( k==start ) {
    /* Complete ball has been travelled without crossing a boundary or finding a
     * sign change: we are in the special case where v1 = v2 = v[istart] = 0 in
     * tria start. In this case, test MG_SMSGN(v1,v2) returns 0 while smsgn is
     * computed to 1, which is non consistent.  */
    assert ( smsgn );
    return 0;
  }

  end1 = k;
  k = start;
  i = MMG5_iprv2[istart];

  /* Second loop: same travel in the opposite sense */
  do{
    adja = &mesh->adja[3*(k-1)+1];
    k = adja[i] / 3;
    i1 = adja[i] % 3;
    i = MMG5_inxt2[i1];

    if ( k==0 ) break;

    pt = &mesh->tria[k];
    ip1 = pt->v[i1];
    ip2 = pt->v[i];

    v1 = sol->m[ip1];
    v2 = sol->m[ip2];

    if ( (fabs(v1) < MMG5_EPS) && (fabs(v2) < MMG5_EPS) ) {
      /* Do not authorize a snap that leads to a triangle with only 0 vertices */
      return 0;
    }

    if ( pt->ref != refstart && pt->edg[i1] != mesh->info.isoref ) {
      smsgn = 0;
      k = 0;
    } else
      smsgn = (fabs(v2) < MMG5_EPS) || ( (fabs(v1) > MMG5_EPS) && MG_SMSGN(v1,v2) ) ? 1 : 0;
    // smsgn = MG_SMSGN(v1,v2) ? 1 : 0;
  }
  while ( smsgn && (k != start) );

  assert ( k!=start );

  /* If first stop was due to an external boundary, the second one must too;
     else, the final triangle for the first travel must be that of the second one */
  if ( k != end1 ) {
    if ( !mmgWarn ) {
      mmgWarn = 1;
      fprintf(stderr,"\n  ## Warning: %s: unsnap at least 1 point "
              "(point %d in tri %d).\n",__func__,MMG2D_indElt(mesh,start),
              MMG2D_indPt(mesh,mesh->tria[start].v[istart]));
    }
    return 0;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set
 *
 * \return 1 if success, 0 if fail
 *
 * Snap values of sol very close to 0 to 0 exactly (to avoid very small
 * triangles in cutting)
 */
int MMG2D_snapval(MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_pTria       pt,pt1;
  MMG5_pPoint      p0;
  double           v1,v2,*tmp;
  int              k,kk,iel,ns,nc,ip,ip1,ip2,npl,nmn,ilist,list[MMG2D_LONMAX+2];
  int8_t           i,j,j1,j2;

  /* Allocate memory for tmp */
  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(double),"temporary table",
                printf("  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(tmp,mesh->npmax+1,double,return 0);

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Snap values of sol that are close to 0 to 0 exactly */
  ns = nc = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) ) continue;
    if ( fabs(sol->m[k]) < MMG5_EPS ) {
      tmp[k] =  sol->m[k];
      p0->flag = 1;
      sol->m[k] = 0.0;
      ns++;
    }
  }

  /* Check that the snapping process has not led to a nonmanifold situation */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      ip1 = pt->v[MMG5_inxt2[i]];
      ip2 = pt->v[MMG5_iprv2[i]];

      p0 = &mesh->point[ip];
      v1 = sol->m[ip1];
      v2 = sol->m[ip2];

      /* Catch a snapped point by a triangle where there is a sign change */
      if ( p0->flag && !(MG_SMSGN(v1,v2)) ) {
        if ( !MMG2D_ismaniball(mesh,sol,k,i) ) {
          sol->m[ip] = tmp[ip];
          nc++;
        }
        p0->flag = 0;
      }
    }
  }

  /* Check that the ls function does not show isolated spots with 0 values (without sign changes) */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      if ( fabs(sol->m[ip]) >= MMG5_EPS ) continue;
      npl = nmn = 0;
      ilist = MMG2D_boulet(mesh,k,i,list);
      for(kk=0; kk<ilist; kk++) {
        iel = list[kk] / 3;
        j = list[kk] % 3;
        j1 = MMG5_inxt2[j];
        j2 = MMG5_iprv2[i];
        pt1 = &mesh->tria[iel];
        ip1 = pt1->v[j1];
        ip2 = pt1->v[j2];
        if ( sol->m[ip1] >= MMG5_EPS ) npl = 1;
        else if ( sol->m[ip1] <= -MMG5_EPS ) nmn = 1;

        if ( sol->m[ip2] >= MMG5_EPS ) npl = 1;
        else if ( sol->m[ip2] <= -MMG5_EPS ) nmn = 1;
      }

      if ( npl == 1 && nmn == 0 )
        sol->m[ip] = 100.0*MMG5_EPS;
      else if ( npl == 0 && nmn == 1 )
        sol->m[ip] = -100.0*MMG5_EPS;
    }
  }

  MMG5_DEL_MEM ( mesh, tmp );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && ns+nc > 0 )
    fprintf(stdout,"     %8d points snapped, %d corrected\n",ns,nc);

  return 1;
}

/* Check whether the ball of vertex i in tria start is manifold;
 by assumption, i inxt[i] is one edge of the implicit boundary */
int MMG2D_chkmaniball(MMG5_pMesh mesh, int start, int8_t istart) {
  MMG5_pTria         pt;
  int                *adja,k,refstart;
  int8_t             i,i1;

  pt = &mesh->tria[start];
  k = start;
  i = istart;
  refstart = pt->ref;

  /* First travel, while another part of the implicit boundary is not met */
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];

    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_inxt2[i];
  }
  while ( k && ( mesh->tria[k].ref == refstart ) );

  /* Case where a boundary is hit: travel in the other sense from start, and make sure
   that a boundary is hit too */
  if ( k == 0 ) {
    k = start;
    i = istart;

    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_iprv2[i];
    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_iprv2[i];

    /* Tested point is connected to two external edges */
    if ( k == 0 ) return 1;

    do {
      adja = &mesh->adja[3*(k-1)+1];
      i1 = MMG5_iprv2[i];

      k = adja[i1] / 3;
      i = adja[i1] % 3;
      i = MMG5_iprv2[i];
    }
    while ( k && ( mesh->tria[k].ref != refstart ) );

    if ( k == 0 ) return 1;
    else          return 0;

  }

  /* General case: go on travelling until another implicit boundary is met */
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];

    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_inxt2[i];
  }
  while ( k && ( mesh->tria[k].ref != refstart ) );

  /* At least 3 boundary segments meeting at p */
  if ( k != start )
    return 0;

  return 1;
}

/* Check whether the resulting two subdomains coming from isovalue
 * discretization are manifold */
int MMG2D_chkmanimesh(MMG5_pMesh mesh) {
  MMG5_pTria      pt,pt1;
  int             *adja,k,cnt,iel;
  int8_t          i,i1;
  static int8_t   mmgWarn=0;

  /* First check: check whether one triangle in the mesh has 3 boundary faces */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];
    cnt = 0;
    for (i=0; i<3; i++) {
      iel = adja[i] / 3;

      if (!iel ) {
        cnt++;
        continue;
      }
      else {
        pt1 = &mesh->tria[iel];
        if ( pt1->ref != pt->ref ) cnt++;
      }
    }
    if( cnt == 3 ) {
      if ( !mmgWarn ) {
        mmgWarn = 1;
        fprintf(stderr,"\n  ## Warning: %s: at least 1 triangle with 3 boundary"
                " edges.\n",__func__);
      }
    }
  }

  /* Second check: check whether the configuration is manifold in the ball of
     each point; each vertex on the implicit boundary is caught in such a way
     that i1 inxt[i1] is one edge of the implicit boundary */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      adja = &mesh->adja[3*(k-1)+1];
      iel = adja[i] / 3;

      if (! iel ) continue;
      pt1 = &mesh->tria[iel];
      if ( pt->ref == pt1->ref || pt->edg[i]!= mesh->info.isoref ) continue;

      i1 = MMG5_inxt2[i];
      if ( !MMG2D_chkmaniball(mesh,k,i1) ) {
        fprintf(stderr,"   *** Topological problem\n");
        fprintf(stderr,"       non manifold curve at point %d %d\n",pt->v[i1], MMG2D_indPt(mesh,pt->v[i1]));
        fprintf(stderr,"       non manifold curve at tria %d (ip %d)\n", MMG2D_indElt(mesh,k),i1);
        return 0;
      }
    }
  }

  if ( mesh->info.imprim > 0 || mesh->info.ddebug )
    fprintf(stdout,"  *** Manifold implicit surface.\n");
  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param sol pointer toward the level-set
 *
 * \return 1 if success, 0 otherwise
 *
 * Removal of small parasitic components (bubbles of material, etc) with volume less than
 * mesh->info.rmc * volume of the mesh.
 *
 */
int MMG2D_rmc(MMG5_pMesh mesh, MMG5_pSol sol){
  MMG5_pTria     pt,pt1,pt2;
  double         volc,voltot,v0,v1,v2;
  int            k,kk,l,ll,ncp,ncm,ip0,ip1,ip2,base,cur,ipile,*pile,*adja;
  int8_t         i,i1,i2;

  ncp = 0;
  ncm = 0;

  /* Erase triangle flags */
  for (k=1; k<=mesh->nt; k++) mesh->tria[k].flag = 0;

  /* Calculate volume of the total mesh */
  voltot = 0.0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    ip0 = pt->v[0];
    ip1 = pt->v[1];
    ip2 = pt->v[2];
    voltot += MMG2D_voltri(mesh,ip0,ip1,ip2);
  }

  /* Memory allocation for pile */
  MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(int),"temporary table",
               printf("  Exit program.\n");
               return 0);
  MMG5_SAFE_CALLOC(pile,mesh->nt+1,int,return 0);

  /* Investigate only positive connected components */
  base = ++mesh->base;

  for (k=1; k<=mesh->nt; k++) {
    ipile = 0;
    volc = 0.0;
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    if ( pt->flag == base ) continue;

    /* Checks signs of the LS function at the 3 vertices of pt */
    ip0 = pt->v[0];
    ip1 = pt->v[1];
    ip2 = pt->v[2];

    v0 = sol->m[ip0];
    v1 = sol->m[ip1];
    v2 = sol->m[ip2];

    if ( v0 <= 0.0 && v1 <= 0.0 && v2 <= 0.0 ) continue;

    /* Add triangle to pile if one vertex is > 0 */
    pt->flag = base;
    pile[ipile] = k;
    ipile++;
    if ( ipile > mesh->nt ) {
      fprintf(stderr,"\n  ## Problem in length of pile; function rmc.\n"
              " Check that the level-set intersect the mesh.\n"
              " Exit program.\n");

      return 0;
    }

    /* Pile up all the positive connected component attached to the first triangle */
    cur = 0;
    do {
      kk = pile[cur];
      pt1 = &mesh->tria[kk];

      /* Add local volume fraction of the positive subdomain to volc */
      volc += MMG2D_vfrac(mesh,sol,kk,1);

      /* Add adjacent triangles to kk via positive vertices to the pile, if need be */
      adja = &mesh->adja[3*(kk-1)+1];
      for (i=0; i<3; i++) {
        ip0 = pt1->v[i];
        if ( sol->m[ip0] <= 0.0 ) continue;

        i1 = MMG5_inxt2[i];
        i2 = MMG5_inxt2[i1];

        /* First neighbor of positive vertex i */
        ll = adja[i1] / 3;
        if ( ll ) {
          pt2 = &mesh->tria[ll];
          if ( pt2->flag != base ) {
            pt2->flag = base;
            pile[ipile] = ll;
            ipile++;
            if ( ipile > mesh->nt ) {
              fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
              return 0;
            }
          }
        }

        /* Second neighbor of positive vertex i */
        ll = adja[i2] / 3;
        if ( ll ) {
          pt2 = &mesh->tria[ll];
          if ( pt2->flag != base ) {
            pt2->flag = base;
            pile[ipile] = ll;
            ipile++;
            if ( ipile > mesh->nt ) {
              fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
              return 0;
            }
          }
        }
      }
    }
    while ( ++cur < ipile );

    /* Remove connected component if its volume is too small */
    if ( volc < mesh->info.rmc*voltot ) {
      for (l=0; l<ipile; l++) {
        pt1 = &mesh->tria[pile[l]];
        for (i=0; i<3; i++) {
          ip0 = pt1->v[i];
          if ( sol->m[ip0] > 0.0 ) sol->m[ip0] = -100*MMG5_EPS;
        }
      }
      ncp++;
    }

  }

  /* Investigate only negative connected components */
  base = ++mesh->base;

  for (k=1; k<=mesh->nt; k++) {
    ipile = 0;
    volc = 0.0;
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    if ( pt->flag == base ) continue;

    /* Checks signs of the LS function at the 3 vertices of pt */
    ip0 = pt->v[0];
    ip1 = pt->v[1];
    ip2 = pt->v[2];

    v0 = sol->m[ip0];
    v1 = sol->m[ip1];
    v2 = sol->m[ip2];

    if ( v0 >= 0.0 && v1 >= 0.0 && v2 >= 0.0 ) continue;

    /* Pile up all the negative connected component attached to the first triangle */
    pt->flag = base;
    pile[ipile] = k;
    ipile++;
    if ( ipile > mesh->nt ) {
      fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
      return 0;
    }

    cur = 0;
    do {
      kk = pile[cur];
      pt1 = &mesh->tria[kk];

      /* Add local volume fraction of the negative subdomain to volc */
      volc += MMG2D_vfrac(mesh,sol,kk,-1);

      /* Add adjacent triangles to kk via negative vertices to the pile, if need be */
      adja = &mesh->adja[3*(kk-1)+1];
      for (i=0; i<3; i++) {
        ip0 = pt1->v[i];
        if ( sol->m[ip0] >= 0.0 ) continue;

        i1= MMG5_inxt2[i];
        i2 = MMG5_inxt2[i1];

        /* First neighbor of negative vertex i */
        ll = adja[i1] / 3;
        if ( ll ) {
          pt2 = &mesh->tria[ll];
          if ( pt2->flag != base ) {
            pt2->flag = base;
            pile[ipile] = ll;
            ipile++;
            if ( ipile > mesh->nt ) {
              fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
              return 0;
            }
          }
        }

        /* Second neighbor of negative vertex i */
        ll = adja[i2] / 3;
        if ( ll ) {
          pt2 = &mesh->tria[ll];
          if ( pt2->flag != base ) {
            pt2->flag = base;
            pile[ipile] = ll;
            ipile++;
            if ( ipile > mesh->nt ) {
              fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
              return 0;
            }
          }
        }

      }
    }
    while ( ++cur < ipile );

    /* Remove connected component if its volume is too small */
    if ( volc < mesh->info.rmc*voltot ) {
      for (l=0; l<ipile; l++) {
        pt1 = &mesh->tria[pile[l]];
        for (i=0; i<3; i++) {
          ip0 = pt1->v[i];
          if ( sol->m[ip0] < 0.0 ) sol->m[ip0] = 100*MMG5_EPS;
        }
      }
      ncm++;
    }
  }

  /* Erase triangle flags */
  for (k=1; k<=mesh->nt; k++) mesh->tria[k].flag = 0;

  /* Release memory */
  MMG5_DEL_MEM(mesh,pile);

  if ( mesh->info.imprim > 0 || mesh->info.ddebug ) {
    printf("\n  *** Removed %d positive parasitic bubbles and %d negative parasitic bubbles\n",ncp,ncm);
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh
 * \param sol pointer toward the level-set
 * \param met pointer toward a metric (non-mandatory)
 *
 * \return 1 if success, 0 otherwise
 *
 * Effective discretization of the 0 level set encoded in sol in the mesh
 *
 */
int MMG2D_cuttri_ls(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSol met){
  MMG5_pTria   pt;
  MMG5_pPoint  p0,p1;
  MMG5_Hash    hash;
  double       v0,v1,s,c[2];
  int          k,ip0,ip1,nb,np,nt,ns,refint,refext,vx[3];
  int8_t       i,i0,i1,ier;

  /* Reset flag field for points */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Evaluate the number of intersected edges by the 0 level set */
  nb = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      i0 = MMG5_inxt2[i];
      i1 = MMG5_inxt2[i0];

      ip0 = pt->v[i0];
      ip1 = pt->v[i1];

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];

      if ( p0->flag && p1->flag ) continue;

      v0 = sol->m[ip0];
      v1 = sol->m[ip1];

      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
        nb++;
        if ( !p0->flag ) p0->flag = nb;
        if ( !p1->flag ) p1->flag = nb;
      }
    }
  }
  if ( !nb ) return 1;

  /* Create the intersection points between the edges in the mesh and the 0
   * level set */
  if ( !MMG5_hashNew(mesh,&hash,nb,2*nb) ) return 0;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      i0 = MMG5_inxt2[i];
      i1 = MMG5_inxt2[i0];

      ip0 = pt->v[i0];
      ip1 = pt->v[i1];

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];

      np = MMG5_hashGet(&hash,ip0,ip1);
      if ( np ) continue;

      if ( !MMG5_isSplit(mesh,pt->ref,&refint,&refext) ) continue;

      v0 = sol->m[ip0];
      v1 = sol->m[ip1];

      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;

      /* Intersection point between edge p0p1 and the 0 level set */
      s = v0/(v0-v1);
      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);

      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);

      np = MMG2D_newPt(mesh,c,0);
      if ( !np ) {
       /* reallocation of point table */
        MMG2D_POINT_REALLOC(mesh,met,np,mesh->gap,
                            fprintf(stderr,"\n  ## Error: %s: unable to"
                                    " allocate a new point.\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            return 0;,
                            c,0);
      }
      sol->m[np] = 0.0;
      /* If there is a metric in the mesh, interpolate it at the new point */
      if ( met && met->m )
        MMG2D_intmet(mesh,met,k,i,np,s);

      MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
    }
  }

  /* Proceed to splitting by calling patterns */
  nt  = mesh->nt;
  ns  = 0;
  ier = 1;
  for (k=1; k<=nt; k++) {

    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    pt->flag = 0;

    for (i=0; i<3; i++) {
      i0 = MMG5_inxt2[i];
      i1 = MMG5_inxt2[i0];

      ip0 = pt->v[i0];
      ip1 = pt->v[i1];

      vx[i] = MMG5_hashGet(&hash,ip0,ip1);

      if ( vx[i] ) MG_SET(pt->flag,i);
    }

    switch( pt->flag ) {
      /* 1 edge split -> 0-+ */
      case 1: case 2: case 4:
        ier = MMG2D_split1(mesh,met,k,vx);
        ns++;
        break;

      /* 2 edge split -> +-- or -++ */
      case 3: case 5: case 6:
        ier = MMG2D_split2(mesh,met,k,vx);
        ns++;
        break;

      default:
        assert(pt->flag==0);
        break;
    }
    if ( !ier ) return 0;
  }

  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",ns);

  MMG5_DEL_MEM(mesh,hash.item);
  return ns;

}

/* Set references to the new triangles */
int MMG2D_setref_ls(MMG5_pMesh mesh, MMG5_pSol sol){
  MMG5_pTria    pt;
  double        v,v1;
  int           k,ip,ip1,ier,ref,refint,refext;
  int8_t        i,i1,i2,nmn,npl,nz;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    ref = pt->ref;
    nmn = npl = nz = 0;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      v = sol->m[ip];

      if ( v > 0.0 )
        npl++;
      else if ( v < 0.0 )
        nmn++;
      else
        nz++;
    }

    assert(nz < 3);
    ier = MMG5_isSplit(mesh,ref,&refint,&refext);

    if ( npl ) {
      if ( ier ) {
        assert ( !nmn );
        pt->ref = refext;
      }
    }
    else {
      if ( ier ) {
        assert ( !npl );
        pt->ref = refint;
      }
    }

    /* Set mesh->info.isoref ref at ls edges and at the points of these edges */
    if ( nz == 2 ) {
      for (i=0; i<3; i++) {
        ip  = pt->v[MMG5_inxt2[i]];
        ip1 = pt->v[MMG5_iprv2[i]];
        v   = sol->m[ip];
        v1  = sol->m[ip1];
        if ( v == 0.0 && v1 == 0.0) {
          pt->edg[i]  = mesh->info.isoref;
          pt->tag[i] |= MG_REF;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].ref = mesh->info.isoref;
          mesh->point[pt->v[i2]].ref = mesh->info.isoref;
        }
      }
    }

  }

  return 1;
}

/* Main function of the -ls mode */
int MMG2D_mmg2d6(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_pSol met) {
  int k;

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION\n");

  if ( mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: Isosurface extraction not available with"
            " hybrid meshes. Exit program.\n");
    return 0;
  }

  /* Work only with the 0 level set */
  for (k=1; k<= sol->np; k++)
    sol->m[k] -= mesh->info.ls;

  /* Transfer the boundary edge references to the triangles */
  if ( !MMG2D_assignEdge(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }

  /* Snap values of the level set function which are very close to 0 to 0 exactly */
  if ( !MMG2D_snapval(mesh,sol) ) {
    fprintf(stderr,"\n  ## Wrong input implicit function. Exit program.\n");
    return 0;
  }

  /* Removal of small parasitic components */
  if ( mesh->info.rmc > 0. && !MMG2D_rmc(mesh,sol) ) {
    fprintf(stderr,"\n  ## Error in removing small parasitic components. Exit program.\n");
    return 0;
  }

  /* No need to keep adjacencies from now on */
  MMG5_DEL_MEM(mesh,mesh->adja);

  /* Reset the mesh->info.isoref field everywhere it appears */
  if ( !MMG2D_resetRef(mesh) ) {
    fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
    return 0;
  }

  /* Effective splitting of the crossed triangles */
  if ( !MMG2D_cuttri_ls(mesh,sol,met) ) {
    fprintf(stderr,"\n  ## Problem in cutting triangles. Exit program.\n");
    return 0;
  }

  /* Set references on the interior / exterior triangles*/
  if ( !MMG2D_setref_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }

  /* Creation of adjacency relations in the mesh */
  if ( !MMG2D_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* Check that the resulting mesh is manifold */
  if ( !MMG2D_chkmanimesh(mesh) ) {
    fprintf(stderr,"\n  ## No manifold resulting situation. Exit program.\n");
    return 0;
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,sol->m);
  sol->np = 0;

  MMG5_DEL_MEM( mesh,mesh->info.mat );

  return 1;
}
