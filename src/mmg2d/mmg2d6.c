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
#include "libmmg2d_private.h"
#include "mmg2dexterns.h"

/* Check whether the ball of vertex i in tria start is manifold;
 by assumption, i inxt[i] is one edge of the implicit boundary */
int MMG2D_chkmaniball(MMG5_pMesh mesh, MMG5_int start, int8_t istart) {
  MMG5_pTria         pt;
  MMG5_int           refstart,*adja,k;
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

/**
 * \param mesh pointer toward the mesh.
 * \return 1 if the mesh is manifold, 0 otherwise.
 *
 * Check whether the resulting two subdomains occupying mesh are manifold.
 *
 */
int MMG2D_chkmanimesh(MMG5_pMesh mesh) {
  MMG5_pTria      pt,pt1;
  MMG5_int        *adja,k;
  MMG5_int        cnt,iel;
  int8_t          i,i1;
  static int8_t   mmgWarn = 0;

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
        fprintf(stderr,"       non manifold curve at point %" MMG5_PRId " %" MMG5_PRId "\n",pt->v[i1], MMG2D_indPt(mesh,pt->v[i1]));
        fprintf(stderr,"       non manifold curve at tria %" MMG5_PRId " (ip %d)\n", MMG2D_indElt(mesh,k),i1);
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
  MMG5_int     k,ip0,ip1,nb,np,nt,ns,vx[3],refint,refext;
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
    fprintf(stdout,"     %7" MMG5_PRId " splitted\n",ns);

  MMG5_DEL_MEM(mesh,hash.item);
  return ns;

}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \return 1.
 *
 * Set references to tris according to the sign of the level set function.
 *
 */
int MMG2D_setref_ls(MMG5_pMesh mesh, MMG5_pSol sol){
  MMG5_pTria    pt;
  double        v,v1;
  int           ier;
  MMG5_int      k,ip,ip1,ref,refint,refext;
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

    if ( ier ) {
      if ( npl ) {
        assert ( !nmn );
        pt->ref = refext;
      }
      else {
        assert ( nmn );
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
  MMG5_int k;

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
  if ( !MMG5_snpval_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Wrong input implicit function. Exit program.\n");
    return 0;
  }

  /* Removal of small parasitic components */
  if ( mesh->info.rmc > 0. && !MMG5_rmc(mesh,sol) ) {
    fprintf(stderr,"\n  ## Error in removing small parasitic components. Exit program.\n");
    return 0;
  }

  /* No need to keep adjacencies from now on */
  MMG5_DEL_MEM(mesh,mesh->adja);

  /* Reset the mesh->info.isoref field everywhere it appears */
  if ( !MMG5_resetRef(mesh) ) {
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
