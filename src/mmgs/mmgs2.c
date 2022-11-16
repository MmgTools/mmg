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
 * \file mmgs/mmgs2.c
 * \brief Create implicit surface in mesh.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmgs_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param start index of starting tria.
 * \param istart local index of point that we check (in tria \a start)
 * \return 1 if the ball is manifold, 0 otherwise.
 *
 * Check whether the ball of vertex i in tria start is manifold;
 *
 * \warning i inxt[i] is one edge of the implicit boundary.
 *
 */
int MMGS_chkmaniball(MMG5_pMesh mesh, MMG5_int start, int8_t istart) {
  MMG5_int           *adja,k;
  int8_t             i,i1;

  k = start;
  i = istart;

  i1 = MMG5_iprv2[i];

#ifndef NDEBUG
  MMG5_pTria pt = &mesh->tria[start];
  assert( MG_EDG(pt->tag[i1]) && (pt->edg[i1]==mesh->info.isoref) );
#endif

  /* First travel, while another part of the implicit boundary is not met */
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];

    k = adja[i1] / 3;
    i = adja[i1] % 3;

    if ( !k || mesh->tria[k].edg[i]==mesh->info.isoref ) break;

    i = MMG5_inxt2[i];
  }
  while ( k!=start );

  assert(k!=start); //unexpected case

  /* Case where a boundary is hit: travel in the other sense from start, and make sure
   that a boundary is hit too */
  if ( k == 0 ) {
    k = start;
    i = istart;

    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_iprv2[i];
    k = adja[i1] / 3;

    /* Check of the way the point is caught (the left-hand edge is not an external edge) */
    assert ( k );

    i = adja[i1] % 3;
    i = MMG5_iprv2[i];

    do {
      adja = &mesh->adja[3*(k-1)+1];
      i1 = MMG5_iprv2[i];

      k = adja[i1] / 3;
      i = adja[i1] % 3;

      if ( (!k) || mesh->tria[k].edg[i]==mesh->info.isoref ) break;

      i = MMG5_iprv2[i];
    }
    while ( k!=start );

    assert(k!=start); //unexpected case

    return !k;
  }

  /* General case: go on travelling until another implicit boundary is met */
  i = MMG5_inxt2[i];
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];

    k = adja[i1] / 3;
    i = adja[i1] % 3;

    if ( (!k) || mesh->tria[k].edg[i]==mesh->info.isoref ) break;

    i = MMG5_inxt2[i];
  }
  while ( k!=start );

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
static
int MMGS_chkmanimesh(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
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
        if ( pt->edg[i] == mesh->info.isoref ) cnt++;
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
      if ( pt->edg[i] != mesh->info.isoref ) continue;

      i1 = MMG5_inxt2[i];
      if ( !MMGS_chkmaniball(mesh,k,i1) ) {
        fprintf(stderr,"   *** Topological problem\n");
        fprintf(stderr,"       non manifold curve at point %" MMG5_PRId " %" MMG5_PRId "\n",pt->v[i1], MMGS_indPt(mesh,pt->v[i1]));
        fprintf(stderr,"       non manifold curve at tria %" MMG5_PRId " (ip %d)\n", MMGS_indElt(mesh,k),i1);
        return 0;
      }
    }
  }

  if ( mesh->info.imprim > 0 || mesh->info.ddebug )
    fprintf(stdout,"  *** Manifold implicit surface.\n");

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \param met pointer toward a metric (non-mandatory).
 * \return 1 if success, 0 otherwise.
 *
 * Proceed to discretization of the implicit function carried by sol into mesh,
 * once values of sol have been snapped/checked
 *
 */
static int MMGS_cuttri_ls(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_pSol met){
  MMG5_pTria   pt;
  MMG5_pPoint  p0,p1;
  MMG5_Hash    hash;
  double       c[3],v0,v1,s;
  MMG5_int     vx[3],k,np;
  MMG5_int     ip0,ip1,ns,nt,ier,nb;
  int8_t       ia;
  /* reset point flags and h */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* compute the number nb of intersection points on edges */
  nb = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    for (ia=0; ia<3; ia++) {
      ip0 = pt->v[MMG5_inxt2[ia]];
      ip1 = pt->v[MMG5_iprv2[ia]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];
      if ( p0->flag && p1->flag )  continue;
      v0  = sol->m[ip0];
      v1  = sol->m[ip1];
      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
        if ( !p0->flag ) {
          p0->flag = ++nb;
        }
        if ( !p1->flag ) {
          p1->flag = ++nb;
        }
      }
    }
  }
  if ( ! nb )  return 1;

  /* Create intersection points at 0 isovalue and set flags to trias */
  if ( !MMG5_hashNew(mesh,&hash,nb,3*nb) ) return 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (ia=0; ia<3; ia++) {
      ip0 = pt->v[MMG5_inxt2[ia]];
      ip1 = pt->v[MMG5_iprv2[ia]];
      np  = MMG5_hashGet(&hash,ip0,ip1);
      if ( np )  continue;

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      v0 = sol->m[ip0];
      v1 = sol->m[ip1];
      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;

      s = v0 / (v0-v1);

      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);
      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
      c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

      np = MMGS_newPt(mesh,c,NULL);
      if ( !np ) {
        MMGS_POINT_REALLOC(mesh,sol,np,0.2,
                            fprintf(stderr,"\n  ## Error: %s: unable to"
                                    " allocate a new point\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            return 0
                            ,c,NULL);
      }
      sol->m[np] = 0.0;

      /* If user provide a metric, interpolate it at the new point */
      if ( met && met->m ) {
        if ( met->size > 1 ) {
          ier = MMGS_intmet33_ani(mesh,met,k,ia,np,s);
        }
        else {
          ier = intmet_iso(mesh,met,k,ia,np,s);
        }
        if ( ier <= 0 ) {
          // Unable to compute the metric
          fprintf(stderr,"\n  ## Error: %s: unable to"
                  " interpolate the metric during the level-set"
                  " discretization\n",__func__);
          return 0;
        }
      }

      MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
    }
  }

  /* Proceed to splitting, according to flags to tris */
  nt  = mesh->nt;
  ns  = 0;
  ier = 1;
  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    pt->flag = 0;
    memset(vx,0,3*sizeof(MMG5_int));
    for (ia=0; ia<3; ia++) {
      vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_inxt2[ia]],pt->v[MMG5_iprv2[ia]]);
      if ( vx[ia] ) {
        MG_SET(pt->flag,ia);
      }
    }
    switch (pt->flag) {
    case 1: /* 1 edge split */
      ier = MMGS_split1(mesh,met,k,0,vx);
      ns++;
      break;

    case 2: /* 1 edge split */
      ier = MMGS_split1(mesh,met,k,1,vx);
      ns++;
      break;

    case 4: /* 1 edge split */
      ier = MMGS_split1(mesh,met,k,2,vx);
      ns++;
      break;

    case 3: case 5: case 6: /* 2 edges split */
      ier = MMGS_split2(mesh,met,k,vx);
      ns++;
      break;

    default :
      assert(pt->flag == 0);
      break;
    }
    if ( !ier ) return 0;
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7" MMG5_PRId " splitted\n",ns);

  /* reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

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
static int MMGS_setref_ls(MMG5_pMesh mesh, MMG5_pSol sol) {
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

    /* Keep the initial triangle references of the mesh if iso==2, set
     * positive and negative ls refs otherwise */
    if ( mesh->info.iso != 2 ) {

      /* find if current reference should be splitted and the new positive and negative refs */
      ier = MMG5_isSplit(mesh,ref,&refint,&refext);
      if ( ier ) {
        if ( npl ) {
          assert( !nmn );
          pt->ref = refext;
        }
        else {
          assert ( nmn );
          pt->ref = refint;
        }
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

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set
 * \param met pointer toward a metric (optionnal)
 *
 * \return 0 if fail, 1 otherwise.
 *
 * Create implicit surface in mesh.
 *
 */
int MMGS_mmgs2(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol met) {
  MMG5_int k;

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION\n");

  /* Work only with the 0 level set */
  for (k=1; k<= sol->np; k++)
    sol->m[k] -= mesh->info.ls;

  /* Snap values of level set function if need be, then discretize it */
  if ( !MMGS_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: hashing problem (1). Exit program.\n",
            __func__);
    return 0;
  }

  if ( !MMG5_snpval_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem with implicit function. Exit program.\n");
    return 0;
  }

  /* Removal of small parasitic components */
  if ( mesh->info.rmc > 0. && !MMG5_rmc(mesh,sol) ) {
    fprintf(stderr,"\n  ## Error in removing small parasitic components. Exit program.\n");
    return 0;
  }

  MMG5_DEL_MEM(mesh,mesh->adja);

  if ( mesh->info.iso != 2 ) {
    /* Reset the mesh->info.isoref field everywhere it appears */
    if ( !MMG5_resetRef(mesh) ) {
      fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
      return 0;
    }
  }

  if ( !MMGS_cuttri_ls(mesh,sol,met) ) {
    fprintf(stderr,"\n  ## Problem in discretizing implicit function. Exit program.\n");
    return 0;
  }

  if ( !MMGS_setref_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }

  /* Creation of adjacency relations in the mesh */
  if ( !MMGS_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* Check that the resulting mesh is manifold */
  if ( !MMGS_chkmanimesh(mesh) ) {
    fprintf(stderr,"\n  ## No manifold resulting situation. Exit program.\n");
    return 0;
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,sol->m);
  sol->np = 0;

  return 1;
}
