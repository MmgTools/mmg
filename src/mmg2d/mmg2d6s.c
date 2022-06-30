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

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set
 *
 * \return 1 if success, 0 if fail
 *
 * Snap values of sol very close to 0 to 0 exactly in the case of surface ls splitting
 * (to avoid very small triangles in cutting)
 */
int MMG2D_snapval_lssurf(MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_pTria       pt;
  MMG5_pPoint      p0,p1;
  double           v0,v1,*tmp;
  int              k,ns,nc,ip0,ip1;
  int8_t           i,i0,i1;

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

  /* Unsnap values that have been put to 0, entailing a component reduced to a point */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;

    for (i=0; i<3; i++) {
      if ( !(pt->tag[i] & MG_BDY) ) continue;
      i0 = MMG5_inxt2[i];
      i1 = MMG5_inxt2[i0];

      ip0 = pt->v[i0];
      ip1 = pt->v[i1];

      v0 = sol->m[ip0];
      v1 = sol->m[ip1];

      p0 = &mesh->point[pt->v[i0]];
      p1 = &mesh->point[pt->v[i1]];

      if ( fabs(v0) < MMG5_EPS && fabs(v1) < MMG5_EPS ) {
        if ( p0->flag ) {
          if ( tmp[ip0] < 0.0 )
            sol->m[ip0] = -100.0*MMG5_EPS;
          else
            sol->m[ip0] = 100.0*MMG5_EPS;
          nc++;
          p0->flag = 0;
        }
        else if ( p1->flag ) {
          if ( tmp[ip1] < 0.0 )
            sol->m[ip1] = -100.0*MMG5_EPS;
          else
            sol->m[ip1] = 100.0*MMG5_EPS;
          nc++;
          p1->flag = 0;
        }
      }
    }
  }

  MMG5_DEL_MEM ( mesh, tmp );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && ns+nc > 0 )
    fprintf(stdout,"     %8d points snapped, %d corrected\n",ns,nc);

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 *
 * Reset mesh->info.isoref vertex references to 0.
 *
 */
int MMG2D_resetRef_lssurf(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0,p1;
  int             k,ref;
  int8_t          i,i0,i1;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;

    for (i=0; i<3; i++) {
      if ( !(pt->tag[i] & MG_BDY) ) continue;

      if( !MMG5_getStartRef(mesh,pt->edg[i],&ref) ) return 0;
      pt->edg[i] = ref;

      i0 = MMG5_inxt2[i];
      i1 = MMG5_inxt2[i0];

      p0 = &mesh->point[pt->v[i0]];
      p1 = &mesh->point[pt->v[i1]];

      if ( p0->ref == mesh->info.isoref ) p0->ref = 0;
      if ( p1->ref == mesh->info.isoref ) p1->ref = 0;
    }
  }

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
int MMG2D_cuttri_lssurf(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSol met){
  MMG5_pTria   pt;
  MMG5_pPoint  p0,p1;
  MMG5_Hash    hash;
  double       v0,v1,s,c[2];
  int          k,ip0,ip1,nb,np,nt,ns,refint,refext,vx[3];
  int8_t       i,i0,i1,ier;

  /* Reset flag field for points */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Estimate the number of intersected edges or surface edges by the 0 level set */
  nb = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {

      if ( !(pt->tag[i] & MG_BDY) ) continue;

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

  /* Create the intersection points between the surface edges in the mesh and the 0
   * level set */
  if ( !MMG5_hashNew(mesh,&hash,nb,2*nb) ) return 0;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {

      if ( !(pt->tag[i] & MG_BDY) ) continue;

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
int MMG2D_setref_lssurf(MMG5_pMesh mesh, MMG5_pSol sol){
  MMG5_pTria     pt;
  double         v,v1;
  int            k,ier,ip1,ref,refint,refext;
  int8_t         i,i1,i2,j,nmn,npl,nz;

  /* Travel all surface edges (via triangles) */
  for(k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      if ( !(pt->tag[i] & MG_BDY) ) continue;
      ref = pt->edg[i];
      nmn = npl = nz = 0;

      i1 = i;
      for (j=0; j<2; j++) {
        i1  = MMG5_inxt2[i1];
        ip1 = pt->v[i1];
        v1  = sol->m[ip1];

        if ( v1 > 0.0 )
          npl++;
        else if ( v1 < 0.0 )
          nmn++;
        else
          nz++;
      }

      assert(nz < 2);
      ier = MMG5_isSplit(mesh,ref,&refint,&refext);
      if ( npl ) {
        if (ier ) {
          assert ( !nmn );
          pt->edg[i] = refext;
        }
      }
      else {
        if ( ier ) {
          assert ( !npl );
          pt->edg[i] = refint;
        }
      }
    }
  }

  /* Set MG_ISO ref at vertices on the ls */
  for (k=1; k<=mesh->np; k++) {
    v = sol->m[k];
    if ( v == 0.0 ) {
      mesh->point[k].ref = MG_ISO;
    }
  }

  return 1;
}

/* Main function of the -ls surf mode */
int MMG2D_mmg2d6s(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_pSol met) {
  int k;

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION (BOUNDARY PART)\n");

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

  /* Creation of tria adjacency relations in the mesh */
  if ( !MMG2D_hashTria(mesh) ) {
     fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* Set tags to triangles from geometric configuration */
  if ( !MMG2D_setadj(mesh) ) {
    fprintf(stderr,"\n  ## Problem in function setadj. Exit program.\n");
    return 0;
  }

  /* Snap values of the level set function which are very close to 0 to 0 exactly */
  if ( !MMG2D_snapval_lssurf(mesh,sol) ) {
    fprintf(stderr,"\n  ## Wrong input implicit function. Exit program.\n");
    return 0;
  }

  /* RMC : on verra */
  if ( mesh->info.rmc > 0 ) {
    fprintf(stdout,"\n  ## Warning: rmc option not implemented for boundary"
            " isosurface extraction.\n");
  }

  /* No need to keep adjacencies from now on */
  MMG5_DEL_MEM(mesh,mesh->adja);

  /* Reset the mesh->info.isoref field everywhere it appears */
  if ( !MMG2D_resetRef_lssurf(mesh) ) {
    fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
    return 0;
  }

  /* Effective splitting of the crossed triangles */
  if ( !MMG2D_cuttri(mesh,sol,met) ) {
    fprintf(stderr,"\n  ## Problem in cutting triangles. Exit program.\n");
    return 0;
  }

  /* Set references on the interior / exterior triangles*/
  if ( !MMG2D_setref_lssurf(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }

  /* Creation of adjacency relations in the mesh */
  if ( !MMG2D_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,sol->m);
  sol->np = 0;

  MMG5_DEL_MEM( mesh,mesh->info.mat );

  return 1;
}
