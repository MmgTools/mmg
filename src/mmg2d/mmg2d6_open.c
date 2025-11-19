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
 * \file mmg2d/mmg2d6_open.c
 * \brief Isosurface discretization for open surface.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "libmmg2d_private.h"
#include "mmgexterns_private.h"
#include "mmg2dexterns_private.h"

/**
 * \param mesh pointer to the mesh
 *
 * Reset mesh->info.isoref vertex and edge references to 0 in the open mode.
 *
 */
int MMG5_resetRef_ls_open(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0;
  MMG5_int        ref,k;
  int8_t          i;
  
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;

    for (i=0; i<3; i++) {
      p0 = &mesh->point[pt->v[i]];
      if ( pt->edg[i] == mesh->info.isoref ) {
        pt->tag[i] &= ~MG_BDY;
        pt->tag[i] &= ~MG_REF;
        pt->edg[i] = 0;
      }
      if ( p0->ref == mesh->info.isoref ) {
        p0->ref = 0;
        p0->tag = 0;
      }
    }
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set values.
 * \return 1.
 *
 * Set references to edges and vertices according to the sign of the level set functions phi and psi.
 *
 */
int MMG5_setref_lsopen(MMG5_pMesh mesh, MMG5_pSol phi, MMG5_pSol psi) {
  MMG5_pTria    pt;
  double        v,v1,vps,vps1;
  int           ier;
  MMG5_int      k,ip,ip1;
  int8_t        i,i1,i2;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    /* Set mesh->info.isoref ref at ls edges and at the points of these edges */
    for (i=0; i<3; i++) {
      ip   = pt->v[MMG5_inxt2[i]];
      ip1  = pt->v[MMG5_iprv2[i]];
      v    = phi->m[ip];
      v1   = phi->m[ip1];
      vps  = psi->m[ip];
      vps1 = psi->m[ip1];
      if ( v == 0.0 && v1 == 0.0 && vps <= 0.0 && vps1 <= 0.0 ) {
        pt->edg[i]  = mesh->info.isoref;
        pt->tag[i] |= MG_REF;
        i1 = MMG5_inxt2[i];
        i2 = MMG5_inxt2[i1];
        mesh->point[pt->v[i1]].ref = mesh->info.isoref;
        mesh->point[pt->v[i2]].ref = mesh->info.isoref;
      }
    }
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param phi pointer to the level set function phi
 * \param psi pointer to the second level set function psi
 *
 * \return 1 if success, 0 otherwise
 *
 * Effective discretization of the part of the 0 level set of psi contained in triangles having at least one vertex psi \leq 0.
 *
 */
int MMG2D_cuttri_lsopen_phi(MMG5_pMesh mesh, MMG5_pSol phi, MMG5_pSol psi){
  MMG5_pTria   pt;
  MMG5_pPoint  p0,p1;
  MMG5_Hash    hash;
  double       v0,v1,s,c[2];
  int          k,ns,nt,nb,np,ip0,ip1,vx[3];
  int8_t       i,ier;
  
  /* Reset flag field for points */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;
    
  /* Estimate the number of intersected edges or surface edges by the 0 level set */
  nb = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    
    /* Consider only triangles having one vertex with psi \leq 0 */
    for (i=0; i<3; i++) {
      ip0 = pt->v[i];
      v0  = psi->m[ip0];
      if ( v0 <= MMG5_EPSD2 ) break;
    }
    if ( i == 3 ) continue;
    
    for (i=0; i<3; i++) {
      ip0 = pt->v[MMG5_inxt2[i]];
      ip1 = pt->v[MMG5_iprv2[i]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];

      if ( p0->flag && p1->flag ) continue;
      
      v0 = phi->m[ip0];
      v1 = phi->m[ip1];
      
      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
        nb++;
        if ( !p0->flag ) p0->flag = nb;
        if ( !p1->flag ) p1->flag = nb;
      }
    }
  }

  if ( !nb ) return 1;
    
  /* Create the intersection points between the edges in the mesh and the 0 level set */
  if ( !MMG5_hashNew(mesh,&hash,nb,2*nb) ) return 0;
  
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    
    for (i=0; i<3; i++) {
      ip0 = pt->v[i];
      v0  = psi->m[ip0];
      if ( v0 <= MMG5_EPSD2 ) break;
    }
    if ( i == 3 ) continue;

    for (i=0; i<3; i++) {
      ip0 = pt->v[MMG5_inxt2[i]];
      ip1 = pt->v[MMG5_iprv2[i]];

      np = MMG5_hashGet(&hash,ip0,ip1);
      if ( np ) continue;
      
      v0 = phi->m[ip0];
      v1 = phi->m[ip1];

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      
      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;
      
      /* Intersection point between edge p0p1 and the 0 level set */
      s = v0 / (v0-v1);
      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);

      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);

      np = MMG2D_newPt(mesh,c,0);
      if ( !np ) {
       /* reallocation of point table */
        MMG2D_POINT_REALLOC(mesh,psi,np,mesh->gap,
                            fprintf(stderr,"\n  ## Error: %s: unable to"
                                    " allocate a new point.\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            return 0;,
                            c,0);
      }
      phi->m[np] = 0.0;
      
      /* Interpolate psi at newly created vertex */
      psi->m[np] = psi->m[ip0] + s*(psi->m[ip1]-psi->m[ip0]);

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
      ip0 = pt->v[MMG5_inxt2[i]];
      ip1 = pt->v[MMG5_iprv2[i]];

      vx[i] = MMG5_hashGet(&hash,ip0,ip1);

      if ( vx[i] ) MG_SET(pt->flag,i);
    }

    switch( pt->flag ) {
      /* 1 edge split -> 0-+ */
      case 1: case 2: case 4:
        ier = MMG2D_split1(mesh,psi,k,vx);
        ns++;
        break;

      /* 2 edge split -> +-- or -++ */
      case 3: case 5: case 6:
        ier = MMG2D_split2(mesh,psi,k,vx);
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
 * \param mesh pointer to the mesh
 * \param phi pointer to the level set function phi
 * \param psi pointer to the second level set function psi
 *
 * \return 1 if success, 0 otherwise
 *
 * Effective discretization of the part of 0 level set of psi intersecting the 0 level set of phi.
 *
 */
int MMG2D_cuttri_lsopen_psi(MMG5_pMesh mesh, MMG5_pSol phi, MMG5_pSol psi){
  MMG5_pTria   pt;
  MMG5_pPoint  p0,p1;
  MMG5_Hash    hash;
  double       v0,v1,s,c[2];
  int          k,ns,nt,nb,np,ip0,ip1,vx[3];
  int8_t       i,ier;
  
  /* Reset flag field for points */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;
  
  /* Estimate the number of concerned edges */
  nb = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    
    for (i=0; i<3; i++) {
      ip0 = pt->v[MMG5_inxt2[i]];
      ip1 = pt->v[MMG5_iprv2[i]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];
      
      v0 = phi->m[ip0];
      v1 = phi->m[ip1];

      if ( fabs(v0) > MMG5_EPSD2 || fabs(v1) > MMG5_EPSD2 ) continue;
      if ( p0->flag && p1->flag ) continue;
      
      v0 = psi->m[ip0];
      v1 = psi->m[ip1];
      
      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
        nb++;
        if ( !p0->flag ) p0->flag = nb;
        if ( !p1->flag ) p1->flag = nb;
      }
    }
  }
  if ( !nb ) return 1;
  
  /* Create the intersection points between the edges in the mesh and the 0 level set of psi */
  if ( !MMG5_hashNew(mesh,&hash,nb,2*nb) ) return 0;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      ip0 = pt->v[MMG5_inxt2[i]];
      ip1 = pt->v[MMG5_iprv2[i]];

      np = MMG5_hashGet(&hash,ip0,ip1);
      if ( np ) continue;
      
      v0 = phi->m[ip0];
      v1 = phi->m[ip1];
      if ( fabs(v0) > MMG5_EPSD2 || fabs(v1) > MMG5_EPSD2 ) continue;
    
      v0 = psi->m[ip0];
      v1 = psi->m[ip1];

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];

      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;

      /* Intersection point between edge p0p1 and the 0 level set */
      s = v0 / (v0-v1);
      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);

      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);

      np = MMG2D_newPt(mesh,c,0);
      if ( !np ) {
       /* reallocation of point table */
        MMG2D_POINT_REALLOC(mesh,phi,np,mesh->gap,
                            fprintf(stderr,"\n  ## Error: %s: unable to"
                                    " allocate a new point.\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            return 0;,
                            c,0);
      }
      psi->m[np] = 0.0;

      /* Interpolate phi at newly created vertex */
      phi->m[np] = phi->m[ip0] + s*(phi->m[ip1]-phi->m[ip0]);

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
      ip0 = pt->v[MMG5_inxt2[i]];
      ip1 = pt->v[MMG5_iprv2[i]];

      vx[i] = MMG5_hashGet(&hash,ip0,ip1);

      if ( vx[i] ) MG_SET(pt->flag,i);
    }

    switch( pt->flag ) {
      /* 1 edge split -> 0-+ */
      case 1: case 2: case 4:
        ier = MMG2D_split1(mesh,psi,k,vx);
        ns++;
        break;

      /* 2 edge split -> +-- or -++ */
      case 3: case 5: case 6:
        ier = MMG2D_split2(mesh,psi,k,vx);
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

/* Main function of the -lsopen mode */
int MMG2D_mmg2d6_open(MMG5_pMesh mesh, MMG5_pSol phi,MMG5_pSol psi) {
  
  /* Set information to preserve internal edges which may not delimit subdomains with different references */
  mesh->info.opnbdy = 1;
  
  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION (OPEN CURVE)\n");
  }

  if ( mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: Isosurface extraction not available with"
            " hybrid meshes. Exit program.\n");
    return 0;
  }
  
  /* Transfer the boundary edge references to the triangles */
  if ( !MMG2D_assignEdge(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }

  /* Snapping */
  if ( !MMG5_snpval_ls(mesh,phi) ) {
    fprintf(stderr,"\n  ## Wrong input implicit function. Exit program.\n");
    return 0;
  }
  
  if ( !MMG5_snpval_ls(mesh,psi) ) {
    fprintf(stderr,"\n  ## Wrong input implicit function. Exit program.\n");
    return 0;
  }
  
  /* Removal of small parasitic components: to do */

  /* No need to keep adjacencies from now on */
  MMG5_DEL_MEM(mesh,mesh->adja);
  
  /* Reset the mesh->info.isoref field everywhere it appears */
  if ( !mesh->info.kiso ) {
    if ( !MMG5_resetRef_ls_open(mesh) ) {
      fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
      return 0;
    }
  }
  
  /* Effective splitting of the crossed triangles */
  if ( !MMG2D_cuttri_lsopen_phi(mesh,phi,psi) ) {
    fprintf(stderr,"\n  ## Problem in cutting triangles. Exit program.\n");
    return 0;
  }
  
  if ( !MMG2D_cuttri_lsopen_psi(mesh,phi,psi) ) {
    fprintf(stderr,"\n  ## Problem in cutting triangles. Exit program.\n");
    return 0;
  }
  
  /* Set references to inserted edges and vertices */
  if ( !MMG5_setref_lsopen(mesh,phi,psi) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }
  
  /* Creation of adjacency relations in the mesh */
  if ( !MMG2D_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }
  
  /* Clean memory */
  MMG5_DEL_MEM(mesh,phi->m);
  phi->np = 0;
  MMG5_DEL_MEM(mesh,psi->m);
  psi->np = 0;
  
  return 1;
}
