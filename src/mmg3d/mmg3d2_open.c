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
 * \file mmg3d/mmg3d2.c
 * \brief Create implicit surface in mesh.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */
 
#include "libmmg3d.h"
#include "libmmg3d_private.h"
#include "mmg3dexterns_private.h"

extern int8_t ddb;

/**
 * \param mesh pointer to the mesh
 *
 * Reset mesh->info.isoref vertex and edge references to 0 in the open mode.
 *
 */
int MMG3D_resetRef_lsopen(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0;
  MMG5_int        ref,k;
  int8_t          i;
  
  printf("MMG3D_resetRef_ls_open: TO DO\n");
  
  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param phi pointer to the level set function phi
 * \param psi pointer to the second level set function psi
 *
 * \return 1 if success, 0 otherwise
 *
 * Effective discretization of the part of the 0 level set of psi contained in tetras having at least one vertex psi \leq 0.
 *
 */
int MMG3D_cuttet_lsopen_phi(MMG5_pMesh mesh, MMG5_pSol phi, MMG5_pSol psi){
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1;
  MMG5_Hash     hash;
  double        v0,v1,s,c[3];
  int           ier;
  MMG5_int      ip0,ip1,k,nb,np,ne,ns,src,refint,refext,vx[6];
  int8_t        i,j,ia,npneg;
  static int8_t mmgWarn = 0;

  
  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;
    
  /* Compute the number nb of intersection points on edges */
  nb = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    
    /* Consider only tets having one vertex with \psi \leq 0 */
    for (i=0; i<4; i++) {
      ip0 = pt->v[i];
      v0  = psi->m[ip0];
      if ( v0 <= MMG5_EPSD2 ) break;
    }
    if ( i == 4 ) continue;
    
    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];
      if ( p0->flag && p1->flag )  continue;
      v0  = phi->m[ip0];
      v1  = phi->m[ip1];
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
  
  /* Create intersection points at 0 isovalue */
  if ( !MMG5_hashNew(mesh,&hash,nb,7*nb) ) return (0);
  
  /* Hash all required edges and put ip = -1 in hash structure */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /* Avoid split of edges belonging to a required tet */
    if ( pt->tag & MG_REQ ) {
      for (ia=0; ia<6; ia++) {
        ip0 = pt->v[MMG5_iare[ia][0]];
        ip1 = pt->v[MMG5_iare[ia][1]];
        np  = -1;
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return (-1);
      }
      continue;
    }
    
    if ( !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];
    for (ia=0; ia<4; ia++) {
      if ( !(pxt->ftag[ia] & MG_BDY) ) continue;

      for (j=0; j<3; j++) {
        if ( !(pxt->tag[ MMG5_iarf[ia][j] ] & MG_REQ) ) continue;

        ip0 = pt->v[MMG5_idir[ia][MMG5_inxt2[j]]];
        ip1 = pt->v[MMG5_idir[ia][MMG5_iprv2[j]]];
        np  = -1;
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
      }
    }
  }
  
  /* Hash remaining split edges */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    
    /* Subtle point: considering only tets cut by psi requires the use of patterns not corresponding to LS disc */
    /* On the other hand, considering a pattern issued from LS disc and deactivating points amounts to setting these to 0 -> leads to already implemented pattern */
    /* for (i=0; i<4; i++) {
      ip0 = pt->v[i];
      v0  = psi->m[ip0];
      if ( v0 <= MMG5_EPSD2 ) break;
    }
    if ( i == 4 ) continue;*/

    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      np  = MMG5_hashGet(&hash,ip0,ip1);

      if ( np > 0 )  continue;

      if ( !MMG5_isSplit(mesh,pt->ref,&refint,&refext) ) continue;

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      v0 = phi->m[ip0];
      v1 = phi->m[ip1];
      if ( fabs(v0) <= MMG5_EPSD2 || fabs(v1) <= MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;
      
      npneg = ( np < 0 );

      s = v0 / (v0-v1);
      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);
      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
      c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

#ifdef USE_POINTMAP
      src = p0->src;
#else
      src = 1;
#endif
      np = MMG3D_newPt(mesh,c,0,src);
      if ( !np ) {
        MMG5_int oldnpmax = mesh->npmax;
        MMG3D_POINT_REALLOC(mesh,phi,np,MMG5_GAP,
                             fprintf(stderr,"\n  ## Error: %s: unable to"
                                     " allocate a new point\n",__func__);
                             MMG5_INCREASE_MEM_MESSAGE();
                             return 0
                             ,c,0,src);
      }
      
      /* Values of phi and psi at newly created vertex */
      phi->m[np] = 0.0;
      psi->m[np] = psi->m[ip0] + s*(psi->m[ip1]-psi->m[ip0]);

      /* Case where a required edge is split */
      if ( npneg ) {
        if ( !mmgWarn ) {
          mmgWarn = 1;
          fprintf(stderr,"  ## Warning: %s: the level-set intersect at least"
                  " one required entity. Required entity ignored.\n\n",__func__);
        }
        MMG5_hashUpdate(&hash,ip0,ip1,np);
      }
      else {
        MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
      }
    }
  }
  
  /* Proceed to splitting, according to flags to tets */
  ne  = mesh->ne;
  ns  = 0;
  ier = 1;
  
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    
    pt->flag = 0;
    memset(vx,0,6*sizeof(MMG5_int));
    
    for (ia=0; ia<6; ia++) {
      vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_iare[ia][0]],pt->v[MMG5_iare[ia][1]]);
      if ( vx[ia] > 0 )  MG_SET(pt->flag,ia);
    }
        
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
      ier = MMG5_split1(mesh,psi,k,vx,1);
      ns++;
      break;

    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      ier = MMG5_split2sf(mesh,psi,k,vx,1);
      ns++;
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      ier = MMG5_split3cone(mesh,psi,k,vx,1);
      ns++;
      break;

    case 30: case 45: case 51:
      ier = MMG5_split4op(mesh,psi,k,vx,1);
      ns++;
      break;

    default :
      if ( pt->flag != 0 ) {
        printf("   Error in function cuttet_lsopen_phi: use of non implemented pattern; exit.\n");
        exit(0);
      }
      break;
    }
    if ( !ier ) return 0;
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7" MMG5_PRId " splitted\n",ns);

  MMG5_DEL_MEM(mesh,hash.item);
  return (ns);
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
int MMG3D_cuttet_lsopen_psi(MMG5_pMesh mesh, MMG5_pSol phi, MMG5_pSol psi){
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1;
  MMG5_Hash     hash;
  double        v0,v1,s,c[3];
  MMG5_int      ip0,ip1,k,nb,ns,ne,np,src,refint,refext,vx[6];
  int           ier;
  int8_t        ia,j,npneg;
  static int8_t mmgWarn = 0;

  /* Reset flag field for points */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;
    
  /* Compute the number nb of intersection points on edges */
  nb = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    
    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];
      if ( p0->flag && p1->flag )  continue;
      
      v0  = phi->m[ip0];
      v1  = phi->m[ip1];
      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 ) continue;
      
      v0  = psi->m[ip0];
      v1  = psi->m[ip1];
      
      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
        nb++;
        if ( !p0->flag ) p0->flag = nb;
        if ( !p1->flag ) p1->flag = nb;
      }
    }
  }
  if ( ! nb )  return 1;
  
  /* Create intersection points at 0 isovalue */
  if ( !MMG5_hashNew(mesh,&hash,nb,7*nb) ) return (0);
  
  /* Hash all required edges and put ip = -1 in hash structure */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /* Avoid split of edges belonging to a required tet */
    if ( pt->tag & MG_REQ ) {
      for (ia=0; ia<6; ia++) {
        ip0 = pt->v[MMG5_iare[ia][0]];
        ip1 = pt->v[MMG5_iare[ia][1]];
        np  = -1;
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return (-1);
      }
      continue;
    }
    
    if ( !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];
    for (ia=0; ia<4; ia++) {
      if ( !(pxt->ftag[ia] & MG_BDY) ) continue;

      for (j=0; j<3; j++) {
        if ( !(pxt->tag[ MMG5_iarf[ia][j] ] & MG_REQ) ) continue;

        ip0 = pt->v[MMG5_idir[ia][MMG5_inxt2[j]]];
        ip1 = pt->v[MMG5_idir[ia][MMG5_iprv2[j]]];
        np  = -1;
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
      }
    }
  }
  
  /* Hash remaining split edges */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      np  = MMG5_hashGet(&hash,ip0,ip1);

      if ( np > 0 )  continue;
      if ( !MMG5_isSplit(mesh,pt->ref,&refint,&refext) ) continue;
      
      v0 = phi->m[ip0];
      v1 = phi->m[ip1];
      if ( fabs(v0) > MMG5_EPSD2 || fabs(v1) > MMG5_EPSD2 ) continue;

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      v0 = psi->m[ip0];
      v1 = psi->m[ip1];
      
      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;

      npneg = ( np < 0 );

      s = v0 / (v0-v1);
      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);
      
      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
      c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

#ifdef USE_POINTMAP
      src = p0->src;
#else
      src = 1;
#endif
      np = MMG3D_newPt(mesh,c,0,src);
      if ( !np ) {
        MMG5_int oldnpmax = mesh->npmax;
        MMG3D_POINT_REALLOC(mesh,psi,np,MMG5_GAP,
                             fprintf(stderr,"\n  ## Error: %s: unable to"
                                     " allocate a new point\n",__func__);
                             MMG5_INCREASE_MEM_MESSAGE();
                             return 0
                             ,c,0,src);
      }
      
      /* Values of phi and psi at newly created vertex */
      psi->m[np] = 0.0;
      phi->m[np] = phi->m[ip0] + s*(phi->m[ip1]-phi->m[ip0]);

      /* Case where a required edge is split */
      if ( npneg ) {
        if ( !mmgWarn ) {
          mmgWarn = 1;
          fprintf(stderr,"  ## Warning: %s: the level-set intersect at least"
                  " one required entity. Required entity ignored.\n\n",__func__);
        }
        MMG5_hashUpdate(&hash,ip0,ip1,np);
      }
      else
        MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
    }
  }
  
  /* Proceed to splitting, according to flags to tets */
  ne  = mesh->ne;
  ns  = 0;
  ier = 1;
  
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    
    pt->flag = 0;
    memset(vx,0,6*sizeof(MMG5_int));
    
    for (ia=0; ia<6; ia++) {
      vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_iare[ia][0]],pt->v[MMG5_iare[ia][1]]);
      if ( vx[ia] > 0 )  MG_SET(pt->flag,ia);
    }
    
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
      ier = MMG5_split1(mesh,psi,k,vx,1);
      ns++;
      break;

    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      ier = MMG5_split2sf(mesh,psi,k,vx,1);
      ns++;
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      ier = MMG5_split3cone(mesh,psi,k,vx,1);
      ns++;
      break;

    case 30: case 45: case 51:
      ier = MMG5_split4op(mesh,psi,k,vx,1);
      ns++;
      break;

    default :
      if ( pt->flag != 0 ) {
        printf("   Error in function cuttet_lsopen_phi: use of non implemented pattern; exit.\n");
        exit(0);
      }
      break;
    }
    if ( !ier ) return 0;
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7" MMG5_PRId " splitted\n",ns);
      
  MMG5_DEL_MEM(mesh,hash.item);
  return (ns);
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set values.
 * \return 1.
 *
 * Set references to faces, edges and vertices according to the sign of the level set functions phi and psi.
 *
 */
int MMG3D_setref_lsopen(MMG5_pMesh mesh, MMG5_pSol phi, MMG5_pSol psi) {
  MMG5_pTetra   pt,pt1;
  MMG5_pxTetra  pxt,pxt1;
  double        vopp,v0,v1,v2,vp0,vp1,vp2;
  MMG5_int      *adja,k,jel,ip,ip0,ip1,ip2;
  int           i,voy;

  /* xtetra structure and adjacencies must be set */
  if ( !mesh->xtetra ) {
    fprintf(stderr,"\n  ## Error: %s: the xtetra array must be allocated.\n",
      __func__);
    return 0;
  }
  if ( !mesh->adja ) {
    fprintf(stderr,"\n  ## Error: %s: the ajda array must be allocated.\n",
      __func__);
    return 0;
  }
      
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      /* Face is already stored */
      if ( !adja[i] ) continue;
      
      jel = adja[i] / 4;
      voy = adja[i] % 4;
      pt1 = &mesh->tetra[jel];
      
      ip  = pt->v[i];
      ip0 = pt->v[MMG5_idir[i][0]];
      ip1 = pt->v[MMG5_idir[i][1]];
      ip2 = pt->v[MMG5_idir[i][2]];
      
      vopp = phi->m[ip];
      v0   = phi->m[ip0];
      v1   = phi->m[ip1];
      v2   = phi->m[ip2];
      
      vp0  = psi->m[ip0];
      vp1  = psi->m[ip1];
      vp2  = psi->m[ip2];
      
      if ( v0 != 0.0 || v1 != 0.0 || v2 != 0.0 || vp0 > 0.0 || vp1 > 0.0 || vp2 > 0.0 ) continue;
            
      /* Create both faces */
      if ( !pt->xt ) {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                            "larger xtetra table",
                            mesh->xt--;
                            fprintf(stderr,"  Exit program.\n"); return 0;);
        }
        pt->xt = mesh->xt;
      }
      
      if ( !pt1->xt ) {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                            "larger xtetra table",
                            mesh->xt--;
                            fprintf(stderr,"  Exit program.\n"); return 0;);
        }
        pt1->xt = mesh->xt;
      }
            
      pxt = &mesh->xtetra[pt->xt];
      if ( !(pxt->ftag[i] & MG_BDY) ) {
        pxt->ftag[i] |= MG_BDY;
        pxt->ref[i]  = mesh->info.isoref;
        if ( vopp < 0.0 ) MG_SET(pxt->ori,i);
      }
      
      pxt1 = &mesh->xtetra[pt1->xt];
      if ( !(pxt1->ftag[voy] & MG_BDY) ) {
        pxt1->ftag[voy] |= MG_BDY;
        pxt1->ref[voy]  = mesh->info.isoref;
        if ( vopp >= 0.0 ) MG_SET(pxt1->ori,voy);
      }
    }
  }

  return (1);
}


/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set.
 * \param met pointer to  a metric (optionnal).
 * \return 0 if fail, 1 otherwise.
 *
 * Create implicit surface in mesh.
 *
 */
int MMG3D_mmg3d2_open(MMG5_pMesh mesh,MMG5_pSol phi,MMG5_pSol psi) {
  char str[16]="";
  
  /* Preserve triangles separating tetras with same ref */
  mesh->info.opnbdy = 1;
    
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION (OPEN SURFACE) %s\n",str);

  if ( mesh->nprism || mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: Isosurface extraction not available with"
            " hybrid meshes. Exit program.\n");
    return 0;
  }
  
  /* Snap values of level set function if need be */
  /* Like in 2d... not sure it is needed */

  /* Create adjacencies (and pack deleted elts) */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* Make triangle orientations compatible with w/r tetras */
  if ( !MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
    return 0;
  }

  /* Correct triangles (remove twins, triangles that do not belong to a bdy if not opnbdy, ...) */
  if ( !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }

  /* Build hash table for initial edges */
  if ( !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    return 0;
  }
  
  /* Set the references from triangles to tetrahedra vis the structure xtetra */
  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }
  
  /* Reset the mesh->info.isoref field everywhere it appears */
  if ( !mesh->info.kiso ) {
    if ( !MMG3D_resetRef_lsopen(mesh) ) {
      fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
      return 0;
    }
  }
  
  /* Removal of small parasitic components: to do */
  
  /* Free memory */
  MMG5_DEL_MEM(mesh,mesh->adja);
  MMG5_DEL_MEM(mesh,mesh->adjt);
  MMG5_DEL_MEM(mesh,mesh->tria);
  mesh->nt = 0;

  /* Effective splitting of the crossed triangles */
  if ( !MMG3D_cuttet_lsopen_phi(mesh,phi,psi) ) {
    fprintf(stderr,"\n  ## Problem in cutting tetras. Exit program.\n");
    return 0;
  }
  
  if ( !MMG3D_cuttet_lsopen_psi(mesh,phi,psi) ) {
    fprintf(stderr,"\n  ## Problem in cutting triangles. Exit program.\n");
    return 0;
  }
  
  /* Re-create adjacencies */
  MMG5_DEL_MEM(mesh,mesh->adja);
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }
    
  /* Set references to inserted edges and vertices */
  if ( !MMG3D_setref_lsopen(mesh,phi,psi) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }
      
  /* Save mesh */
  // MMG3D_packMesh(mesh,phi,psi);
  // MMG3D_saveMesh(mesh,mesh->nameout);
  // exit(0);
  
  /* Clean memory */
  MMG5_DEL_MEM(mesh,phi->m);
  phi->np = 0;
  MMG5_DEL_MEM(mesh,psi->m);
  psi->np = 0;
  
  return 1;
}
