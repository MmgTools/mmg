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

extern int8_t ddb;

/**
 * \param mesh pointer toward the mesh.
 *
 * Reset mesh->info.isoref vertex and tetra references to 0.
 *
 * \warning to improve: for now, entities linked to the old ls (corners,required
 * points, normals/tangents, triangles and edges) are deleted in loadMesh. It
 * would be better to analyze wich entities must be keeped and which one must be
 * deleted depending on the split/nosplit infos.
 */
int MMG3D_resetRef_lssurf(MMG5_pMesh mesh) {
 
  printf("TO DO !!!\n");
  
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set function.
 * \return 1 if success, 0 if fail.
 *
 * Snap values of the level set function very close to 0 to exactly 0,
 * and prevent nonmanifold patterns from being generated.
 *
 */
int MMG3D_snpval_lssurf(MMG5_pMesh mesh,MMG5_pSol sol) {
  
  printf("TO DO !!!!\n");
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \param met pointer toward a metric (non-mandatory).
 * \return 1 if success, 0 otherwise.
 *
 * Proceed to discretization of the trace of the implicit function
 * carried by sol into the surface part of mesh,
 * once values of sol have been snapped/checked
 *
 */
int MMG3D_cuttet_lssurf(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_pSol met){
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1;
  MMG5_Hash     hash;
  double        c[3],v0,v1,s;
  int           vx[6],nb,k,ip0,ip1,np,ns,ne,ier,src,refext,refint;
  int8_t        ia,iface,j,npneg;
  static int8_t mmgWarn = 0;

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Compute the number nb of intersection points on edges */
  nb = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    for (iface=0; iface<4; iface++) {
      if ( !(pxt->ftag[iface] & MG_BDY) ) continue;
      for (j=0; j<3; j++) {
        ia = MMG5_iarf[iface][j];
        ip0 = pt->v[MMG5_iare[ia][0]];
        ip1 = pt->v[MMG5_iare[ia][1]];
        
        p0  = &mesh->point[ip0];
        p1  = &mesh->point[ip1];
        
        if ( p0->flag && p1->flag )  continue;
        
        v0  = sol->m[ip0]-mesh->info.ls;
        v1  = sol->m[ip1]-mesh->info.ls;
        
        if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
          nb ++;
          if ( !p0->flag ) p0->flag = nb;
          if ( !p1->flag ) p1->flag = nb;
        }
      }
    }
  }
  if ( ! nb )  return 1;

  /* Create intersection points at 0 isovalue and set flags to tetras */
  if ( !MMG5_hashNew(mesh,&hash,nb,7*nb) ) return 0;
  
  /* Hash all required edges, and put ip = -1 in hash structure */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /* avoid split of edges belonging to a required tet */
    if ( pt->tag & MG_REQ ) {
      for (ia=0; ia<6; ia++) {
        ip0 = pt->v[MMG5_iare[ia][0]];
        ip1 = pt->v[MMG5_iare[ia][1]];
        np  = -1;
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
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

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    for (iface=0; iface<4; iface++) {
      for (j=0; j<3; j++) {
        ia = MMG5_iarf[iface][j];
        ip0 = pt->v[MMG5_iare[ia][0]];
        ip1 = pt->v[MMG5_iare[ia][1]];

        np  = MMG5_hashGet(&hash,ip0,ip1);
        if ( np>0 )  continue;
        if ( !MMG5_isSplit(mesh,pt->ref,&refint,&refext) ) continue;
        
        p0 = &mesh->point[ip0];
        p1 = &mesh->point[ip1];
        v0 = sol->m[ip0]-mesh->info.ls;
        v1 = sol->m[ip1]-mesh->info.ls;
        if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
        else if ( MG_SMSGN(v0,v1) )  continue;
        else if ( !p0->flag || !p1->flag )  continue;
        
        npneg = (np<0);

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
          int oldnpmax = mesh->npmax;
          MMG3D_POINT_REALLOC(mesh,sol,np,MMG5_GAP,
                             fprintf(stderr,"\n  ## Error: %s: unable to"
                                     " allocate a new point\n",__func__);
                             MMG5_INCREASE_MEM_MESSAGE();
                             return 0
                             ,c,0,src);
          if( met ) {
            if( met->m ) {
              MMG5_ADD_MEM(mesh,(met->size*(mesh->npmax-met->npmax))*sizeof(double),
                         "larger solution",
                         MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,);
                         mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point);
                         mesh->npmax = oldnpmax;
                         mesh->np = mesh->npmax-1;
                         mesh->npnil = 0;
                         return 0);
              MMG5_SAFE_REALLOC(met->m,met->size*(met->npmax+1),
                              met->size*(mesh->npmax+1),
                              double,"larger solution",
                              MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,);
                              mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point);
                              mesh->npmax = oldnpmax;
                              mesh->np = mesh->npmax-1;
                              mesh->npnil = 0;
                              return 0);
            }
            met->npmax = mesh->npmax;
          }
        }
        
        sol->m[np] = mesh->info.ls;

        /* If user provide a metric, interpolate it at the new point */
        if ( met && met->m ) {
          if ( met->size > 1 ) {
            ier = MMG3D_intmet33_ani(mesh,met,k,ia,np,s);
          }
          else {
            ier = MMG5_intmet_iso(mesh,met,k,ia,np,s);
          }
          if ( ier <= 0 ) {
            // Unable to compute the metric
            fprintf(stderr,"\n  ## Error: %s: unable to"
                    " interpolate the metric during the level-set"
                    " discretization\n",__func__);
            return 0;
          }
        }
        
        if ( npneg ) {
          /* We split a required edge */
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
  }

  /* Proceed to splitting, according to flags to tets */
  ne  = mesh->ne;
  ns  = 0;
  ier = 1;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    pt->flag = 0;
    memset(vx,0,6*sizeof(int));
    
    for (ia=0; ia<6; ia++) {
      vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_iare[ia][0]],pt->v[MMG5_iare[ia][1]]);
      if ( vx[ia] > 0 )  MG_SET(pt->flag,ia);
    }
    
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
      ier = MMG5_split1(mesh,met,k,vx,1);
      ns++;
      break;

    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      ier = MMG5_split2sf(mesh,met,k,vx,1);
      ns++;
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      ier = MMG5_split3cone(mesh,met,k,vx,1);
      ns++;
      break;

    case 30: case 45: case 51:
      ier = MMG5_split4op(mesh,met,k,vx,1);
      ns++;
      break;

    default :
      assert(pt->flag == 0);
      break;
    }
    if ( !ier ) return 0;
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",ns);

  MMG5_DEL_MEM(mesh,hash.item);
  return ns;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \return 1.
 *
 * Set references to surface triangles according to the sign of the level set function.
 *
 */
int MMG3D_setref_lssurf(MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  double        v;
  int           k,ip,ref,refint,refext,ier;
  int8_t        nmns,npls,nz,i,j;
  
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    
    if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++) {
      if ( !(pxt->ftag[i] & MG_BDY) ) continue;
      ref = pxt->ref[i];
      nmns = npls = nz = 0;
      
      for (j=0; j<3; j++) {
        ip = pt->v[MMG5_idir[i][j]];
        v  = sol->m[ip]-mesh->info.ls;
        if ( v > 0.0 )
          npls++;
        else if ( v < 0.0 )
          nmns++;
        else
          nz ++;
      }
      
      assert(nz < 3);
      ier = MMG5_isSplit(mesh,ref,&refint,&refext);
      
      if ( npls ) {
        if ( ier ) {
          assert(!nmns);
          pxt->ref[i] = refext;
        }
      }
      else {
        if ( ier ) {
          assert(nmns);
          pxt->ref[i] = refint;
        }
      }
    }
  }
  return 1;
}
