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
 * \file mmg3d/mmg3d1_delone.c
 * \brief Perform volume and surface mesh adaptation in delaunay mode.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 * Perform volume and surface mesh adaptation in delaunay mode (\a
 * MMG_PATTERN preprocessor flag set to OFF).
 *
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"
#include "mmg3dexterns.h"

#ifndef MMG_PATTERN

int8_t  ddb;

#define MMG3D_THRES_DEL     1.6
#define MMG3D_LOPTL_DEL     1.41
#define MMG3D_LFILTS_DEL    0.7
#define MMG3D_LFILTL_DEL    0.2

/**
 * \param mesh pointer toward mesh
 * \param k index of input tetra
 * \param imax index of edge in tetra \a k
 * \param i index of boundary face of tetra from which we will work
 * \param j index of edge in face \a i
 * \param pxt boundary tetra associated to \a k
 * \param ip1 first vertex of edge \a i
 * \param ip2 second vertex of edge \a i
 * \param p0 point \a ip1
 * \param p1 point \a ip2
 * \param ref edge ref (to fill)
 * \param tag edge tag (to fill)
 * \param o coordinates of new point along bezier edge (to fill)
 * \param to tangent at new point \a o (to fill if needed)
 * \param no1 first normal at new point \a o (to fill if needed)
 * \param no2 second normal at new point (to fill if needed)
 * \param list pointer toward edge shell (to fill)
 * \param ilist 2x edge shell size (+1 for a bdy edge)
 *
 * \return -1 for strong failure
 * \return 0 if we fail to compute new point and want to pass to next elt of the main loop
 * \return 1 if we fail to compute new point and want to try to collapse to short edge
 * \return 2 if we can compute new point.
 *
 * Build Bezier edge from the boundary face of a boundary tetra and compute
 * position and feature of new point along this edge.
 *
 */
static inline
int8_t MMG3D_build_bezierEdge(MMG5_pMesh mesh,MMG5_int k,
                              int8_t imax,int8_t i, int8_t j,
                              MMG5_pxTetra pxt,
                              MMG5_int ip1,MMG5_int ip2,
                              MMG5_pPoint p0, MMG5_pPoint p1,
                              MMG5_int *ref,int16_t *tag,
                              double o[3],double to[3],double no1[3],
                              double no2[3],int64_t *list,int *ilist) {
  MMG5_Tria   ptt;
  double      v[3];

  if ( (p0->tag & MG_PARBDY) && (p1->tag & MG_PARBDY) ) {
    /* Skip edge with extremities on parallel interfaces */
    return 0;
  }
  if ( !(MG_GET(pxt->ori,i)) ) {
    /* Treat triangles at interface of 2 subdomains from well oriented face */
    return 0;
  }

  *ref = pxt->edg[MMG5_iarf[i][j]];
  *tag = pxt->tag[MMG5_iarf[i][j]];
  if ( (*tag) & MG_REQ ) {
    /* No need to split required edges */
    return 0;
  }

  (*tag) |= MG_BDY;
  *ilist = MMG5_coquil(mesh,k,imax,list);
  if ( !(*ilist) ) {
    /* On of the tetra of the edge shell is required: we cannot split the edge */
    return 0;
  }
  else if ( (*ilist)<0 ) {
    /* Shell computation has failed */
    return -1;
  }

  /** a/ computation of bezier edge */
  if ( (*tag) & MG_NOM ){
    /* Edge is non-manifold */
    if( !MMG5_BezierNom(mesh,ip1,ip2,0.5,o,no1,to) ) {
      /* Unable to treat edge: pass to next elt */
      return 0;
    }
    else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
      assert( 0<=i && i<4 && "unexpected local face idx");
      MMG5_tet2tri(mesh,k,i,&ptt);
      MMG5_nortri(mesh,&ptt,no1);
    }
  }
  else if ( (*tag) & MG_GEO ) {
    /* Edge is ridge */
    if ( !MMG5_BezierRidge(mesh,ip1,ip2,0.5,o,no1,no2,to) ) {
      /* Unable to treat edge: pass to next elt */
//#warning why a continue here?
      return 0;
    }

    if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
      if ( !MMG3D_normalAndTangent_at_sinRidge(mesh,k,i,j,pxt,no1,no2,to) ) {
        return -1;
      }
    }
  }
  else if ( (*tag) & MG_REF ) {
    /* Edge is ref */
    if ( !MMG5_BezierRef(mesh,ip1,ip2,0.5,o,no1,to) ) {
      /* Unable to treat long edge: try to collapse short one */
      return 1;
    }
    if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
//#warning creation of sin-sin ref edge to see if it works without normal realloc
      assert( 0<=i && i<4 && "unexpected local face idx");
      MMG5_tet2tri(mesh,k,i,&ptt);
      MMG5_nortri(mesh,&ptt,no1);
    }
  }
  else {
    /* Longest edge is regular */
    if ( !MMG5_norface(mesh,k,i,v) ) {
      /* Unable to treat long edge: try to collapse short one */
      return 1;
    }
    if ( !MMG5_BezierReg(mesh,ip1,ip2,0.5,v,o,no1) ) {
      /* Unable to treat longest edge: try to collapse the shorter one */
      return 1;
    }
  }
  return 2;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param k index of tetra in which we work.
 * \param imax index in \a k of edge that we consider for split.
 * \param lmax length of edge \a imax.
 * \param lmaxtet length of largest edge of tetra \a k.
 * \param 1 if we want to check tetra with 4 ridge metrics.
 * \param ifilt pointer to store the number of vertices filtered by the PROctree.
 * \param ns pointer toward count of splits (has to be updated)
 * \param warn pointer to store a flag that warn the user in case of
 * reallocation error.
 * \param countMemFailure number of memory errors (to update)
 *
 * \return -2 for low failure (mesh has to be saved).
 *
 * \return -1 for strong failure.
 *
 * \return 0 if edge cannot be splitted and if we want to pass to next loop step
 * (next element or next tetra edge)
 *
 * \return 1 if edge cannot be splitted and we want to try to collapse too long
 * edge.
 *
 * \return 2 if edge has been splitted and we want to treat next element.
 *
 * \return 3 if nothing has been done (no error but no split either).
 *
 * Try to split \a imax if too large.
 *
 */
static inline
int MMG3D_mmg3d1_delone_split(MMG5_pMesh mesh, MMG5_pSol met,
                              MMG3D_pPROctree *PROctree,MMG5_int k,
                              int8_t imax,double lmax,double lmaxtet,
                              int8_t chkRidTet,MMG5_int *ifilt,MMG5_int *ns,
                              int *warn,int8_t *countMemFailure ) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1,ppt;
  MMG5_pxPoint  pxp;
  double        o[3],to[3],no1[3],no2[3],lfilt;
  int64_t       list[MMG3D_LMAX+2];
  MMG5_int      ip1,ip2;
  MMG5_int      src,ip,ref;
  int           ilist;
  int16_t       tag;
  int8_t        j,i,i1,i2,ifa0,ifa1;

  if ( lmax < MMG3D_LOPTL_DEL )  {
    /* Edge is small enough: nothing to do */
    return 3;
  }

  /** Edge is too long: try to split it */
  pt = &mesh->tetra[k];
  pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

  /* Try to treat the edge from a bdy face if possible */
  ifa0 = MMG5_ifar[imax][0];
  ifa1 = MMG5_ifar[imax][1];
  i  = (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
  j  = MMG5_iarfinv[i][imax];
  i1 = MMG5_idir[i][MMG5_inxt2[j]];
  i2 = MMG5_idir[i][MMG5_iprv2[j]];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  p0  = &mesh->point[ip1];
  p1  = &mesh->point[ip2];

  if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
    /** Edge belongs to a boundary face: try to split using patterns */
    /* Construction of bezier edge */
    int8_t ier = MMG3D_build_bezierEdge(mesh,k,imax,i,j,pxt,ip1,ip2,p0,p1,
                                        &ref,&tag,o,to,no1,no2,list,&ilist);
    if ( ier < 0 ) {
      /* Strong failure */
      return -1;
    }
    else if ( !ier ) {
      /* Unable to split edge: pass to next elt */
      return 0;
    }
    else if ( ier == 1 ) {
      /* Unable to split edge: try to collapse shortest edge */
      return 1;
    }

    /** b/ Edge splitting */
#ifdef USE_POINTMAP
    src = mesh->point[ip1].src;
#else
    src = 1;
#endif
    ip = MMG3D_newPt(mesh,o,tag,src);
    if ( !ip ){
      /* reallocation of point table */
      MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                          *warn=1;++(*countMemFailure);
                          return 1,
                          o,tag,src);
    }

    if ( met && met->m ) {
      if ( MMG5_intmet(mesh,met,k,imax,ip,0.5)<=0 ) {
        MMG3D_delPt(mesh,ip);
        return 1;
      }
    }
    ier = MMG3D_simbulgept(mesh,met,list,ilist,ip);
    assert ( (!mesh->info.ddebug) || (mesh->info.ddebug && ier != -1) );
    if ( ier == 2 || ier < 0 ) {
      /* sharp angle failure */
      MMG3D_delPt(mesh,ip);
      return 1;
    }
    else if ( ier == 0 ) {
      /* very bad quality failure */
      ier = MMG3D_dichoto1b(mesh,met,list,ilist,ip);
    }
    if ( ier == 1 ) {
      ier = MMG5_split1b(mesh,met,list,ilist,ip,1,1,chkRidTet);
    }

    /* if we realloc memory in MMG5_split1b pt and pxt pointers are not valid */
    pt = &mesh->tetra[k];
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    if ( ier < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to split.\n",__func__);
      MMG3D_delPt(mesh,ip);
      return -1;
    }
    else if ( ier == 0 || ier == 2 ) {
      MMG3D_delPt(mesh,ip);
      return 1;
    } else {
      (*ns)++;
      ppt = &mesh->point[ip];

      if ( MG_EDG(tag) || (tag & MG_NOM) )
        ppt->ref = ref;
      else
        ppt->ref = pxt->ref[i];
      ppt->tag = tag;

      pxp = &mesh->xpoint[ppt->xp];
      if ( tag & MG_NOM ){
        memcpy(pxp->n1,no1,3*sizeof(double));
        memcpy(ppt->n,to,3*sizeof(double));
      }
      else if ( tag & MG_GEO ) {
        memcpy(pxp->n1,no1,3*sizeof(double));
        memcpy(pxp->n2,no2,3*sizeof(double));
        memcpy(ppt->n,to,3*sizeof(double));
      }
      else if ( tag & MG_REF ) {
        memcpy(pxp->n1,no1,3*sizeof(double));
        memcpy(ppt->n,to,3*sizeof(double));
      }
      else
        memcpy(pxp->n1,no1,3*sizeof(double));
    }
    return 2;
    /* End of case of a bdy face */
  }
  else if(pt->xt){
    /** Tetra has a xtetra but the longest edge do not belong to a bdy face:
     * do nothing to avoid splitting of a bdy edge from a non bdy face (due
     * to collapses, a tetra with no bdy faces may have a xtetra and
     * boundary tags or no tags on boundary edge). */
    return 0;
  } else {
    /** Case of a tetra without xtetra (no boundary faces): split non-bdy
     * edges with Delauney kernel. */
    /* Note that it is possible that non bdy tetra contains a bdy edge, here
     * only non bdy edge are considered */
    ilist = MMG5_coquil(mesh,k,imax,list);
    if ( !ilist ){
      /* Unable to compute edge shell: treat next element */
      return 0;
    }
    else if ( ilist<0 ) {
      return -1;
    }
    else if(ilist%2) {
      /* Edge is bdy: we want to treat it from a bdy face */
      return 1;
    }
    o[0] = 0.5*(p0->c[0] + p1->c[0]);
    o[1] = 0.5*(p0->c[1] + p1->c[1]);
    o[2] = 0.5*(p0->c[2] + p1->c[2]);
#ifdef USE_POINTMAP
    src = mesh->point[ip1].src;
#else
    src = 1;
#endif
    ip = MMG3D_newPt(mesh,o,MG_NOTAG,src);

    if ( !ip )  {
      /* reallocation of point table */
      MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                          *warn=1;++(*countMemFailure);
                          return 1,
                          o,MG_NOTAG,src);
    }
    if ( met && met->m ) {
      if ( MMG5_intmet(mesh,met,k,imax,ip,0.5)<=0 ) {
        MMG3D_delPt(mesh,ip);
        return 1;
      }
    }

    /* Delaunay */
    if ( lmaxtet< MMG3D_THRES_DEL ) {
      lfilt = MMG3D_LFILTS_DEL;
    }
    else {
      lfilt = MMG3D_LFILTL_DEL;
    }

    int ier = 1;
    if ( *PROctree ) {
      ier = MMG3D_PROctreein(mesh,met,*PROctree,ip,lfilt);
    }

    if ( ier == 0 ) {
      /* PROctree allocated and PROctreein refuse the insertion */
      MMG3D_delPt(mesh,ip);
      (*ifilt)++;
      return 1;
    }
    else if ( ier < 0 ) {
      /* PROctree allocated but PROctreein fail due to lack of memory */
      MMG3D_freePROctree ( mesh,PROctree );
      MMG3D_delPt(mesh,ip);
      (*ifilt)++;
      return 1;
    } else {
      int lon = MMG5_cavity(mesh,met,k,ip,list,ilist/2,MMG5_EPSOK);
      if ( lon < 1 ) {
        MMG3D_delPt(mesh,ip);
        return 1;
      } else {
        int ret = MMG5_delone(mesh,met,ip,list,lon);
        if ( ret > 0 ) {
          if ( *PROctree )
            MMG3D_addPROctree(mesh,*PROctree,ip);
          (*ns)++;
          return 2;
        }
        else if ( ret == 0 ) {
          MMG3D_delPt(mesh,ip);
          return 1;
        }
        else {
          /* Allocation problem ==> savemesh */
          MMG3D_delPt(mesh,ip);
          return -2;
        }
      }
    }
  }

  return 3;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param k index of tetra in which we work.
 * \param imin index in \a k of edge that we consider for collapse.
 * \param lmin length of edge \a imin.
 * \param nc pointer toward count of collapses (has to be updated)
 *
 * \return -1 for strong failure.
 *
 * \return 0 if edge cannot be collapsed and if we want to pass to next loop
 * step (next element or next tetra edge)
 *
 * \return 1 if edge has been collapsed.
 *
 * \return 2 if nothing has been done (no error but no collapse either).
 *
 * Try to collapse edge \a imin it too small.
 *
 */
static inline
int MMG3D_mmg3d1_delone_collapse(MMG5_pMesh mesh, MMG5_pSol met,
                                 MMG3D_pPROctree *PROctree,MMG5_int k,
                                 int8_t imin,double lmin,MMG5_int* nc) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1;
  int64_t       list[MMG3D_LMAX+2];
  MMG5_int      lists[MMG3D_LMAX+2],ip1,ip2;
  int           ilist,ilists;
  int8_t        j,i,i1,i2,ifa0,ifa1;

  if(lmin > MMG3D_LOPTS) {
    /* Edge is large enough: nothing to do */
    return 2;
  }

  // Case of an internal tetra with 4 ridges vertices.
  if ( lmin == 0 ) {
    /* Case of an internal tetra with 4 ridges vertices */
//#warning is it possible to merge this edge ??
    return 0;
  }

  pt = &mesh->tetra[k];
  pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

  ifa0 = MMG5_ifar[imin][0];
  ifa1 = MMG5_ifar[imin][1];
  i  =  (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
  j  = MMG5_iarfinv[i][imin];
  i1 = MMG5_idir[i][MMG5_inxt2[j]];
  i2 = MMG5_idir[i][MMG5_iprv2[j]];

  assert( 0<=i1 && i1<4 && "unexpected local index for vertex");
  assert( 0<=i2 && i2<4 && "unexpected local index for vertex");

  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  p0  = &mesh->point[ip1];
  p1  = &mesh->point[ip2];

  /* Ignore OLDPARBDY tag of p0 */
  int16_t tag = p0->tag;
  tag &= ~MG_OLDPARBDY;
  if ( (tag > p1->tag) || (tag & MG_REQ) ) {
    /* Unable to merge edge: pass to next element */
    return 0;
  }

  ilist = 0;
  if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
    /* Case of a boundary face */
    tag = pxt->tag[MMG5_iarf[i][j]];
    if ( tag & MG_REQ ) {
      return 0;
    }
    tag |= MG_BDY;
    if ( p0->tag > tag ) {
      return 0;
    }
    if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) {
      return 0;
    }

    int16_t isnm = (p0->tag & MG_NOM);
    if (MMG5_boulesurfvolp(mesh,k,i1,i, list,&ilist,lists,&ilists,isnm) < 0 ) {
      return -1;
    }

    ilist = MMG5_chkcol_bdy(mesh,met,k,i,j,list,ilist,lists,ilists,0,0,2,0,0);
    if ( ilist > 0 ) {
      int ier = MMG5_colver(mesh,met,list,ilist,i2,2);
      if ( ier < 0 ) {
        return -1;
      }
      else if(ier) {
        MMG3D_delPt(mesh,ier);
        (*nc)++;
        return 1;
      }
    }
    else if (ilist < 0 ) {
      return -1;
    }
  }
  else {
    /* Case of an internal face */
    if ( p0->tag & MG_BDY ) {
      return 0;
    }

    ilist = MMG5_boulevolp(mesh,k,i1,list);
    ilist = MMG5_chkcol_int(mesh,met,k,i,j,list,ilist,2);
    if ( ilist > 0 ) {
      int ier = MMG5_colver(mesh,met,list,ilist,i2,2);
      if ( ilist < 0 ) {
        return 0;
      }

      if ( ier < 0 ) {
        return -1;
      }
      else if(ier) {
        if ( *PROctree ) {
          MMG3D_delPROctree(mesh,*PROctree,ier);
        }
        MMG3D_delPt(mesh,ier);
        (*nc)++;
        return 1;
      }
    }
    else if (ilist < 0 ) {
      return -1;
    }
  }
  return 2;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param ne number of elements.
 * \param ifilt pointer to store the number of vertices filtered by the PROctree.
 * \param ns pointer to store the number of vertices insertions.
 * \param nc pointer to store the number of collapse.
 * \param warn pointer to store a flag that warn the user in case of
 * reallocation difficulty.
 * \return -1 if fail and we don't save the mesh, 0 if fail but we try to save
 * the mesh, 1 otherwise.
 *
 * \a adpsplcol loop: split edges longer than \ref MMG3D_LOPTL_DEL and
 * collapse edges shorter than \ref MMG3D_LOPTS.
 *
 */
static inline int
MMG5_adpsplcol(MMG5_pMesh mesh, MMG5_pSol met,MMG3D_pPROctree *PROctree,
                MMG5_int ne,MMG5_int* ifilt,MMG5_int* ns,MMG5_int* nc,int* warn) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  double        len,lmax;
  MMG5_int      k,base;
  double        lmin;
  double        lmaxtet,lmintet;
  int           ier,imaxtet,imintet;
  int8_t        imin,imax,chkRidTet,countMemFailure;
  static int8_t mmgWarn0 = 0;

  countMemFailure = 0;

  base = ++mesh->mark;

  if ( met->size==6 )  chkRidTet=1;
  else chkRidTet=0;

  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt)  || (pt->tag & MG_REQ) )   continue;
    else if ( pt->mark < base-2 )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /** Step 1: find longest and shortest edge  (and try to manage them) */
    imax = -1; lmax = 0.0;
    imin = -1; lmin = DBL_MAX;
    int ii;
    for (ii=0; ii<6; ii++) {
      if ( pt->xt && (pxt->tag[ii] & MG_REQ) )  continue;
      len = MMG5_lenedg(mesh,met,ii,pt);

      if ( len > lmax ) {
        lmax = len;
        imax = ii;
      }
      if ( len < lmin ) {
        lmin = len;
        imin = ii;
      }
    }
    /* Check that we have found valid edges */
    if ( imax==-1 ) {
      if ( (mesh->info.ddebug || mesh->info.imprim > 5 ) ) {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  # Warning: %s: all edges of tetra %" MMG5_PRId " are"
                  " boundary and required.\n",
                  __func__,k);
        }
      }
      continue;
    }
    if ( imin==-1 ) {
      if ( (mesh->info.ddebug || mesh->info.imprim > 5 ) ) {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  # Warning: %s: all edges of tetra %" MMG5_PRId " are"
                  " boundary and required.\n",
                  __func__,k);
        }
      }
      continue;
    }

    ier = MMG3D_mmg3d1_delone_split(mesh,met,PROctree,k,imax,lmax,lmax,
                                    chkRidTet,ifilt,ns,warn,&countMemFailure);
    if ( ier == -2 ) {
      /* Low failure: try to save mesh and exit lib */
      return 0;
    }
    if ( ier == -1 ) {
      /* Strong failure: exit lib without saving mesh */
      return -1;
    }
    else if ( !ier ) {
      /* Unable to treat largest edge: pass to next element */
      continue;
    }
    else if ( ier == 2 ) {
      /* Edge has been splitted: pass to next element */
      continue;
    }
    assert ( (ier==1 || ier==3) && "Check return val of delone_split");

    if ( countMemFailure > 10 ) {
      printf("  ## Error:%s: too much reallocation errors. Exit program.\n",__func__);
      return -1;
    }

    /** If unable to treat largest edge with ier==1 return value or if edge has
     * not been splitted but slpit_delone has not raised any error: try to
     * collapse shortest edge. */

    /** 2. Try to merge smallest edge: if collapse is not possible, pass to next
     * element */
    ier = MMG3D_mmg3d1_delone_collapse(mesh,met,PROctree,k,imin,lmin,nc);
    if ( ier < 0 ) {
      /* Strong failure */
      return -1;
    }
    else if ( !ier ) {
      /* Unable to treat smallest edge: pass to next element */
      continue;
    }
    else if ( ier == 1 ) {
      /* Smallest edge has been collapsed: pass to next element */
      continue;
    }

    /** Step 2: longest and shortest edges are stucked => try the other edges */
    imaxtet = imax;
    imintet = imin;
    lmaxtet = lmax;
    lmintet = lmin;
    assert(lmin);

    for (ii=0; ii<6; ii++) {
      if ( pt->xt && (pxt->tag[ii] & MG_REQ) )  continue;
      if ( (ii==imintet) && (lmintet < MMG3D_LOPTS)) continue;
      if ( (ii==imaxtet) && (lmaxtet > MMG3D_LOPTL_DEL) ) continue;

      len = MMG5_lenedg(mesh,met,ii,pt);

      imax = ii;
      lmax = len;
      imin = ii;
      lmin = len;

      /** 1. Try to split too long edge */
      ier = MMG3D_mmg3d1_delone_split(mesh,met,PROctree,k,imax,lmax,lmaxtet,
                                      chkRidTet,ifilt,ns,warn,&countMemFailure);
      if ( ier == -2 ) {
        /* Low failure: try to save mesh and exit lib */
        return 0;
      }
      if ( ier == -1 ) {
        /* Strong failure: exit lib without saving mesh */
        return -1;
      }
      else if ( !ier ) {
        /* Unable to treat too large edge: pass to next edge of element */
        continue;
      }
      else if ( ier == 2 ) {
        /* Edge has been splitted: pass to next element */
        break;
      }

      if ( countMemFailure > 10 ) {
        fprintf(stderr,"  ## Error:%s: too much reallocation errors."
                " Exit program.\n",__func__);
        return -1;
      }

      assert ( (ier==1 || ier==3) && "Check return val of delone_split");

      /** If unable to treat too large edge with ier==1 return value or if edge
       * has not been splitted but slpit_delone has not raised any error: try to
       * collapse too short edge. */

      /** 2. Try to merge smallest edge: if collapse is not possible, pass to
       * next element */
      ier = MMG3D_mmg3d1_delone_collapse(mesh,met,PROctree,k,imin,lmin,nc);
      if ( ier < 0 ) {
        /* Strong failure */
        return -1;
      }
      else if ( !ier ) {
        /* Unable to treat too small edge: pass to next edge of element */
        continue;
      }
      else if ( ier == 1 ) {
        /* Edge has been collapsed: pass to next element */
        break;
      }
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Mesh optimization during insertion phase.
 *
 */
static int
MMG5_optbad(MMG5_pMesh mesh, MMG5_pSol met,MMG3D_pPROctree PROctree) {
  int           it,maxit;
  MMG5_int      nf,nnf,nnm,nm,nw;
  double        crit;

  /* shape optim */
  it = nnm = nnf = 0;
  maxit = 3;
  crit = 1.053;

  do {
    /* treatment of bad elements*/
    nw = MMG3D_opttyp(mesh,met,PROctree,-1);
    /* badly shaped process */
    if ( !mesh->info.noswap ) {
      nf = MMG5_swpmsh(mesh,met,PROctree,2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
                __func__);
        return 0;
      }
      nnf += nf;

      nf += MMG5_swptet(mesh,met,crit,MMG3D_SWAP06,PROctree,2,mesh->mark-1);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      /* move for tria with qual<1., tetra with qual<1, internal move
       * allowed, surface degradation forbidden, volume degradation during the
       * surface move forbidden and volume degradation during volumic move
       * forbidden. Perform 1 iter max (0). */
      nm = MMG5_movtet(mesh,met,PROctree,MMG3D_MAXKAL,MMG3D_MAXKAL,1,1,1,1,0,mesh->mark-1);
      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh.\n",
                __func__);
        return 0;
      }
    }
    else  nm = 0;
    nnm += nm;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nw+nf+nm > 0 ){
      fprintf(stdout,"                                          ");
      fprintf(stdout,"  %8" MMG5_PRId " improved, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",nw,nf,nm);
    }
  }
  while( ++it < maxit && nw+nm+nf > 0 );

  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnf > 0 || nnm > 0) )
      fprintf(stdout,"                                                 "
              "        "
              "      %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved, %d iter. \n",nnf,nnm,it);
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param warn set to 1 if we can't insert point due to lack of memory.
 * \return -1 if fail and we dont try to end the remesh process,
 * 0 if fail but we try to end the remesh process and 1 if success.
 *
 * Split edges longer than \ref MMG3D_LOPTL_DEL and collapse edges shorter
 * than \ref MMG3D_LOPTS.
 *
 */
static int
MMG5_adpdel(MMG5_pMesh mesh,MMG5_pSol met,MMG3D_pPROctree *PROctree, int* warn) {
  int        ier;
  int        it,maxit,noptim;
  MMG5_int   ns,nc,ne,nnm,nm,nnf,nf,nnc,nns,nfilt,ifilt;
  double     maxgap,dd,declic,declicsurf;

  /* Iterative mesh modifications */
  it = nnc = nns = nnf = nnm = nfilt = 0;
  noptim = 0;
  maxit = 10;
  mesh->gap = maxgap = 0.5;
  declic = 0.5/MMG3D_ALPHAD;
  declicsurf = 1./3.46;

  do {
    if ( !mesh->info.noinsert ) {
      *warn=0;
      ns = nc = 0;
      ifilt = 0;
      ne = mesh->ne;
      ier = MMG5_adpsplcol(mesh,met,PROctree,ne,&ifilt,&ns,&nc,warn);
      if ( ier<=0 ) return -1;
    } /* End conditional loop on mesh->info.noinsert */
    else  ns = nc = ifilt = 0;

    if ( !mesh->info.noswap ) {
      nf = MMG5_swpmsh(mesh,met,*PROctree,2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
      nnf += nf;
      nf += MMG5_swptet(mesh,met,MMG3D_SSWAPIMPROVE,declic,*PROctree,2,mesh->mark-2);

      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
    } else {
      nf = 0;
    }


    if ( !mesh->info.nomove ) {
      /* move for tria with qual<declicsurf, tetra with qual<declic, internal
       * move allowed, surface degradation forbidden, volume degradation during
       * the surface move authorized and volume degradation during volumic move
       * forbidden. Perform 2 iter max (1). */
      nm = MMG5_movtet(mesh,met,*PROctree,declicsurf,declic,1,1,0,1,1,mesh->mark-2);

      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: Unable to improve mesh.\n",__func__);
        return 0;
      }
    }
    else  nm = 0;

    nnm += nm;
    nnc += nc;
    nns += ns;
    nnf += nf;
    nfilt += ifilt;

    /* decrease size of gap for reallocation */

    if ( mesh->gap > maxgap/(double)maxit )
      mesh->gap -= maxgap/(double)maxit;
    else
      mesh->gap -= mesh->gap/(double)maxit;


    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc+nm+nf > 0)
      fprintf(stdout,"     %8"MMG5_PRId" filtered, %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed,"
              " %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",ifilt,ns,nc,nf,nm);

    /*optimization*/
    dd = MMG5_abs(nc-ns);
    if ( !noptim && (it==5 || ((dd < 5) || (dd < 0.05*MG_MAX(nc,ns)) || !(ns+nc))) ) {
      MMG5_optbad(mesh,met,*PROctree);
      noptim = 1;
    }

    if( it > 5 ) {
      //  if ( ns < 10 && MMG5_abs(nc-ns) < 3 )  break;
      //else if ( it > 3 && MMG5_abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
      dd = MMG5_abs(nc-ns);
      if ( dd < 5 || dd < 0.05*MG_MAX(nc,ns) )   break;
      //else if ( it > 12 && nc >= ns )  break;
    }
  }
  while( ++it < maxit && (noptim || nc+ns > 0) );

  if ( mesh->info.imprim > 0 ) {
    if ( (abs(mesh->info.imprim) < 5) && ( nnc || nns ) ) {
      fprintf(stdout,"     %8"MMG5_PRId" filtered, %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed,"
              " %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved, %d iter.\n",nfilt,nns,nnc,nnf,nnm, it);
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Mesh optimization for LES computation (improve the element skewness).
 *
 */
static int
MMG5_optetLES(MMG5_pMesh mesh, MMG5_pSol met,MMG3D_pPROctree PROctree) {
  int       it,maxit;
  MMG5_int  nnf,nf,nw,nm,nnm;
  double    declic;

  it = nnm = nnf = 0;
  maxit = 10;
  declic = 1.01;
  ++mesh->mark;
  do {
    /* treatment of bad elements*/
    if(it < 5) {
      nw = MMG3D_opttyp(mesh,met,PROctree,mesh->mark-2);
    }
    else {
      nw = 0;
    }

    /* badly shaped process */
    if ( !mesh->info.noswap ) {
      nf = MMG5_swptet(mesh,met,declic,MMG3D_SWAP06,PROctree,2,mesh->mark-2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      /* move for tria with qual<1, tetra with qual<1, internal
       * move allowed, surface degradation forbidden, volume degradation during
       * the surface move forbidden and volume degradation during volumic move
       * forbidden. Perform 4 iter max (3). */
      nm = MMG5_movtet(mesh,met,PROctree,MMG3D_MAXKAL,MMG3D_MAXKAL,1,1,1,1,3,mesh->mark-2);
      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh.\n",
          __func__);
        return 0;
      }
    }
    else  nm = 0;
    nnm += nm;

//be careful, this procedure can degrade the worst elt
    if ( !mesh->info.nomove && (it==2)) {
      MMG3D_optlap(mesh,met);
    }

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nw+nf+nm > 0 ){
      fprintf(stdout,"                                          ");
      fprintf(stdout,"  %8" MMG5_PRId " improved, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",nw,nf,nm);
    }
  }
  while( ++it < maxit && nw+nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    /* move for tria with qual<declicsurf, tetra with qual<declic, internal
     * move allowed, surface degradation forbidden, volume degradation during
       * the surface move authorized and volume degradation during volumic move
       * forbidden. Perform 4 iter max (3). */
    nm = MMG5_movtet(mesh,met,PROctree,MMG3D_MAXKAL,MMG3D_MAXKAL,1,1,1,1,3,mesh->mark-2);
    if ( nm < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to improve mesh.\n",__func__);
      return 0;
    }
  }
  else  nm = 0;
  nnm += nm;
  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 ) {
    fprintf(stdout,"                                            "
            "                                ");
    fprintf(stdout,"     %8" MMG5_PRId " moved\n",nm);
  }


  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnf > 0 || nnm > 0) )
      fprintf(stdout,"                                                 "
              "        "
              "      %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved, %d iter. \n",nnf,nnm,it);
  }
  return 1;
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Mesh optimization using egde swapping and point relocation.
 *
 */
static int
MMG5_optet(MMG5_pMesh mesh, MMG5_pSol met,MMG3D_pPROctree PROctree) {
  MMG5_pTetra   pt;
  int           it,maxit;
  MMG5_int      nnf,nf,nw,k,nnm,nm;
  double        crit,declic;

  /* shape optim */
  it = nnm = nnf = 0;
  maxit  = 10;
  crit   = MMG3D_SSWAPIMPROVE;
  declic = 0.7/MMG3D_ALPHAD;
  /* mark reinitialization in order to see at least one time each tetra */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    pt->mark = mesh->mark;
  }

  do {
  ++mesh->mark;

    /* treatment of bad elements */
    if(it < 5) {
      nw = MMG3D_opttyp(mesh,met,PROctree,mesh->mark-1);
    }
    else
      nw = 0;
    /* badly shaped process */
    if ( !mesh->info.noswap ) {
      nf = MMG5_swpmsh(mesh,met,PROctree,2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
      nnf += nf;

      nf += MMG5_swptet(mesh,met,crit,declic,PROctree,2,mesh->mark-1);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      /* move for tria with qual<1., tetra with qual<declic, internal move
       * allowed, surface degradation forbidden, volume degradation during the
       * surface move forbidden and volume degradation during volumic move
       * authorized. Perform 1 iter max (0). */
      nm = MMG5_movtet(mesh,met,PROctree,MMG3D_MAXKAL,declic,1,1,1,1,0,mesh->mark-1);
      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh.\n",__func__);
        return 0;
      }
    }
    else  nm = 0;
    nnm += nm;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nw+nf+nm > 0 ){
      fprintf(stdout,"                                          ");
      fprintf(stdout,"  %8" MMG5_PRId " improved, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",nw,nf,nm);
    }

    if ( it > 3 ) {
      if ( !nw && (!nm || !nf) )   break;
    }
  }
  while( ++it < maxit && nw+nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    /* move for tria with qual<1., tetra with qual<1., internal move allowed,
     * surface degradation forbidden, volume degradation during the surface and
     * volume move forbidden. Perform 4 iter max. */
    nm = MMG5_movtet(mesh,met,PROctree,MMG3D_MAXKAL,MMG3D_MAXKAL,1,1,1,1,3,mesh->mark-2);
    if ( nm < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: Unable to improve mesh.\n",__func__);
      return 0;
    }
  }
  else  nm = 0;
  nnm += nm;
  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 ) {
    fprintf(stdout,"                                            "
            "                                ");
    fprintf(stdout,"     %8" MMG5_PRId " moved\n",nm);
  }

  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnf > 0 || nnm > 0) )
      fprintf(stdout,"                                                 "
              "        "
              "      %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved, %d iter. \n",nnf,nnm,it);
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param permNodGlob if provided, strore the global permutation of nodes
 * \return 0 if failed, 1 otherwise.
 *
 * Analyze tetrahedra and split long / collapse short, according to
 * prescribed metric.
 *
 */
static int
MMG5_adptet_delone(MMG5_pMesh mesh,MMG5_pSol met,MMG3D_pPROctree *PROctree,
                   MMG5_int * permNodGlob) {
  MMG5_int  nnf,nf;
  int       warn,ns;

  /** Step 1: few iters of swaps */
  if ( !mesh->info.noswap ) {
    nnf = MMG5_swpmsh(mesh,met,*PROctree,2);
    if ( nnf < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
              __func__);
      return 0;
    }
    nf = MMG5_swptet(mesh,met,MMG3D_SSWAPIMPROVE,MMG3D_SWAP06,*PROctree,2,mesh->mark-2);
    if ( nf < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: Unable to improve mesh. Exiting.\n",
              __func__);
      return 0;
    }
    nnf+=nf;
  } else {
    nnf = nf = 0;
  }

  if ( mesh->info.ddebug ) {
    fprintf(stdout," ------------- Delaunay: INITIAL SWAP %7"MMG5_PRId"\n",nnf);
    MMG3D_outqua(mesh,met);
  }

  /** Step 2: few iters of splits, collapses, swaps and moves */
  warn = 0;

  ns = MMG5_adpdel(mesh,met,PROctree,&warn);

  if ( ns < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to complete mesh. Exit program.\n",
      __func__);
    return 0;
  }

  if ( warn ) {
    fprintf(stderr,"\n  ## Error: %s:",__func__);
    fprintf(stderr," unable to allocate a new point in last call of adpspl.\n");
    fprintf(stderr,"  ## Check the mesh size or ");
    fprintf(stderr,"increase the maximal authorized memory with the -m option.\n");
    fprintf(stderr,"  ## Uncomplete mesh. Exiting\n" );
    return 0;
  }

  /* renumerotation if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,permNodGlob) )
    return 0;

  /** Step 3: Last wave of improvements: few iters of bad elts treatment, swaps
   * and moves */
  if(mesh->info.optimLES) {
    if(!MMG5_optetLES(mesh,met,*PROctree)) return 0;
  }
  else {
    if(!MMG5_optet(mesh,met,*PROctree)) return 0;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param permNodGlob if provided, strore the global permutation of nodes
 * \return 0 if failed, 1 if success.
 *
 * Main adaptation routine.
 *
 */
int MMG5_mmg3d1_delone(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *permNodGlob) {
  MMG3D_pPROctree PROctree = NULL;

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  if ( mesh->info.iso && !MMG5_chkmani(mesh) ) {
    fprintf(stderr,"\n  ## Non orientable implicit surface. Exit program.\n");
    return 0;
  }

  /**--- stage 1: geometric mesh  */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !MMG5_anatet(mesh,met,1,0) ) {
    fprintf(stderr,"\n  ## Unable to split mesh. Exiting.\n");
    return 0;
  }

  /* Debug: export variable MMG_SAVE_ANATET1 to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_ANATET1") ) {
    printf("  ## WARNING: EXIT AFTER ANATET-1."
           " (MMG_SAVE_ANATET1 env variable is exported).\n");
    return 1;
  }

  /**--- stage 2: computational mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");


  /* define metric map */
  if ( !MMG3D_defsiz(mesh,met) ) {
    fprintf(stderr,"\n  ## Metric undefined. Exit program.\n");
    return 0;
  }

  /* Debug: export variable MMG_SAVE_DEFSIZ to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_DEFSIZ") ) {
    printf("  ## WARNING: EXIT AFTER DEFSIZ."
           " (MMG_SAVE_DEFSIZ env variable is exported).\n");
    return 1;
  }

  MMG5_gradation_info(mesh);

  if ( mesh->info.hgrad > 0. ) {
    if ( !MMG3D_gradsiz(mesh,met) ) {
      fprintf(stderr,"\n  ## Gradation problem. Exit program.\n");
      return 0;
    }
  }
  if ( mesh->info.hgradreq > 0. ) {
    MMG3D_gradsizreq(mesh,met);
  }

  /* Debug: export variable MMG_SAVE_GRADSIZ to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_GRADSIZ") ) {
    printf("  ## WARNING: EXIT AFTER GRADSIZ."
           " (MMG_SAVE_GRADSIZ env variable is exported).\n");
    return 1;
  }

  /*update quality*/
  if ( !MMG3D_tetraQual(mesh,met,1) ) return 0;

  if ( !MMG5_anatet(mesh,met,2,0) ) {
    fprintf(stderr,"\n  ## Unable to split mesh. Exiting.\n");
    return 0;
  }

  /* Debug: export variable MMG_SAVE_ANATET2 to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_ANATET2") ) {
    printf("  ## WARNING: EXIT AFTER ANATET-2."
           " (MMG_SAVE_ANATET2 env variable is exported).\n");
    return 1;
  }

  /* renumerotation if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,permNodGlob) ) {
    return 0;
  }

  if ( mesh->info.PROctree > 0 ) {
    if ( !MMG3D_initPROctree(mesh,&PROctree,mesh->info.PROctree) ) {
      if ( PROctree ) {
        /*free PROctree*/
        MMG3D_freePROctree(mesh,&PROctree);
      }
    }
  }

  if ( !MMG5_adptet_delone(mesh,met,&PROctree,permNodGlob) ) {
    fprintf(stderr,"\n  ## Unable to adapt. Exit program.\n");
    if ( PROctree ) {
      /*free PROctree*/
      MMG3D_freePROctree(mesh,&PROctree);
    }
    return 0;
  }

  /* in test phase: check if no element with 2 bdry faces */
  if ( !MMG5_chkfemtopo(mesh) ) {
    fprintf(stderr,"\n  ## Topology of mesh unsuited for fem computations. Exit program.\n");
    if ( PROctree ) {
      /*free PROctree*/
      MMG3D_freePROctree(mesh,&PROctree);
    }
    return 0;
  }

  int ier = 1;

  if ( mesh->info.iso && !MMG5_chkmani(mesh) ) {
    fprintf(stderr,"\n  ## Non orientable implicit surface. Exit program.\n");
    ier = 0;
  }

  if ( PROctree ) {
    /*free PROctree*/
    MMG3D_freePROctree(mesh,&PROctree);
  }

  return ier;
}

#endif
