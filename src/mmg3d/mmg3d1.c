/* =============================================================================
**  This file is part of the MMG3D 5 software package for the tetrahedral
**  mesh modification.
**  Copyright (c) 2014 Inria / Université de Bordeaux, IMB / UPMC, LJLL.
**
**  MMG3D 5 is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  MMG3D 5 is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with MMG3D 5 (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the MMG3D 5 distribution only if you accept them.
** =============================================================================
*/

#include "mmg3d.h"

char  ddb;

/** set triangle corresponding to face ie of tetra k */
void tet2tri(pMesh mesh,int k,char ie,Tria *ptt) {
  pTetra  pt;
  pxTetra pxt;
  char    i;

  pt = &mesh->tetra[k];
  memset(ptt,0,sizeof(Tria));
  ptt->v[0] = pt->v[idir[ie][0]];
  ptt->v[1] = pt->v[idir[ie][1]];
  ptt->v[2] = pt->v[idir[ie][2]];
  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    ptt->ref = pxt->ref[ie];
    for (i=0; i<3; i++) {
      ptt->edg[i] = pxt->edg[iarf[ie][i]];
      ptt->tag[i] = pxt->tag[iarf[ie][i]];
    }
  }
  else {
    for (i=0; i<3; i++) {
      ptt->edg[i] = 0;
      ptt->tag[i] = 0;
    }
  }
}

/** find acceptable position for splitting */
int dichoto(pMesh mesh,pSol met,int k,int *vx) {
  pTetra  pt;
  pPoint  pa,pb,ps;
  double  o[6][3],p[6][3];
  float   to,tp,t;
  int     ia,ib,ier,it,maxit;
  char    i;

  pt = &mesh->tetra[k];
  /* get point on surface and along segment for edge split */
  for (i=0; i<6; i++) {
    memset(p[i],0,3*sizeof(double));
    memset(o[i],0,3*sizeof(double));
    if ( vx[i] > 0 ) {
      ia = pt->v[iare[i][0]];
      ib = pt->v[iare[i][1]];
      pa = &mesh->point[ia];
      pb = &mesh->point[ib];
      ps = &mesh->point[vx[i]];
      o[i][0] = 0.5 * (pa->c[0] + pb->c[0]);
      o[i][1] = 0.5 * (pa->c[1] + pb->c[1]);
      o[i][2] = 0.5 * (pa->c[2] + pb->c[2]);
      p[i][0] = ps->c[0];
      p[i][1] = ps->c[1];
      p[i][2] = ps->c[2];
    }
  }
  maxit = 4;
  it = 0;
  tp = 1.0;
  to = 0.0;
  do {
    /* compute new position */
    t = 0.5 * (tp + to);
    for (i=0; i<6; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
        ps->c[2] = o[i][2] + t*(p[i][2] - o[i][2]);
      }
    }
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32:
      ier = split1_sim(mesh,met,k,vx);
      break;
    case 11: case 21: case 38: case 56:
      ier = split3_sim(mesh,met,k,vx);
      break;
    default:
      ier = split2sf_sim(mesh,met,k,vx);
      break;
    }
    if ( ier )
      to = t;
    else
      tp = t;
    /* if we realloc mem in the split function, pt is not valid anymore */
    pt = &mesh->tetra[k];
  }
  while ( ++it < maxit );
  /* restore coords of last valid pos. */
  if ( !ier ) {
    t = to;
    for (i=0; i<6; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
        ps->c[2] = o[i][2] + t*(p[i][2] - o[i][2]);
      }
    }
  }
  return(1);
}

/** Find acceptable position for split1b, passing the shell of considered edge, starting from o */
int dichoto1b(pMesh mesh,int *list,int ret,double o[3],double ro[3]) {
  pTetra  pt;
  pPoint  p0,p1;
  int     iel,np,nq,it,maxit;
  double  m[3],c[3],tp,to,t;
  char    ia,ier;

  iel = list[0] / 6;
  ia  = list[0] % 6;
  pt  = &mesh->tetra[iel];

  np = pt->v[iare[ia][0]];
  nq = pt->v[iare[ia][1]];
  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];

  /* midpoint along edge */
  m[0] = 0.5*(p0->c[0] + p1->c[0]);
  m[1] = 0.5*(p0->c[1] + p1->c[1]);
  m[2] = 0.5*(p0->c[2] + p1->c[2]);

  maxit = 4;
  it    = 0;
  ier   = 0;
  tp    = 1.0;
  to    = 0.0;
  do {
    t = 0.5*(to + tp);
    c[0] = m[0] + t*(o[0]-m[0]);
    c[1] = m[1] + t*(o[1]-m[1]);
    c[2] = m[2] + t*(o[2]-m[2]);

    ier = simbulgept(mesh,list,ret,c);
    if ( ier )
      to = t;
    else
      tp = t;
  }
  while ( ++it < maxit );
  if ( !ier )  t = to;
  ro[0] = m[0] + t*(o[0]-m[0]);
  ro[1] = m[1] + t*(o[1]-m[1]);
  ro[2] = m[2] + t*(o[2]-m[2]);

  return(1);
}

/** return edges of (virtual) triangle pt that need to be split w/r Hausdorff criterion */
char chkedg(pMesh mesh,Tria *pt,char ori) {
  pPoint   p[3];
  xPoint  *pxp;
  double   n[3][3],t[3][3],nt[3],*n1,*n2,t1[3],t2[3];
  double   ps,ps2,ux,uy,uz,ll,il,alpha,dis,hma2,hausd;
  int      ia,ib,ic;
  char     i,i1,i2;
  pPar     par;

  ia   = pt->v[0];
  ib   = pt->v[1];
  ic   = pt->v[2];
  p[0] = &mesh->point[ia];
  p[1] = &mesh->point[ib];
  p[2] = &mesh->point[ic];
  pt->flag = 0;
  hma2 = LLONG*LLONG*mesh->info.hmax*mesh->info.hmax;
  
 /* for (i=0; i<3; i++) { */
 /*   for (i1=0; i1<3; i1++)  */
 /*     t[i][i1] = 10000000; */
 /* } */
  /* normal recovery */
  for (i=0; i<3; i++) {
    if ( MG_SIN(p[i]->tag) ) {
      nortri(mesh,pt,n[i]);
      if(!ori) {
        n[i][0] *= -1.0;
        n[i][1] *= -1.0;
        n[i][2] *= -1.0;
      }
    }
    else if (p[i]->tag & MG_NOM){
      nortri(mesh,pt,n[i]);
      if(!ori) {
        n[i][0] *= -1.0;
        n[i][1] *= -1.0;
        n[i][2] *= -1.0;
      }
      assert(p[i]->xp);
      pxp = &mesh->xpoint[p[i]->xp];
      memcpy(&t[i],pxp->t,3*sizeof(double));
    }
    else {
      assert(p[i]->xp);
      pxp = &mesh->xpoint[p[i]->xp];
      if ( MG_EDG(p[i]->tag) ) {
        memcpy(&t[i],pxp->t,3*sizeof(double));
        nortri(mesh,pt,nt);
        if(!ori) {
          nt[0] *= -1.0;
          nt[1] *= -1.0;
          nt[2] *= -1.0;
        }
        ps  = pxp->n1[0]*nt[0] + pxp->n1[1]*nt[1] + pxp->n1[2]*nt[2];
        ps2 = pxp->n2[0]*nt[0] + pxp->n2[1]*nt[1] + pxp->n2[2]*nt[2];
        if ( fabs(ps) > fabs(ps2) )
          memcpy(&n[i],pxp->n1,3*sizeof(double));
        else
          memcpy(&n[i],pxp->n2,3*sizeof(double));
      }
      else
        memcpy(&n[i],pxp->n1,3*sizeof(double));
    }
  }

  /* local hausdorff for triangle */
  hausd = mesh->info.hausd;
  for (i=0; i<mesh->info.npar; i++) {
    par = &mesh->info.par[i];
    if ( (par->elt == MMG5_Triangle) && (pt->ref == par->ref ) )
      hausd = par->hausd;
  }

  /* analyze edges */
  for (i=0; i<3; i++) {
    i1 = inxt2[i];
    i2 = iprv2[i];

    /* check length */
    ux = p[i2]->c[0] - p[i1]->c[0];
    uy = p[i2]->c[1] - p[i1]->c[1];
    uz = p[i2]->c[2] - p[i1]->c[2];
    ll = ux*ux + uy*uy + uz*uz;
    if ( ll < EPSD )  continue;
    else if ( ll > hma2 ) {
      MG_SET(pt->flag,i);
      continue;
    }
    il = 1.0 / sqrt(ll);

    /* Hausdorff w/r tangent direction */
    if ( MG_EDG(pt->tag[i]) || ( pt->tag[i] & MG_NOM )) {
      if ( MG_SIN(p[i1]->tag) ) {
        t1[0] = il * ux;
        t1[1] = il * uy;
        t1[2] = il * uz;
      }
      else {
	if(!((p[i1]->tag & MG_NOM) ||  MG_EDG(p[i1]->tag) ) ) {
	  // 	if(t[i1][0] > 10) {
	  fprintf(stdout,"warning geometrical problem %d -- %d %d -- %e\n",p[i1]->tag,
		  MG_SIN(p[i1]->tag ),p[i1]->tag & MG_NOM,t[i1][0]);
	    return(0);
	  }
        memcpy(t1,t[i1],3*sizeof(double));
        ps = t1[0]*ux + t1[1]*uy + t1[2]*uz;
        if ( ps < 0.0 ) {
          t1[0] *= -1.0;
          t1[1] *= -1.0;
          t1[2] *= -1.0;
        }
      }
      if ( MG_SIN(p[i2]->tag) ) {
        t2[0] = -il * ux;
        t2[1] = -il * uy;
        t2[2] = -il * uz;
      }
      else {
  	if(!((p[i2]->tag & MG_NOM) || MG_EDG(p[i2]->tag) ) ) {
	    fprintf(stdout,"2. warning geometrical problem\n");
	    return(0);
	  }
        memcpy(t2,t[i2],3*sizeof(double));
        ps = - ( t2[0]*ux + t2[1]*uy + t2[2]*uz );
        if ( ps < 0.0 ) {
          t2[0] *= -1.0;
          t2[1] *= -1.0;
          t2[2] *= -1.0;
        }
      }
    }
    else {
      n1 = n[i1];
      n2 = n[i2];
      if ( !BezierTgt(p[i1]->c,p[i2]->c,n1,n2,t1,t2) ) {
        t1[0] = ux * il;
        t1[1] = uy * il;
        t1[2] = uz * il;

        t2[0] = -ux * il;
        t2[1] = -uy * il;
        t2[2] = -uz * il;
      }
    }
    alpha = BezierGeod(p[i1]->c,p[i2]->c,t1,t2);
    ps  = t1[0]*ux + t1[1]*uy + t1[2]*uz;
    ps *= il;
    dis = alpha*alpha*fabs(1.0 - ps*ps);
    if ( dis > hausd*hausd ) {
      MG_SET(pt->flag,i);
      continue;
    }
    ps  = -( t2[0]*ux + t2[1]*uy + t2[2]*uz );
    ps *= il;
    dis = alpha*alpha*fabs(1.0 - ps*ps);

    if ( dis > hausd*hausd ) {
      MG_SET(pt->flag,i);
      continue;
    }
  }
  return(pt->flag);
}

/** Search for boundary edges that could be swapped for geometric approximation */
/*static*/ int swpmsh(pMesh mesh,pSol met
#ifndef PATTERN
		      ,pBucket bucket
#endif
		      ) {
  pTetra   pt;
  pxTetra  pxt;
  int      k,it,list[LMAX+2],ilist,ret,it1,it2,ns,nns,maxit;
  char     i,j,ia,ier;

  it = nns = 0;
  maxit = 2;
  do {
    ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( (!MG_EOK(pt)) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
      else if ( !pt->xt ) continue;
      pxt = &mesh->xtetra[pt->xt];

      for (i=0; i<4; i++) {
        ier = 0;
        if ( !(pxt->ftag[i] & MG_BDY) ) continue;
        for (j=0; j<3; j++) {
          ia  = iarf[i][j];
          if ( (pxt->tag[ia] & MG_REQ) ) continue;
          ret = coquilface(mesh,k,ia,list,&it1,&it2);
          ilist = ret / 2;
          if ( ret < 0 )  return(-1);
          /* CAUTION: trigger collapse with 2 elements */
          if ( ilist <= 1 )  continue;
          ier = chkswpbdy(mesh,list,ilist,it1,it2);
          if ( ier ) {
#ifdef PATTERN
            ier = swpbdy(mesh,met,list,ret,it1);
#else
            ier = swpbdy(mesh,met,list,ret,it1,bucket);
#endif
            if ( ier > 0 )  ns++;
            else if ( ier < 0 )  return(-1);
            break;
          }
        }
        if ( ier )  break;
      }
    }
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8d edge swapped\n",nns);

  return(nns);
}

/** Internal edge flipping */
static int swptet(pMesh mesh,pSol met,double crit) {
  pTetra   pt;
  pxTetra  pxt;
  int      list[LMAX+2],ilist,k,it,nconf,maxit,ns,nns,ier;
  char     i;

  maxit = 2;
  it = nns = 0;

  do {
    ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
      if ( pt->qual > 0.0288675 /*0.6/ALPHAD*/ )  continue;

      for (i=0; i<6; i++) {
        /* Prevent swap of a ref or tagged edge */
        if ( pt->xt ) {
          pxt = &mesh->xtetra[pt->xt];
          if ( pxt->edg[i] || pxt->tag[i] ) continue;
        }

        nconf = chkswpgen(mesh,k,i,&ilist,list,crit);
        if ( nconf ) {
#ifdef PATTERN
          ier = swpgen(mesh,met,nconf,ilist,list);
#else
          ier = swpgen(mesh,met,nconf,ilist,list,NULL);
#endif
          if ( ier > 0 )  ns++;
          else if ( ier < 0 ) return(-1);
          break;
        }
      }
    }
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8d edge swapped\n",nns);

  return(nns);
}

/** Analyze tetrahedra and move points so as to make mesh more uniform */
static int movtet(pMesh mesh,pSol met,int maxit) {
  pTetra        pt;
  pPoint        ppt;
  pxTetra       pxt;
  double        *n;
  int           i,k,ier,nm,nnm,ns,lists[LMAX+2],listv[LMAX+2],ilists,ilistv,it;
  int           improve;
  unsigned char j,i0,base;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** OPTIMIZING MESH\n");

  base = 1;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = base;

  it = nnm = 0;
  do {
    base++;
    nm = ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;

      /* point j on face i */
      for (i=0; i<4; i++) {
        for (j=0; j<3; j++) {
          if ( pt->xt ) {
            pxt = &mesh->xtetra[pt->xt];
            if ( pxt->tag[iarf[i][j]] & MG_REQ )  continue;
          }
          else  pxt = 0;
          i0  = idir[i][j];
          ppt = &mesh->point[pt->v[i0]];
          if ( ppt->flag == base )  continue;
          else if ( MG_SIN(ppt->tag) )  continue;
#ifdef SINGUL
          else if ( ppt->tag & MG_SGL )  continue;
          else if ( mesh->info.sing && pt->xt && (pxt->tag[iarf[i][j]] & MG_SGL) )
            continue;
#endif
          if ( maxit != 1 ) {
            ppt->flag = base;
            improve   = 1;
          }
          else {
            improve = 0;
          }
          ier = 0;
          if ( ppt->tag & MG_BDY ) {
            /* Catch a boundary point by a boundary face */
            if ( !pt->xt || !(MG_BDY & pxt->ftag[i]) )  continue;
            else if( ppt->tag & MG_NOM ){
              if( mesh->adja[4*(k-1)+1+i] ) continue;
              if( !(ier=bouleext(mesh,k,i0,i,listv,&ilistv,lists,&ilists)) )  continue;
              else if ( ier>0 )
                ier = movbdynompt(mesh,listv,ilistv,lists,ilists);
              else
                return(-1);
            }
            else if ( ppt->tag & MG_GEO ) {
              if ( !(ier=boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists)) )
                continue;
              else if ( ier>0 )
                ier = movbdyridpt(mesh,listv,ilistv,lists,ilists);
              else
                return(-1);
            }
            else if ( ppt->tag & MG_REF ) {
              if ( !(ier=boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists)) )
                continue;
              else if ( ier>0 )
                ier = movbdyrefpt(mesh,listv,ilistv,lists,ilists);
              else
                return(-1);
            }
            else {
              if ( !(ier=boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists)) )
                continue;
              else if ( ier<0 )
                return(-1);

              n = &(mesh->xpoint[ppt->xp].n1[0]);
              if ( !directsurfball(mesh, pt->v[i0],lists,ilists,n) )  continue;
              ier = movbdyregpt(mesh,listv,ilistv,lists,ilists);
              if ( ier )  ns++;
            }
          }
          else {
            ilistv = boulevolp(mesh,k,i0,listv);
            if ( !ilistv )  continue;
            ier = movintpt(mesh,listv,ilistv,improve);
          }
          if ( ier ) {
            nm++;
            if(maxit==1){
              ppt->flag = base;
            }
          }
        }
      }
    }
    nnm += nm;
    if ( mesh->info.ddebug )  fprintf(stdout,"     %8d moved, %d geometry\n",nm,ns);
  }
  while( ++it < maxit && nm > 0 );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nnm )
    fprintf(stdout,"     %8d vertices moved, %d iter.\n",nnm,it);

  return(nnm);
}

/** attempt to collapse small edges */
/*static*/ int coltet(pMesh mesh,pSol met,char typchk) {
  pTetra     pt;
  pxTetra    pxt;
  pPoint     p0,p1;
  double     ll,ux,uy,uz,hmi2;
  int        k,nc,list[LMAX+2],ilist,base,nnm;
  char       i,j,tag,ip,iq,isnm;
  int        ier;

  nc = nnm = 0;
  hmi2 = mesh->info.hmin*mesh->info.hmin;

  /* init of point flags, otherwise it can be uninitialized */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  for (k=1; k<=mesh->ne; k++) {
    base = ++mesh->base;
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )   continue;

    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;
    for (i=0; i<4; i++) {
      ier = 0;
      for (j=0; j<3; j++) {
        if ( pt->xt && (pxt->tag[iarf[i][j]] & MG_REQ) )  continue;
        ip = idir[i][inxt2[j]];
        iq = idir[i][iprv2[j]];

        p0 = &mesh->point[pt->v[ip]];
        p1 = &mesh->point[pt->v[iq]];
        if ( p0->flag == base )  continue;
        else if ( (p0->tag & MG_REQ) || (p0->tag > p1->tag) )  continue;
#ifdef SINGUL
        else if ( mesh->info.sing && (p0->tag & MG_SGL) )  continue;
#endif

        /* check length */
        if ( typchk == 1 ) {
          ux = p1->c[0] - p0->c[0];
          uy = p1->c[1] - p0->c[1];
          uz = p1->c[2] - p0->c[2];
          ll = ux*ux + uy*uy + uz*uz;
          if ( ll > hmi2 )  continue;
        }
        else if ( typchk == 2 ) {
          ll = lenedg(mesh,met,pt->v[ip],pt->v[iq]);
          if ( ll > LSHRT )  continue;
        }

        /* boundary face: collapse ip on iq */
        if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
          tag = pxt->tag[iarf[i][j]];
          tag |= MG_BDY;

          isnm = ( tag & MG_NOM );
          if ( isnm ) {
            if ( mesh->adja[4*(k-1)+1+i] )  continue;
          }
          if ( (tag & MG_REQ) || p0->tag > tag )  continue;
          ilist = chkcol_bdy(mesh,k,i,j,list);
        }
        /* internal face */
        else {
          isnm = 0;
          if ( p0->tag & MG_BDY )  continue;
          ilist = chkcol_int(mesh,met,k,i,j,list,typchk);
        }

        if ( ilist > 0 ) {
          ier = colver(mesh,list,ilist,iq);
          if ( ier < 0 ) return(-1);
          else if ( ier ) {
            delPt(mesh,ier);
            break;
          }
        }
        else if (ilist < 0 ) return(-1);
      }
      if ( ier ) {
        p1->flag = base;
        if ( isnm )  nnm++;
        nc++;
        break;
      }
    }
  }
  if ( nc > 0 && (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) )
    fprintf(stdout,"     %8d vertices removed, %8d non manifold,\n",nc,nnm);

  return(nc);
}

/** analyze volume tetra and split if needed */
static int anatetv(pMesh mesh,pSol met,char typchk) {
  pTetra   pt;
  pPoint   p1,p2;
  xTetra  *pxt;
  Hash     hash;
  double   ll,o[3],ux,uy,uz,hma2;
  int      vx[6],k,ip,ip1,ip2,nap,ns,ne,memlack;
  char     i,j,ia;

  /** 1. analysis */
  if ( !hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return(-1);
  memlack = ns = nap = 0;
  hma2 = LLONG*LLONG*mesh->info.hmax*mesh->info.hmax;

  /* Hash all boundary and required edges, and put ip = -1 in hash structure */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

   /* avoid split of edges belonging to a required tet */
    if ( pt->tag & MG_REQ ) {
      for (i=0; i<6; i++) {
        ip1 = pt->v[iare[i][0]];
        ip2 = pt->v[iare[i][1]];
        ip  = -1;
        if ( !hashEdge(mesh,&hash,ip1,ip2,ip) )  return(-1);
      }
      continue;
    }

    if ( !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<4; i++) {
      if ( pxt->ftag[i] & MG_BDY ) {
        for (j=0; j<3; j++) {
          ip1 = pt->v[idir[i][inxt2[j]]];
          ip2 = pt->v[idir[i][iprv2[j]]];
          ip  = -1;
          if ( !hashEdge(mesh,&hash,ip1,ip2,ip) )  return(-1);
        }
        break;
      }
    }
  }

  /** 2. Set flags and split internal edges */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    pt->flag = 0;
    for (i=0; i<6; i++) {
      ip  = -1;
      ip1 = pt->v[iare[i][0]];
      ip2 = pt->v[iare[i][1]];
      p1  = &mesh->point[ip1];
      p2  = &mesh->point[ip2];
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        if ( pxt->tag[i] & MG_REQ ) continue;
      }
      else  pxt = 0;
      if ( (p1->tag & MG_BDY) && (p2->tag & MG_BDY) ) {
        ip = hashGet(&hash,ip1,ip2);
      }
      else {
        if (typchk == 1) {
          ux = p2->c[0] - p1->c[0];
          uy = p2->c[1] - p1->c[1];
          uz = p2->c[2] - p1->c[2];
          ll = ux*ux + uy*uy + uz*uz;
          if ( ll > hma2 )
            ip = hashGet(&hash,ip1,ip2);
        }
        else if ( typchk == 2 ) {
          ll = lenedg(mesh,met,ip1,ip2);
          if ( ll > LLONG )
            ip = hashGet(&hash,ip1,ip2);
        }
      }
      if ( ip < 0 ) continue;
      else if ( !ip ) {
        /* new midpoint */
        o[0] = 0.5 * (p1->c[0]+p2->c[0]);
        o[1] = 0.5 * (p1->c[1]+p2->c[1]);
        o[2] = 0.5 * (p1->c[2]+p2->c[2]);
#ifdef SINGUL
        if ( mesh->info.sing && pt->xt && (pxt->tag[i] & MG_SGL) )
          ip = newPt(mesh,o,MG_SGL);
        else
#endif
          ip  = newPt(mesh,o,0);
        if ( !ip ) {
          /* reallocation of point table */
#ifdef SINGUL
          if ( mesh->info.sing && pt->xt && (pxt->tag[i] & MG_SGL) )
            POINT_REALLOC(mesh,met,ip,mesh->gap,
                          printf("  ## Error: unable to allocate a new point\n");
                          INCREASE_MEM_MESSAGE();
                          memlack=1;
                          goto split
                          ,o,MG_SGL);
          else
#endif
            POINT_REALLOC(mesh,met,ip,mesh->gap,
                          printf("  ## Error: unable to allocate a new point\n");
                          INCREASE_MEM_MESSAGE();
                          memlack=1;
                          goto split
                          ,o,0);
          p1  = &mesh->point[ip1];
          p2  = &mesh->point[ip2];
        }

        if ( met->m )
          met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
        if ( !hashEdge(mesh,&hash,ip1,ip2,ip) )  return(-1);
        MG_SET(pt->flag,i);
        nap++;
      }
#ifdef SINGUL
      /* check that we create a point tag MG_SGL but not MG_BDY */
      if ( mesh->info.sing && pt->xt && (pxt->tag[i] & MG_SGL) ) {
        assert(mesh->point[ip].tag & MG_SGL);
        assert( !(mesh->point[ip].tag & MG_BDY) );
      }
#endif
    }
  }
  if ( !nap )  {
    DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
    return(0);
  }

  /** 3. check and split */
 split:
  ns = 0;
  ne = mesh->ne;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    memset(vx,0,6*sizeof(int));
    pt->flag = 0;
    for (ia=0,i=0; i<3; i++) {
      for (j=i+1; j<4; j++,ia++) {
        if ( pt->xt && (mesh->xtetra[pt->xt].tag[ia] & MG_REQ) ) continue;
        vx[ia] = hashGet(&hash,pt->v[i],pt->v[j]);
        if ( vx[ia] > 0 )  MG_SET(pt->flag,ia);
      }
    }

    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
      split1(mesh,met,k,vx);
      ns++;
      break;
    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      split2sf(mesh,met,k,vx);
      ns++;
      break;

    case 33: case 18: case 12: /* 2 opposite edges split */
      split2(mesh,met,k,vx);
      ns++;
      break;

    case 11: case 21: case 38: case 56: /* 3 edges on the same faces splitted */
      split3(mesh,met,k,vx);
      ns++;
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      split3cone(mesh,met,k,vx);
      ns++;
      break;

    case 35: case 19: case 13: case 37: case 22: case 28: case 26:
    case 14: case 49: case 50: case 44: case 41: /* 3 edges on opposite configuration splitted */
      split3op(mesh,met,k,vx);
      ns++;
      break;

    case 23: case 29: case 53: case 60: case 57: case 58:
    case 27: case 15: case 43: case 39: case 54: case 46: /* 4 edges with 3 lying on the same face splitted */
      split4sf(mesh,met,k,vx);
      ns++;
      break;

      /* 4 edges with no 3 lying on the same face splitted */
    case 30: case 45: case 51:
      split4op(mesh,met,k,vx);
      ns++;
      break;

    case 62: case 61: case 59: case 55: case 47: case 31: /* 5 edges split */
      split5(mesh,met,k,vx);
      ns++;
      break;

    case 63: /* 6 edges split */
      split6(mesh,met,k,vx);
      ns++;
      break;
    }
  }

  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",nap);

  DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
  if ( memlack )  return(-1);
  return(nap);
}

/** analyze tetra and split on geometric criterion */
/*static*/ int anatets(pMesh mesh,pSol met,char typchk) {
  pTetra   pt;
  pPoint   ppt,p1,p2;
  Tria     ptt;
  xTetra  *pxt;
  xPoint  *pxp;
  Bezier   pb;
  Hash     hash;
  double   o[3],no[3],to[3],dd,len;
  int      vx[6],k,ip,ic,it,nap,nc,ni,ne,npinit,ns,ip1,ip2,ier;
  char     i,j,ia,i1,i2;
  static double uv[3][2] = { {0.5,0.5}, {0.,0.5}, {0.5,0.} };

  /** 1. analysis of boundary elements */
  if ( !hashNew(mesh,&hash,mesh->np,7*mesh->np) ) return(-1);
  ns = nap = 0;
  npinit=mesh->np;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) || !pt->xt )  continue;

    /* check boundary face cut w/r Hausdorff or hmax */
    pt->flag = 0;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++){
      if ( pxt->ftag[i] & MG_REQ )  continue;
      if ( pxt->ftag[i] & MG_BDY )  break;
    }
    if ( i == 4 )  continue;

    /* virtual triangle */
    tet2tri(mesh,k,i,&ptt);
    if ( typchk == 1 ) {
      if ( !MG_GET(pxt->ori,i) ) continue;
      if ( !chkedg(mesh,&ptt,MG_GET(pxt->ori,i)) )  continue;
      /* put back flag on tetra */
      for (j=0; j<3; j++){
        if ( pxt->tag[iarf[i][j]] & MG_REQ )  continue;
        if ( MG_GET(ptt.flag,j) )  MG_SET(pt->flag,iarf[i][j]);
      }
    }
    else if ( typchk == 2 ) {
      for (j=0; j<3; j++) {
        ia = iarf[i][j];
        if ( pxt->tag[ia] & MG_REQ )  continue;
        i1  = iare[ia][0];
        i2  = iare[ia][1];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        len = lenedg(mesh,met,ip1,ip2);
        if ( len > LLONG )  MG_SET(pt->flag,ia);
      }
    }
    if ( !pt->flag )  continue;
    ns++;

    /* geometric support */
    ier = bezierCP(mesh,&ptt,&pb,MG_GET(pxt->ori,i));
    assert(ier);

    /* scan edges in face to split */
    for (j=0; j<3; j++) {
      ia = iarf[i][j];
      if ( !MG_GET(pt->flag,ia) )  continue;
      if ( pxt->tag[ia] & MG_REQ )  continue;
      i1  = iare[ia][0];
      i2  = iare[ia][1];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      ip  = hashGet(&hash,ip1,ip2);

      if ( ip > 0 && !(ptt.tag[j] & MG_GEO) )  continue;

      ier = bezierInt(&pb,&uv[j][0],o,no,to);
      /* new point along edge */
      if ( !ip ) {
        ip = newPt(mesh,o,MG_BDY);
        if ( !ip ) {
          /* reallocation of point table */
          POINT_REALLOC(mesh,met,ip,mesh->gap,
                        printf("  ## Error: unable to allocate a new point.\n");
                        INCREASE_MEM_MESSAGE();
                        do {
                          delPt(mesh,mesh->np);
                        } while ( mesh->np>npinit );
                        return(-1)
                        ,o,MG_BDY);
        }
        if ( !hashEdge(mesh,&hash,ip1,ip2,ip) )  return(-1);
        ppt = &mesh->point[ip];
        p1  = &mesh->point[ip1];
        p2  = &mesh->point[ip2];
        if ( met->m )
          met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
        if ( MG_EDG(ptt.tag[j]) || (ptt.tag[j] & MG_NOM) )
          ppt->ref = ptt.edg[j] ? ptt.edg[j] : ptt.ref;
        else
          ppt->ref = ptt.ref;
        ppt->tag |= ptt.tag[j];
        pxp = &mesh->xpoint[ppt->xp];
        memcpy(pxp->n1,no,3*sizeof(double));
        memcpy(pxp->t,to,3*sizeof(double));
        nap++;
      }
      else if ( MG_EDG(ptt.tag[j]) && !(ptt.tag[j] & MG_NOM) ) {
        ppt = &mesh->point[ip];
        assert(ppt->xp);
        pxp = &mesh->xpoint[ppt->xp];

        dd = no[0]*pxp->n1[0]+no[1]*pxp->n1[1]+no[2]*pxp->n1[2];
        if ( dd > 1.0-EPS ) continue;

        memcpy(pxp->n2,no,3*sizeof(double));
        /* a computation of the tangent with respect to these two normals is possible */
        pxp->t[0] = pxp->n1[1]*pxp->n2[2] - pxp->n1[2]*pxp->n2[1];
        pxp->t[1] = pxp->n1[2]*pxp->n2[0] - pxp->n1[0]*pxp->n2[2];
        pxp->t[2] = pxp->n1[0]*pxp->n2[1] - pxp->n1[1]*pxp->n2[0];
        dd = pxp->t[0]*pxp->t[0] + pxp->t[1]*pxp->t[1] + pxp->t[2]*pxp->t[2];
        if ( dd > EPSD2 ) {
          dd = 1.0 / sqrt(dd);
          pxp->t[0] *= dd;
          pxp->t[1] *= dd;
          pxp->t[2] *= dd;
        }
      }
    }
  }
  if ( !ns ) {
    DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
    return(ns);
  }

  /** 2. check if split by adjacent; besides, a triangle may have been splitted and not its adjacent
      (thus, the associated n2 may not exist) : update this normal if need be */
  nc = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /* update face-edge flag */
    for (i=0; i<4; i++) {
      /* virtual triangle */
      memset(&ptt,0,sizeof(Tria));
      if ( pt->xt && pxt->ftag[i] )
        tet2tri(mesh,k,i,&ptt);

      for (j=0; j<3; j++) {
        ia  = iarf[i][j];
        if ( MG_GET(pt->flag,ia) )                continue;
        if ( pt->xt && (pxt->tag[ia] & MG_REQ) )  continue;
        else if ( ptt.tag[j] & MG_REQ )           continue;
        ip1 = pt->v[iare[ia][0]];
        ip2 = pt->v[iare[ia][1]];
        ip  = hashGet(&hash,ip1,ip2);
        if ( ip > 0 ) {
#ifdef SINGUL
          if ( mesh->info.sing && pt->xt && (pxt->tag[ia] & MG_SGL) ) {
            assert(mesh->point[ip].tag & MG_SGL);
            assert( !(mesh->point[ip].tag & MG_BDY) );
          }
#endif
          MG_SET(pt->flag,ia);
          nc++;
          /* ridge on a boundary face */
          if ( !(ptt.tag[j] & MG_GEO) && !(ptt.tag[j] & MG_NOM) )  continue;
          ppt = &mesh->point[ip];
          assert(ppt->xp);
          pxp = &mesh->xpoint[ppt->xp];
          if ( pt->xt )  ier = bezierCP(mesh,&ptt,&pb,MG_GET(pxt->ori,i));
          else  ier = bezierCP(mesh,&ptt,&pb,1);
          ier = bezierInt(&pb,&uv[j][0],o,no,to);

          dd = no[0]*pxp->n1[0]+no[1]*pxp->n1[1]+no[2]*pxp->n1[2];
          if ( dd > 1.0-EPS ) continue;

          memcpy(pxp->n2,no,3*sizeof(double));
        }
      }
    }
  }
  if ( mesh->info.ddebug && nc ) {
    fprintf(stdout,"     %d added\n",nc);
    fflush(stdout);
  }

  /** 3. Simulate splitting and delete points leading to invalid configurations */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  it = 1;
  do {
    ni = nc = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || (pt->tag & MG_REQ) || !pt->flag )  continue;
      memset(vx,0,6*sizeof(int));
      pt->flag = ic = 0;
      for (ia=0,i=0; i<3; i++) {
        for (j=i+1; j<4; j++,ia++) {
          if ( pt->xt && (mesh->xtetra[pt->xt].tag[ia] & MG_REQ) )  continue;
          vx[ia] = hashGet(&hash,pt->v[i],pt->v[j]);
          if ( vx[ia] > 0 ) {
            MG_SET(pt->flag,ia);
            if ( mesh->point[vx[ia]].flag > 2 )  ic = 1;
          }
        }
      }
      if ( !pt->flag )  continue;
      switch (pt->flag) {
      case 1: case 2: case 4: case 8: case 16: case 32:
        ier = split1_sim(mesh,met,k,vx);
        break;
      case 11: case 21: case 38: case 56:
        ier = split3_sim(mesh,met,k,vx);
        break;
      default:
        ier = split2sf_sim(mesh,met,k,vx);
        break;
      }
      if ( ier )  continue;

      nc++;
      if ( ic == 0 && dichoto(mesh,met,k,vx) ) {
        for (ia=0; ia<6; ia++)
          if ( vx[ia] > 0 )  mesh->point[vx[ia]].flag++;
      }
      else {
        for (ia=0,i=0; i<3; i++) {
          for (j=i+1; j<4; j++,ia++) {
            if ( vx[ia] > 0 ) {
              p1 = &mesh->point[pt->v[iare[ia][0]]];
              p2 = &mesh->point[pt->v[iare[ia][1]]];
              ppt = &mesh->point[vx[ia]];
              ppt->c[0] = 0.5 * (p1->c[0] + p2->c[0]);
              ppt->c[1] = 0.5 * (p1->c[1] + p2->c[1]);
              ppt->c[2] = 0.5 * (p1->c[2] + p2->c[2]);
            }
          }
        }
      }
    }
  }
  while( nc > 0 && ++it < 20 );
  if ( mesh->info.ddebug && nc ) {
    fprintf(stdout,"     %d corrected, %d invalid\n",nc,ni);
    fflush(stdout);
  }

  /** 4. splitting */
  ns = 0;
  ne = mesh->ne;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || !pt->flag || (pt->tag & MG_REQ) )  continue;
    memset(vx,0,6*sizeof(int));
    for (ia=0,i=0; i<3; i++) {
      for (j=i+1; j<4; j++,ia++) {
        if ( MG_GET(pt->flag,ia) )  {
          vx[ia] = hashGet(&hash,pt->v[i],pt->v[j]);
          assert(vx[ia]);
        }
      }
    }
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32:  /* 1 edge split */
      split1(mesh,met,k,vx);
      ns++;
      break;
    case 11: case 21: case 38: case 56: /* 1 face (3 edges) subdivided */
      split3(mesh,met,k,vx);
      ns++;
      break;
    default:
      split2sf(mesh,met,k,vx);
      ns++;
      break;
    }
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"       %7d elements splitted\n",nap);

  DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
  return(nap);
}

/** Split edges of length bigger than LOPTL */
static int adpspl(pMesh mesh,pSol met, int* warn) {
  pTetra     pt;
  pxTetra    pxt;
  Tria       ptt;
  pPoint     p0,p1,ppt;
  pxPoint    pxp;
  double     dd,len,lmax,o[3],to[3],ro[3],no1[3],no2[3],v[3];
  int        k,ip,ip1,ip2,list[LMAX+2],ilist,ns,ref,ier;
  char       imax,tag,j,i,i1,i2,ifa0,ifa1;

  *warn=0;
  ns = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )   continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /* find longest edge */
    imax = -1; lmax = 0.0;
    for (i=0; i<6; i++) {
      if ( pt->xt && (pxt->tag[i] & MG_REQ) )  continue;
      ip1  = iare[i][0];
      ip2  = iare[i][1];
      len = lenedg(mesh,met,pt->v[ip1],pt->v[ip2]);
      if ( len > lmax ) {
        lmax = len;
        imax = i;
      }
    }
    if ( imax==-1 )
      fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are required or of length null.\n",
              __FILE__,__LINE__,k);
    if ( lmax < LOPTL )  continue;

    /* proceed edges according to lengths */
    ifa0 = ifar[imax][0];
    ifa1 = ifar[imax][1];
    i  = (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
    j  = iarfinv[i][imax];
    i1 = idir[i][inxt2[j]];
    i2 = idir[i][iprv2[j]];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];
    p0  = &mesh->point[ip1];
    p1  = &mesh->point[ip2];

    /* Case of a boundary face */
    if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
      if ( !(MG_GET(pxt->ori,i)) ) continue;
      ref = pxt->edg[iarf[i][j]];
      tag = pxt->tag[iarf[i][j]];
      if ( tag & MG_REQ )  continue;
      tag |= MG_BDY;
      ilist = coquil(mesh,k,imax,list);
      if ( !ilist )  continue;
      else if ( ilist < 0 )
        return(-1);
      if ( tag & MG_NOM ){
        if( !BezierNom(mesh,ip1,ip2,0.5,o,no1,to) )
          continue;
        else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
          tet2tri(mesh,k,i,&ptt);
          nortri(mesh,&ptt,no1);
          if ( !MG_GET(pxt->ori,i) ) {
            no1[0] *= -1.0;
            no1[1] *= -1.0;
            no1[2] *= -1.0;
          }
        }
      }
      else if ( tag & MG_GEO ) {
        if ( !BezierRidge(mesh,ip1,ip2,0.5,o,no1,no2,to) )
          continue;
        if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
          tet2tri(mesh,k,i,&ptt);
          nortri(mesh,&ptt,no1);
          no2[0] = to[1]*no1[2] - to[2]*no1[1];
          no2[1] = to[2]*no1[0] - to[0]*no1[2];
          no2[2] = to[0]*no1[1] - to[1]*no1[0];
          dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
          if ( dd > EPSD2 ) {
            dd = 1.0 / sqrt(dd);
            no2[0] *= dd;
            no2[1] *= dd;
            no2[2] *= dd;
          }
        }
      }
      else if ( tag & MG_REF ) {
        if ( !BezierRef(mesh,ip1,ip2,0.5,o,no1,to) )
          continue;
      }
      else {
        if ( !norface(mesh,k,i,v) )  continue;
        if ( !BezierReg(mesh,ip1,ip2,0.5,v,o,no1) )
          continue;
      }
      ier = simbulgept(mesh,list,ilist,o);
      if ( !ier ) {
        ier = dichoto1b(mesh,list,ilist,o,ro);
        memcpy(o,ro,3*sizeof(double));
      }
      ip = newPt(mesh,o,tag);
      if ( !ip ) {
        /* reallocation of point table */
        POINT_REALLOC(mesh,met,ip,mesh->gap,
                      *warn=1;
                      break
                      ,o,tag);
      }
      //CECILE
      if ( met->m )
        met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
      //CECILE
      ier = split1b(mesh,met,list,ilist,ip,1);
      /* if we realloc memory in split1b pt and pxt pointers are not valid */
      pt = &mesh->tetra[k];
      pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

      if ( ier < 0 ) {
        fprintf(stdout," ## Error: unable to split.\n");
        return(-1);
      }
      else if ( !ier ) {
        delPt(mesh,ip);
        continue;
      }
      ns++;
      ppt = &mesh->point[ip];
      if ( MG_EDG(tag) || (tag & MG_NOM) )
        ppt->ref = ref;
      else
        ppt->ref = pxt->ref[i];
      if ( met->m )
        met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

      pxp = &mesh->xpoint[ppt->xp];
      if ( tag & MG_NOM ){
        memcpy(pxp->n1,no1,3*sizeof(double));
        memcpy(pxp->t,to,3*sizeof(double));
      }
      else if ( tag & MG_GEO ) {
        memcpy(pxp->n1,no1,3*sizeof(double));
        memcpy(pxp->n2,no2,3*sizeof(double));
        memcpy(pxp->t,to,3*sizeof(double));
      }
      else if ( tag & MG_REF ) {
        memcpy(pxp->n1,no1,3*sizeof(double));
        memcpy(pxp->t,to,3*sizeof(double));
      }
      else
        memcpy(pxp->n1,no1,3*sizeof(double));
    }

    /* Case of an internal face */
    else {
      if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) continue;
      ilist = coquil(mesh,k,imax,list);
      if ( !ilist ) continue;
      else if ( ilist<0 ) return(-1);
      o[0] = 0.5*(p0->c[0] + p1->c[0]);
      o[1] = 0.5*(p0->c[1] + p1->c[1]);
      o[2] = 0.5*(p0->c[2] + p1->c[2]);
#ifdef SINGUL
      if ( mesh->info.sing && pt->xt && (pxt->tag[imax] & MG_SGL) ) {
        ip = newPt(mesh,o,MG_SGL);
      }
      else
#endif
        ip = newPt(mesh,o,MG_NOTAG);
      if ( !ip )  {
        /* reallocation of point table */
        POINT_REALLOC(mesh,met,ip,mesh->gap,
                      *warn=1;
                      break
                      ,o,MG_NOTAG);
      }
      //CECILE
      if ( met->m )
        met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
      //CECILE
      ier = split1b(mesh,met,list,ilist,ip,1);
      if ( ier < 0 ) {
        fprintf(stdout,"  ## Error: unable to split.\n");
        return(-1);
      }
      else if ( !ier ) {
        delPt(mesh,ip);
      }
      else {
        ppt = &mesh->point[ip];
        met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
        ns++;
      }
    }
  }

  return(ns);
}

/** Collapse edges of length smaller than LOPTS */
static int adpcol(pMesh mesh,pSol met) {
  pTetra     pt;
  pxTetra    pxt;
  pPoint     p0,p1;
  double     len,lmin;
  int        k,ip,iq,list[LMAX+2],ilist,nc;
  char       imin,tag,j,i,i1,i2,ifa0,ifa1;
  int        ier;

  nc = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;
    ier = 0;

    /* find shortest edge */
    imin = -1; lmin = DBL_MAX;
    for (i=0; i<6; i++) {
      if ( pt->xt && (pxt->tag[i] & MG_REQ) )  continue;
      i1  = iare[i][0];
      i2  = iare[i][1];
      len = lenedg(mesh,met,pt->v[i1],pt->v[i2]);
      if ( len < lmin ) {
        lmin = len;
        imin = i;
      }
    }
    if ( imin==-1 )
      fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are boundary and required\n",
              __FILE__,__LINE__,k);
    if ( lmin > LOPTS )  continue;

    ifa0 = ifar[imin][0];
    ifa1 = ifar[imin][1];
    i  =  (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
    j  = iarfinv[i][imin];
    i1 = idir[i][inxt2[j]];
    i2 = idir[i][iprv2[j]];
    ip = pt->v[i1];
    iq = pt->v[i2];
    p0 = &mesh->point[ip];
    p1 = &mesh->point[iq];
    if ( (p0->tag > p1->tag) || (p0->tag & MG_REQ) )  continue;
#ifdef SINGUL
    else if ( mesh->info.sing && (p0->tag & MG_SGL) ) {
      if ( !( pt->xt && (pxt->tag[imin] & MG_SGL) ) )  continue;
    }
#endif

    /* Case of a boundary face */
    ilist = 0;
    if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
      tag = pxt->tag[iarf[i][j]];
      if ( tag & MG_REQ )  continue;
      tag |= MG_BDY;
      if ( p0->tag > tag )   continue;
      if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) continue;
      ilist = chkcol_bdy(mesh,k,i,j,list);
    }
    /* Case of an internal face */
    else {
      if ( p0->tag & MG_BDY )  continue;
      ilist = chkcol_int(mesh,met,k,i,j,list,2);
    }
    if ( ilist > 0 ) {
      ier = colver(mesh,list,ilist,i2);
      if ( ier < 0 )  return(-1);
      else if ( ier ) {
        delPt(mesh,ier);
        nc++;
      }
    }
    else if (ilist < 0 )  return(-1);
  }

  return(nc);
}

/** Analyze tetrahedra and split long / collapse short, according to prescribed metric */
static int adptet(pMesh mesh,pSol met) {
  int      it,nnc,nns,nnf,nnm,maxit,nc,ns,nf,nm;
  int      warn;
  double   maxgap;

  /* Iterative mesh modifications */
  it = nnc = nns = nnf = nnm = warn = 0;
  maxit = 10;
  mesh->gap = maxgap = 0.5;
  do {
    if ( !mesh->info.noinsert ) {
      ns = adpspl(mesh,met,&warn);
      if ( ns < 0 ) {
        fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
        return(0);
      }
    }
    else  ns = 0;

#ifdef USE_SCOTCH
    /*check enough vertex to renum*/
    if ( mesh->info.renum && (it == 1) && (mesh->np/2. > BOXSIZE) && mesh->np>100000 ) {
      /* renumbering begin */
      if ( mesh->info.imprim > 5 )
        fprintf(stdout,"  -- RENUMBERING. \n");
      if ( !renumbering(BOXSIZE,mesh, met) ) {
        fprintf(stdout,"  ## Unable to renumbering mesh. \n");
        fprintf(stdout,"  ## Try to run without renumbering option (-rn 0)\n");
        return(0);
      }

      if ( mesh->info.imprim > 5) {
        fprintf(stdout,"  -- PHASE RENUMBERING COMPLETED. \n");
      }

      if ( mesh->info.ddebug )  chkmsh(mesh,1,0);
      /* renumbering end */
    }
#endif

    if ( !mesh->info.noinsert ) {
      nc = adpcol(mesh,met);
      if ( nc < 0 ) {
        fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
        return(0);
      }
#ifdef DEBUG
      if ( nc ) { printf("APS ADPCOL == %d\n",nc);
        prilen(mesh,met);
      }
#endif
    }
    else  nc = 0;

    if ( !mesh->info.nomove ) {
      nm = movtet(mesh,met,1);
      if ( nm < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nm = 0;

    if ( !mesh->info.noswap ) {
#ifdef PATTERN
      nf = swpmsh(mesh,met);
#else
      nf = swpmsh(mesh,met,NULL);
#endif
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = swptet(mesh,met,1.053);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    nnc += nc;
    nns += ns;
    nnf += nf;
    nnm += nm;
    /* decrease size of gap for reallocation */
    if ( mesh->gap > maxgap/(double)maxit )
      mesh->gap -= maxgap/(double)maxit;
    else
      mesh->gap -= mesh->gap/(double)maxit;

    if ( (abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && ns+nc > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved\n",ns,nc,nf,nm);
    if ( ns < 10 && abs(nc-ns) < 3 )  break;
    else if ( it > 3 && abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
  }
  while( ++it < maxit && nc+ns > 0 );

  if ( warn ) {
    fprintf(stdout,"  ## Error:");
    fprintf(stdout," unable to allocate a new point in last call of adpspl.\n");
    INCREASE_MEM_MESSAGE();
    fprintf(stdout,"  ## Uncomplete mesh. Exiting\n" );
    return(0);
  }

#ifdef USE_SCOTCH
  /*check enough vertex to renum*/
  if ( mesh->info.renum && (mesh->np/2. > BOXSIZE) && mesh->np>100000 ) {
    /* renumbering begin */
    if ( mesh->info.imprim > 5 )
      fprintf(stdout,"  -- RENUMBERING. \n");
    renumbering(BOXSIZE,mesh, met);

    if ( mesh->info.imprim > 5) {
      fprintf(stdout,"  -- PHASE RENUMBERING COMPLETED. \n");
    }
    if ( mesh->info.ddebug )  chkmsh(mesh,1,0);
    /* renumbering end */
  }
#endif

  /*shape optim*/
  it = 0;
  maxit = 2;
  do {
    /* badly shaped process */
    /*ier = badelt(mesh,met);
      if ( ier < 0 ) {
      fprintf(stdout,"  ## Unable to remove bad elements.\n");
      return(0);
      }*/

    if ( !mesh->info.nomove ) {
      nm = movtet(mesh,met,0);
      if ( nm < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh.\n");
        return(0);
      }
      nnm += nm;
    }
    else  nm = 0;

    if ( !mesh->info.noswap ) {
#ifdef PATTERN
      nf = swpmsh(mesh,met);
#else
      nf = swpmsh(mesh,met,NULL);
#endif
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = swptet(mesh,met,1.053);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    if ( (abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && nf+nm > 0 ){
      fprintf(stdout,"                                            ");
      fprintf(stdout,"%8d swapped, %8d moved\n",nf,nm);
    }
  }
  while( ++it < maxit && nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    nm = movtet(mesh,met,3);
    if ( nm < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh.\n");
      return(0);
    }
    nnm += nm;
  }
  else  nm = 0;

  if ( (abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && nm > 0 ){
    fprintf(stdout,"                                            ");
    fprintf(stdout,"                  %8d moved\n",nm);
  }


  if ( abs(mesh->info.imprim) < 4 && (nnc > 0 || nns > 0) )
    fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved, %d iter. \n",
            nns,nnc,nnf,nnm,it);

  return(1);
}

/** split tetra into 4 when more than 1 boundary face */
/*static*/ int anatet4(pMesh mesh, pSol met) {
  pTetra      pt;
  pPoint      ppt;
  pxTetra     pxt;
  int         k,ns;
  char        nf,j;

  ns = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
    nf = 0;
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      for (j=0; j<4; j++)
        if ( pxt->ftag[j] & MG_BDY )  nf++;
    }
    if ( nf > 1 ) {
      if ( !split4bar(mesh,met,k) ) return(-1);
      ns++;
    }
    else {
      nf = 0;
      for (j=0; j<4; j++) {
        ppt = &mesh->point[pt->v[j]];
        if ( ppt->tag & MG_BDY )  nf++;
      }
      if ( nf == 4 ) {
        if ( !split4bar(mesh,met,k) ) return(-1);
        ns++;
      }
    }
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d boundary elements splitted\n",ns);
  return(ns);
}


/** analyze tetrahedra and split if needed */
/*static*/ int anatet(pMesh mesh,pSol met,char typchk) {
  int     ier,nc,ns,nf,nnc,nns,nnf,it,maxit;

  /* analyze tetras : initial splitting */
  nns = nnc = nnf = it = 0;
  maxit = 5;
  mesh->gap = 0.5;
  do {
    /* memory free */
    DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));
    if ( !mesh->info.noinsert ) {

      /* split tetra with more than 2 bdry faces */
      ier = anatet4(mesh,met);
      if ( ier < 0 )  return(0);
      ns = ier;

      /* analyze surface tetras */
      ier = anatets(mesh,met,typchk);

      if ( ier < 0 ) {
        fprintf(stdout,"  ## Unable to complete surface mesh. Exit program.\n");
        return(0);
      }
      ns += ier;
      /* analyze internal tetras */
      ier = anatetv(mesh,met,typchk);
      if ( ier < 0 ) {
        fprintf(stdout,"  ## Unable to complete volume mesh. Exit program.\n");
        return(0);
      }
      ns += ier;
    }
    else  ns = 0;

    if ( !hashTetra(mesh,1) ) {
      fprintf(stdout,"  ## Hashing problem. Exit program.\n");
      return(0);
    }
    if ( typchk == 2 && it == maxit-1 )  mesh->info.fem = 1;

    /* collapse short edges */
    if ( !mesh->info.noinsert ) {
      nc = coltet(mesh,met,typchk);
      if ( nc < 0 ) {
        fprintf(stdout,"  ## Unable to collapse mesh. Exiting.\n");
        return(0);
      }
    }
    else  nc = 0;

    /* attempt to swap */
    if ( !mesh->info.noswap ) {
#ifdef PATTERN
      nf = swpmsh(mesh,met);
#else
      nf = swpmsh(mesh,met,NULL);
#endif
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = swptet(mesh,met,1.1);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    nnc += nc;
    nns += ns;
    nnf += nf;
    if ( (abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && ns+nc+nf > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped\n",ns,nc,nf);
    if ( it > 3 && abs(nc-ns) < 0.1 * MG_MAX(nc,ns) )  break;
  }
  while ( ++it < maxit && ns+nc+nf > 0 );

  if ( (abs(mesh->info.imprim) < 4 || mesh->info.ddebug ) && nns+nnc > 0 )
    fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %d iter.\n",nns,nnc,nnf,it);

#ifdef USE_SCOTCH
    /*check enough vertex to renum*/
    if ( mesh->info.renum && (mesh->np/2. > BOXSIZE) && mesh->np>100000 ) {
      /* renumbering begin */
      if ( mesh->info.imprim > 5 )
        fprintf(stdout,"  -- RENUMBERING. \n");
      if ( !renumbering(BOXSIZE,mesh, met) ) {
        fprintf(stdout,"  ## Unable to renumbering mesh. \n");
        fprintf(stdout,"  ## Try to run without renumbering option (-rn 0)\n");
        return(0);
      }

      if ( mesh->info.imprim > 5) {
        fprintf(stdout,"  -- PHASE RENUMBERING COMPLETED. \n");
      }

      if ( mesh->info.ddebug )  chkmsh(mesh,1,0);
      /* renumbering end */
    }
#endif

  return(1);
}

/** main adaptation routine */
int mmg3d1(pMesh mesh,pSol met) {

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  if ( mesh->info.iso && !chkmani(mesh) ) {
    fprintf(stdout,"  ## Non orientable implicit surface. Exit program.\n");
    return(0);
  }

  /**--- stage 1: geometric mesh */
  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !anatet(mesh,met,1) ) {
    fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }

  /**--- stage 2: computational mesh */
  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");

  /* define metric map */
  if ( !defsiz(mesh,met) ) {
    fprintf(stdout,"  ## Metric undefined. Exit program.\n");
    return(0);
  }

  if ( mesh->info.hgrad > 0. && !gradsiz(mesh,met) ) {
    fprintf(stdout,"  ## Gradation problem. Exit program.\n");
    return(0);
  }

  if ( !anatet(mesh,met,2) ) {
    fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }

#ifdef DEBUG
  puts("---------------------------Fin anatet---------------------");
  outqua(mesh,met);
#endif
  if ( !adptet(mesh,met) ) {
    fprintf(stdout,"  ## Unable to adapt. Exit program.\n");
    return(0);
  }

#ifdef DEBUG
  puts("---------------------Fin adptet-----------------");
  outqua(mesh,met);
#endif
  /* in test phase: check if no element with 2 bdry faces */
  if ( !chkfemtopo(mesh) ) {
    fprintf(stdout,"  ## Topology of mesh unsuited for fem computations. Exit program.\n");
    return(0);
  }

  if ( mesh->info.iso && !chkmani(mesh) ) {
    fprintf(stdout,"  ## Non orientable implicit surface. Exit program.\n");
    return(0);
  }

  return(1);
}
