#include "mmg3d.h"

extern Info  info;
char  ddb;

/* set triangle corresponding to face ie of tetra k */
void tet2tri(pMesh mesh,int k,char ie,Tria *ptt) {
  pTetra  pt;
  pxTetra pxt;
  char    i,i1,i2;

  pt = &mesh->tetra[k];
  memset(ptt,0,sizeof(Tria));
  ptt->v[0] = pt->v[idir[ie][0]];
  ptt->v[1] = pt->v[idir[ie][1]];
  ptt->v[2] = pt->v[idir[ie][2]];
  for (i=0; i<3; i++) {
    i1 = inxt2[i];
    i2 = iprv2[i];
    hGet(&mesh->htab,ptt->v[i1],ptt->v[i2],&ptt->edg[i],&ptt->tag[i]);
  }
  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    ptt->ref = pxt->ref[ie];
  }
}

/* find acceptable position for splitting */
static int dichoto(pMesh mesh,pSol met,int k,int *vx) {
  pTetra  pt;
  pPoint  pa,pb,ps;
  double  o[6][3],p[6][3];
  float   to,tp,t,oo,op;
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
  tp = oo = 1.0;
  to = op = 0.0;
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
    if ( ier ) {
      to = t;
      oo = to;
      op = tp;
    }
    else
      tp = t;
  }
  while ( ++it < maxit );
  /* restore coords of last valid pos. */
  if ( !ier ) {
    if ( op > oo )
      t = 0.5 * (op + oo);
    else
      t = 0.0;
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

/* Find acceptable position for split1b, passing the shell of considered edge, starting from o */
int dichoto1b(pMesh mesh,int *list,int ret,double o[3],double ro[3]) {
  pTetra  pt;
  pPoint  p0,p1;
  int     iel,np,nq,it,maxit;
  double  m[3],c[3],tp,to,op,oo,t;
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
  it = 0;
  ier = 0;
  tp = oo = 1.0;
  to = op = 0.0;
  do {
    t = 0.5*(to + tp);
    c[0] = m[0] + t*(o[0]-m[0]);
    c[1] = m[1] + t*(o[1]-m[1]);
    c[2] = m[2] + t*(o[2]-m[2]);

    ier = simbulgept(mesh,list,ret,c);
    if ( ier ) {
      to = t;
      oo = to;
      op = tp;
    }
    else
      tp = t;
  }
  while ( ++it < maxit );
  if ( !ier ) {
    if ( op > oo )
      t = 0.5*(op+oo);
    else
      t = 0.0;
  }
  ro[0] = m[0] + t*(o[0]-m[0]);
  ro[1] = m[1] + t*(o[1]-m[1]);
  ro[2] = m[2] + t*(o[2]-m[2]);

  return(1);
}

/* return edges of (virtual) triangle pt that need to be split w/r Hausdorff criterion */
char chkedg(pMesh mesh,Tria *pt) {
  pPoint   p[3];
  xPoint  *pxp;
  double   n[3][3],t[3][3],nt[3],*n1,*n2,t1[3],t2[3];
  double   ps,ps2,ux,uy,uz,ll,il,alpha,dis;
  int      ia,ib,ic;
  char     i,i1,i2;

  ia   = pt->v[0];
  ib   = pt->v[1];
  ic   = pt->v[2];
  p[0] = &mesh->point[ia];
  p[1] = &mesh->point[ib];
  p[2] = &mesh->point[ic];
  pt->flag = 0;

  /* normal recovery */
  for (i=0; i<3; i++) {
    if ( MG_SIN(p[i]->tag) ) {
      nortri(mesh,pt,n[i]);
    }
    else if (p[i]->tag & MG_NOM){
      nortri(mesh,pt,n[i]);
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
    else if ( ll > info.hmax*info.hmax ) {
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
    if ( dis > info.hausd*info.hausd ) {
      MG_SET(pt->flag,i);
      continue;
    }
    ps  = -( t2[0]*ux + t2[1]*uy + t2[2]*uz );
    ps *= il;
    dis = alpha*alpha*fabs(1.0 - ps*ps);

    if ( dis > info.hausd*info.hausd ) {
      MG_SET(pt->flag,i);
      continue;
    }
  }
  return(pt->flag);
}

/* Search for boundary edges that could be swapped for geometric approximation */
static int swpmsh(pMesh mesh) {
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
      if ( !MG_EOK(pt) || pt->ref < 0 )   continue;
      else if ( !pt->xt ) continue;
      pxt = &mesh->xtetra[pt->xt];

      for (i=0; i<4; i++) {
	ier = 0;
	if ( !(pxt->ftag[i] & MG_BDY) ) continue;
	for (j=0; j<3; j++) {
	  ia  = iarf[i][j];
	  ret = coquilface(mesh,k,ia,list,&it1,&it2);
	  ilist = ret / 2;

	  /* CAUTION: trigger collapse with 2 elements */
	  if ( ilist <= 1 )  continue;
	  ier = chkswpbdy(mesh,list,ilist,it1,it2);
	  if ( ier ) {
	    ier = swpbdy(mesh,list,ret,it1);
	    ns++;
	    break;
	  }
	}
	if ( ier )  break;
      }
    }
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );
  if ( (abs(info.imprim) > 5 || info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8d edge swapped\n",nns);

  return(nns);
}

/* Internal edge flipping */
static int swptet(pMesh mesh,pSol met) {
  pTetra   pt;
  int      list[LMAX+2],ilist,k,it,nconf,maxit,ns,nns;
  char     i;

  maxit = 2;
  it = nns = 0;

  do {
    ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<6; i++) {
	nconf = chkswpgen(mesh,k,i,&ilist,list);
	if ( nconf ) {
	  ns++;
	  swpgen(mesh,nconf,ilist,list);
	  break;
	}
      }
    }
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );
  if ( (abs(info.imprim) > 5 || info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8d edge swapped\n",nns);

  return(nns);
}

/* Analyze tetrahedra and move points so as to make mesh more uniform */
static int movtet(pMesh mesh,pSol met,int maxit) {
  pTetra        pt;
  pPoint        ppt;
  pxTetra       pxt;
  double        *n;
  int           i,k,ier,nm,nnm,ns,lists[LMAX+2],listv[LMAX+2],ilists,ilistv,it;
  unsigned char j,i0,base;

  if ( abs(info.imprim) > 5 || info.ddebug )
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
      if ( !MG_EOK(pt) || pt->ref < 0 )   continue;

      /* point j on face i */
      for (i=0; i<4; i++) {
	for (j=0; j<3; j++) {
	  i0  = idir[i][j];
	  ppt = &mesh->point[pt->v[i0]];
	  if ( ppt->flag == base )  continue;
	  else if ( MG_SIN(ppt->tag) )  continue;

	  ier = 0;
	  if ( ppt->tag & MG_BDY ) {
	    pxt = &mesh->xtetra[pt->xt];
	    if ( !(MG_BDY & pxt->ftag[i]) )  continue; // Catch a boundary point by a boundary face
	    else if( ppt->tag & MG_NOM ){
	      if( mesh->adja[4*(k-1)+1+i] ) continue;
	      if( !bouleext(mesh,k,i0,i,listv,&ilistv,lists,&ilists) )  continue;
	      ier = movbdynompt(mesh,listv,ilistv,lists,ilists);
	    }
	    else if ( ppt->tag & MG_GEO ) {
	      if ( !boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists) )  continue;
	      ier = movbdyridpt(mesh,listv,ilistv,lists,ilists);
	    }
	    else if ( ppt->tag & MG_REF ) {
	      if ( !boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists) )  continue;
	      ier = movbdyrefpt(mesh,listv,ilistv,lists,ilists);
	    }
	    else {
	      if ( !boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists) )  continue;
	      n = &(mesh->xpoint[ppt->xp].n1[0]);
	      if ( !directsurfball(mesh, pt->v[i0],lists,ilists,n) )  continue;
	      ier = movbdyregpt(mesh,listv,ilistv,lists,ilists);
	      if ( ier )  ns++;
	    }
	  }
	  else {
	    ilistv = boulevolp(mesh,k,i0,listv);
	    if ( !ilistv )  continue;
	    ier = movintpt(mesh,listv,ilistv);
	  }
	  if ( ier ) {
	    nm++;
	    ppt->flag = base;
	  }
	}
      }
    }
    nnm += nm;
    if ( info.ddebug )  fprintf(stdout,"     %8d moved, %d geometry\n",nm,ns);
  }
  while( ++it < maxit && nm > 0 );

  if ( (abs(info.imprim) > 5 || info.ddebug) && nnm )
    fprintf(stdout,"     %8d vertices moved, %d iter.\n",nnm,it);

  return(nnm);
}

/* attempt to collapse small edges */
static int coltet(pMesh mesh,pSol met,char typchk) {
  pTetra     pt;
  pxTetra    pxt;
  pPoint     p0,p1;
  double     ll,ux,uy,uz;
  int        k,nc,ref,list[LMAX+2],ilist,base,nnm;
  char       i,j,tag,ip,iq,ier,isnm;

  nc = nnm = 0;
  for (k=1; k<=mesh->ne; k++) {
    base = ++mesh->base;
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || MG_SIN(pt->tag))   continue;

    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;
    for (i=0; i<4; i++) {
      for (j=0; j<3; j++) {
	ier = 0;
	ip = idir[i][inxt2[j]];
	iq = idir[i][iprv2[j]];

	p0 = &mesh->point[pt->v[ip]];
	p1 = &mesh->point[pt->v[iq]];
	if ( p0->flag == base )  continue;
	else if ( p0->tag > p1->tag )  continue;

	/* check length */
	if ( typchk == 1 ) {
	  ux = p1->c[0] - p0->c[0];
	  uy = p1->c[1] - p0->c[1];
	  uz = p1->c[2] - p0->c[2];
	  ll = ux*ux + uy*uy + uz*uz;
	  if ( ll > info.hmin*info.hmin )  continue;
	}
	else if ( typchk == 2 ) {
	  ll = lenedg(mesh,pt->v[ip],pt->v[iq]);
	  if ( ll > LSHRT )  continue;
	}

	/* boundary face: collapse ip on iq */
	if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
	  hGet(&mesh->htab,pt->v[ip],pt->v[iq],&ref,&tag);
	  tag |= MG_BDY;

	  isnm = (tag & MG_NOM);
	  if(isnm){
	    if(mesh->adja[4*(k-1)+1+i]) continue;
	  }

	  if ( MG_SIN(tag) || p0->tag > tag )  continue;
	  ilist = chkcol_bdy(mesh,k,i,j,list);
	}
	/* internal face */
	else {
	  isnm = 0;
	  if ( p0->tag & MG_BDY )  continue;
	  ilist = chkcol_int(mesh,k,i,j,list);
	}
	if ( ilist ) {
	  ier = colver(mesh,list,ilist,iq);
	  if ( ier )  break;
	}
      }
      if ( ier ) {
	p1->flag = base;
	if(isnm) nnm++;
	nc++;
	break;
      }
    }
  }
  if ( nc > 0 && (abs(info.imprim) > 5 || info.ddebug) )
    fprintf(stdout,"     %8d vertices removed, %8d non manifold,\n",nc,nnm);

  return(nc);
}

/* analyze volume tetra and split if needed */
static int anatetv(pMesh mesh,pSol met,char typchk) {
  pTetra   pt;
  pPoint   ppt,p1,p2;
  xTetra  *pxt;
  Hash     hash;
  double   ll,o[3],ux,uy,uz;
  int      vx[6],k,ip,ip1,ip2,nap,ns,ne;
  char     i,j,ia;

  /* 1. analysis */
  hashNew(&hash,mesh->np,7*mesh->np);
  ns = nap = 0;

  /* Hash all boundary edges, and put ip = -1 in hash structure */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || MG_SIN(pt->tag) || !pt->xt )  continue;

    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<4; i++) {
      if ( pxt->ftag[i] & MG_BDY ) {
	for (j=0; j<3; j++) {
	  ip1 = pt->v[idir[i][inxt2[j]]];
	  ip2 = pt->v[idir[i][iprv2[j]]];
	  ip  = -1;
	  assert(hashEdge(&hash,ip1,ip2,ip));
	}
	break;
      }
    }
  }

  /* 2. Set flags and split internal edges */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || MG_SIN(pt->tag) )  continue;
    pt->flag = 0;
    for (i=0; i<6; i++) {
      ip  = -1;
      ip1 = pt->v[iare[i][0]];
      ip2 = pt->v[iare[i][1]];
      p1  = &mesh->point[ip1];
      p2  = &mesh->point[ip2];
      if ( (p1->tag & MG_BDY) && (p2->tag & MG_BDY) ) {
	ip = hashGet(&hash,ip1,ip2);
      }
      else {
	ux = p2->c[0] - p1->c[0];
	uy = p2->c[1] - p1->c[1];
	uz = p2->c[2] - p1->c[2];
	ll = ux*ux + uy*uy + uz*uz;
	if ( ll > info.hmax*info.hmax )
	  ip = hashGet(&hash,ip1,ip2);

	else if ( typchk == 2 ) {
	  ll = lenedg(mesh,ip1,ip2);
	  if ( ll > LLONG )
	    ip = hashGet(&hash,ip1,ip2);
	}
      }
      /* new midpoint */
      if ( ip == 0 ) {
	o[0] = 0.5 * (p1->c[0]+p2->c[0]);
	o[1] = 0.5 * (p1->c[1]+p2->c[1]);
	o[2] = 0.5 * (p1->c[2]+p2->c[2]);
	ip = newPt(mesh,o,0);
	ppt = &mesh->point[ip];
	ppt->h = 0.5 * (p1->h + p2->h);
	hashEdge(&hash,ip1,ip2,ip);
	MG_SET(pt->flag,i);
	nap++;
      }
    }
  }

  /* 3. check and split */
  ns = 0;
  ne = mesh->ne;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || MG_SIN(pt->tag) )  continue;
    memset(vx,0,6*sizeof(int));
    pt->flag = 0;
    for (ia=0,i=0; i<3; i++) {
      for (j=i+1; j<4; j++,ia++) {
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

  if ( (info.ddebug || abs(info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",nap);

  free(hash.item);
  return(nap);
}

/* analyze tetra and split on geometric criterion */
static int anatets(pMesh mesh,pSol met,char typchk) {
  pTetra   pt;
  pPoint   ppt,p1,p2;
  Tria     ptt;
  xTetra  *pxt;
  xPoint  *pxp;
  Bezier   pb;
  Hash     hash;
  double   o[3],no[3],to[3],dd,len;
  int      vx[6],k,ip,ic,it,nap,nc,ni,ne,ns,ip1,ip2,ier;
  char     i,j,ia,i1,i2;
  static double uv[3][2] = { {0.5,0.5}, {0.,0.5}, {0.5,0.} };

  /* 1. analysis of boundary elements */
  hashNew(&hash,mesh->np,7*mesh->np);
  ns = nap = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || MG_SIN(pt->tag) || !pt->xt )  continue;

    /* check boundary face cut w/r Hausdorff or hmax */
    pt->flag = 0;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++)
      if ( pxt->ftag[i] & MG_BDY )  break;
    if ( i == 4 )  continue;

    /* virtual triangle */
    tet2tri(mesh,k,i,&ptt);
    if ( typchk == 1 ) {
      if ( !chkedg(mesh,&ptt) )  continue;
      /* put back flag on tetra */
      for (j=0; j<3; j++)
	if ( MG_GET(ptt.flag,j) )  MG_SET(pt->flag,iarf[i][j]);
    }
    else if ( typchk == 2 ) {
      for (j=0; j<3; j++) {
	ia = iarf[i][j];
	i1  = iare[ia][0];
	i2  = iare[ia][1];
	ip1 = pt->v[i1];
	ip2 = pt->v[i2];
	len = lenedg(mesh,ip1,ip2);
	if ( len > LOPTL )  MG_SET(pt->flag,ia);
      }
    }
    if ( !pt->flag )  continue;
    ns++;

    /* geometric support */
    ier = bezierCP(mesh,&ptt,&pb);
    assert(ier);

    /* scan edges in face to split */
    for (j=0; j<3; j++) {
      ia = iarf[i][j];
      if ( !MG_GET(pt->flag,ia) )  continue;
      i1  = iare[ia][0];
      i2  = iare[ia][1];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      ip  = hashGet(&hash,ip1,ip2);
      if ( ip > 0 && !(ptt.tag[j] & MG_GEO) )  continue;

      /* new point along edge */
      ier = bezierInt(&pb,&uv[j][0],o,no,to);
      if ( !ip ) {
	ip = newPt(mesh,o,MG_BDY);
	assert(ip);
	hashEdge(&hash,ip1,ip2,ip);
	ppt = &mesh->point[ip];
	p1  = &mesh->point[ip1];
	p2  = &mesh->point[ip2];
	ppt->h    = 0.5  *(p1->h + p2->h);
	if ( MG_EDG(ptt.tag[j]) || (ptt.tag[j] & MG_NOM) )
	  ppt->ref = ptt.edg[j];
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
    free(hash.item);
    return(ns);
  }

  /* 2. check if split by adjacent; besides, a triangle may have been splitted and not its adjacent
     (thus, the associated n2 may not exist) : update this normal if need be */
  nc = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || MG_SIN(pt->tag) )  continue;
    pxt = &mesh->xtetra[pt->xt];

    /* update face-edge flag */
    for (i=0; i<4; i++) {
      /* virtual triangle */
      memset(&ptt,0,sizeof(Tria));
      if ( pt->xt && pxt->ftag[i] )
	tet2tri(mesh,k,i,&ptt);

      for (j=0; j<3; j++) {
	ia  = iarf[i][j];
	if ( MG_GET(pt->flag,ia) )  continue;
	else if ( MG_SIN(ptt.tag[j]) )  continue;
	ip1 = pt->v[iare[ia][0]];
	ip2 = pt->v[iare[ia][1]];
	ip  = hashGet(&hash,ip1,ip2);
	if ( ip > 0 ) {
	  MG_SET(pt->flag,ia);
	  nc++;
	  /* ridge on a boundary face */
	  if ( !(ptt.tag[j] & MG_GEO) && !(ptt.tag[j] & MG_NOM) )  continue;
	  ppt = &mesh->point[ip];
	  assert(ppt->xp);
	  pxp = &mesh->xpoint[ppt->xp];
	  ier = bezierCP(mesh,&ptt,&pb);
	  ier = bezierInt(&pb,&uv[j][0],o,no,to);
	  memcpy(pxp->n2,no,3*sizeof(double));
	}
      }
    }
  }
  if ( info.ddebug && nc ) {
    fprintf(stdout,"     %d added\n",nc);
    fflush(stdout);
  }

  /* 3. Simulate splitting and delete points leading to invalid configurations */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  it = 1;
  do {
    ni = nc = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || MG_SIN(pt->tag) || !pt->flag )  continue;
      memset(vx,0,6*sizeof(int));
      pt->flag = ic = 0;
      for (ia=0,i=0; i<3; i++) {
	for (j=i+1; j<4; j++,ia++) {
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
  if ( info.ddebug && nc ) {
    fprintf(stdout,"     %d corrected, %d invalid\n",nc,ni);
    fflush(stdout);
  }

  /* 4. splitting */
  ns = 0;
  ne = mesh->ne;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || !pt->flag )  continue;
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
  if ( (info.ddebug || abs(info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"       %7d elements splitted\n",nap);

  free(hash.item);
  return(nap);
}

static int adpspl(pMesh mesh,pSol met) {
  pTetra     pt;
  pxTetra    pxt;
  Tria       ptt;
  pPoint     p0,p1,ppt;
  pxPoint    pxp;
  double     dd,len,lmax,o[3],to[3],ro[3],no1[3],no2[3],v[3];
  int        k,ip,iq,list[LMAX+2],ilist,ns,ref;
  char       imax,tag,j,i,i1,i2,ier,ifa0,ifa1;

  ns = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )   continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /* find longest edge */
    imax = -1; lmax = 0.0;
    for (i=0; i<6; i++) {
      ip  = iare[i][0];
      iq  = iare[i][1];
      len = lenedg(mesh,pt->v[ip],pt->v[iq]);
      if ( len > lmax ) {
	lmax = len;
	imax = i;
      }
    }
    if ( lmax < LOPTL )  continue;

    /* proceed edges according to lengths */
    ifa0 = ifar[imax][0];
    ifa1 = ifar[imax][1];
    i  = (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
    j  = iarfinv[i][imax];
    i1 = idir[i][inxt2[j]];
    i2 = idir[i][iprv2[j]];
    ip = pt->v[i1];
    iq = pt->v[i2];
    p0 = &mesh->point[ip];
    p1 = &mesh->point[iq];

    /* Case of a boundary face */
    if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
      hGet(&mesh->htab,ip,iq,&ref,&tag);
      if ( MG_SIN(tag) )  continue;
      tag |= MG_BDY;
      //ilist = coquil(mesh,k,iarf[i][j],list);
      ilist = coquil(mesh,k,imax,list);
      if ( !ilist )  continue;

      if ( tag & MG_NOM ){
	if( !BezierNom(mesh,ip,iq,0.5,o,no1,to) )
	  continue;
	else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
	  tet2tri(mesh,k,i,&ptt);
	  nortri(mesh,&ptt,no1);
	}
      }
      else if ( tag & MG_GEO ) {
	if ( !BezierRidge(mesh,ip,iq,0.5,o,no1,no2,to) )
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
	else if ( tag & MG_REF ) {
	  if ( !BezierRef(mesh,ip,iq,0.5,o,no1,to) )
	    continue;
	  else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
	    tet2tri(mesh,k,i,&ptt);
	    nortri(mesh,&ptt,no1);
	  }
	}
	else {
	  if ( !norface(mesh,k,i,v) )  continue;
	  else if ( !BezierReg(mesh,ip,iq,0.5,v,o,no1) )
	    continue;
	}
	ier = simbulgept(mesh,list,ilist,o);
	if ( !ier ) {
	  ier = dichoto1b(mesh,list,ilist,o,ro);
	  memcpy(o,ro,3*sizeof(double));
	}
	ip = newPt(mesh,o,MG_NOTAG);
	if ( !ip )  break;
	split1b(mesh,list,ilist,ip);
	ns++;
	ppt = &mesh->point[ip];
	if ( MG_EDG(tag) || (tag & MG_NOM) )
	  ppt->ref = ref;
	else
	  ppt->ref = pxt->ref[i];
	ppt->tag = tag;
	ppt->h   = 0.5 * (p0->h + p1->h);
	mesh->xp++;
	assert(mesh->xp < mesh->xpmax);
	ppt->xp = mesh->xp;
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
    }
    /* Case of an internal face */
    else {
      if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) )  break;     // C'est pas plutot continue ici ? On va arreter super vite non ?
      ilist = coquil(mesh,k,imax,list);
      if ( !ilist )  continue;
      o[0] = 0.5*(p0->c[0] + p1->c[0]);
      o[1] = 0.5*(p0->c[1] + p1->c[1]);
      o[2] = 0.5*(p0->c[2] + p1->c[2]);
      ip = newPt(mesh,o,MG_NOTAG);
      if ( ! ip )  break;
      split1b(mesh,list,ilist,ip);               //Et on teste pas du tout les qualités ici ?
      ppt = &mesh->point[ip];
      ppt->h = 0.5 * (p0->h + p1->h);
      ns++;
    }
  }

  return(ns);
}

static int adpcol(pMesh mesh,pSol met) {
  pTetra     pt;
  pxTetra    pxt;
  pPoint     p0,p1;
  double     len,lmin;
  int        k,ip,iq,list[LMAX+2],ilist,nc,ref;
  char       imin,tag,j,i,i1,i2,ier,ifa0,ifa1;

  nc = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;
    ier = 0;

    /* find shortest edge */
    imin = -1; lmin = DBL_MAX;
    for (i=0; i<6; i++) {
      i1  = iare[i][0];
      i2  = iare[i][1];
      len = lenedg(mesh,pt->v[i1],pt->v[i2]);
      if ( len < lmin ) {
	lmin = len;
	imin = i;
      }
    }
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
    if ( (p0->tag > p1->tag) )  continue;

    /* Case of a boundary face */
    ilist = 0;
    if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
      hGet(&mesh->htab,ip,iq,&ref,&tag);
      if ( MG_SIN(tag) )  continue;
      tag |= MG_BDY;
      if ( p0->tag > tag )   continue;
      if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) continue;
      ilist = chkcol_bdy(mesh,k,i,j,list);
    }
    /* Case of an internal face */
    else {
      if ( p0->tag & MG_BDY )  continue;
      ilist = chkcol_int(mesh,k,i,j,list);
    }
    if ( ilist ) {
      ier = colver(mesh,list,ilist,i2);
      nc += ier;
    }
  }

  return(nc);
}

/* Analyze tetrahedra and split long / collapse short, according to prescribed metric */
static int adptet(pMesh mesh,pSol met) {
  int        ier,it,nnc,nns,nnf,nnm,maxit,nc,ns,nf,nm;

  /* Iterative mesh modifications */
  it = nnc = nns = nnf = nnm = 0;
  maxit = 5;
  do {
    ns = adpspl(mesh,met);
    if ( ns < 0 ) {
      fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
      return(0);
    }

    nc = adpcol(mesh,met);
    if ( nc < 0 ) {
      fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
      return(0);
    }

    nm = movtet(mesh,met,1);
    if ( nm < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }

    nf = swpmsh(mesh);
    if ( nf < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
    nnf += nf;

    nf = swptet(mesh,met);
    if ( nf < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
    nnc += nc;
    nns += ns;
    nnf += nf;
    nnm += nm;
    if ( (abs(info.imprim) > 4 || info.ddebug) && ns+nc > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved\n",ns,nc,nf,nm);
    if ( ns < 10 && abs(nc-ns) < 3 )  break;
    else if ( it > 3 && abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
  }
  while( ++it < maxit && nc+ns > 0 );

  /* badly shaped process */
  ier = badelt(mesh,met);
  if ( ier < 0 ) {
    fprintf(stdout,"  ## Unable to remove bad elements.\n");
    return(0);
  }

  nm = movtet(mesh,met,3);
  if ( nm < 0 ) {
    fprintf(stdout,"  ## Unable to improve mesh.\n");
    return(0);
  }
  nnm += nm;
  if ( abs(info.imprim) < 5 && (nnc > 0 || nns > 0) )
    fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved, %d iter. \n",nns,nnc,nnf,nnm,it);

  return(1);
}

/* split tetra into 4 when more than 1 boundary face */
static int anatet4(pMesh mesh, pSol met) {
  pTetra      pt;
  pPoint      ppt;
  pxTetra     pxt;
  int         k,ns;
  char        nf,j;

  ns = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )   continue;
    nf = 0;
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      for (j=0; j<4; j++)
	if ( pxt->ftag[j] & MG_BDY )  nf++;
    }
    if ( nf > 1 ) {
      split4bar(mesh,met,k);
      ns++;
    }
    else {
      nf = 0;
      for (j=0; j<4; j++) {
	ppt = &mesh->point[pt->v[j]];
	if ( ppt->tag & MG_BDY )  nf++;
      }
      if ( nf == 4 ) {
	split4bar(mesh,met,k);
	ns++;
      }
    }
  }
  if ( (info.ddebug || abs(info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d boundary elements splitted\n",ns);
  return(ns);
}


/* analyze tetrahedra and split if needed */
static int anatet(pMesh mesh,pSol met,char typchk) {
  int     ier,nc,ns,nf,nnc,nns,nnf,it,maxit;

  /* analyze tetras : initial splitting */
  nns = nnc = nnf = it = 0;
  maxit = 5;

  do {
    /* memory free */
    free(mesh->adja);
    mesh->adja = 0;

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

    if ( !hashTetra(mesh) ) {
      fprintf(stdout,"  ## Hashing problem. Exit program.\n");
      return(0);
    }
    if ( typchk == 2 && it == maxit-1 )  info.fem = 1;

    /* collapse short edges */
    nc = coltet(mesh,met,typchk);
    if ( nc < 0 ) {
      fprintf(stdout,"  ## Unable to collapse mesh. Exiting.\n");
      return(0);
    }

    /* attempt to swap */
    nf = swpmsh(mesh);
    if ( nf < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
    nnf += nf;

    nf = swptet(mesh,met);
    if ( nf < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }

    nnc += nc;
    nns += ns;
    nnf += nf;
    if ( (abs(info.imprim) > 4 || info.ddebug) && ns+nc > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped\n",ns,nc,nf);
    if ( it > 3 && abs(nc-ns) < 0.1 * MG_MAX(nc,ns) )  break;
  }
  while ( ++it < maxit && ns+nc+nf > 0 );

  if ( (abs(info.imprim) < 5 || info.ddebug ) && nns+nnc > 0 )
    fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %d iter.\n",nns,nnc,nnf,it);

  return(1);
}

/* main adaptation routine */
int mmg3d1(pMesh mesh,pSol met) {
  if ( abs(info.imprim) > 4 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  /*--- stage 1: geometric mesh */
  if ( abs(info.imprim) > 4 || info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !anatet(mesh,met,1) ) {
    fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }

  /*--- stage 2: computational mesh */
  if ( abs(info.imprim) > 4 || info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");

  /* define metric map */
  if ( !defsiz(mesh,met) ) {
    fprintf(stdout,"  ## Metric undefined. Exit program.\n");
    return(0);
  }

  if ( info.hgrad > 0. && !gradsiz(mesh) ) {
    fprintf(stdout,"  ## Gradation problem. Exit program.\n");
    return(0);
  }

  if ( !anatet(mesh,met,2) ) {
    fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }

  if ( !adptet(mesh,met) ) {
    fprintf(stdout,"  ## Unable to adapt. Exit program.\n");
    return(0);
  }

  return(1);
}
