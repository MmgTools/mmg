#include "mmg3d.h"

extern Info  info;
char ddb;

/* set triangle corresponding to face ie of tetra k */
inline void tet2tri(pMesh mesh,int k,char ie,Tria *ptt) {
  pTetra  pt;
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
int dichoto1b(pMesh mesh,int *list,int ret,double o[3],double ro[3]){
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
    if ( MG_EDG(pt->tag[i]) ) {
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
      if ( !BezierTgt(p[i1]->c,p[i2]->c,&n[i1][0],&n[i2][0],t1,t2) ) {
        t1[0] = ux * il;
        t1[1] = uy * il;
        t1[2] = uz * il;

        t2[0] = - ux * il;
        t2[1] = - uy * il;
        t2[2] = - uz * il;
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

/* Travel mesh and search for boundary edges that should be swapped for
 geometric approximation improvement */
static int swpmsh(pMesh mesh) {
  pTetra   pt;
  pxTetra  pxt;
  int      k,it,maxit,list[LMAX+2],ilist,ret,it1,it2,ns,nns;
  char     i,j,ia,ier;

  if ( abs(info.imprim) > 5 || info.ddebug )
    fprintf(stdout,"  ** IMPROVING SURFACE MESH\n");

  nns = it = 0;
  maxit = 5;
  pxt = 0;
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

          ddb = ilist == 1;
          /* CAUTION: trigger collapse with 2 elements */
          if ( ilist == 1 )  continue; 
          ier = chkswpbdy(mesh,list,ilist,it1,it2);
          if ( ier ) {
            swpbdy(mesh,list,ret,it1);
            ns++;
            break;
          }
        }
        if ( ier )  break; 
      }
    }
    nns += ns;
  }
  while( ++it < maxit && ns > 0 ); 

  if ( abs(info.imprim) > 4 || info.ddebug )
    fprintf(stdout,"     %8d edge swapped, %d iter.\n",nns,it);

  return(nns);
}

/* Analyze tetrahedra and move points so as to make mesh more uniform */
static int movtet(pMesh mesh,pSol met,int maxit){
  pTetra          pt;
  pPoint          ppt;
  pxTetra         pxt;
  double          *n;
  int             i,k,ier,nm,nnm,lists[LMAX+2],listv[LMAX+2],ilists,ilistv,it;
  unsigned char   j,i0,base;
	int ns;

  if ( abs(info.imprim) > 5 || info.ddebug )
    fprintf(stdout,"  ** OPTIMIZING MESH\n");

  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  it = nnm = 0;
  do {
    base = ++mesh->base;
		nm = ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || pt->ref < 0 )   continue;

      /* point j on face i */
      for (i=0; i<4; i++) {
        for (j=0; j<3; j++) {
          i0  = idir[i][j];  
          ppt = &mesh->point[pt->v[i0]];
          if ( MG_SIN(ppt->tag) )  continue;
          //if ( ppt->flag == base )  continue;

					ier = 0;
          if ( ppt->tag & MG_BDY ) {
            pxt = &mesh->xtetra[pt->xt];
            if ( !(MG_BDY & pxt->ftag[i]) )  continue; // Catch a boundary point by a boundary face            
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
		if ( info.ddebug )  fprintf(stdout,"     %8d moved dont %d surface\n",nm,ns);
  }
  while( ++it < maxit && nm > 0 );

  if ( abs(info.imprim) > 4 || info.ddebug )
    fprintf(stdout,"     %8d points moved, %d iter.\n",nnm,it);

  return(1);
}

/* attempt to collapse small edges */
static int coltet(pMesh mesh,pSol met){
  pTetra     pt;
  pxTetra    pxt;
  pPoint     p0,p1;
  double     ll,ux,uy,uz;
  int        it,maxit,k,nc,nnc,np,ref,ne,list[LMAX+2],ilist,base;
  char       i,j,tag,ia,ip,iq,ier;

  nnc = it = 0;
  maxit = 20;
  do {
    nc = 0;
    ne = mesh->ne;
    for (k=1; k<=mesh->ne; k++) {
      base = ++mesh->base;
      pt = &mesh->tetra[k];   
      if ( !MG_EOK(pt) || MG_SIN(pt->tag))   continue;

      pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;    
      for (i=0; i<4; i++) {
        for (j=0; j<3; j++) {
          ier = 0;
          ia = iarf[i][j];
          ip = idir[i][inxt2[j]];
          iq = idir[i][iprv2[j]];
        
          p0 = &mesh->point[pt->v[ip]];
          p1 = &mesh->point[pt->v[iq]];
          if ( p0->flag == base )  continue;
          else if ( p0->tag > p1->tag )  continue;

          /* check length */
          ux = p1->c[0] - p0->c[0];
          uy = p1->c[1] - p0->c[1];
          uz = p1->c[2] - p0->c[2];
          ll = ux*ux + uy*uy + uz*uz;

          if ( ll >= info.hmin*info.hmin )  continue;

          /* boundary face: collapse ip on iq */
          if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
            hGet(&mesh->htab,pt->v[ip],pt->v[iq],&ref,&tag);
            tag |= MG_BDY;        
            if ( MG_SIN(tag) || p0->tag > tag )  continue;
            ilist = chkcol_bdy(mesh,k,i,j,list);
          }
          /* internal face */
          else {
            if ( p0->tag & MG_BDY )  continue;
            ilist = chkcol_int(mesh,k,i,j,list);
          }
          if ( ilist ) {
            np++;          
            ier = colver(mesh,list,ilist,iq);
            if ( ier )  break;
          }
        }
        if ( ier ) {
          p1->flag = base;
          nc++;
          break;
        }
      }
    }
    nnc += nc;
    if ( info.ddebug && nc > 0 )
      fprintf(stdout,"     %d vertices removed\n",nc);
  }
  while ( ++ it < maxit && nc > 0 );
  if ( ( abs(info.imprim) > 4 || info.ddebug ) && nnc > 0 )
    fprintf(stdout,"     %8d vertices removed,  %d iter.\n",nnc,it);

  return(1);
}

/* analyze volume tetra and split if needed */
static int anatetv(pMesh mesh,pSol met) {
  pTetra   pt;
  pPoint   p1,p2; 
  xTetra  *pxt;
  Hash     hash;
  double   ll,o[3],ux,uy,uz;
  int      vx[6],k,ip,ip1,ip2,ns,ne,n3cone,n3sf,n5,n2,n2op,n1,n6,n3op,n4sf,n4op;
  char     i,j,ia;

  /* 1. analysis */
  hashNew(&hash,mesh->np,7*mesh->np);
  ns = 0;
 
  n3cone = 0 ; n3sf = 0 ; n5 = 0 ; n2 = 0 ; n2op = 0;
  n1 =0 ; n6 =0; n3op = 0; n4sf = 0; n4op = 0;  
  
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
      }
      /* new midpoint */
      if ( ip == 0 ) {
        o[0] = 0.5 * (p1->c[0]+p2->c[0]);
        o[1] = 0.5 * (p1->c[1]+p2->c[1]);
        o[2] = 0.5 * (p1->c[2]+p2->c[2]);
        ip = newPt(mesh,o,0);
        hashEdge(&hash,ip1,ip2,ip);
        MG_SET(pt->flag,i);
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
n1++;
      break;
    case 48: case 24: case 40: case 6: case 34: case 36: 
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      split2sf(mesh,met,k,vx);
      ns++;
n2++;
      break;
     
    case 33: case 18: case 12: /* 2 opposite edges split */ 
      split2(mesh,met,k,vx);
      ns++;
n2op++;
      break;
    
    case 11: case 21: case 38: case 56: /* 3 edges on the same faces splitted */
      split3(mesh,met,k,vx);
      ns++;
n3sf++;
      break;
     
    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      split3cone(mesh,met,k,vx);
      ns++;
n3cone++;
      break;
      
    case 35: case 19: case 13: case 37: case 22: case 28: case 26: 
    case 14: case 49: case 50: case 44: case 41: /* 3 edges on opposite configuration splitted */
      split3op(mesh,met,k,vx);
      ns++;
n3op++;
      break;
    
    case 23: case 29: case 53: case 60: case 57: case 58: 
    case 27: case 15: case 43: case 39: case 54: case 46: /* 4 edges with 3 lying on the same face splitted */
      split4sf(mesh,met,k,vx);
      ns++;
n4sf++;
      break;
    
    /* 4 edges with no 3 lying on the same face splitted */
    case 30: case 45: case 51: 
      split4op(mesh,met,k,vx);
    ns++;
    n4op++;
    break;
    case 62: case 61: case 59: case 55: case 47: case 31: /* 5 edges split */
      split5(mesh,met,k,vx);
      ns++;
n5++;
      break;
    case 63: /* 6 edges split */
      split6(mesh,met,k,vx);
      ns++;
n6++;
      break;
    }
  }
  
  if ( info.ddebug ) {
    printf("nombre de splits a 1 arete : %d \n",n1);
    printf("nombre de splits a 2 aretes sur la meme face : %d, opposees %d\n",n2,n2op);
    printf("nombre de splits a 3 aretes sur la meme face : %d, en cone %d, opposees : %d\n",n3sf,n3cone,n3op);
    printf("nombre de splits a 4 aretes, avec 3 sur la meme face : %d, opposees : %d \n",n4sf,n4op);
    printf("nombre de splits a 5 aretes : %d \n",n5);
    printf("nombre de splits a 6 arete : %d \n",n6);
  }
  if ( (info.ddebug || abs(info.imprim) > 5) && ns > 0 )  
    fprintf(stdout,"     %7d splitted\n",ns);

  free(hash.item);
  return(ns);
}

/* analyze tetra and split on geometric criterion */
static int anatets(pMesh mesh,pSol met) {
  pTetra   pt;
  pPoint   ppt,p1,p2;
  Tria     ptt;
  xTetra  *pxt;
  xPoint  *pxp;
  Bezier   pb;
  Hash     hash;
  double   o[3],no[3],to[3],dd;
  int      vx[6],k,ip,nc,ni,ne,ns,ip1,ip2,ier;
  char     i,j,ia,i1,i2;
  static double uv[3][2] = { {0.5,0.5}, {0.,0.5}, {0.5,0.} }; 

  /* 1. analysis of boundary elements */
  hashNew(&hash,mesh->np,7*mesh->np);
  ns = 0;
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

    if ( !chkedg(mesh,&ptt) )  continue;
    /* put back flag on tetra */
    for (j=0; j<3; j++)
      if ( MG_GET(ptt.flag,j) )  MG_SET(pt->flag,iarf[i][j]);
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
      if ( !MG_EDG(ptt.tag[j]) && ip > 0 )  continue;

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
        ppt->ref  = ptt.edg[j];
        ppt->tag |= ptt.tag[j];

        pxp = &mesh->xpoint[ppt->xp];
        memcpy(pxp->n1,no,3*sizeof(double));
        memcpy(pxp->t,to,3*sizeof(double));
      }
      else if ( MG_EDG(ptt.tag[j]) ) {
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
          if ( !(ptt.tag[j] & MG_GEO) )  continue;
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
  do {
    ni = nc = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || MG_SIN(pt->tag) || !pt->flag )  continue;
      memset(vx,0,6*sizeof(int));
      pt->flag = 0;
      for (ia=0,i=0; i<3; i++) {
        for (j=i+1; j<4; j++,ia++) {
          vx[ia] = hashGet(&hash,pt->v[i],pt->v[j]);
          if ( vx[ia] > 0 )  MG_SET(pt->flag,ia);
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
      /* check for badly shaped element: split on straight edge */
      if ( orvol(mesh->point,pt->v) < NULKAL ) {
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
      else {
        /* attempt to find valid position */
        if ( !dichoto(mesh,met,k,vx) ) {
          pt->flag = 0;
          for (ia=0,i=0; i<3; i++) {
            for (j=i+1; j<4; j++,ia++) {
              if ( vx[ia] > 0 ) {
                ni++;
                delPt(mesh,vx[ia]);
                hashPop(&hash,pt->v[i],pt->v[j]);
              }
            }
          }
        }
      }   
    }
  }
  while( nc > 0 );
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
    fprintf(stdout,"       %7d elements splitted\n",ns);

  free(hash.item);  
  return(ns);
}

/* Analyze tetrahedra and split long / collapse short, according to prescribed metric */
static int adptet(pMesh mesh,pSol met) {
  pTetra     pt;
  pxTetra    pxt;
  pPoint     p0,p1,ppt;
  pxPoint    pxp;
  double     len,o[3],to[3],ro[3],no1[3],no2[3],v[3],edglen[6],dist1[6];
  int        k,it,nnc,maxit,list[LMAX+2],ilist,nc,ns,nns,ref,ip;
	char       tag,j,i,l,ia,indp,indq,ier,ifa0,ifa1,alert,sort[6];

  if ( abs(info.imprim) > 5 || info.ddebug )
    fprintf(stdout,"  ** ADAPTING MESH\n"); 

  /* Iterative mesh modifications */  
  it = nnc = nns = alert = 0;
  maxit = 20;
  do {
    nc = ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )   continue;
      pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;
      ier = 0;

      /* Store the 6 lengths of edges of face i, then sort it by priority order */
      for(ia=0; ia<6; ia++){
        indp = iare[ia][0];
        indq = iare[ia][1];
        edglen[ia] = lenedg(mesh,pt->v[indp],pt->v[indq]);
        dist1[ia]  = edglen[ia] > 1.0 ? (1.0/edglen[ia]) : edglen[ia]; 
      }
      nsort(6,dist1,sort);

      /* proceed edges according to lengths */
      for (l=0; l<6; l++) {
				ia   = sort[l];
        ifa0 = ifar[ia][0];
        ifa1 = ifar[ia][1];
        i = ifa0;
        if ( pt->xt && (pxt->ftag[ifa1] & MG_BDY) )
          i = ifa1;
        j = iarfinv[i][ia];      
        indp = idir[i][inxt2[j]];
        indq = idir[i][iprv2[j]];
        p0 = &mesh->point[pt->v[indp]];
        p1 = &mesh->point[pt->v[indq]];

        len = lenedg(mesh,pt->v[indp],pt->v[indq]);
				if ( len >= LSHRT && len <= LLONG )  continue;
  
        /* Case of a boundary face */
        if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
          hGet(&mesh->htab,pt->v[indp],pt->v[indq],&ref,&tag);
          if ( MG_SIN(tag) )  continue;
          tag |= MG_BDY;  
        
          /* short edge collapse */
          if ( len < LSHRT ) {
            if ( (p0->tag > p1->tag) )  continue; 
            else if ( p0->tag > tag )   continue;

            ilist = chkcol_bdy(mesh,k,i,j,list);
            if ( ilist ) {
              ier = colver(mesh,list,ilist,indq);
              nc += ier;
              if ( ier )  break; 
            }
          }
          /* long edge split */ 
          else {
            ilist = coquil(mesh,k,ia,list);
            assert(ilist);           
            if ( tag & MG_GEO ) {
              if ( !BezierRidge(mesh,pt->v[indp],pt->v[indq],0.5,o,no1,no2,to) )
                continue;
            }
            else if ( tag & MG_REF ) {
              if ( !BezierRef(mesh,pt->v[indp],pt->v[indq],0.5,o,no1,to) )
                continue;
            }
            else {
              if ( !norface(mesh,k,i,v) )  continue;
              else if ( !BezierReg(mesh,pt->v[indp],pt->v[indq],0.5,v,o,no1) )
                continue;  
            }    
            ier = simbulgept(mesh,list,ilist,o);               
            if ( !ier ) {
              ier = dichoto1b(mesh,list,ilist,o,ro);
              memcpy(o,ro,3*sizeof(double));
            } 
            ip = newPt(mesh,o,MG_NOTAG);
            if ( !ip ) {
							alert = 1;
							break;
	          }
            //assert(ip);
            split1b(mesh,list,ilist,ip);
            ns++;

            ppt = &mesh->point[ip];
            ppt->tag = tag;
            ppt->ref = ref; 
            mesh->xp++;
            assert(mesh->xp < mesh->xpmax);
            ppt->xp = mesh->xp;
            pxp = &mesh->xpoint[ppt->xp];
            ppt->h = 0.5*(p0->h + p1->h);
                  
            if ( tag & MG_GEO ) {
              memcpy(pxp->n1,no1,3*sizeof(double));
              memcpy(pxp->n2,no2,3*sizeof(double));
              memcpy(pxp->t,to,3*sizeof(double));
            }
            else if ( tag & MG_REF ) {
              memcpy(pxp->n1,no1,3*sizeof(double));
              memcpy(pxp->t,to,3*sizeof(double));
            }
            else {
              memcpy(pxp->n1,no1,3*sizeof(double));
            }
            break;
          }              
        }
   
        /* Case of an internal face */
        else {
	        /* short edge collapse */
          if ( len < LSHRT ) {
            if ( p0->tag & MG_BDY )  continue;
            else if ( p0->tag > p1->tag ) continue;
           
            ilist = chkcol_int(mesh,k,i,j,list);
            if ( ilist ) {
              ier = colver(mesh,list,ilist,indq);
              nc += ier;
              if ( ier )  break;
            }
          }
          else {
            if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) )  break;                             
            ilist = coquil(mesh,k,ia,list);
            assert(ilist);

            o[0] = 0.5*(p0->c[0] + p1->c[0]);
            o[1] = 0.5*(p0->c[1] + p1->c[1]);
            o[2] = 0.5*(p0->c[2] + p1->c[2]);
                
            ier = simbulgept(mesh,list,ilist,o); 
            if ( !ier )  continue;

            ip = newPt(mesh,o,MG_NOTAG);
						if ( ! ip ) {
							alert = 1;
							break;
						}
            //assert(ip);
            split1b(mesh,list,ilist,ip);
            ns++;
            ppt = &mesh->point[ip];
            ppt->h = 0.5*(p0->h + p1->h);
            break;
          }
        }
      }
			if ( alert )  break;   
    }
    nnc += nc;
    nns += ns;
  }
  while( ++it < maxit && nc+ns > 0 && ! alert );

  if ( abs(info.imprim) > 4 && (nnc > 0 || nns > 0) )
    fprintf(stdout,"     %8d splitted, %d collapsed, %d iter. \n",nns,nnc,it);

  return(1);
}


/* split tetra into 4 when more than 1 boundary face */
static int anatet4(pMesh mesh, pSol met) {
  pTetra      pt;
  pPoint      ppt;
  pxTetra     pxt;
  int         k,ne,ns;
  char        nf,j;

  ns = 0;
  ne = mesh->ne;  
  for (k=1; k<=ne; k++) {
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
    fprintf(stdout,"     %d elements splitted on boundary\n",ns);
  return(ns);
}


/* analyze tetrahedra and split if needed */
static int anatet(pMesh mesh,pSol met) {
  int     is,ns,nns,it,maxit;

  /* memory free */
  free(mesh->adja);
  mesh->adja = 0;

  /* analyze tetras : initial splitting */
  nns = it = 0;
  maxit = 5;
  do {
    /* split tetra with more than 2 bdry faces */
    is = anatet4(mesh,met);
    if ( is < 0 )  return(0);
    ns = is;

    /* analyze surface tetras */
    is = anatets(mesh,met);
    if ( is < 0 ) {
      fprintf(stdout,"  ## Unable to complete surface mesh. Exit program.\n");
      return(0);
    }
    ns += is;

    /* analyze internal tetras */
    is = anatetv(mesh,met);
    if ( is < 0 ) {
      fprintf(stdout,"  ## Unable to complete volume mesh. Exit program.\n");
      return(0);
    }
    ns  += is;
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );

  if ( !hashTetra(mesh) ) {
    fprintf(stdout,"  ## Hashing problem. Exit program.\n");
    return(0);
  }

  if ( ( abs(info.imprim) > 4 || info.ddebug ) && nns > 0 )
    fprintf(stdout,"     %8d elements splitted, %d iter.\n",nns,it);

  return(1);
}

/* main adaptation routine */
int mmg3d1(pMesh mesh,pSol met) {
  int        is,nc,ns,nnc,nns,it,itt,maxit;

  if ( abs(info.imprim) > 4 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  /*--- stage 1: geometric mesh */
  if ( abs(info.imprim) > 4 || info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");  

  if ( !anatet(mesh,met) ) {
    fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }
  if ( !coltet(mesh,met) ) {
    fprintf(stdout,"  ## Unable to collapse mesh. Exiting.\n");
    return(0);
  }
  if ( !swpmsh(mesh) ) {   
		saveMesh(mesh);
    fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
    return(0);
	}
  if ( !movtet(mesh,met,2) ) {
    fprintf(stdout,"  ## Unable to optimize mesh. Exiting.\n");
    return(0);
  }

  /*--- stage 2: computational mesh */
  if ( abs(info.imprim) > 4 || info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");  

  /* define metric map */
  if ( !defsiz_iso(mesh,met) ) {
   fprintf(stdout,"  ## Metric undefined. Exit program.\n");
   return(0);
  }
  if ( info.hgrad > 0. && !gradsiz_iso(mesh) ) {
   fprintf(stdout,"  ## Gradation problem. Exit program.\n");
   return(0);
  }
  if ( !adptet(mesh,met) ) {
    fprintf(stdout,"  ## Unable to adapt. Exit program.\n");
    return(0);
	}
  if ( !swpmsh(mesh) ) {
    fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
    return(0);
	}
  if ( !movtet(mesh,met,5) ) {
    fprintf(stdout,"  ## Unable to optimize mesh. Exiting.\n");
    return(0);
  }

  return(1);



 goto skip;
  
  maxit = 20;
  
  nnc= it = 0;
  do {
    nc   = coltet(mesh,met);
    nnc += nc;
  }
  while ( ++it < maxit && nc > 0 );
  
  chkmsh(mesh,0,0);
  return(1);
  
  /*--- stage 1: adapt mesh on geometric criteria */
  itt = 0;
  do {
    free(mesh->adja);
    mesh->adja = 0;

    /* analyze tetras : initial splitting */
    it  = 0;
    nnc = nns = 0;  
    do {
      is = anatet4(mesh,met);
      if ( is < 0 ) {
        fprintf(stdout,"  ## Unable to modify surface mesh. Exit program.\n");
        return(0);
      }
      nns = is;
      /* analyze surface tetras */
      is = anatets(mesh,met);
      if ( is < 0 ) {
        fprintf(stdout,"  ## Unable to complete surface mesh. Exit program.\n");
        return(0);
      }
      ns = is;
      /* analyze internal tetras */
      is = anatetv(mesh,met);
      if ( is < 0 ) {
        fprintf(stdout,"  ## Unable to complete volume mesh. Exit program.\n");
        return(0);
      }
      ns  += is;
      nns += ns;
    }
    while ( ++it < maxit && ns > 0 );

    if ( !hashTetra(mesh) ) {
      fprintf(stdout,"  ## Hashing problem. Exit program.\n");
      return(0);
    }
    
    nnc= it = 0;
    
   /* do {
      nc   = coltet(mesh,met);
      nnc += nc;
    }
    while ( ++it < maxit && nnc > 0 );*/
    
 }
  while ( ++itt < maxit && nns+nnc > 0 );

  if ( (info.ddebug || abs(info.imprim) > 4) && nns+nnc ) {
    fprintf(stdout,"     %d elements splitted, %d collapsed\n",nns,nnc);
  }

  return(1);

skip:

  /* analyze triangles: adapt */
  if ( !adptetra(mesh,met) ) {
   fprintf(stdout,"  ## Unable to adapt. Exit program.\n");
   return(0);
  }
    
  /* analyze triangles: move */
/*  if ( !movetetra(mesh) ) {
    fprintf(stdout,"  ## Unable to proceed adaptation. Exit program.\n");
    return(0);
  } 
  */
  chkmsh(mesh,0,0);
  
  return(1);
   
  /*--- stage 2: adapt mesh on size */
  nns = nnc = it = 0;
  maxit = 10;  
  do {
    /* collapse */
    nc = coltet(mesh,met);
    nnc += nc;
  }
  while ( ++it < maxit && ns+nc > 0 );
  if ( (info.ddebug || abs(info.imprim) > 4) && nnc ) {
    fprintf(stdout,"     %d collapsed\n",nnc);
  }   
     
  /* analyze triangles: move */
/*  if ( !movetetra(mesh) ) {
    fprintf(stdout,"  ## Unable to proceed adaptation. Exit program.\n");
    return(0);
  } 
*/    
  
  return(1);
}
   