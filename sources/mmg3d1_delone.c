#include "mmg3d.h"

char  ddb;

#define LOPTLDEL     1.3//1.41
#define LOPTSDEL     0.6
int MMG_npuiss,MMG_nvol,MMG_npres;
/** set triangle corresponding to face ie of tetra k */
void tet2tri(pMesh mesh,int k,char ie,Tria *ptt);

/** find acceptable position for splitting */
static int dichotocpy(pMesh mesh,pSol met,int k,int *vx) {
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
    /* if we realloc mem in the split function, pt is not valid anymore */
    pt = &mesh->tetra[k];
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

/** Find acceptable position for split1b, passing the shell of considered edge, starting from o */
int dichoto1b(pMesh mesh,int *list,int ret,double o[3],double ro[3]);

/** return edges of (virtual) triangle pt that need to be split w/r Hausdorff criterion */
char chkedg(pMesh mesh,Tria *pt,char ori);


/** Search for boundary edges that could be swapped for geometric approximation */
static int swpmshcpy(pMesh mesh,pSol met) {
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
      if ( (!MG_EOK(pt)) || pt->ref < 0 /*|| (pt->tag & MG_REQ)*/ )   continue;
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
            ier = swpbdy(mesh,met,list,ret,it1);
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
/*static*/ int swptetdel(pMesh mesh,pSol met,double crit) {
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
          ier = swpgen(mesh,met,nconf,ilist,list);
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
/*static*/ int movtetdel(pMesh mesh,pSol met,int maxitin) {
  pTetra        pt;
  pPoint        ppt;
  pxTetra       pxt;
  double        *n;
  int           i,k,ier,nm,nnm,ns,lists[LMAX+2],listv[LMAX+2],ilists,ilistv,it;
  int           improve;
  unsigned char j,i0,base;
  int internal,maxit;

  if(maxitin<0) {
    internal = 0;
    maxit = abs(maxitin);
  } else {
    maxit = maxitin;
    internal=1;
  }
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
      if ( !MG_EOK(pt) || pt->ref < 0 /*|| (pt->tag & MG_REQ)*/ )   continue;

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
              if( !(ier=bouleext(mesh,k,i0,i,listv,&ilistv,lists,&ilists)) )
                continue;
              else if ( ier>0 )
                ier = movbdynompt(mesh,listv,ilistv,lists,ilists);
              else  return(-1);
            }
            else if ( ppt->tag & MG_GEO ) {
              if ( !(ier=boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists)) )
                continue;
              else if ( ier>0)
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
          else if(internal){
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
/*static*/ int coltetdel(pMesh mesh,pSol met,char typchk) {
  pTetra     pt;
  pxTetra    pxt;
  pPoint     p0,p1;
  double     ll,ux,uy,uz,hmi2;
  int        k,nc,list[LMAX+2],ilist,base,nnm;
  char       i,j,tag,ip,iq,isnm;
  int        ier;

  nc = nnm = 0;
  hmi2 = mesh->info.hmin*mesh->info.hmin;

  for (k=1; k<=mesh->ne; k++) {
    base = ++mesh->base;
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )   continue;

    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;
    for (i=0; i<4; i++) {
      for (j=0; j<3; j++) {
        if ( pt->xt && (pxt->tag[iarf[i][j]] & MG_REQ) )  continue;
        ier = 0;
        ip = idir[i][inxt2[j]];
        iq = idir[i][iprv2[j]];

        p0 = &mesh->point[pt->v[ip]];
        p1 = &mesh->point[pt->v[iq]];
        if ( p0->flag == base )  continue;
        else if ( (p0->tag & MG_REQ) || (p0->tag > p1->tag) )  continue;

        /* check length */
        if ( typchk == 1 ) {
          ux = p1->c[0] - p0->c[0];
          uy = p1->c[1] - p0->c[1];
          uz = p1->c[2] - p0->c[2];
          ll = ux*ux + uy*uy + uz*uz;
          if ( ll > hmi2 )  continue;
        }
        else if ( typchk == 2 ) {
          ll = MMG5_lenedg(mesh,met,pt->v[ip],pt->v[iq]);
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
        else if (ilist < 0 )  return(-1);
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
/*static*/ int anatetvdel(pMesh mesh,pSol met,char typchk) {
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
        if ( !hashEdge(mesh,&hash,ip1,ip2,ip) ) return(-1);
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
          if ( !hashEdge(mesh,&hash,ip1,ip2,ip) ) return(-1);
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
          ll = MMG5_lenedg(mesh,met,ip1,ip2);
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
        ip  = newPt(mesh,o,0);
        if ( !ip ) {
          /* reallocation of point table */
          POINT_REALLOC(mesh,met,ip,0.5,
                        printf("  ## Error: unable to allocate a new point\n");
                        printf("  ## Check the mesh size or increase");
                        printf(" the allocated memory with the -m option.\n");
                        memlack=1;
                        goto split
                        ,o,0);
          p1  = &mesh->point[ip1];
          p2  = &mesh->point[ip2];
        }

        if ( met->m )
          met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
        if ( !hashEdge(mesh,&hash,ip1,ip2,ip) ) return(-1);
        MG_SET(pt->flag,i);
        nap++;
      }
    }
  }
  if ( !nap )  {
    DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
    return(0);
  }

  /** 3. check and split */
#ifdef DEBUG
  for(k=0;k<12;k++){
    for(j=0;j<7;j++){
      tabtmp[k][j]=0;
    }
  }
#endif
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
#ifdef DEBUG
      tabtmp[0][0]++;
      tabtmp[0][1]/=2;
      tabtmp[0][2]/=2;
#endif
      break;
    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      split2sf(mesh,met,k,vx);
      ns++;
#ifdef DEBUG
      tabtmp[1][0]++;
      tabtmp[1][1]/=2;
      tabtmp[1][2]/=2;
#endif
      break;

    case 33: case 18: case 12: /* 2 opposite edges split */
      split2(mesh,met,k,vx);
      ns++;
#ifdef DEBUG
      tabtmp[2][0]++;
#endif
      break;

    case 11: case 21: case 38: case 56: /* 3 edges on the same faces splitted */
      split3(mesh,met,k,vx);
      ns++;
#ifdef DEBUG
      tabtmp[3][0]++;
      tabtmp[3][1]/=2;
      tabtmp[3][2]/=2;
#endif
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      split3cone(mesh,met,k,vx);
      ns++;
#ifdef DEBUG
      tabtmp[4][0]++;
      tabtmp[4][1]/=2;
      tabtmp[4][2]/=2;
#endif
      break;

    case 35: case 19: case 13: case 37: case 22: case 28: case 26:
    case 14: case 49: case 50: case 44: case 41: /* 3 edges on opposite configuration splitted */
      split3op(mesh,met,k,vx);
      ns++;
#ifdef DEBUG
      tabtmp[5][0]++;
#endif
      break;

    case 23: case 29: case 53: case 60: case 57: case 58:
    case 27: case 15: case 43: case 39: case 54: case 46: /* 4 edges with 3 lying on the same face splitted */
      split4sf(mesh,met,k,vx);
      ns++;
#ifdef DEBUG
      tabtmp[8][0]++;
#endif
      break;

      /* 4 edges with no 3 lying on the same face splitted */
    case 30: case 45: case 51:
      split4op(mesh,met,k,vx);
      ns++;
#ifdef DEBUG
      tabtmp[9][0]++;
#endif
      break;
    case 62: case 61: case 59: case 55: case 47: case 31: /* 5 edges split */
      split5(mesh,met,k,vx);
      ns++;
#ifdef DEBUG
      tabtmp[10][0]++;
#endif
      break;
    case 63: /* 6 edges split */
      split6(mesh,met,k,vx);
#ifdef DEBUG
      tabtmp[11][0]++;
#endif
      ns++;
      break;
    }
  }

  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",nap);

  DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
#ifdef DEBUG
  if(ns) {
    printf("RESULTATS ns %d\n",ns);
    ns = 0;
    for(k=0;k<12;k++){
      printf(" k %d : %5.1f %5.1f %5.1f -init- %e %e \n",k+1,tabtmp[k][0],tabtmp[k][1],
             tabtmp[k][2],
             tabtmp[k][3],tabtmp[k][4]/*,tabtmp[k][5],tabtmp[k][6]*/);
      ns+=tabtmp[k][0];
    }
    printf("on trouve %d split \n",ns);
  }
#endif
  if(memlack)  return(-1);
  return(nap);
}

/** analyze tetra and split on geometric criterion */
/*static*/ int anatetsdel(pMesh mesh,pSol met,char typchk) {
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
        len = MMG5_lenedg(mesh,met,ip1,ip2);
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
          POINT_REALLOC(mesh,met,ip,0.5,
                        printf("  ## Error: unable to allocate a new point.\n");
                        printf("  ## Check the mesh size or increase ");
                        printf("the allocated memory with the -m option.\n");
                        do {
                          delPt(mesh,mesh->np);
                        } while ( mesh->np>npinit );
                        return(-1)
                        ,o,MG_BDY);
        }
        if ( !hashEdge(mesh,&hash,ip1,ip2,ip) ) return(-1);
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
      if ( ic == 0 && dichotocpy(mesh,met,k,vx) ) {
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
    if ( !MG_EOK(pt) || !pt->flag /*|| (pt->tag & MG_REQ)*/ )  continue;
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
/*static*/ int adpspl_delone(pMesh mesh,pSol met,pBucket bucket, int* warn) {
  pTetra     pt;
  pxTetra    pxt;
  Tria       ptt;
  pPoint     p0,p1,ppt;
  pxPoint    pxp;
  double     dd,len,lmax,o[3],to[3],ro[3],no1[3],no2[3],v[3];
  int        k,ip,ip1,ip2,list[LMAX+2],ilist,ns,ref;
  char       imax,tag,j,i,i1,i2,ifa0,ifa1;
  int        ifilt,lon,ret/*,ne*/,ier;

  *warn=0;
  ns = 0;
  ifilt = 0;
  /*ne = mesh->ne;*/
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) /*|| (pt->tag & MG_REQ)*/ )   continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /* find longest edge */
    imax = -1; lmax = 0.0;
    for (i=0; i<6; i++) {
      if ( pt->xt && (pxt->tag[i] & MG_REQ) )  continue;
      ip1  = iare[i][0];
      ip2  = iare[i][1];
      len = MMG5_lenedg(mesh,met,pt->v[ip1],pt->v[ip2]);
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
      else if ( ilist<0 ) return(-1);
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
        POINT_REALLOC(mesh,met,ip,0.2,
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
        fprintf(stdout,"  ## Error: unable to split.\n");
        return(-1);
      }
      else if ( !ier ) {
        delPt(mesh,ip);
        continue;
      }
      ns++;
      addBucket(mesh,bucket,ip);

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
    else if(pt->xt){
      if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) continue;
      ilist = coquil(mesh,k,imax,list);
      if ( !ilist )  continue;
      else if ( ilist<0 ) return(-1);
      o[0] = 0.5*(p0->c[0] + p1->c[0]);
      o[1] = 0.5*(p0->c[1] + p1->c[1]);
      o[2] = 0.5*(p0->c[2] + p1->c[2]);
      ip = newPt(mesh,o,MG_NOTAG);
      if ( !ip )  {
        /* reallocation of point table */
        POINT_REALLOC(mesh,met,ip,0.2,
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
      else if ( !ier ) { //Et on teste pas du tout les qualités ici ?
        delPt(mesh,ip);
      }
      else {
        ppt = &mesh->point[ip];
        met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
        ns++;
      }
    }
    else {    /* Case of an internal face */
      if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) continue;
      ilist = coquil(mesh,k,imax,list);
      if ( !ilist )  continue;
      else if ( ilist<0 ) return(-1);
      o[0] = 0.5*(p0->c[0] + p1->c[0]);
      o[1] = 0.5*(p0->c[1] + p1->c[1]);
      o[2] = 0.5*(p0->c[2] + p1->c[2]);
      ip = newPt(mesh,o,MG_NOTAG);
      if ( !ip )  {
        /* reallocation of point table */
        POINT_REALLOC(mesh,met,ip,0.2,
                      *warn=1;
                      break
                      ,o,MG_NOTAG);
      }
      //CECILE
      if ( met->m )
        met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
      //CECILE

      if ( !buckin_iso(mesh,met,bucket,ip) ) {
        delPt(mesh,ip);
        ifilt++;
        continue;
      }
      lon = cavity(mesh,met,k,ip,list,ilist/2);
      if ( lon < 1 ) {
        delPt(mesh,ip);
        continue;
      } else {
        ret = delone(mesh,met,ip,list,lon);

        if ( ret > 0 ) {
          ppt = &mesh->point[ip];
          met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
          //chkmsh(mesh,0,0);
          addBucket(mesh,bucket,ip);
          ns++;
          continue;
        }
        else if ( ret == 0 ) {
          delPt(mesh,ip);
          continue;
        }
        else {
          delPt(mesh,ip);
          continue;
        }
      }
      /* ier = split1b(mesh,met,list,ilist,ip,1); */
      /* if ( ier<0 ) { */
      /*   fprintf(stdout,"%s:%d: Error: unable to split\n" */
      /*           ,__FILE__,__LINE__); */
      /*   return(-1); */
      /* } */
      /* else if ( !ier ) { //Et on teste pas du tout les qualités ici ? */
      /*   delPt(mesh,ip); */
      /* } */
      /* else { */
      /*   ppt = &mesh->point[ip]; */
      /*   met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]); */
      /*   ns++; */
      /* } */
    }
  }

  printf("on a filtre %7d\n",ifilt);
  return(ns);
}
/** Split edges of length bigger than LOPTL */
int adpsplcol(pMesh mesh,pSol met,pBucket bucket, int* warn) {
  pTetra     pt;
  pxTetra    pxt;
  Tria       ptt;
  pPoint     p0,p1,ppt;
  pxPoint    pxp;
  double     dd,dd2,len,lmax,o[3],to[3],ro[3],no1[3],no2[3],v[3];
  int        k,ip,ip1,ip2,list[LMAX+2],ilist,ns,ref;
  char       imax,tag,j,i,i1,i2,ifa0,ifa1;
  int        ifilt,lon,ret,ne,ier;
  double     lmin;
  int        imin,iq,nc,it,nnc,nns,nnf,nnm,maxit,nf,nm;
  int ii,MMG_npd;
  /* Iterative mesh modifications */
  it = nnc = nns = nnf = nnm = 0;
  maxit = 10;
  MMG_npuiss=MMG_nvol=MMG_npres =MMG_npd=0 ;
  do {
    if ( !mesh->info.noinsert ) {
      *warn=0;
      ns = nc = 0;
      nf = nm = 0;
      ifilt = 0;
      ne = mesh->ne;
      for (k=1; k<=ne; k++) {
        pt = &mesh->tetra[k];
        if ( !MG_EOK(pt)  || (pt->tag & MG_REQ) )   continue;
        pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

        /* find longest and shortest edge */
        imax = -1; lmax = 0.0;
        imin = -1; lmin = DBL_MAX;
        for (ii=0; ii<6; ii++) {
          if ( pt->xt && (pxt->tag[i] & MG_REQ) )  continue;
          ip1  = iare[ii][0];
          ip2  = iare[ii][1];
          len = MMG5_lenedg(mesh,met,pt->v[ip1],pt->v[ip2]);
          if ( len > lmax ) {
            lmax = len;
            imax = ii;
          }
          if ( len < lmin ) {
            lmin = len;
            imin = ii;
          }
        }
        if ( imax==-1 )
          fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are boundary and required\n",
                  __FILE__,__LINE__,k);
        if ( imin==-1 )
          fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are boundary and required\n",
                  __FILE__,__LINE__,k);
        /* imax = ii; */
        /* lmax = len; */
        /* imin = ii; */
        /* lmin = len; */
        if ( lmax >= LOPTLDEL )  {

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
            else if ( ilist<0 ) return(-1);
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
                goto collapse;//continue;
            }
            else {
              //CECILE : je comprend pas pourquoi la normale est mauvaise a la fin
              goto collapse;
              if ( !norface(mesh,k,i,v) )  goto collapse;//continue;
              else {
                dd  = v[0]*no1[0]+v[1]*no1[1]+v[2]*no1[2];
                dd2 = v[0]*no2[0]+v[1]*no2[1]+v[2]*no2[2];

                if ( dd>=dd2 ) {
                  if ( !BezierReg(mesh,ip1,ip2,0.5,v,o,no1) )
                    goto collapse;//continue;
                }
                else {
                  if ( !BezierReg(mesh,ip1,ip2,0.5,v,o,no2) )
                    goto collapse;//continue;
                }
              }
            }
            ier = simbulgept(mesh,list,ilist,o);
            if ( !ier ) {
              ier = dichoto1b(mesh,list,ilist,o,ro);
              memcpy(o,ro,3*sizeof(double));
            }
            ip = newPt(mesh,o,tag);

            if ( !ip ){
              /* reallocation of point table */
              POINT_REALLOC(mesh,met,ip,0.2,
                            *warn=1;
                            goto collapse//break
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
              fprintf(stdout,"  ## Error: unable to split.\n");
              return(-1);
            }
            else if ( !ier ) {
              delPt(mesh,ip);
              goto collapse;//continue;
            } else {
              ns++;
              //addBucket(mesh,bucket,ip);

              ppt = &mesh->point[ip];
              if ( MG_EDG(tag) || (tag & MG_NOM) )
                ppt->ref = ref;
              else
                ppt->ref = pxt->ref[i];
              ppt->tag = tag;
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
            continue;//break;//imax continue;
          }
          else if(pt->xt){
            if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) continue;
            ilist = coquil(mesh,k,imax,list);
            if ( !ilist )    continue;
            else if ( ilist<0 ) return(-1);
            o[0] = 0.5*(p0->c[0] + p1->c[0]);
            o[1] = 0.5*(p0->c[1] + p1->c[1]);
            o[2] = 0.5*(p0->c[2] + p1->c[2]);
            ip = newPt(mesh,o,MG_NOTAG);

            if ( !ip )  {
              /* reallocation of point table */
              POINT_REALLOC(mesh,met,ip,0.2,
                            *warn=1;
                            goto collapse//break
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
            else if ( !ier ) { //Et on teste pas du tout les qualités ici ?
              delPt(mesh,ip);
              goto collapse;//continue;
            }
            else {
              ppt = &mesh->point[ip];
              met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
              addBucket(mesh,bucket,ip);
              ns++;
              continue;//break;//imax continue;
            }
            printf("on doit pas passer la\n");
            /* Case of an internal face */
          } else {
            if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) continue;
            ilist = coquil(mesh,k,imax,list);
            if ( !ilist )    continue;
            else if ( ilist<0 ) return(-1);
            o[0] = 0.5*(p0->c[0] + p1->c[0]);
            o[1] = 0.5*(p0->c[1] + p1->c[1]);
            o[2] = 0.5*(p0->c[2] + p1->c[2]);
            ip = newPt(mesh,o,MG_NOTAG);

            if ( !ip )  {
              /* reallocation of point table */
              POINT_REALLOC(mesh,met,ip,0.2,
                            *warn=1;
                            goto collapse//break
                            ,o,MG_NOTAG);
            }
            //CECILE
            if ( met->m )
              met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
            //CECILE
            //LA DELONE

            if ( !buckin_iso(mesh,met,bucket,ip) ) {
              delPt(mesh,ip);
              ifilt++;
              goto collapse;////continue;
            } else {
              lon = cavity(mesh,met,k,ip,list,ilist/2);
              if ( lon < 1 ) {
                MMG_npd++;
                delPt(mesh,ip);
                goto collapse;//continue;
              } else {
                ret = delone(mesh,met,ip,list,lon);
                if ( ret > 0 ) {
                  ppt = &mesh->point[ip];
                  met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
                  //chkmsh(mesh,0,0);
                  addBucket(mesh,bucket,ip);
                  ns++;
                  continue;//break;//imax continue;
                }
                else if ( ret == 0 ) {
                  MMG_npd++;
                  delPt(mesh,ip);
                  goto collapse;//continue;
                }
                else {
                  MMG_npd++;
                  delPt(mesh,ip);
                  goto collapse;//continue;
                }
                printf("on passe pas la1\n");
              }
              printf("on passe pas la2\n");
            }
            printf("on passe pas la3\n");
          }
          printf("on passe pas la3\n");
        }
      collapse:
        if(lmin > LOPTSDEL) continue;
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

        /* Case of a boundary face */
        ilist = 0;
        if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
          tag = pxt->tag[iarf[i][j]];
          if ( tag & MG_REQ )  continue;
          tag |= MG_BDY;
          if ( p0->tag > tag )   continue;
          if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) continue;
          ilist = chkcol_bdy(mesh,k,i,j,list);
          if ( ilist > 0 ) {
            ier = colver(mesh,list,ilist,i2);
            //nc += ier;
            if ( ier < 0 ) return(-1);
            else if(ier) {
              //delBucket(mesh,bucket,ier);
              delPt(mesh,ier);
              nc++;
              continue;//break;//imax continue;
            }
          }
          else if (ilist < 0 )  return(-1);
        }
        /* Case of an internal face */
        else {
          if ( p0->tag & MG_BDY )  continue;
          ilist = chkcol_int(mesh,met,k,i,j,list,2);
          if ( ilist > 0 ) {
            ier = colver(mesh,list,ilist,i2);
            if ( ilist < 0 ) continue;
            //nc += ier;
            if ( ier < 0 ) return(-1);
            else if(ier) {
              delBucket(mesh,bucket,ier);
              delPt(mesh,ier);
              nc++;
              continue;//break;//imax continue;
            }
          }
          else if (ilist < 0 )  return(-1);
        }

        // }//end for ii

      }
    } /* End conditional loop on mesh->info.noinsert */
    else  ns = nc = 0;

    /* prilen(mesh,met); */
    /*     fprintf(stdout,"    REJECTED : %5d\n",MMG_npd); */
    /* fprintf(stdout,"          VOL      : %6.2f %%    %5d \n", */
    /*         100*(MMG_nvol/(float) */
    /*              MMG_npd),MMG_nvol); */
    /* fprintf(stdout,"          PUISS    : %6.2f %%    %5d \n", */
    /*         100*(MMG_npuiss/(float) MMG_npd),MMG_npuiss); */
    /* fprintf(stdout,"         PROCHE    : %6.2f %%    %5d \n", */
    /*         100*(MMG_npres/(float) MMG_npuiss),MMG_npres); */
    /* MMG_npd=0; */
    /* MMG_npuiss=0; */
    /* MMG_nvol=0; */
    /* MMG_npres=0; */
    if ( !mesh->info.noswap ) {
      nf = swpmshcpy(mesh,met);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;
      nf = swptetdel(mesh,met,1.053);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;
    nnf+=nf;
    fprintf(stdout,"$$$$$$$$$$$$$$$$$$ ITER SWAP %7d\n",nnf);

    if ( !mesh->info.nomove ) {
      nm = movtetdel(mesh,met,-1);
      if ( nm < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh.\n");
        return(0);
      }
    }
    else  nm = 0;
    nnm += nm;
    nnc += nc;
    nns += ns;
    if ( 1 || ((abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc > 0) )
      fprintf(stdout,"     %8d filtered %8d splitted, %8d collapsed, %8d swapped, %8d moved\n",ifilt,ns,nc,nf,nm);
    if ( ns < 10 && abs(nc-ns) < 3 )  break;
    else if ( it > 3 && abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
  }
  while( ++it < maxit && nc+ns > 0 );

  return(1);
}
/** Collapse edges of length smaller than LOPTS */
/*static*/ int adpcoldel(pMesh mesh,pSol met,pBucket bucket) {
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
    if ( !MG_EOK(pt) /*|| (pt->tag & MG_REQ)*/ )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;
    ier = 0;

    /* find shortest edge */
    imin = -1; lmin = DBL_MAX;
    for (i=0; i<6; i++) {
      if ( pt->xt && (pxt->tag[i] & MG_REQ) )  continue;
      i1  = iare[i][0];
      i2  = iare[i][1];
      len = MMG5_lenedg(mesh,met,pt->v[i1],pt->v[i2]);
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
    if ( ilist ) {
      ier = colver(mesh,list,ilist,i2);
      if ( ier < 0 ) return(-1);
      else if ( ier ) {
        delPt(mesh,ier);
        nc++;
      }
    }
  }

  return(nc);
}

/** Analyze tetrahedra and split long / collapse short, according to prescribed metric */
/*static*/ int adptet1(pMesh mesh,pSol met,pBucket bucket) {
  int      it,nnf,nnm,maxit,ns,nf,nm;
  int      warn;

  /*initial swap*/
  if ( !mesh->info.noswap ) {
    nf = swpmshcpy(mesh,met);
    if ( nf < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
    nnf = nf;
    nf = swptetdel(mesh,met,1.053);
    if ( nf < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
    nnf+=nf;
  } else  nnf = nf = 0;
  fprintf(stdout,"$$$$$$$$$$$$$$$$$$ INITIAL SWAP %7d\n",nnf);

  /* Iterative mesh modifications */
  warn = 0;
#ifdef DEBUG
  tabtmp[0][0]=0;
  tabtmp[0][1]=0;     tabtmp[0][2]=0;
#endif
  ns = adpsplcol(mesh,met,bucket,&warn);

#ifdef DEBUG
  if ( ns ) { printf("APS ADPSPLCOL == %d\n",ns);
    prilen(mesh,met);
    printf(" histo %5.1f  %5.1f %5.1f\n", tabtmp[0][0],tabtmp[0][1],tabtmp[0][2]);}
#endif
  if ( ns < 0 ) {
    fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
    return(0);
  }

  if ( warn ) {
    fprintf(stdout,"  ## Error:");
    fprintf(stdout," unable to allocate a new point in last call of adpspl.\n");
    fprintf(stdout,"  ## Check the mesh size or ");
    fprintf(stdout,"increase the allocated memory with the -m option.\n");
    fprintf(stdout,"  ## Uncomplete mesh. Exiting\n" );
    return(0);
  }

#ifdef USE_SCOTCH
  /*check enough vertex to renum*/
  if ( mesh->info.renum && (mesh->np/2. > BOXSIZE) ) {
    /* renumbering begin */
    if ( mesh->info.imprim > 5 )
      fprintf(stdout,"renumbering");
    renumbering(BOXSIZE,mesh, met);

    if ( mesh->info.imprim > 5) {
      fprintf(stdout,"  -- PHASE RENUMBERING COMPLETED. \n");
    }
    if ( mesh->info.ddebug )  chkmsh(mesh,1,0);
    /* renumbering end */
  }
#endif

  /*shape optim*/
  it = nnm = nnf = 0;
  maxit = 2;
  do {
    /* badly shaped process */
    /*ier = badelt(mesh,met);
      if ( ier < 0 ) {
      fprintf(stdout,"  ## Unable to remove bad elements.\n");
      return(0);
      }*/
    if ( !mesh->info.nomove ) {
      nm = movtetdel(mesh,met,0);
      if ( nm < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh.\n");
        return(0);
      }
    }
    else  nm = 0;
    nnm += nm;

    if ( !mesh->info.noswap ) {
      nf = swpmshcpy(mesh,met);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = swptetdel(mesh,met,1.053);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nf+nm > 0 ){
      fprintf(stdout,"                                            ");
      fprintf(stdout,"%8d swapped, %8d moved\n",nf,nm);
    }
  }
  while( ++it < maxit && nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    nm = movtetdel(mesh,met,3);
    if ( nm < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh.\n");
      return(0);
    }
  }
  else  nm = 0;
  nnm += nm;
  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 )
    fprintf(stdout,"                                            ");
  fprintf(stdout,"                  %8d moved\n",nm);

  return(1);
}

  /** Analyze tetrahedra and split long / collapse short, according to prescribed metric */
/*static*/ int adptetdel(pMesh mesh,pSol met,pBucket bucket) {
  int      it,nnc,nns,nnf,nnm,maxit,nc,ns,nf,nm;
  int      warn;

  //ATTENTION MARCHE PAS................................
  if ( !mesh->info.noswap ) {
    nf = swpmshcpy(mesh,met);
    if ( nf < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
    nnf = nf;
    nf = swptetdel(mesh,met,1.053);
    if ( nf < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
  }
  else nf = nnf = 0;
  nnf+=nf;
  fprintf(stdout,"INTIA SWAP %7d\n",nnf);

  /* Iterative mesh modifications */
  it = nnc = nns = nnm = warn = 0;
  maxit = 10;
  nf = nm = 0;
  do {
#ifdef DEBUG
    tabtmp[0][0]=0;
    tabtmp[0][1]=0;     tabtmp[0][2]=0;
#endif
    if ( !mesh->info.noinsert ) {
      ns = adpspl_delone(mesh,met,bucket,&warn);

#ifdef DEBUG
      if ( ns ) { printf("APS ADPSPL == %d\n",ns);
        prilen(mesh,met);
        printf(" histo %5.1f  %5.1f %5.1f\n", tabtmp[0][0],tabtmp[0][1],tabtmp[0][2]);}
#endif
      if ( ns < 0 ) {
        fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
        return(0);
      }

      nc = adpcoldel(mesh,met,bucket);
#ifdef DEBUG
      if(nc){ printf("APS ADPCOL == %d\n",nc);
        prilen(mesh,met);}
#endif
      if ( nc < 0 ) {
        fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
        return(0);
      }
    }
    else  ns = nc = 0;

    if ( !mesh->info.noswap ) {
      nf = swpmshcpy(mesh,met);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;
    nnf += nf;

    if ( !mesh->info.nomove ) {
      nm = movtetdel(mesh,met,-1);
      if ( nm < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh.\n");
        return(0);
      }
    }
    else  nm = 0;

    nnm += nm;
    nnc += nc;
    nns += ns;
    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved\n",ns,nc,nf,nm);
    // if ( ns < 10 && abs(nc-ns) < 3 )  break;
    // else if ( it > 3 && abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
    if ( it > 3 && ns < 0  )  break;
  }
  while( ++it < maxit && nc+ns > 0 );

  if ( warn ) {
    fprintf(stdout,"  ## Error:");
    fprintf(stdout," unable to allocate a new point in last call of adpspl.\n");
    fprintf(stdout,"  ## Check the mesh size or ");
    fprintf(stdout,"increase the allocated memory with the -m option.\n");
    fprintf(stdout," ## Uncomplete mesh. Exiting\n" );
    return(0);
  }

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
      nm = movtetdel(mesh,met,0);
      if ( nm < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh.\n");
        return(0);
      }
      nnm += nm;
    }
    else  nm = 0;

    if ( !mesh->info.noswap ) {
      nf = swpmshcpy(mesh,met);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = swptetdel(mesh,met,1.053);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nf+nm > 0 ){
      fprintf(stdout,"                                            ");
      fprintf(stdout,"%8d swapped, %8d moved\n",nf,nm);
    }
  }
  while( ++it < maxit && nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    nm = movtetdel(mesh,met,3);
    if ( nm < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh.\n");
      return(0);
    }
    nnm += nm;
  }
  else  nm = 0;

  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 ){
    fprintf(stdout,"                                            ");
    fprintf(stdout,"                  %8d moved\n",nm);
  }


  if ( abs(mesh->info.imprim) < 5 && (nnc > 0 || nns > 0) )
    fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved, %d iter. \n",
            nns,nnc,nnf,nnm,it);

  return(1);
}

/** split tetra into 4 when more than 1 boundary face */
/*static*/ int anatet4del(pMesh mesh, pSol met) {
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
/*static*/ int anatetdel(pMesh mesh,pSol met,char typchk) {
  int     ier,nc,ns,nf,nnc,nns,nnf,it,maxit;

  /* analyze tetras : initial splitting */
  nns = nnc = nnf = it = 0;
  maxit = 5;
  do {
    /* memory free */
    DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));
#ifdef DEBUG
    puts("AVT ANATET4");
    prilen(mesh,met);
#endif
    if ( !mesh->info.noinsert ) {
      /* split tetra with more than 2 bdry faces */
      ier = anatet4del(mesh,met);
#ifdef DEBUG
      if ( ier ) { printf("APS ANATET4 == %d\n",ier);
        prilen(mesh,met);}
#endif
      if ( ier < 0 )  return(0);
      ns = ier;

      /* analyze surface tetras */
      ier = anatetsdel(mesh,met,typchk);
#ifdef DEBUG
      if ( ier ) { printf("APS ANATETS == %d\n",ier);
        prilen(mesh,met);}
#endif
      if ( ier < 0 ) {
        fprintf(stdout,"  ## Unable to complete surface mesh. Exit program.\n");
        return(0);
      }
      ns += ier;

      /* analyze internal tetras */
      ier = anatetvdel(mesh,met,typchk);
#ifdef DEBUG
      if(ier){ printf("APS ANATETV == %d\n",ier);
        prilen(mesh,met);}
#endif
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
      nc = coltetdel(mesh,met,typchk);
      if ( nc < 0 ) {
        fprintf(stdout,"  ## Unable to collapse mesh. Exiting.\n");
        return(0);
      }
    }
    else  nc = 0;

    /* attempt to swap */
    if ( !mesh->info.noswap ) {
      nf = swpmshcpy(mesh,met);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = swptetdel(mesh,met,1.1);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    nnc += nc;
    nns += ns;
    nnf += nf;
    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc+nf > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped\n",ns,nc,nf);
    if ( it > 3 && abs(nc-ns) < 0.1 * MG_MAX(nc,ns) )  break;
  }
  while ( ++it < maxit && ns+nc+nf > 0 );

  if ( (abs(mesh->info.imprim) < 5 || mesh->info.ddebug ) && nns+nnc > 0 )
    fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %d iter.\n",nns,nnc,nnf,it);
#ifdef DEBUG
  puts("FIN ANATET");
  prilen(mesh,met);
#endif
  return(1);
}

/** main adaptation routine */
int mmg3d1_delone(pMesh mesh,pSol met) {
  pBucket bucket;

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  if ( mesh->info.iso && !chkmani(mesh) ) {
    fprintf(stdout,"  ## Non orientable implicit surface. Exit program.\n");
    return(0);
  }

  /**--- stage 1: geometric mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !anatetdel(mesh,met,1) ) {
    fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }
#ifdef DEBUG
  outqua(mesh,met);
#endif

  /**--- stage 2: computational mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
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

  if ( !anatetdel(mesh,met,2) ) {
    fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }
#ifdef DEBUG
  puts("---------------------------Fin anatet---------------------");
  outqua(mesh,met);
#endif
#ifdef USE_SCOTCH
    /*check enough vertex to renum*/
    if ( mesh->info.renum  && (mesh->np/2. > BOXSIZE) ) {
      /* renumbering begin */
      if ( mesh->info.imprim > 5 )
        fprintf(stdout,"renumbering");
      renumbering(BOXSIZE,mesh, met);

      if ( mesh->info.imprim > 5) {
        fprintf(stdout,"  -- PHASE RENUMBERING COMPLETED. \n");
      }

      if ( mesh->info.ddebug )  chkmsh(mesh,1,0);
      /* renumbering end */
    }
#endif

  /* CEC : create filter */
  bucket = newBucket(mesh,64);//M_MAX(mesh->mesh->info.bucksiz,BUCKSIZ));
  if ( !bucket )  return(0);

  if ( !adptet1(mesh,met,bucket) ) {
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

  /*free bucket*/
  DEL_MEM(mesh,bucket->head,(bucket->size*bucket->size*bucket->size+1)*sizeof(int));
  DEL_MEM(mesh,bucket->link,(mesh->npmax+1)*sizeof(int));
  DEL_MEM(mesh,bucket,sizeof(Bucket));

  return(1);
}
