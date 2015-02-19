#include "mmg3d.h"

extern char  ddb;

/** Check whether collapse ip -> iq could be performed, ip internal ;
 *  'mechanical' tests (positive jacobian) are not performed here */
int chkcol_int(pMesh mesh,pSol met,int k,char iface,char iedg,int *list,char typchk) {
  pTetra   pt,pt0;
  pPoint   p0;
  double   calold,calnew,caltmp,lon;
  int      j,iel,ilist,nq;
  char     i,jj,ip,iq;

  ip  = idir[iface][inxt2[iedg]];
  iq  = idir[iface][iprv2[iedg]];
  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  nq  = pt->v[iq];
  ilist = boulevolp(mesh,k,ip,list);
  lon = 1.e20;
  if ( typchk == 2 && met->m ) {
    lon = lenedg(mesh,met,pt->v[ip],nq);
    lon = MG_MIN(lon,LSHRT);
    lon = MG_MAX(1.0/lon,LLONG);
  }
  calold = calnew = DBL_MAX;
  for (j=0; j<ilist; j++) {
    iel = list[j] / 4;
    ip  = list[j] % 4;
    pt  = &mesh->tetra[iel];
    /* exclude elements from shell */
    for (jj=0; jj<4; jj++)  if ( pt->v[jj] == nq )  break;
    if ( jj < 4 )  continue;
    memcpy(pt0,pt,sizeof(Tetra));

    /* prevent from recreating internal edge between boundaries */
    if ( mesh->info.fem ) {
      p0 = &mesh->point[nq];
      if ( p0->tag & MG_BDY ) {
        i = ip;
        for (jj=0; jj<3; jj++) {
          i = inxt3[i];
          p0 = &mesh->point[pt->v[i]];
          if ( p0->tag & MG_BDY )  return(0);
        }
      }
    }

    pt0->v[ip] = nq;
    calold = MG_MIN(calold,pt->qual);
    caltmp = orcal(mesh,0);
    if ( caltmp < EPSD )  return(0);
    calnew = MG_MIN(calnew,caltmp);
    /* check length */
    if ( typchk == 2 && met->m ) {
      for (jj=0; jj<6; jj++) {
        if ( lenedg(mesh,met,pt0->v[iare[jj][0]],pt0->v[iare[jj][1]]) > lon )  return(0);
      }
    }
  }
  if ( calold < NULKAL && calnew <= calold )  return(0);
  else if ( calnew < NULKAL || calnew < 0.3*calold )  return(0);

  return(ilist);
}

/** Topological check on the surface ball of np and nq in collapsing np->nq ;
 *  iface = boundary face on which lie edge iedg - in local face num.
 *  (pq, or ia in local tet notation) */
static int topchkcol_bdy(pMesh mesh,int k,int iface,char iedg,int *lists,int ilists) {
  pTetra   pt;
  pxTetra  pxt;
  int      nump,numq,piv0,piv,iel,jel,nap,nbp,naq,nbq,nro,adj,*adja;
  char     ip,iq,ipiv,iopp,i,j,jface,ipa,ipb,isface;

  pt = &mesh->tetra[k];
  ip = idir[iface][inxt2[iedg]];
  iq = idir[iface][iprv2[iedg]];
  nump = pt->v[ip];
  numq = pt->v[iq];

  /* Pivot in enumeration of the surface ball of np */
  ipiv = idir[iface][inxt2[idirinv[iface][ip]]];
  piv0 = pt->v[ipiv];

  /* Surface ball has been enumerated as f1,...,f2 - f1,f2 = both triangles of surface shell */
  if ( piv0 == numq ) {
    /*  Point nap, facing the first vanishing face in surface ball of p */
    nro = pt->v[idir[iface][iprv2[idirinv[iface][ip]]]];

    jel = lists[1] / 4;
    jface = lists[1] % 4;

    pt = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != nro ) break;
    }
    assert(j<3);

    nap = pt->v[i];

    /* Unfold shell of (nq,nro), starting from (k,iface), with pivot np */
    adj = k;
    piv = nump;
    do {
      iel = adj;
      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      /* Identification of edge number in tetra iel */
      for (i=0; i<6; i++) {
        ipa = iare[i][0];
        ipb = iare[i][1];
        if ( ((pt->v[ipa] == numq) && (pt->v[ipb] == nro)) || ((pt->v[ipa] == nro)  && (pt->v[ipb] == numq))  ) break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ ifar[i][0] ] == piv ) {
        adj  = adja[ ifar[i][0] ] / 4;
        ipiv = ifar[i][1];
        iopp = ifar[i][0];
        piv  = pt->v[ipiv];
      }
      else {
        adj  = adja[ ifar[i][1] ] / 4;
        ipiv = ifar[i][0];
        iopp = ifar[i][1];
        piv  = pt->v[ipiv];
      }

      isface = 0;
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && ( adj != k ) && !isface );

    naq = piv;
    if ( nap == naq ) {
      /*printf("%s: %d: On devrait rarement passer ici:",__FILE__,__LINE__);
        printf(" k=%d (%d in saveMesh), nap=%d (%d in saveMesh)\n",
        k,indElt(mesh,k),nap,indPt(mesh,nap));*/
      return(0);
    }

    /*  Point nbp, facing the second vanishing face in surface ball of p */
    jel   = lists[ilists-1] / 4;
    jface = lists[ilists-1] % 4;
    pt    = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != numq )  break;
    }
    assert(j<3);

    nro   = pt->v[i];
    jel   = lists[ilists-2] / 4;
    jface = lists[ilists-2] % 4;
    pt    = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != nro )  break;
    }
    assert(j<3);

    nbp = pt->v[i];

    /* Unfold shell of (nq,nro), starting from (jel,jface), with pivot np */
    adj = lists[ilists-1] / 4;
    piv=  nump;
    do {
      iel  = adj;
      pt   = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      /* Identification of edge number in tetra iel */
      for (i=0; i<6; i++) {
        ipa = iare[i][0];
        ipb = iare[i][1];
        if ( ((pt->v[ipa] == numq) && (pt->v[ipb] == nro)) || ((pt->v[ipa] == nro) && (pt->v[ipb] == numq))  ) break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ ifar[i][0] ] == piv ) {
        adj  = adja[ ifar[i][0] ] / 4;
        ipiv = ifar[i][1];
        iopp = ifar[i][0];
        piv  = pt->v[ipiv];
      }
      else {
        adj  = adja[ ifar[i][1] ] / 4;
        ipiv = ifar[i][0];
        iopp = ifar[i][1];
        piv  = pt->v[ipiv];
      }

      isface = 0;
      if ( pt->xt ) {
        pxt    = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && ( adj != k ) && !isface );

    nbq = piv;
    if ( nbp == nbq ) {
      /*printf("%s: %d: On devrait rarement passer ici:",__FILE__,__LINE__);
        printf(" k=%d (%d in saveMesh), nbp=%d (%d in saveMesh)\n",
        k,indElt(mesh,k),nbp,indPt(mesh,nbp));*/
      return(0);
    }
  }
  /* Surface ball has been enumerated as f1,f2,... */
  else {
    /*  Point nap, facing the fist vanishing face in surface ball of p */
    nro   = piv0;
    jel   = lists[ilists-1] / 4;
    jface = lists[ilists-1] % 4;
    pt    = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != nro )  break;
    }
    assert(j<3);

    nap = pt->v[i];

    /* Unfold shell of (nq,nro), starting from (k,iface), with pivot np */
    adj = k;
    piv = nump;
    do {
      iel = adj;
      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      /* Identification of edge number in tetra iel */
      for (i=0; i<6; i++) {
        ipa = iare[i][0];
        ipb = iare[i][1];
        if ( ((pt->v[ipa] == numq) && (pt->v[ipb] == nro)) || ((pt->v[ipa] == nro) && (pt->v[ipb] == numq))  ) break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ ifar[i][0] ] == piv ) {
        adj  = adja[ ifar[i][0] ] / 4;
        ipiv = ifar[i][1];
        iopp = ifar[i][0];
        piv  = pt->v[ipiv];
      }
      else {
        adj  = adja[ ifar[i][1] ] / 4;
        ipiv = ifar[i][0];
        iopp = ifar[i][1];
        piv  = pt->v[ipiv];
      }

      isface = 0;
      if ( pt->xt ) {
        pxt    = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && ( adj != k ) && !isface );

    naq = piv;
    if ( nap == naq ) {
      /*printf("%s: %d: On devrait rarement passer ici:",__FILE__,__LINE__);
        printf(" k=%d (%d in saveMesh), nap=%d (%d in saveMesh)\n",
        k,indElt(mesh,k),nap,indPt(mesh,nap));*/
      return(0);
    }

    /*  Point nbp, facing the second vanishing face in surface ball of p */
    jel   = lists[1] / 4;
    jface = lists[1] % 4;
    pt    = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != numq )  break;
    }
    assert(j<3);

    nro   = pt->v[i];
    jel   = lists[2] / 4;
    jface = lists[2] % 4;
    pt    = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != nro )  break;
    }
    assert(j<3);

    nbp = pt->v[i];

    /* Unfold shell of (nq,nro), starting from lists[1], with pivot np */
    adj = lists[1] / 4;
    piv =  nump;
    do {
      iel  = adj;
      pt   = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      /* Identification of edge number in tetra iel */
      for (i=0; i<6; i++) {
        ipa = iare[i][0];
        ipb = iare[i][1];
        if ( ((pt->v[ipa] == numq) && (pt->v[ipb] == nro)) || ((pt->v[ipa] == nro) && (pt->v[ipb] == numq))  ) break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ ifar[i][0] ] == piv ) {
        adj  = adja[ ifar[i][0] ] / 4;
        ipiv = ifar[i][1];
        iopp = ifar[i][0];
        piv  = pt->v[ipiv];
      }
      else {
        adj  = adja[ ifar[i][1] ] / 4;
        ipiv = ifar[i][0];
        iopp = ifar[i][1];
        piv  = pt->v[ipiv];
      }

      isface = 0;
      if ( pt->xt ) {
        pxt    = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && ( adj != k ) && !isface );

    nbq = piv;
    if ( nbp == nbq ) {
      /*printf("%s: %d: On devrait rarement passer ici:",__FILE__,__LINE__);
        printf(" k=%d (%d in saveMesh), nap=%d (%d in saveMesh)\n",
        k,indElt(mesh,k),nap,indPt(mesh,nap));*/
      return(0);
    }
  }

  return(1);
}

/** Check whether collapse ip -> iq could be performed, ip boundary point ;
 *  'mechanical' tests (positive jacobian) are not performed here ;
 *  iface = boundary face on which lie edge iedg - in local face num.
 *  (pq, or ia in local tet notation) */
int chkcol_bdy(pMesh mesh,int k,char iface,char iedg,int *listv) {
  pTetra        pt,pt0;
  pxTetra       pxt;
  pPoint        p0;
  Tria          tt;
  double        calold,calnew,caltmp,nprvold[3],nprvnew[3],ncurold[3],ncurnew[3],ps,devold,devnew;
  int           ipp,ilistv,nump,numq,ilists,lists[LMAX+2],l,iel,nbbdy,ndepmin,ndepplus;
  char          iopp,ia,ip,tag,i,iq,i0,i1,ier,isminp,isplp;

  pt   = &mesh->tetra[k];
  pxt  = 0;
  pt0  = &mesh->tetra[0];
  ia   = iarf[iface][iedg];
  ip   = idir[iface][inxt2[iedg]];
  nump = pt->v[ip];
  numq = pt->v[idir[iface][iprv2[iedg]]];
  p0   = &mesh->point[nump];
  assert(p0->tag & MG_BDY);
  assert(p0->xp);

  ndepmin = ndepplus = 0;
  isminp  = isplp = 0;

  /* collect triangles and tetras around ip */
  if ( p0->tag & MG_NOM ) {
    if ( bouleext(mesh,k,ip,iface,listv,&ilistv,lists,&ilists) < 0 )
      return(-1);
  }
  else {
    if ( boulesurfvolp(mesh,k,ip,iface,listv,&ilistv,lists,&ilists) < 0 )
      return(-1);
  }

  /* prevent collapse in case surface ball has 3 triangles */
  if ( ilists <= 2 )  return(0);  // ATTENTION, Normalement, avec 2 c est bon !

  /* Surfacic ball is enumerated with first tet having (pq) as edge n° iprv2[ip] on face iopp */
  startedgsurfball(mesh,nump,numq,lists,ilists);

  /* check tetra quality */
  calold = calnew = DBL_MAX;
  for (l=0; l<ilistv; l++) {
    iel = listv[l] / 4;
    ipp = listv[l] % 4;
    pt  = &mesh->tetra[iel];

    if ( pt->ref == MG_MINUS ) isminp = 1;
    else if ( pt->ref == MG_PLUS ) isplp = 1;

    /* Topological test for tetras of the shell */
    for (iq=0; iq<4; iq++)
      if ( pt->v[iq] == numq )  break;

    if ( iq < 4 ) {
      nbbdy = 0;
      if ( pt->xt )  pxt = &mesh->xtetra[pt->xt];
      for (i=0; i<4; i++) {
        if ( pt->xt && (pxt->ftag[i] & MG_BDY) )  nbbdy++;
      }

      /* Topological problem triggered when one of the two faces of collapsed edge is the only
         internal one : closing a part of the domain */
      if ( nbbdy == 4 )
        return(0);
      else if ( nbbdy == 3 ) {
        for (ia=0; ia<6; ia++) {
          i0 = iare[ia][0];
          i1 = iare[ia][1];
          if ( ((pt->v[i0] == nump) && (pt->v[i1] == numq)) ||
               ((pt->v[i0] == numq) && (pt->v[i1] == nump)) )
            break;
        }
        assert(ia < 6);
        i0 = ifar[ia][0];
        i1 = ifar[ia][1];
        if ( pt->xt && (!(pxt->ftag[i0] & MG_BDY) || !(pxt->ftag[i1] & MG_BDY)) )
          return(0);
      }

      /* Now check that the 2 faces identified by collapse are not boundary */
      if ( pt->xt && (pxt->ftag[ipp] & MG_BDY) && (pxt->ftag[iq] & MG_BDY) )
        return(0);

      continue;
    }

    /* Volume test for tetras outside the shell */
    if ( mesh->info.iso ) {
      if ( !ndepmin && pt->ref == MG_MINUS )
        ndepmin = iel;
      else if ( !ndepplus && pt->ref == MG_PLUS )
        ndepplus = iel;
    }

    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[ipp] = numq;

    calold = MG_MIN(calold, pt->qual);
    caltmp = orcal(mesh,0);

    if ( caltmp < EPSD )  return(0);
    calnew = MG_MIN(calnew,caltmp);
  }
  if ( calold < NULKAL && calnew <= calold )  return(0);
  else if ( calnew < NULKAL || calnew < 0.3*calold )  return(0);

  /* analyze surfacic ball of p */
  for (l=1; l<ilists-1; l++) {
    iel  = lists[l] / 4;
    iopp = lists[l] % 4;
    pt   = &mesh->tetra[iel];
    pxt = &mesh->xtetra[pt->xt];
    assert(pt->xt);

    /* retrieve vertex in tetra */
    for (ip=0; ip<4; ip++)
      if ( pt->v[ip] == nump )  break;
    assert(ip<4);

    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[ip] = numq;

    if ( !norface(mesh,iel,iopp,ncurold) )  return(0);
    if ( !norface(mesh,0,iopp,ncurnew) )    return(0);

    /* check normal flipping */
    ps = ncurold[0]*ncurnew[0] + ncurold[1]*ncurnew[1] + ncurold[2]*ncurnew[2];
    if ( ps < 0.0 )  return(0);

    /* check normal deviation */
    if ( l > 1 ) {
      ia = idirinv[iopp][ip]; /* index of p in tria iopp */
      ia = iprv2[ia];         /* edge between l-1 and l, in local num of tria */
      ia = iarf[iopp][ia];    /* edge between l-1 and l in local num of tetra */

      if ( (!pt->xt) || (!(mesh->xtetra[pt->xt].tag[ia] & MG_GEO)) ) {

        devold = nprvold[0]*ncurold[0] + nprvold[1]*ncurold[1] + nprvold[2]*ncurold[2];
        devnew = nprvnew[0]*ncurnew[0] + nprvnew[1]*ncurnew[1] + nprvnew[2]*ncurnew[2];
        if ( devold < ANGEDG ) {
          if ( devnew < devold )  return(0);
        }
        else if ( devnew < ANGEDG )  return(0);
      }
    }

    /* check Hausdorff distance to geometric support */
    tet2tri(mesh,iel,iopp,&tt);
    if ( l == 1 ) {
      for (i=0; i<3; i++) {
        if ( tt.v[i] == nump )  break;
      }
      assert(i<3);
      /* Index of the third point of the first collapsed triangle */
      i  = inxt2[i];
      ia = inxt2[i];
      tag = pxt->tag[iarf[iopp][ia]];
      tt.tag[ia] = MG_MAX(tt.tag[ia],tag);
    }
    else if ( l == ilists-2 ) {
      for (i=0; i<3; i++) {
        if ( tt.v[i] == nump )  break;
      }
      assert(i<3);
      /* Index of the third point of the first collapsed triangle */
      i  = iprv2[i];
      ia = iprv2[i];
      tag = pxt->tag[iarf[iopp][ia]];
      tt.tag[ia] = MG_MAX(tt.tag[ia],tag);
    }

    for (i=0; i<3; i++) {
      if ( tt.v[i] == nump )  break;
    }
    assert(i<3);
    tt.v[i] = numq;
    if ( chkedg(mesh,&tt,MG_GET(pxt->ori,iopp)) )  return(0);

    memcpy(nprvold,ncurold,3*sizeof(double));
    memcpy(nprvnew,ncurnew,3*sizeof(double));
  }

  /* Ensure collapse does not lead to a non manifold configuration (case of implicit surface)*/
  if ( mesh->info.iso ) {
    ier = chkmanicoll(mesh,k,iface,iedg,ndepmin,ndepplus,isminp,isplp);
    if ( !ier )  return(0);
  }
  /* Topological check for surface ball */
  else {
    ier = topchkcol_bdy(mesh,k,iface,iedg,lists,ilists);
    if ( !ier )  return(0);
  }

  return(ilistv);
}

/** Collapse vertex p = list[0]%4 of tetra list[0]/4 over vertex indq of tetra list[0]/4.
 *  Only physical tests (positive jacobian) are done (i.e. approximation of the surface,
 *  etc... must be performed outside). */
int colver(pMesh mesh,int *list,int ilist,char indq) {
  pTetra          pt,pt1;
  pxTetra         pxt,pxt1;
  xTetra          xt,xts;
  int             i,iel,jel,pel,qel,k,np,nq,*adja,p0,p1;
  unsigned char   ip,iq,j,voy,voyp,voyq,ia,iav;
  unsigned char   ind[ilist][2];
  int             p0_c[ilist],p1_c[ilist];
  char            indar[4][4][2] = {
    /* indar[ip][iq][0/1]: indices of edges which have iq for extremity but not ip*/
    { {-1,-1}, { 3, 4}, { 3, 5}, { 4, 5} },
    { { 1, 2}, {-1,-1}, { 1, 5}, { 2, 5} },
    { { 0, 2}, { 0, 4}, {-1,-1}, { 2, 4} },
    { { 0, 1}, { 0, 3}, { 1, 3}, {-1,-1} } };
#ifdef SINGUL
  int             warn_sing,need_xt;
#endif

  iel = list[0] / 4;
  ip  = list[0] % 4;
  pt  = &mesh->tetra[iel];
  np  = pt->v[ip];
  nq  = pt->v[indq];

#ifdef SINGUL
  warn_sing = 0;
#endif
  memset(p0_c,0,ilist*sizeof(int));
  memset(p1_c,0,ilist*sizeof(int));
  /* Mark elements of the shell of edge (pq) */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 4;
    i   = list[k] % 4;
    pt  = &mesh->tetra[iel];
#ifdef SINGUL
    if ( mesh->info.sing )  pt->flag = 0;
#endif

    for (j=0; j<3; j++) {
      i = inxt3[i];
      if ( pt->v[i] == nq ) {
        /* list edges that we need to update */
        if ( pt->xt ) {
          pxt = &mesh->xtetra[pt->xt];
          ip  = list[k]%4;
          ind[k][0] = indar[ip][i][0];
          if ( pxt->tag[ind[k][0]] || pxt->edg[ind[k][0]] ) {
#ifdef SINGUL
            /* Check if we need to take care about singularities */
            if ( mesh->info.sing && (pxt->tag[ind[k][0]] & MG_SGL) )  warn_sing = 1;
#endif
            if ( iare[ind[k][0]][0]==i )  p0_c[k] = pt->v[iare[ind[k][0]][1]];
            else  p0_c[k] = pt->v[iare[ind[k][0]][0]];
          }
          ind[k][1] = indar[ip][i][1];
          if ( pxt->tag[ind[k][1]] || pxt->edg[ind[k][1]] ) {
#ifdef SINGUL
            /* Check if we need to take care about singularities */
            if ( mesh->info.sing && (pxt->tag[ind[k][1]] & MG_SGL) )
              warn_sing = 1;
#endif
            if ( iare[ind[k][1]][0]==i )  p1_c[k] = pt->v[iare[ind[k][1]][1]];
            else  p1_c[k] = pt->v[iare[ind[k][1]][0]];
          }
        }
        list[k] *= -1;
        break;
      }
    }
  }

  /* avoid recreating existing elt */
  for (k=0; k<ilist; k++) {
    if ( list[k] < 0 )  continue;
    iel = list[k] / 4;
    ip  = list[k] % 4;
    pt  = &mesh->tetra[iel];

    /* update edges of elements that do not belong to the shell of pq */
    if ( !pt->xt ) {
#ifndef SINGUL
      continue;
#else
      if ( !mesh->info.sing || !warn_sing ) continue;

      /* we need a xtetra for all singular edges so we may need to create it. */
      /* First: check that xtetra will not be created by adjacency at next step */
      for ( i=0; i<ilist; i++ ) {
        if ( list[i] > 0 ) continue;
        pt1  = &mesh->tetra[(-list[i])/4];

        iq  = (-list[i]) % 4;
        for (j=0; j<3; j++) {
          iq = inxt3[iq];
          if ( pt1->v[iq] == nq )  break;
        }
        assert(j<3);

        adja = &mesh->adja[4*(iel-1)+1];
        if ( adja[iq]/4 == iel ) break;
      }
      if ( i < ilist ) {
        /* our element is adjacent to an element of the shell */
        continue;
      }
      else {
        /* we need to create the xtetra */
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,xTetra,"larger xtetra table",
                       mesh->xt--;
                       return(-1));
        }
        pt->xt = mesh->xt;
        memset(&mesh->xtetra[pt->xt],0,sizeof(xTetra));
        mesh->xtetra[pt->xt].ori = 15;
        /* mark tetra for which we have created xTetras */
        pt->flag = 1;
      }
#endif
    }
#ifdef SINGUL
    need_xt = 0;
#endif
    pxt = &mesh->xtetra[pt->xt];
    for ( i=0; i<ilist; i++ ) {
      if ( (list[i]>0) || (!(mesh->tetra[-list[i]/4].xt)) )  continue;
      pt1  = &mesh->tetra[-list[i]/4];
      pxt1 = &mesh->xtetra[pt1->xt];
      if ( p0_c[i] ) {
        for ( j=0; j<3; j++) {
          ia = idir[ip][j];
          if ( pt->v[ia]==p0_c[i] ) {
            pxt->tag[arpt[ip][j]] |= pxt1->tag[ind[i][0]];
#ifdef SINGUL
            /* we need the xtetra? */
            if ( mesh->info.sing && (pxt1->tag[ind[i][0]] & MG_SGL) )  need_xt=1;
#endif
            if ( !pxt->edg[arpt[ip][j]] )
              pxt->edg[arpt[ip][j]] = pxt1->edg[ind[i][0]];
            else if ( pxt1->edg[arpt[ip][j]] )
              pxt->edg[arpt[ip][j]] =
                MG_MAX(pxt->edg[arpt[ip][j]],pxt1->edg[ind[i][0]]);
            break;
          }
        }
      }
      if ( p1_c[i] ) {
        for ( j=0; j<3; j++) {
          ia = idir[ip][j];
          if ( pt->v[ia]==p1_c[i] ) {
            pxt->tag[arpt[ip][j]] |= pxt1->tag[ind[i][1]];
#ifdef SINGUL
            /* we need the xtetra? */
            if ( mesh->info.sing && (pxt1->tag[ind[i][1]] & MG_SGL) )  need_xt=1;
#endif
            if ( !pxt->edg[arpt[ip][j]] )
              pxt->edg[arpt[ip][j]] = pxt1->edg[ind[i][1]];
            else if ( pxt1->edg[arpt[ip][j]] )
              pxt->edg[arpt[ip][j]] =
                MG_MAX(pxt->edg[arpt[ip][j]],pxt1->edg[ind[i][1]]);
            break;
          }
        }
      }
    }
#ifdef SINGUL
    /* delete useless created xTetra */
    if ( mesh->info.sing && pt->flag && (!need_xt) ) {
      pt->xt = 0;
      mesh->xt--;
    }
    else
      need_xt = 0;
#endif
    adja = &mesh->adja[4*(iel-1)+1];
    jel  = adja[ip] / 4;
    voy  = adja[ip] % 4;
    if ( !jel )  continue;
    pt = &mesh->tetra[jel];
    if ( pt->v[voy] == nq )  return(0);
  }

  /* deal with the shell of edge (pq) and the implied updates */
  for (k=0; k<ilist; k++) {
    if ( list[k] > 0 )  continue;
    iel = (-list[k]) / 4;
    ip  = (-list[k]) % 4;
    pt  = &mesh->tetra[iel];

    iq  = ip;
    for (j=0; j<3; j++) {
      iq = inxt3[iq];
      if ( pt->v[iq] == nq )  break;
    }
    assert(j<3);

    adja = &mesh->adja[4*(iel-1)+1];

    /* pel = neighbour of iel that belongs to ball of p \setminus shell, same for qel */
    pel  = adja[iq] / 4;
    voyp = adja[iq] % 4;
    qel  = adja[ip] / 4;
    voyq = adja[ip] % 4;
    /*op = 0;
      if ( pel ) {
      pt1 = &mesh->tetra[pel];
      op  = pt1->v[voyp];
      }
      if ( qel ) {
      pt1 = &mesh->tetra[qel];
      oq  = pt1->v[voyq];
      //assert(op != oq);
      }*/

    /* Update adjacency relations */
    if ( pel ) {
      adja = &mesh->adja[4*(pel-1)+1];
      adja[voyp] = 4*qel+voyq;
    }
    if ( qel ) {
      adja = &mesh->adja[4*(qel-1)+1];
      adja[voyq] = 4*pel+voyp;
    }

    /* Update references for faces (one in pel) ;
       possibly, creation of a new field pxt for pel must be carried out */
    if ( pel ) {
      pt1 = &mesh->tetra[pel];
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        memcpy(&xts,pxt,sizeof(xTetra));
        if ( pt1->xt > 0 ) {
          pxt1 = &mesh->xtetra[pt1->xt];
          pxt1->ref[voyp] = MG_MAX(pxt1->ref[voyp],pxt->ref[ip]);
          pxt1->ftag[voyp] = pxt1->ftag[voyp] | pxt->ftag[ip];

          if ( qel && (pt1->ref < mesh->tetra[qel].ref) )  MG_CLR( pxt1->ori,voyp );
          else   MG_SET(pxt1->ori,voyp);


          /* update tags for edges */
          for ( j=0; j<3; j++ ) {
            ia = iarf[ip][j];
            p0 = pt->v[iare[ia][0]];
            p1 = pt->v[iare[ia][1]];

            for ( i=0; i<3; i++ ) {
              iav=iarf[voyp][i];
              if ( p0==nq ) {
                if ( ((pt1->v[iare[iav][0]]==np) && (pt1->v[iare[iav][1]]==p1)) ||
                     ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==np)) )
                  break;
              }
              else if ( p1==nq ) {
                if ( ((pt1->v[iare[iav][0]]==np) && (pt1->v[iare[iav][1]]==p0)) ||
                     ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==np)) )
                  break;
              }
              else {
                if ( ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==p1)) ||
                     ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==p0)) )
                  break;
              }
            }
            assert(i!=3);
            pxt1->tag[iav] = pxt1->tag[iav] | pxt->tag[ia];
          }
        }
        else {
          pxt1 = &xt;
          memset(pxt1,0,sizeof(xTetra));
          pxt1->ref[voyp] = pxt->ref[ip];
          pxt1->ftag[voyp] = pxt->ftag[ip];
          pxt1->ori = 15;
          if ( !MG_GET(pxt->ori,ip) )  MG_CLR(pxt1->ori,voyp);

          /* update tags for edges */
          for ( j=0; j<3; j++ ) {
            ia = iarf[ip][j];
            p0 = pt->v[iare[ia][0]];
            p1 = pt->v[iare[ia][1]];
            if ( pxt->tag[ia] ) {
              for ( i=0; i<3; i++ ) {
                iav=iarf[voyp][i];
                if ( p0==nq ) {
                  if ( ((pt1->v[iare[iav][0]]==np) && (pt1->v[iare[iav][1]]==p1)) ||
                       ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==np)) )
                    break;
                }
                else if ( p1==nq ) {
                  if ( ((pt1->v[iare[iav][0]]==np ) && (pt1->v[iare[iav][1]]==p0)) ||
                       ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==np )) )
                    break;
                }
                else {
                  if ( ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==p1)) ||
                       ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==p0)) )
                    break;
                }
              }
              assert(i!=3);
              pxt1->tag[iav] = pxt->tag[ia];
            }
          }
          /* Recover the already used place by pxt */
          pt1->xt = pt->xt;
          memcpy(pxt,pxt1,sizeof(xTetra));
        }
      }
      else {
        /* Only the values corresponding to pt become 0 */
        if ( pt1->xt > 0 ) {
          pxt1 = &mesh->xtetra[pt1->xt];
          pxt1->ref[voyp]  = 0;
          pxt1->ftag[voyp] = 0;
          MG_SET(pxt1->ori,voyp);
        }
      }

      if ( qel ) {
        pt1 = &mesh->tetra[qel];
        if ( pt->xt ) {
          pxt = &xts;
          if ( pt1->xt > 0 ) {
            pxt1 = &mesh->xtetra[pt1->xt];
            pxt1->ref[voyq]  = MG_MAX(pxt1->ref[voyq],pxt->ref[iq]);
            pxt1->ftag[voyq] = (pxt1->ftag[voyq] | pxt->ftag[iq]);

            if ( pel && (pt1->ref < mesh->tetra[pel].ref) )  MG_CLR( pxt1->ori,voyq );
            else   MG_SET(pxt1->ori,voyq);

            /* update tags for edges */
            for ( j=0; j<3; j++ ) {
              ia = iarf[ip][j];
              p0 = pt->v[iare[ia][0]];
              p1 = pt->v[iare[ia][1]];
              for ( i=0; i<3; i++ ) {
                iav=iarf[voyq][i];
                if ( p0==np ) {
                  if ( ((pt1->v[iare[iav][0]]==nq) && (pt1->v[iare[iav][1]]==p1)) ||
                       ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==nq)) )
                    break;
                }
                else if ( p1==np ) {
                  if ( ((pt1->v[iare[iav][0]]==nq ) && (pt1->v[iare[iav][1]]==p0)) ||
                       ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==nq )) )
                    break;
                }
                else {
                  if ( ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==p1)) ||
                       ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==p0)) )
                    break;
                }
              }
              assert(i!=3);
              pxt1->tag[iav] = pxt1->tag[iav] | pxt->tag[ia];
            }
          }
          else {
            pxt1 = &xt;
            memset(pxt1,0,sizeof(xTetra));
            pxt1->ref[voyq] = pxt->ref[iq];
            pxt1->ftag[voyq] = pxt->ftag[iq];
            pxt1->ori = 15;
            if ( !MG_GET(pxt->ori,iq) )  MG_CLR(pxt1->ori,voyq);
            /* update tags for edges */
            for ( j=0; j<3; j++ ) {
              ia = iarf[iq][j];
              p0 = pt->v[iare[ia][0]];
              p1 = pt->v[iare[ia][1]];
              if ( pxt->tag[ia] ) {
                for ( i=0; i<3; i++ ) {
                  iav=iarf[voyq][i];
                  if ( p0==np ) {
                    if ( ((pt1->v[iare[iav][0]]==nq) && (pt1->v[iare[iav][1]]==p1)) ||
                         ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==nq)) )
                      break;
                  }
                  else if ( p1==np ) {
                    if ( ((pt1->v[iare[iav][0]]==nq ) && (pt1->v[iare[iav][1]]==p0)) ||
                         ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==nq )) )
                      break;
                  }
                  else {
                    if ( ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==p1)) ||
                         ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==p0)) )
                      break;
                  }
                }
                assert(i!=3);
                pxt1->tag[iav] = pxt->tag[ia];
              }
            }
            /* Create new field xt */
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,xTetra,
                           "larger xtetra table",
                           mesh->xt--;
                           return(-1));
            }
            pt1->xt = mesh->xt;
            pxt = &mesh->xtetra[pt1->xt];
            memcpy(pxt,pxt1,sizeof(xTetra));
          }
        }
        else {
          /* Only the values corresponding to pt become 0 */
          if ( pt1->xt > 0 ) {
            pxt1 = &mesh->xtetra[pt1->xt];
            pxt1->ref[voyq]  = 0;
            pxt1->ftag[voyq] = 0;
            MG_SET(pxt1->ori,voyq);
          }
        }
      }
    }
    else {
      assert(pt->xt);
      pxt = &mesh->xtetra[pt->xt];
      if ( qel ) {
        pt1 = &mesh->tetra[qel];
        if ( pt1->xt > 0 ) {
          pxt1 = &mesh->xtetra[pt1->xt];
          pxt1->ref[voyq]  = pxt->ref[iq];
          pxt1->ftag[voyq] = pxt->ftag[iq];

          MG_SET(pxt1->ori,voyq);

          /* update tags for edges */
          for ( j=0; j<3; j++ ) {
            ia = iarf[iq][j];
            p0 = pt->v[iare[ia][0]];
            p1 = pt->v[iare[ia][1]];
            if ( pxt->tag[ia] ) {
              for ( i=0; i<3; i++ ) {
                iav=iarf[voyq][i];
                if ( p0==np ) {
                  if ( ((pt1->v[iare[iav][0]]==nq) && (pt1->v[iare[iav][1]]==p1)) ||
                       ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==nq)) )
                    break;
                }
                else if ( p1==np ) {
                  if ( ((pt1->v[iare[iav][0]]==nq ) && (pt1->v[iare[iav][1]]==p0)) ||
                       ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==nq )) )
                    break;
                }
                else {
                  if ( ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==p1)) ||
                       ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==p0)) )
                    break;
                }
              }
              assert(i!=3);
              pxt1->tag[iav] = pxt->tag[ia];
            }
          }
        }
        else {
          pxt1 = &xt;
          memset(pxt1,0,sizeof(xTetra));
          pxt1->ref[voyq]  = pxt->ref[iq];
          pxt1->ftag[voyq] = pxt->ftag[iq];
          pxt1->ori = 15;

          /* update tags for edges */
          for ( j=0; j<3; j++ ) {
            ia = iarf[iq][j];
            p0 = pt->v[iare[ia][0]];
            p1 = pt->v[iare[ia][1]];
            if ( pxt->tag[ia] ) {
              for ( i=0; i<3; i++ ) {
                iav = iarf[voyq][i];
                if ( p0==np ) {
                  if ( ((pt1->v[iare[iav][0]]==nq) && (pt1->v[iare[iav][1]]==p1)) ||
                       ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==nq)) )
                    break;
                }
                else if ( p1==np ) {
                  if ( ((pt1->v[iare[iav][0]]==nq ) && (pt1->v[iare[iav][1]]==p0)) ||
                       ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==nq )) )
                    break;
                }
                else {
                  if ( ((pt1->v[iare[iav][0]]==p0) && (pt1->v[iare[iav][1]]==p1)) ||
                       ((pt1->v[iare[iav][0]]==p1) && (pt1->v[iare[iav][1]]==p0)) )
                    break;
                }
              }
              assert(i!=3);
              pxt1->tag[iav] = pxt->tag[ia];
            }
          }
          /* Recover the already used place by pxt */
          pt1->xt = pt->xt;
          memcpy(pxt,pxt1,sizeof(xTetra));
        }
      }
    }
    delElt(mesh,iel);
  }

  /* Update vertices coordinates for elements that do not belong to the shell of (pq) */
  for (k=0; k<ilist;  k++) {
    if ( list[k] < 0 )  continue;
    iel = list[k] / 4;
    ip  = list[k] % 4;
    pt  = &mesh->tetra[iel];
    pt->v[ip] = nq;
    pt->qual=orcal(mesh,iel);
  }
  return(np);
}
