#include "mmg3d.h"

extern Info  info;
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
    if ( info.fem ) {
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
  else if ( calnew < 0.3*calold )  return(0);

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
    /*  Point nap, facing the fist vanishing face in surface ball of p */
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
      /* printf("On devrait rarement passer ici\n"); */
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
      /* printf("On devrait rarement passer ici\n"); */
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
      /* printf("On devrait rarement passer ici\n"); */
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
      /* printf("On devrait rarement passer ici\n"); */
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
  int           ipp,ilistv,nump,numq,ilists,lists[LMAX+2],l,iel,ref,nbbdy,ndepmin,ndepplus;
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
  if ( p0->tag & MG_NOM )
    bouleext(mesh,k,ip,iface,listv,&ilistv,lists,&ilists);
  else
    boulesurfvolp(mesh,k,ip,iface,listv,&ilistv,lists,&ilists);

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
    if ( info.iso ) {
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
  else if ( calnew < 0.3*calold )  return(0);

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

      hGet(&mesh->htab,pt->v[iare[ia][0]],pt->v[iare[ia][1]],&ref,&tag);
      if ( !(tag & MG_GEO) ) {
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
      hGet(&mesh->htab,tt.v[i],numq,&ref,&tag);
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
      hGet(&mesh->htab,tt.v[i],numq,&ref,&tag);
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
  if ( info.iso ) {
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
  int             iel,jel,pel,qel,k,np,nq,*adja,ref,p0,p1,up;
  unsigned char   ip,iq,i,j,voy,voyp,voyq,ia,iav;
  char            tag;

  iel = list[0] / 4;
  ip  = list[0] % 4;
  pt  = &mesh->tetra[iel];
  np  = pt->v[ip];
  nq  = pt->v[indq];
  /* flag to know if we need to update tag in xtetra */
  up = 0;
  if ( pt->xt ) {
    if ( mesh->xtetra[pt->xt].tag[iarf[ip][0]] ||
         mesh->xtetra[pt->xt].tag[iarf[ip][1]] ||
         mesh->xtetra[pt->xt].tag[iarf[ip][2]] ) {
      up = 1;
    }
  }

  /* Mark elements of the shell of edge (pq) */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 4;
    i   = list[k] % 4;
    pt  = &mesh->tetra[iel];

    for (j=0; j<3; j++) {
      i = inxt3[i];
      if ( pt->v[i] == nq ) {
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

    /* Update references for edges (pa)->(qa) when pqa is a face of the mesh */
    for (j=0; j<4; j++) {
      if ( j == ip || j == iq )  continue;
      hPop(&mesh->htab,np,pt->v[j],&ref,&tag);
      if ( tag || ref )
        hEdge(&mesh->htab,nq,pt->v[j],ref,tag);
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
          if ( !MG_GET(pxt->ori,ip) && MG_GET(pxt1->ori,voyp) )  MG_CLR(pxt1->ori,voyp);
          else if ( pxt->ftag[ip] && !MG_GET(pxt1->ori,voyp) )   MG_SET(pxt1->ori,voyp);

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
          pxt1->ori = 15;
        }
      }

      if ( qel ) {
        pt1 = &mesh->tetra[qel];
        if ( pt->xt ) {
          pxt = &xts;
          if ( pt1->xt > 0 ) {
            pxt1 = &mesh->xtetra[pt1->xt];
            pxt1->ref[voyq] = MG_MAX(pxt1->ref[voyq],pxt->ref[iq]);
            pxt1->ftag[voyq] = (pxt1->ftag[voyq] | pxt->ftag[iq]);
            if ( !pxt1->ftag[voyq] && !MG_GET(pxt->ori,iq) )  MG_CLR(pxt1->ori,voyq);

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
            pxt1->ori = 15;
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
          pxt1->ref[voyq] = pxt->ref[iq];
          pxt1->ftag[voyq] = pxt->ftag[iq];

          if ( !pxt1->ftag[voyq] && !MG_GET(pxt->ori,iq) )  MG_CLR(pxt1->ori,voyq);

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
    for (j=0; j<4; j++) {
      if ( j==ip )  continue;
      hPop(&mesh->htab,np,pt->v[j],&ref,&tag);
      if ( tag || ref )
        hEdge(&mesh->htab,nq,pt->v[j],ref,tag);
      /* if needed update tetra tag */
      if ( up && pt->xt ) {
        if ( hGet(&mesh->htab,nq,pt->v[j],&ref,&tag) && tag ) {
          for ( ia=0; ia<6; ia++ ) {
            p0 = pt->v[iare[ia][0]];
            p1 = pt->v[iare[ia][1]];
            if ( (p0==np && p1==pt->v[j]) || (p0==pt->v[j] && p1==np) )
              break;
          }
          assert(ia<6);
          mesh->xtetra[pt->xt].tag[ia] |= tag;
        }
      }
    }
    pt->v[ip] = nq;
    pt->qual=orcal(mesh,iel);
  }

  if ( mesh->point[np].tag & MG_BDY )
    hPop(&mesh->htab,np,nq,&ref,&tag);

  delPt(mesh,np);

  return(1);
}
