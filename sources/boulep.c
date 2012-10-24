#include "mmg3d.h"

extern char ddb;

/* return average normal of triangles sharing P without crossing ridge */
int boulen(pMesh mesh,int start,int ip,double *nn) {
  pTria    pt;
  double   n[3],dd;
  int     *adja,k;
  char     i,i1,i2;

  pt = &mesh->tria[start];
  if ( !MG_EOK(pt) )  return(0);
  nn[0] = nn[1] = nn[2] = 0.0;

  /* store neighbors */
  k  = start;
  i  = ip;
  i1 = inxt2[i];
  do {
    pt = &mesh->tria[k];
    nortri(mesh,pt,n);
    nn[0] += n[0];  nn[1] += n[1];  nn[2] += n[2];

    if ( pt->tag[i1] & MG_GEO ) {
      k = 0;
      break;
    }
    adja = &mesh->adjt[3*(k-1)+1];
    k  = adja[i1] / 3;
    i2 = adja[i1] % 3;
    i1 = iprv2[i2];
  }
  while ( k && k != start );

  if ( k == 0 ) {
    k  = start;
    i  = ip;
    i2 = iprv2[i];
    pt = &mesh->tria[k];
    do {
      if ( pt->tag[i2] & MG_GEO )  break;

      adja = &mesh->adjt[3*(k-1)+1];
      k  = adja[i2] / 3;
      if ( k == 0 )  break;
      i1 = adja[i2] % 3;
      i2 = inxt2[i1];
      pt = &mesh->tria[k];

      nortri(mesh,pt,n);

      nn[0] += n[0];  nn[1] += n[1];  nn[2] += n[2];
    }
    while ( k && k != start );
  }

  /* normalize */
  dd = nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2];
  if ( dd > EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    nn[0] *= dd;
    nn[1] *= dd;
    nn[2] *= dd;
    return(1);
  }

  return(0);
}


/* return tangent to curve at ip */
int boulec(pMesh mesh,int start,int ip,double *tt) {
  pTria    pt;
  pPoint   p0,p1,p2;
  double   dd;
  int     *adja,k;
  char     i,i1,i2;

  pt = &mesh->tria[start];
  if ( !MG_EOK(pt) )       return(0);
  p0 = &mesh->point[pt->v[ip]];
  if ( !MG_EDG(p0->tag) )  return(0);

  /* check other triangle vertices */
  k  = start;
  i  = ip;
  i1 = inxt2[i];
  i2 = iprv2[i];
  p1 = p2 = 0;
  do {
    pt = &mesh->tria[k];
    if ( MG_EDG(pt->tag[i1]) ) {
      p1 = &mesh->point[pt->v[i2]];
      k  = 0;
      break;
    }
    adja = &mesh->adjt[3*(k-1)+1];
    k  = adja[i1] / 3;
    i2 = adja[i1] % 3;
    i1 = iprv2[i2];
  }
  while ( k && k != start );

  /* check if open boundary hit */
  if ( k == 0 ) {
    k  = start;
    i  = ip;
    i1 = inxt2[i];
    i2 = iprv2[i];
    do {
      pt = &mesh->tria[k];
      if ( MG_EDG(pt->tag[i2]) ) {
        p2 = &mesh->point[pt->v[i1]];
        break;
      }
      adja = &mesh->adjt[3*(k-1)+1];
      k  = adja[i2] / 3;
      i1 = adja[i2] % 3;
      i2 = inxt2[i1];
    }
    while ( k );
  }

  if ( !p1 || !p2 )
    return(0);
  else if ( p1 == p2 )
    p2 = p0;

  /* tangent approx */
  tt[0] = p2->c[0] - p1->c[0];
  tt[1] = p2->c[1] - p1->c[1];
  tt[2] = p2->c[2] - p1->c[2];
  dd = tt[0]*tt[0] + tt[1]*tt[1] + tt[2]*tt[2];
  if ( dd > EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    tt[0] *= dd;
    tt[1] *= dd;
    tt[2] *= dd;
  }

  return(1);
}


/* store edges and return number (ref+geo) incident to ip */
int bouler(pMesh mesh,int start,int ip,int *list) {
  pTria    pt;
  int     *adja,k,nr;
  char     i,i1,i2;

  pt  = &mesh->tria[start];
  if ( !MG_EOK(pt) )  return(0);

  /* check other triangle vertices */
  k  = start;
  i  = ip;
  nr = 0;
  do {
    i1 = inxt2[i];
    if ( MG_EDG(pt->tag[i1])) {
      i2 = iprv2[i];
      nr++;
      list[nr] = pt->v[i2];
      if ( nr > LMAX-2 )  return(-nr);
    }
    adja = &mesh->adjt[3*(k-1)+1];
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = inxt2[i];
    pt = &mesh->tria[k];
  }
  while ( k && k != start );

  /* reverse loop */
  if ( k != start ) {
    k = start;
    i = ip;
    do {
      pt = &mesh->tria[k];
      i2 = iprv2[i];
      if ( MG_EDG(pt->tag[i2]) ) {
        i1 = inxt2[i];
        nr++;
        list[nr] = pt->v[i1];
        if ( nr > LMAX-2 )  return(-nr);
      }
      adja = &mesh->adjt[3*(k-1)+1];
      k  = adja[i2] / 3;
      i  = adja[i2] % 3;
      i  = iprv2[i];
    }
    while ( k && k != start );
  }

  return(nr);
}

/* Return volumic ball (i.e. filled with tetrahedra) of point ip in tetra start.
   Results are stored under the form 4*kel + jel , kel = number of the tetra, jel = local
   index of p within kel */
int boulevolp (pMesh mesh, int start, int ip, int * list){
  pTetra  pt,pt1;
  int    *adja,nump,ilist,base,cur,k,k1;
  char    j,l,i;

  base = ++mesh->base;
  pt   = &mesh->tetra[start];
  nump = pt->v[ip];
  ilist = 0;

  /* Store initial tetrahedron */
  pt->flag = base;
  list[ilist] = 4*start + ip;
  ilist++;

  /* Explore list and travel by adjacency through elements sharing p */
  cur = 0;
  while ( cur < ilist ) {
    k = list[cur] / 4;
    i = list[cur] % 4; // index of point p in tetra k
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = inxt3[i];
      k1 = adja[i] / 4;
      if ( !k1 )  continue;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;
      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);
      /* overflow */
      if ( ilist > LMAX-3 )  return(0);
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }
  return(ilist);
}

/* Define normal and tangent vectors at a non manifold point (ip in start, supported by
   face iface), enumerating its (outer)surfacic ball ; return sng = whether point is singular
   or not */
int boulenm(pMesh mesh,int start,int ip,int iface,double n[3],double t[3]) {
  pTetra   pt;
  pPoint   p0,p1,ppt;
  double   dd,nt[3],l0,l1;
  int      base,nump,nr,nnm,k,piv,na,nb,adj,nvstart,fstart,aux,ref,ip0,ip1;
  int     *adja;
  char     iopp,ipiv,i,ipa,ipb,isface,tag;

  base = ++mesh->base;
  nr  = nnm = 0;
  ip0 = ip1 = 0;

  memset(n,0.0,3*sizeof(double));
  memset(t,0.0,3*sizeof(double));

  pt   = &mesh->tetra[start];
  nump = pt->v[ip];
  k    = start;

  na  = pt->v[ip];
  nb  = pt->v[idir[iface][inxt2[idirinv[iface][ip]]]];
  piv = pt->v[idir[iface][iprv2[idirinv[iface][ip]]]];

  iopp   = iface;
  fstart = 4*k+iopp;
  do {
    /* computation of normal and tangent at nump */
    if ( norface(mesh,k,iopp,nt) ) {
      n[0] += nt[0];
      n[1] += nt[1];
      n[2] += nt[2];
    }
    hGet(&mesh->htab,na,nb,&ref,&tag);
    if ( MG_EDG(tag) && !(tag & MG_NOM) )
      nr++;
    else if ( tag & MG_NOM ) {
      nnm++;
      if ( !ip0 )
        ip0 = nb;
      else
        ip1 = nb;
    }

    /* A boundary face has been hit : change travel edge */
    aux = nb;
    nb = piv;
    piv = aux;
    nvstart = k;
    adj = k;

    /* Now unfold shell of edge (na,nb) starting from k (included) */
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      if ( pt->flag != base ) {
        for (i=0; i<4; i++)
          if ( pt->v[i] == nump )  break;
        assert(i<4);
        pt->flag = base;
      }

      /* identification of edge number in tetra k */
      for (i=0; i<6; i++) {
        ipa = iare[i][0];
        ipb = iare[i][1];
        if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
             (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ ifar[i][0] ] == piv ) {
        adj = adja[ ifar[i][0] ] / 4;
        ipiv = ifar[i][1];
        iopp = ifar[i][0];
        piv = pt->v[ipiv];
      }
      else {
        adj = adja[ ifar[i][1] ] / 4;
        ipiv = ifar[i][0];
        iopp = ifar[i][1];
        piv = pt->v[ipiv];
      }
      isface = (adja[iopp] == 0);
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  while ( 4*k+iopp != fstart );

  if ( (nr > 0 && nnm > 0) || nnm != 2 )  return(0);

  dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( dd > EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    n[0] *= dd;
    n[1] *= dd;
    n[2] *= dd;
  }
  assert( ip0 && ip1 && (ip0 != ip1) );

  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];
  ppt = &mesh->point[nump];

  l0 = (ppt->c[0] - p0->c[0])*(ppt->c[0] - p0->c[0]) \
    + (ppt->c[1] - p0->c[1])*(ppt->c[1] - p0->c[1]) + (ppt->c[2] - p0->c[2])*(ppt->c[2] - p0->c[2]);
  l1 = (ppt->c[0] - p1->c[0])*(ppt->c[0] - p1->c[0]) \
    + (ppt->c[1] - p1->c[1])*(ppt->c[1] - p1->c[1]) + (ppt->c[2] - p1->c[2])*(ppt->c[2] - p1->c[2]);
  l0 = sqrt(l0);
  l1 = sqrt(l1);

  if ( (l0 < EPSD2) || (l1 < EPSD2) ) {
    t[0] = p1->c[0] - p0->c[0];
    t[1] = p1->c[1] - p0->c[1];
    t[2] = p1->c[2] - p0->c[2];
  }
  else if ( l0 < l1 ) {
    dd = l0 / l1;
    t[0] = dd*(p1->c[0] - ppt->c[0]) + ppt->c[0] - p0->c[0];
    t[1] = dd*(p1->c[1] - ppt->c[1]) + ppt->c[1] - p0->c[1];
    t[2] = dd*(p1->c[2] - ppt->c[2]) + ppt->c[2] - p0->c[2];
  }
  else {
    dd = l1 / l0;
    t[0] = dd*(p0->c[0] - ppt->c[0]) + ppt->c[0] - p1->c[0];
    t[1] = dd*(p0->c[1] - ppt->c[1]) + ppt->c[1] - p1->c[1];
    t[2] = dd*(p0->c[2] - ppt->c[2]) + ppt->c[2] - p1->c[2];
  }
  dd = t[0]*n[0] + t[1]*n[1] + t[2]*n[2];
  t[0] -= dd*n[0];
  t[1] -= dd*n[1];
  t[2] -= dd*n[2];

  dd = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
  if ( dd > EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    t[0] *= dd;
    t[1] *= dd;
    t[2] *= dd;
  }

  return(1);
}

/* Return volumic ball of a surfacic point p, as well as the part of its surfacic ball
   supported in the outer boundary starting from tet start, with point ip, and face if in tetra
   volumic ball ; list[k] = 4*number of tet + index of point
   surfacic ball : list[k] = 4*number of tet + index of FACE */
int bouleext(pMesh mesh, int start, int ip, int iface, int *listv, int *ilistv, int *lists, int*ilists){
  pTetra pt,pt1;
  int base,nump,k,k1,*adja,piv,na,nb,adj,cur,nvstart,fstart,aux;
  char iopp,ipiv,i,j,l,ipa,ipb,isface;

  base = ++mesh->base;
  *ilists = 0;
  *ilistv = 0;

  pt = &mesh->tetra[start];
  nump = pt->v[ip];
  k = start;

  na = pt->v[ip];
  nb = pt->v[idir[iface][inxt2[idirinv[iface][ip]]]];
  piv = pt->v[idir[iface][iprv2[idirinv[iface][ip]]]];

  iopp = iface;
  fstart = 4*k+iopp;

  do {
    /* A boundary face has been hit : change travel edge */
    lists[(*ilists)] = 4*k+iopp;
    (*ilists)++;
    assert(*ilists < LMAX);

    aux = nb;
    nb = piv;
    piv = aux;
    nvstart = k;
    adj = k;

    /* Now unfold shell of edge (na,nb) starting from k (included)*/
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      if ( pt->flag != base ) {
        for (i=0; i<4; i++)
          if ( pt->v[i] == nump )  break;
        assert(i<4);
        listv[(*ilistv)] = 4*k+i;
        (*ilistv)++;
        pt->flag = base;
      }

      /* identification of edge number in tetra k */
      for (i=0; i<6; i++) {
        ipa = iare[i][0];
        ipb = iare[i][1];
        if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
             (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ ifar[i][0] ] == piv ) {
        adj = adja[ ifar[i][0] ] / 4;
        ipiv = ifar[i][1];
        iopp = ifar[i][0];
        piv = pt->v[ipiv];
      }
      else {
        adj = adja[ ifar[i][1] ] / 4;
        ipiv = ifar[i][0];
        iopp = ifar[i][1];
        piv = pt->v[ipiv];
      }
      isface = (adja[iopp] == 0);
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  while ( 4*k+iopp != fstart );

  /* Now, surfacic ball is complete ; finish travel of volumic ball */
  cur = 0;
  while ( cur < (*ilistv) ) {
    k = listv[cur]/4;
    i = listv[cur]%4;
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = inxt3[i];
      k1 = adja[i]/4;
      if ( !k1 )  continue;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;

      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);
      /* overflow */
      assert ( (*ilistv) <= LMAX-3 );
      listv[(*ilistv)] = 4*k1+j;
      (*ilistv)++;
    }
    cur++;
  }

  return(1);
}

/* Return volumic ball of a SURFACE point p, as well as its surfacic ball, starting from tetra
   start, with point ip, and face if in tetra
   volumic ball ; list[k] = 4*number of tet + index of point
   surfacic ball : list[k] = 4*number of tet + index of FACE */
int boulesurfvolp(pMesh mesh,int start,int ip,int iface,int *listv,int *ilistv,int *lists,int*ilists) {
  pTetra pt,pt1;
  pxTetra pxt;
  int base,nump,k,k1,*adja,piv,na,nb,adj,cur,nvstart,fstart,aux;
  char iopp,ipiv,i,j,l,ipa,ipb,isface;

  base = ++mesh->base;
  *ilists = 0;
  *ilistv = 0;

  pt = &mesh->tetra[start];
  nump = pt->v[ip];
  k = start;

  na = pt->v[ip];
  nb = pt->v[idir[iface][inxt2[idirinv[iface][ip]]]];
  piv = pt->v[idir[iface][iprv2[idirinv[iface][ip]]]];

  iopp = iface;
  fstart = 4*k+iopp;

  do {
    /* A boundary face has been hit : change travel edge */
    lists[(*ilists)] = 4*k+iopp;
    (*ilists)++;
    assert(*ilists < LMAX);

    aux = nb;
    nb = piv;
    piv = aux;
    nvstart = k;
    adj = k;

    /* Now unfold shell of edge (na,nb) starting from k (included)*/
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      if ( pt->flag != base ) {
        for (i=0; i<4; i++)
          if ( pt->v[i] == nump )  break;
        assert(i<4);
        listv[(*ilistv)] = 4*k+i;
        (*ilistv)++;
        pt->flag = base;
      }

      /* identification of edge number in tetra k */
      for (i=0; i<6; i++) {
        ipa = iare[i][0];
        ipb = iare[i][1];
        if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
             (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ ifar[i][0] ] == piv ) {
        adj = adja[ ifar[i][0] ] / 4;
        ipiv = ifar[i][1];
        iopp = ifar[i][0];
        piv = pt->v[ipiv];
      }
      else {
        adj = adja[ ifar[i][1] ] / 4;
        ipiv = ifar[i][0];
        iopp = ifar[i][1];
        piv = pt->v[ipiv];
      }
      isface = 0;
      if(pt->xt){
        pxt = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  while ( 4*k+iopp != fstart );

  /* Now, surfacic ball is complete ; finish travel of volumic ball */
  cur = 0;  // Check numerotation
  while ( cur < (*ilistv) ) {
    k = listv[cur]/4;
    i = listv[cur]%4; // index of point p in tetra k
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = inxt3[i];
      k1 = adja[i]/4;
      if ( !k1 )  continue;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;

      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);
      /* overflow */
      assert ( (*ilistv) <= LMAX-3 );
      listv[(*ilistv)] = 4*k1+j;
      (*ilistv)++;
    }
    cur++;
  }

  return(1);
}

/* Get tag of edge ia in tetra start by travelling its shell until meeting a boundary face */
inline int gettag(pMesh mesh,int start,int ia,int *tag,int *edg) {
  pTetra        pt;
  pxTetra       pxt;
  int           na,nb,*adja,adj,piv;
  unsigned char i,ipa,ipb;

  if ( start < 1 )  return(0);
  pt = &mesh->tetra[start];
  if ( !MG_EOK(pt) ) return(0);

  na   = pt->v[ iare[ia][0] ];
  nb   = pt->v[ iare[ia][1] ];

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[ifar[ia][0]] / 4;
  piv = pt->v[ifar[ia][1]];

  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    if ( (pxt->ftag[ifar[ia][0]] & MG_BDY) || (pxt->ftag[ifar[ia][1]] & MG_BDY) ) {
      *tag = pxt->tag[ia];
      *edg = pxt->edg[ia];
      return(1);
    }
  }

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = iare[i][0];
      ipb = iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      if ( (pxt->ftag[ifar[i][0]] & MG_BDY) || (pxt->ftag[ifar[i][1]] & MG_BDY) ) {
        *tag = pxt->tag[i];
        *edg = pxt->edg[i];
        return(1);
      }
    }

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ ifar[i][0] ] == piv ) {
      adj = adja[ ifar[i][0] ] / 4;
      piv = pt->v[ ifar[i][1] ];
    }
    else {
      adj = adja[ ifar[i][1] ] /4;
      piv = pt->v[ ifar[i][0] ];
    }
  }
  return(0);
}

/* Set tag and edg of edge ia (if need be) in tetra start by travelling its shell */
inline int settag(pMesh mesh,int start,int ia,int tag,int edg) {
  pTetra        pt;
  pxTetra       pxt;
  int           na,nb,*adja,adj,piv;
  unsigned char i,ipa,ipb;

  if ( start < 1 )  return(0);
  pt = &mesh->tetra[start];
  if ( !MG_EOK(pt) ) return(0);

  na   = pt->v[ iare[ia][0] ];
  nb   = pt->v[ iare[ia][1] ];

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[ifar[ia][0]] / 4;
  piv = pt->v[ifar[ia][1]];

  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    pxt->tag[ia] = tag;
    pxt->edg[ia] = edg;
  }
  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = iare[i][0];
      ipb = iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      pxt->tag[i] = tag;
      pxt->edg[i] = edg;
    }
    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ ifar[i][0] ] == piv ) {
      adj = adja[ ifar[i][0] ] / 4;
      piv = pt->v[ ifar[i][1] ];
    }
    else {
      adj = adja[ ifar[i][1] ] /4;
      piv = pt->v[ ifar[i][0] ];
    }
  }

  /* If all shell has been travelled, stop, else, travel it the other sense */
  if ( adj == start )  return(1);
  assert(!adj);

  pt = &mesh->tetra[start];
  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[ifar[ia][1]] / 4;
  piv = pt->v[ifar[ia][0]];

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = iare[i][0];
      ipb = iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      pxt->tag[i] = tag;
      pxt->edg[i] = edg;
    }
    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ ifar[i][0] ] == piv ) {
      adj = adja[ ifar[i][0] ] / 4;
      piv = pt->v[ ifar[i][1] ];
    }
    else {
      adj = adja[ ifar[i][1] ] /4;
      piv = pt->v[ ifar[i][0] ];
    }
  }
  return(1);
}

/* Find all tets sharing edge ia of tetra start
   return 2*ilist if shell is closed, 2*ilist +1 otherwise */
int coquil(pMesh mesh,int start,int ia,int * list) {
  pTetra  pt;
  int     ilist,*adja,piv,adj,na,nb,ipa,ipb;
  char    i;

  if ( start < 1 )  return(0);
  pt = &mesh->tetra[start];
  if ( !MG_EOK(pt) ) return(0);

  na   = pt->v[ iare[ia][0] ];
  nb   = pt->v[ iare[ia][1] ];
  ilist = 0;
  list[ilist] = 6*start+ia;
  ilist++;

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[ifar[ia][0]] / 4; // start travelling by face (ia,0)
  piv = pt->v[ifar[ia][1]];

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = iare[i][0];
      ipb = iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    list[ilist] = 6*adj +i;
    ilist++;
    /* overflow */
    assert( ilist <= LMAX-3 );

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ ifar[i][0] ] == piv ) {
      adj = adja[ ifar[i][0] ] / 4;
      piv = pt->v[ ifar[i][1] ];
    }
    else {
      assert(pt->v[ ifar[i][1] ] == piv );
      adj = adja[ ifar[i][1] ] /4;
      piv = pt->v[ ifar[i][0] ];
    }
  }

  /* At this point, the first travel, in one direction, of the shell is complete. Now, analyze why
     the travel ended. */
  if ( adj == start )  return(2*ilist);
  assert(!adj); // a boundary has been detected

  adj = list[ilist-1] / 6;
  i   = list[ilist-1] % 6;
  ilist = 0;

  /* Start back everything from this tetra adj */
  list[ilist] = 6*adj + i;
  ilist++;
  /* overflow */
  assert( ilist <= LMAX-3 );

  adja = &mesh->adja[4*(adj-1)+1];
  if ( pt->v[ ifar[i][0] ] == piv ) {
    adj = adja[ ifar[i][0] ] / 4;
    piv = pt->v[ ifar[i][1] ];
  }
  else {
    adj = adja[ ifar[i][1] ] /4;
    piv = pt->v[ ifar[i][0] ];
  }

  while ( adj ) {
    pt = &mesh->tetra[adj];
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = iare[i][0];
      ipb = iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    list[ilist] = 6*adj +i;
    ilist++;
    /* overflow */
    assert( ilist <= LMAX-2 );

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ ifar[i][0] ] == piv ) {
      adj = adja[ ifar[i][0] ] / 4;
      piv = pt->v[ ifar[i][1] ];
    }
    else {
      adj = adja[ ifar[i][1] ] /4;
      piv = pt->v[ ifar[i][0] ];
    }
  }
  assert(!adj);
  return(2*ilist+1);
}

/* Identify whether edge ia in start is a boundary edge by unfolding its shell */
int srcbdy(pMesh mesh,int start,int ia) {
  pTetra      pt;
  pxTetra     pxt;
  int         na,nb,adj,piv,*adja;
  char        ipa,ipb,iadj,i;

  pt = &mesh->tetra[start];
  na = pt->v[iare[ia][0]];
  nb = pt->v[iare[ia][1]];

  adja = &mesh->adja[4*(start-1)+1];
  iadj = ifar[ia][0];

  if(pt->xt){
    pxt = &mesh->xtetra[pt->xt];
    if( pxt->ftag[iadj] & MG_BDY )
      return(1);
  }

  adj = adja[iadj] / 4;
  piv = pt->v[ifar[ia][1]];

  while( adj && ( adj != start ) ) {
    pt = &mesh->tetra[adj];

    /* identification of edge number in tetra adj */
    for(i=0; i<6; i++) {
      ipa = iare[i][0];
      ipb = iare[i][1];
      if( ( pt->v[ipa] == na && pt->v[ipb] == nb ) || ( pt->v[ipa] == nb && pt->v[ipb] == na ))
        break;

    }
    assert(i<6);

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ ifar[i][0] ] == piv ) {
      iadj = ifar[i][0];
      adj = adja[ iadj ] / 4;
      piv = pt->v[ ifar[i][1] ];
    }
    else {
      iadj = ifar[i][1];
      adj = adja[ iadj ] /4;
      piv = pt->v[ ifar[i][0] ];
    }

    if(pt->xt){
      pxt = &mesh->xtetra[pt->xt];
      if( pxt->ftag[iadj] & MG_BDY )
        return(1);
    }
  }

  return(0);
}

/* Find all tets sharing edge ia of tetra start, and stores boundary faces when met
   it1 & it2 = 6*iel + iface, iel = index of tetra, iface = index of face in tetra
   return 2*ilist if shell is closed, 2*ilist +1 otherwise */
int coquilface(pMesh mesh,int start,int ia,int *list,int *it1,int *it2) {
  pTetra   pt;
  pxTetra  pxt;
  int     *adja,piv,adj,na,nb,ipa,ipb,ilist,pradj;
  char     i,iface,isbdy;

  if ( start < 1 )  return(0);
  pt = &mesh->tetra[start];
  if ( !MG_EOK(pt) ) return(0);

  na   = pt->v[ iare[ia][0] ];
  nb   = pt->v[ iare[ia][1] ];

  ilist = 0;
  list[ilist] = 6*start+ia;
  ilist++;

  *it1 = 0;
  *it2 = 0;

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[ifar[ia][0]] / 4;
  piv = pt->v[ifar[ia][1]];

  pxt = 0;
  if ( pt->xt )
    pxt = &mesh->xtetra[pt->xt];

  iface = ifar[ia][1];
  isbdy = pt->xt ? pxt->ftag[iface] : 0;
  if ( isbdy )
    *it1 = 4*start + iface;

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    pxt = 0;
    if ( pt->xt )
      pxt = &mesh->xtetra[pt->xt];

    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = iare[i][0];
      ipb = iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    list[ilist] = 6*adj +i;
    ilist++;
    /* overflow */
    assert ( ilist <= LMAX-2 );

    /* set new tetra for travel */
    pradj = adj;
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ ifar[i][0] ] == piv ) {
      adj = adja[ ifar[i][0] ] / 4;
      piv = pt->v[ ifar[i][1] ];
      iface = ifar[i][1];
    }
    else {
      assert(pt->v[ ifar[i][1] ] == piv );
      adj = adja[ ifar[i][1] ] /4;
      piv = pt->v[ ifar[i][0] ];
      iface = ifar[i][0];
    }
    isbdy = pt->xt ? pxt->ftag[iface] : 0;

    if ( isbdy ) {
      if ( *it1 == 0 )
        *it1 = 4*pradj+iface;
      else {
        if ( *it2 != 0 ) {
          return(0);
          //saveMesh(mesh);
        }
        assert( *it2 == 0 );
        *it2 = 4*pradj+iface;
      }
    }
  }

  /* At this point, the first travel, in one direction, of the shell is complete. Now, analyze why
     the travel ended. */
  if ( adj == start ) {
    assert(*it1 && *it2);
    assert(*it1 != *it2);
    return(2*ilist);
  }

  /* A boundary has been detected : slightly different configuration */
  assert(!adj);
  adj = list[ilist-1] / 6;
  i   = list[ilist-1] % 6;
  ilist = 0;

  /* Start back everything from this tetra adj */
  pradj = adj;
  /* overflow */
  assert(ilist <= LMAX-2);
  pxt = 0;
  if ( pt->xt )
    pxt = &mesh->xtetra[pt->xt];

  adja = &mesh->adja[4*(adj-1)+1];
  if ( pt->v[ ifar[i][0] ] == piv ) {
    iface = ifar[i][1];
  }
  else {
    iface = ifar[i][0];
  }
  isbdy = pt->xt ? pxt->ftag[iface] : 0;
  if ( isbdy ) {
    *it1 = 4*pradj + iface;
  }

  while ( adj ) {
    pt = &mesh->tetra[adj];

    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = iare[i][0];
      ipb = iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na) )  break;
    }
    assert(i<6);
    list[ilist] = 6*adj +i;
    ilist++;
    /* overflow */
    assert ( ilist <= LMAX-2 );

    /* set new tetra for travel */
    pradj = adj;
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ ifar[i][0] ] == piv ) {
      adj = adja[ ifar[i][0] ] / 4;
      piv = pt->v[ ifar[i][1] ];
      iface = ifar[i][0];
    }
    else {
      adj = adja[ ifar[i][1] ] /4;
      piv = pt->v[ ifar[i][0] ];
      iface = ifar[i][1];
    }
  }

  assert(!adj);
  *it2 = 4*pradj + iface;

  assert( *it1 && *it2);
  assert( *it1 != *it2);
  return(2*ilist+1);
}
