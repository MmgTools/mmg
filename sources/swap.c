#include "mmg3d.h"

extern Info info;
extern char ddb;

/** Check whether edge whose shell is provided should be swapped for
    geometric approximation purposes (the 2 surface triangles are also provided) */
int chkswpbdy(pMesh mesh,int *list,int ilist,int it1,int it2) {
  pTetra   pt,pt0;
  pxTetra  pxt;
  pPoint   p0,p1,ppt0;
  Tria     tt1,tt2;
  double   b0[3],b1[3],v[3],c[3],ux,uy,uz,ps,disnat,dischg,cal1,cal2,calnat,calchg,calold,calnew,caltmp;
  int      iel,iel1,iel2,np,nq,na1,na2,k,nminus,nplus;
  char     ifa1,ifa2,ia,ip,iq,ia1,ia2,j,isshell;

  iel = list[0] / 6;
  ia  = list[0] % 6;
  pt  = &mesh->tetra[iel];
  pt0 = &mesh->tetra[0];
  ppt0= &mesh->point[0];

  np = pt->v[iare[ia][0]];
  nq = pt->v[iare[ia][1]];

  /* No swap of geometric edge */
  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    if ( (pxt->edg[ia]>0) || MG_EDG(pxt->tag[ia]) || MG_SIN(pxt->tag[ia]) ||
         (pxt->tag[ia] & MG_NOM) )  return(0);
  }

  /* No swap when either internal or external component has only 1 element */
  if ( info.iso ) {
    nminus = nplus = 0;
    for (k=0; k<ilist; k++) {
      iel = list[k] / 6;
      pt = &mesh->tetra[iel];
      if ( pt->ref == MG_MINUS )
        nminus++;
      else
        nplus++;
    }
    if ( nplus == 1 || nminus == 1 )  return(0);
  }
  iel1 = it1 / 4;
  ifa1 = it1 % 4;
  iel2 = it2 / 4;
  ifa2 = it2 % 4;
  tet2tri(mesh,iel1,ifa1,&tt1);
  tet2tri(mesh,iel2,ifa2,&tt2);

  for (ia1=0; ia1<3; ia1++) {
    if ( (tt1.v[ia1] != np) && (tt1.v[ia1] != nq) )  break;
  }
  assert( ia1 < 3 );
  assert( (tt1.v[inxt2[ia1]] == np && tt1.v[iprv2[ia1]] == nq) ||
          (tt1.v[inxt2[ia1]] == nq && tt1.v[iprv2[ia1]] == np) );
  na1 = tt1.v[ia1];
  /*if ( !((tt1.v[inxt2[ia1]] == np && tt1.v[iprv2[ia1]] == nq)
         || (tt1.v[inxt2[ia1]] == nq && tt1.v[iprv2[ia1]] == np))) {
    return(0);
  } remplace par l'assert: impossible non?*/

  for (ia2=0; ia2<3; ia2++) {
    if ( (tt2.v[ia2] != np) && (tt2.v[ia2] != nq) )  break;
  }

  assert ( ia2 < 3 );
  assert ( (tt2.v[inxt2[ia2]] == np && tt2.v[iprv2[ia2]] == nq) ||
           (tt2.v[inxt2[ia2]] == nq && tt2.v[iprv2[ia2]] == np) );
  na2 = tt2.v[ia2];
  /*if ( !((tt2.v[inxt2[ia2]] == np && tt2.v[iprv2[ia2]] == nq)
         || (tt2.v[inxt2[ia2]] == nq && tt2.v[iprv2[ia2]] == np))) {
    return(0);
  } remplace par l'assert: impossible non?*/

  /* Check non convexity (temporarily use b0,b1)*/
  norpts(mesh,tt1.v[ia1],tt1.v[inxt2[ia1]],tt2.v[ia2],b0);
  norpts(mesh,tt2.v[ia2],tt2.v[inxt2[ia2]],tt1.v[ia1],b1);
  ps = b0[0]*b1[0] + b0[1]*b1[1] + b0[2]*b1[2];
  if ( ps < ANGEDG ) return(0);

  /* Compare contributions to Hausdorff distance in both configurations */
  norface(mesh,iel1,ifa1,v);

  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  BezierEdge(mesh,np,nq,b0,b1,0,v);
  c[0] = b0[0] - (p0->c[0] + ATHIRD*ux);
  c[1] = b0[1] - (p0->c[1] + ATHIRD*uy);
  c[2] = b0[2] - (p0->c[2] + ATHIRD*uz);

  disnat = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];

  c[0] = b1[0] - (p1->c[0] - ATHIRD*ux);
  c[1] = b1[1] - (p1->c[1] - ATHIRD*uy);
  c[2] = b1[2] - (p1->c[2] - ATHIRD*uz);

  disnat = MG_MAX(disnat, c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  disnat = MG_MAX(disnat,info.hausd * info.hausd);

  p0 = &mesh->point[na1];
  p1 = &mesh->point[na2];
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  BezierEdge(mesh,na1,na2,b0,b1,0,v);
  c[0] = b0[0] - (p0->c[0] + ATHIRD*ux);
  c[1] = b0[1] - (p0->c[1] + ATHIRD*uy);
  c[2] = b0[2] - (p0->c[2] + ATHIRD*uz);

  dischg = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];

  c[0] = b1[0] - (p1->c[0] - ATHIRD*ux);
  c[1] = b1[1] - (p1->c[1] - ATHIRD*uy);
  c[2] = b1[2] - (p1->c[2] - ATHIRD*uz);

  dischg = MG_MAX(dischg,c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  dischg = MG_MAX(dischg,info.hausd * info.hausd);

  if ( dischg > disnat )   return(0);
  cal1 = caltri(mesh,&tt1);
  cal2 = caltri(mesh,&tt2);
  calnat = MG_MIN(cal1,cal2);
  for (j=0; j<3; j++) {
    if ( tt1.v[j] == nq )  tt1.v[j] = na2;
    if ( tt2.v[j] == np )  tt2.v[j] = na1;
  }
  cal1 = caltri(mesh,&tt1);
  cal2 = caltri(mesh,&tt2);
  calchg = MG_MIN(cal1,cal2);
  if ( calchg < 1.01 * calnat )  return(0);

  /* Check mechanical validity of forthcoming operations */
  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];
  ppt0->c[0] = 0.5*(p0->c[0] + p1->c[0]);
  ppt0->c[1] = 0.5*(p0->c[1] + p1->c[1]);
  ppt0->c[2] = 0.5*(p0->c[2] + p1->c[2]);

  /* Check validity of insertion of midpoint on edge (pq), then collapse of m on a1 */
  calold = calnew = DBL_MAX;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 6;
    pt  = &mesh->tetra[iel];
    memcpy(pt0,pt,sizeof(Tetra));
    calold = MG_MIN(calold, pt->qual);

    ia1 = ia2 = ip = iq = -1;
    for (j=0; j< 4; j++) {
      if (pt->v[j] == np)  ip = j;
      else if (pt->v[j] == nq) iq = j;
      else if ( ia1 < 0 ) ia1 = j;
      else ia2 = j;
    }
    assert((ip >= 0) && (iq >= 0) && (ia1 >= 0) && (ia2 >= 0));
    isshell = (pt->v[ia1] == na1 || pt->v[ia2] == na1);

    /* 2 elts resulting from split and collapse */
    pt0->v[ip] = 0;
    if ( orcal(mesh,0) < NULKAL )  return(0);
    if ( !isshell ) {
      pt0->v[ip] = na1;
      caltmp = orcal(mesh,0);
      calnew = MG_MIN(calnew,caltmp);
    }
    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[iq] = 0;
    if ( orcal(mesh,0) < NULKAL )  return(0);

    if ( !isshell ) {
      pt0->v[iq] = na1;
      caltmp = orcal(mesh,0);
      calnew = MG_MIN(calnew,caltmp);
    }
  }
  if ( calold < NULKAL && calnew <= calold )  return(0);
  else if ( calnew < 0.3 * calold )  return(0);

  return(1);
}

/** Swap boundary edge whose shell is provided ; it1 = boundary face
    carrying the beforehand tested terminal point for collapse */
int swpbdy(pMesh mesh,pSol met,int *list,int ret,int it1) {
  pTetra   pt,pt1;
  pPoint   p0,p1;
  int      iel,iel1,ilist,np,nq,na,nm;
  double   c[3];
  char     ia,iface1,j,ipa,im;
	int      ier;

  iel = list[0] / 6;
  ia  = list[0] % 6;
  pt  = &mesh->tetra[iel];

  np = pt->v[iare[ia][0]];
  nq = pt->v[iare[ia][1]];
  na = 0;

  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];

  /* search for na = the point on quadrangle surfacic configuration on which collapse
     validity has been checked in chkswpbdy */
  iel1 = it1 / 4;
  iface1 = it1 % 4;
  pt1 = &mesh->tetra[iel1];

  for (j=0; j<3;j++) {
    ipa = idir[iface1][j];
    if ( (pt1->v[ipa] != np)&&(pt1->v[ipa] != nq) ) {
      na = pt1->v[ipa];
      break;
    }
  }
  assert(na);

  /* Create midpoint m on edge (pq), then split edge */
  c[0] = 0.5*( p0->c[0] + p1->c[0]);
  c[1] = 0.5*( p0->c[1] + p1->c[1]);
  c[2] = 0.5*( p0->c[2] + p1->c[2]);
  nm = newPt(mesh,c,MG_BDY);
  if ( !nm ) {
    fprintf(stdout,"  ## Warning: unable to allocate a new point.\n");
    fprintf(stdout,"  ## Check the mesh size or ");
    fprintf(stdout,"increase the allocated memory with the -m option.\n");
    return(0);
  }
  if ( met->m )  met->m[nm] = 0.5 *(met->m[np]+met->m[nq]);
  split1b(mesh,met,list,ret,nm,0);

  /* Collapse m on na after taking (new) ball of m */
  memset(list,0,(LMAX+2)*sizeof(int));
  for (j=0; j<3; j++) {
    im = idir[iface1][j];
    if ( pt1->v[im] == nm )  break;
  }
  if ( pt1->v[im] != nm ){
    delPt(mesh,nm);
    printf("%s:%d: Warning pt1->v[im] != nm\n",__FILE__,__LINE__);
    return(0);
  }
  ilist = boulevolp(mesh,iel1,im,list);

  assert(list[0]/4 == iel1);
  assert(pt1->v[ipa] == na);

  ier = colver(mesh,list,ilist,ipa);
	if ( ier ) {
		delPt(mesh,ier);
		ier = 1;
	}

  return(ier);
}
