#include "mmg3d.h"

extern Info info;

/** Check whether swap of edge ia in start should be performed, and return 4*k+i =
    index of point corresponding to the swapped configuration ; shell of edge is
    built during the process */
int chkswpgen(pMesh mesh,int start,int ia,int *ilist,int *list) {
  pTetra    pt,pt0;
  pPoint    p0;
  double    calold,calnew,caltmp;
  int       na,nb,np,adj,piv,npol,ref,refdom,k,l,iel;
  int       *adja,pol[LMAX+2];
  char      i,ipa,ipb,tag,ip,ier;

  pt  = &mesh->tetra[start];
  if ( !MG_EOK(pt) ) return(0);
  refdom = pt->ref;

  pt0 = &mesh->tetra[0];
  na  = pt->v[iare[ia][0]];
  nb  = pt->v[iare[ia][1]];
  if((&mesh->tetra[start])->qual==-10)
    (&mesh->tetra[start])->qual=orcal(mesh,start);
  calold = (&mesh->tetra[start])->qual;
  if ( ALPHAD*calold > 0.6 )  return(0);

  /* Prevent swap of a ref or tagged edge */
  hGet(&mesh->htab,na,nb,&ref,&tag);
  if ( ref || tag )  return(0);

  /* Store shell of ia in list, and associated pseudo polygon in pol */
  (*ilist) = 0;
  npol = 0;
  list[(*ilist)] = 6*start+ia;
  (*ilist)++;
  adja = &mesh->adja[4*(start-1)+1];
  adj  = adja[ifar[ia][0]] / 4;      // start travelling by face (ia,0)
  piv  = pt->v[ifar[ia][1]];
  pol[npol] = 4*start + ifar[ia][1];
  npol++;

  while ( adj && adj != start ) {
    pt = &mesh->tetra[adj];

    /* Edge is on a boundary between two different domains */
    if ( pt->ref != refdom )  return(0);
    if( (&mesh->tetra[adj])->qual ==-10)
      (&mesh->tetra[adj])->qual=orcal(mesh,adj);
    calold = MG_MIN(calold, (&mesh->tetra[adj])->qual);
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = iare[i][0];
      ipb = iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    list[(*ilist)] = 6*adj +i;
    (*ilist)++;
    /* overflow */
    if ( (*ilist) > LMAX-3 )  return(0);

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ ifar[i][0] ] == piv ) {
      pol[npol] = 4*adj + ifar[i][1];
      npol++;
      adj = adja[ ifar[i][0] ] / 4;
      piv = pt->v[ ifar[i][1] ];
    }
    else {
      assert(pt->v[ ifar[i][1] ] == piv);
      pol[npol] = 4*adj + ifar[i][0];
      npol++;
      adj = adja[ ifar[i][1] ] /4;
      piv = pt->v[ ifar[i][0] ];
    }
  }
  if ( calold*ALPHAD > 0.5 )  return(0);

  /* Prevent swap of an external boundary edge */
  if ( !adj )  return(0);
  assert(npol == (*ilist)); // du coup, apres on pourra virer npol

  /* Find a configuration that enhances the worst quality within the shell */
  for (k=0; k<npol; k++) {
    iel = pol[k] / 4;
    ip  = pol[k] % 4;
    np  = mesh->tetra[iel].v[ip];
    calnew = 1.0;
    ier = 1;

    if ( info.fem ) {
      p0 = &mesh->point[np];
      if ( p0->tag & MG_BDY ) {
        for (l=0; l<npol;l++) {
          if ( k < npol-1 ) {
            if ( l == k || l == k+1 )  continue;
          }
          else {
            if ( l == npol-1 || l == 0 )  continue;
          }
          iel = pol[l] / 4;
          ip  = pol[l] % 4;
          pt = &mesh->tetra[iel];
          p0 = &mesh->point[pt->v[ip]];
          if ( p0->tag & MG_BDY ) {
            ier = 0;
            break;
          }
        }
      }
      if ( !ier )  continue;
      ier = 1;
    }

    for (l=0; l<(*ilist); l++) {
      /* Do not consider tets of the shell of collapsed edge */
      if ( k < npol-1 ) {
        if ( l == k || l == k+1 )  continue;
      }
      else {
        if ( l == npol-1 || l == 0 )  continue;
      }
      iel = list[l] / 6;
      i   = list[l] % 6;
      pt  = &mesh->tetra[iel];

      /* First tetra obtained from iel */
      memcpy(pt0,pt,sizeof(Tetra));
      pt0->v[iare[i][0]] = np;
      caltmp = orcal(mesh,0);
      calnew = MG_MIN(calnew,caltmp);
      /* Second tetra obtained from iel */
      memcpy(pt0,pt,sizeof(Tetra));
      pt0->v[iare[i][1]] = np;
      caltmp = orcal(mesh,0);
      calnew = MG_MIN(calnew,caltmp);
      ier = (calnew > 1.1*calold);
      if ( !ier )  break;
    }
    if ( ier )  return(pol[k]);
  }
  return(0);
}

/** Perform swap of edge whose shell is passed according to configuration nconf */
void swpgen(pMesh mesh,pSol met,int nconf,int ilist,int *list) {
  pTetra    pt;
  pPoint    p0,p1;
  int       iel,na,nb,np,nball,ret,start;
  double    m[3];
  char      ia,ip,iq;

  iel = list[0] / 6;
  ia  = list[0] % 6;

  pt = &mesh->tetra[iel];
  na = pt->v[iare[ia][0]];
  nb = pt->v[iare[ia][1]];
  p0 = &mesh->point[na];
  p1 = &mesh->point[nb];

  /* Temporarily create midpoint at swapped edge */
  m[0] = 0.5*(p0->c[0] + p1->c[0]);
  m[1] = 0.5*(p0->c[1] + p1->c[1]);
  m[2] = 0.5*(p0->c[2] + p1->c[2]);

  np  = newPt(mesh,m,0);
  assert(np);
  if ( met->m )  met->m[np] = 0.5*(met->m[na]+met->m[nb]);

  /** First step : split of edge (na,nb) */
  ret = 2*ilist + 0;
  split1b(mesh,met,list,ret,np,0);

  /** Second step : collapse of np towards enhancing configuration */
  start = nconf / 4;
  iq = nconf % 4;

  pt = &mesh->tetra[start];
  for (ip=0; ip<4; ip++) {
    if ( pt->v[ip] == np )  break;
  }
  assert(ip<4);

  memset(list,0,(LMAX+2)*sizeof(int));
  nball = boulevolp(mesh,start,ip,list);

  colver(mesh,list,nball,iq);
}
