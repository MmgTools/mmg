#include "mmg3d.h"

extern Info  info;
extern char ddb;

inline double lenedg_ani(pMesh mesh,int ip1,int ip2) {
  return(0.0);
}

/* Compute length of edge [ip1,ip2] according to the size prescription */
inline double lenedg_iso(pMesh mesh,int ip1,int ip2) {
  pPoint   p1,p2;
  double   h1,h2,l,r,len;

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  h1 = p1->h;
  h2 = p2->h;
  l = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]) \
    + (p2->c[2]-p1->c[2])*(p2->c[2]-p1->c[2]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;
  len = fabs(r) < EPS ? l / h1 : l / (h2-h1) * log(r+1.0);

  return(len);
}

/* Return quality of surface triangle */
inline double caltri(pMesh mesh,pTria ptt) {
  double   *a,*b,*c,cal,abx,aby,abz,acx,acy,acz,bcx,bcy,bcz,rap;

  a = &mesh->point[ptt->v[0]].c[0];
  b = &mesh->point[ptt->v[1]].c[0];
  c = &mesh->point[ptt->v[2]].c[0];

  /* area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  cal  = (aby*acz - abz*acy) * (aby*acz - abz*acy);
  cal += (abz*acx - abx*acz) * (abz*acx - abx*acz);
  cal += (abx*acy - aby*acx) * (abx*acy - aby*acx);
  if ( cal < EPSD2 )  return(0.0);

  /* qual = 2.*surf / length */
  rap  = abx*abx + aby*aby + abz*abz;
  rap += acx*acx + acy*acy + acz*acz;
  rap += bcx*bcx + bcy*bcy + bcz*bcz;
  if ( rap < EPSD2 )  return(0.0);

  return(sqrt(cal) / rap);
}

/* compute tetra oriented quality of iel (return 0.0 when element is inverted) */
inline double orcal(pMesh mesh,int iel) {
  pTetra     pt;
  double     abx,aby,abz,acx,acy,acz,adx,ady,adz,bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     vol,v1,v2,v3,rap;
  double     *a,*b,*c,*d;

  pt = &mesh->tetra[iel];
  a = mesh->point[pt->v[0]].c;
  b = mesh->point[pt->v[1]].c;
  c = mesh->point[pt->v[2]].c;
  d = mesh->point[pt->v[3]].c;

  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  rap = abx*abx + aby*aby + abz*abz;

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  rap += acx*acx + acy*acy + acz*acz;

  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];
  rap += adx*adx + ady*ady + adz*adz;

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;

  if ( vol < EPSD2 )  return(0.0);

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];
  rap += bcx*bcx + bcy*bcy + bcz*bcz;

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];
  rap += bdx*bdx + bdy*bdy + bdz*bdz;

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];
  rap += cdx*cdx + cdy*cdy + cdz*cdz;
  if ( rap < EPSD2 )  return(0.0);

  /* quality = vol / len^3/2 */
  rap = rap * sqrt(rap);

  return(vol / rap);
}


/* compute tetra quality iso */
inline double caltet_iso(pMesh mesh,pSol met,int ia,int ib,int ic,int id) {
  double     abx,aby,abz,acx,acy,acz,adx,ady,adz,bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     vol,v1,v2,v3,rap;
  double    *a,*b,*c,*d;

  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;
  d = mesh->point[id].c;

  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  rap = abx*abx + aby*aby + abz*abz;

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  rap += acx*acx + acy*acy + acz*acz;

  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];
  rap += adx*adx + ady*ady + adz*adz;

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;
  if ( vol < EPSD2 )  return(0.0);

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];
  rap += bcx*bcx + bcy*bcy + bcz*bcz;

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];
  rap += bdx*bdx + bdy*bdy + bdz*bdz;

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];
  rap += cdx*cdx + cdy*cdy + cdz*cdz;
  if ( rap < EPSD2 )  return(0.0);

  /* quality = vol / len^3/2 */
  rap = rap * sqrt(rap);
  return(vol / rap);
}


inline double caltet_ani(pMesh mesh,pSol met,int ia,int ib,int ic,int id) {
  return(0.0);
}


/* compute face normal */
inline int nortri(pMesh mesh,pTria pt,double *n) {
  double   *a,*b,*c,dd,abx,aby,abz,acx,acy,acz,det;

  a = mesh->point[pt->v[0]].c;
  b = mesh->point[pt->v[1]].c;
  c = mesh->point[pt->v[2]].c;

  /* area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];

  n[0] = aby*acz - abz*acy;
  n[1] = abz*acx - abx*acz;
  n[2] = abx*acy - aby*acx;
  det  = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( det < EPSD2 )  return(0);

  dd = 1.0 / sqrt(det);
  n[0] *= dd;
  n[1] *= dd;
  n[2] *= dd;
  return(1);
}


int badelt(pMesh mesh,pSol met) {
  pTetra   pt;
  double   kal;
  int      k,ier,it,maxit,nd;

  it = 0;
  maxit = 1;
  do {
    nd = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      kal = ALPHAD * orcal(mesh,k);
      if ( kal > BADKAL )  continue;

      nd++;
    }
  }
  while ( ++it < maxit && nd > 0 );
  printf("nd = %d\n",nd);
  return(nd);
}

/* print mesh quality histo */
int outqua(pMesh mesh,pSol met) {
  pTetra    pt;
  double   rap,rapmin,rapmax,rapavg,med;
  int      i,k,iel,ok,ir,imax,nex,his[5];

  rapmin  = 2.0;
  rapmax  = 0.0;
  rapavg  = med = 0.0;
  iel     = 0;

  for (k=0; k<5; k++)  his[k] = 0;

  nex = ok = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) ) {
      nex++;
      continue;
    }
    ok++;
    rap = ALPHAD * caltet(mesh,met,pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
    if ( rap < rapmin ) {
      rapmin = rap;
      iel    = ok;
    }
    if ( rap > 0.5 )  med++;
    if ( rap < BADKAL )  info.badkal = 1;
    rapavg += rap;
    rapmax  = MG_MAX(rapmax,rap);
    ir = MG_MIN(4,(int)(5.0*rap));
    his[ir] += 1;
  }

  fprintf(stdout,"\n  -- MESH QUALITY   %d\n",mesh->ne - nex);
  fprintf(stdout,"     BEST   %8.6f  AVRG.   %8.6f  WRST.   %8.6f (%d)\n",
	  rapmax,rapavg / (mesh->ne-nex),rapmin,iel);
  if ( abs(info.imprim) < 5 )  return;

  /* print histo */
  fprintf(stdout,"     HISTOGRAMM:  %6.2f %% > 0.5\n",100.0*(med/(float)(mesh->ne-nex)));
  imax = MG_MIN(4,(int)(5.*rapmax));
  for (i=imax; i>=(int)(5*rapmin); i--) {
    fprintf(stdout,"     %5.1f < Q < %5.1f   %7d   %6.2f %%\n",
	    i/5.,i/5.+0.2,his[i],100.*(his[i]/(float)(mesh->ne-nex)));
  }
}
