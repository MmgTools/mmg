#include "mmg3d.h"

extern char ddb;

inline double lenedg_ani(pMesh mesh,pSol met,int ip1,int ip2) {
  return(0.0);
}

/** Compute length of edge [ip1,ip2] according to the size prescription */
inline double lenedg_iso(pMesh mesh,pSol met,int ip1,int ip2) {
  pPoint   p1,p2;
  double   h1,h2,l,r,len;

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  h1 = met->m[ip1];
  h2 = met->m[ip2];
  l = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]) \
    + (p2->c[2]-p1->c[2])*(p2->c[2]-p1->c[2]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;
  len = fabs(r) < EPS ? l / h1 : l / (h2-h1) * log(r+1.0);

  return(len);
}


/** Return quality of surface triangle */
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

/** compute tetra oriented quality of iel (return 0.0 when element is inverted) */
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


/** compute tetra quality iso */
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


/** compute face normal */
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

/* identify type of element :
  ityp= 0: 4 faces bonnes          (elt ok)
        1: 4 faces bonnes, vol nul (sliver) ou "quasi sliver" ie 4 faces ok
        2: 4 faces ok, vol nul+sommet proche face   (chapeau)
        3: 3 faces bonnes, 1 obtuse    (aileron)
        4: 2 faces bonnes, 2 faces aigu => 1 petite arete
        5: 1 face bonne, 3 petites aretes
        6: 2 faces grandes aretes, 2 faces petites iaretes
        7: 4 faces grandes aretes
        8: 2 faces obtus, 1 faces aigu et une face OK
   item: bad entity
*/


/* nb face obtuse :    nb faces aigu :
ityp :  0: 0		0
        1: 0		0
        2: 0		0
        3: 1		0
        4: 0		2
        5: 0		3
        6: 2		2
        7: 0		4
*/
/* nb gde arete :    nb petite arete :
ityp :  0: 0		0
        1: 0		0
        2: 0		0
        3: 1		0
        4: 0		1
        5: 0		3
        6: 1		1
        7: 0		2
*/
#define EPSVOL 0.001
#define RAPMAX    0.4//0.3//0.25
unsigned char inxt[7]    = { 1,2,0,1,2,0,1 };

int typelt(pMesh mesh,int iel,int *item) {
  pTetra    pt;
  pPoint    pa,pb,pc,pd;
  double    abx,aby,abz,acx,acy,acz,adx,ady,adz,v1,v2,v3,vol;
  double    bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz,h[6],volchk,ssmall;
  double    s[4],dd,rapmin,rapmax,surmin,surmax;
  int       i,k,ia,ib,ic,id,isur,isurmax,isurmin,iarmax,iarmin;
  int       nobtus,naigu,aigu;
  short     i0,i1,i2;
  double lmoy;

  pt = &mesh->tetra[iel];
  if ( !pt->v[0] )  return(-1);

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];
  pa = &mesh->point[ia];
  pb = &mesh->point[ib];
  pc = &mesh->point[ic];
  pd = &mesh->point[id];

  /* volume */
  abx = pb->c[0] - pa->c[0];
  aby = pb->c[1] - pa->c[1];
  abz = pb->c[2] - pa->c[2];

  acx = pc->c[0] - pa->c[0];
  acy = pc->c[1] - pa->c[1];
  acz = pc->c[2] - pa->c[2];

  adx = pd->c[0] - pa->c[0];
  ady = pd->c[1] - pa->c[1];
  adz = pd->c[2] - pa->c[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;

  /* max edge */
  h[0] = abx*abx + aby*aby + abz*abz;
  h[1] = acx*acx + acy*acy + acz*acz;
  h[2] = adx*adx + ady*ady + adz*adz;

  bcx = pc->c[0] - pb->c[0];
  bcy = pc->c[1] - pb->c[1];
  bcz = pc->c[2] - pb->c[2];

  bdx = pd->c[0] - pb->c[0];
  bdy = pd->c[1] - pb->c[1];
  bdz = pd->c[2] - pb->c[2];

  cdx = pd->c[0] - pc->c[0];
  cdy = pd->c[1] - pc->c[1];
  cdz = pd->c[2] - pc->c[2];

  h[3] = bcx*bcx + bcy*bcy + bcz*bcz;
  h[4] = bdx*bdx + bdy*bdy + bdz*bdz;
  h[5] = cdx*cdx + cdy*cdy + cdz*cdz;

  /* face areas */
  dd = cdy*bdz - cdz*bdy;
  s[0] = dd * dd;
  dd = cdz*bdx - cdx*bdz;
  s[0] = s[0] + dd * dd;
  dd = cdx*bdy - cdy*bdx;
  s[0] = s[0] + dd * dd;
  s[0] = sqrt(s[0]);

  s[1] = sqrt(v1*v1 + v2*v2 + v3*v3);

  dd = bdy*adz - bdz*ady;
  s[2] = dd * dd;
  dd = bdz*adx - bdx*adz;
  s[2] = s[2] + dd * dd;
  dd = bdx*ady - bdy*adx;
  s[2] = s[2] + dd * dd;
  s[2] = sqrt(s[2]);

  dd = aby*acz - abz*acy;
  s[3] = dd * dd;
  dd = abz*acx - abx*acz;
  s[3] = s[3] + dd * dd;
  dd = abx*acy - aby*acx;
  s[3] = s[3] + dd * dd;
  s[3] = sqrt(s[3]);

  /* classification */
  rapmin = h[0];
  rapmax = h[0];
  iarmin = 0;
  iarmax = 0;
  for (i=1; i<6; i++) {
    if ( h[i] < rapmin ) {
      rapmin = h[i];
      iarmin = i;
    }
    else if ( h[i] > rapmax ) {
      rapmax = h[i];
      iarmax = i;
    }
  }
  rapmin = sqrt(rapmin);
  rapmax = sqrt(rapmax);
  volchk = EPSVOL * rapmin*rapmin*rapmin;
  //printf("$$$$$$$$$$$$$$$$$$$$$ iel %d %e %e\n",iel,volchk,vol);
  /* small volume: types 1,2,3,4 */
  if ( vol < volchk ) {
    puts("volume nul : type 1,2,3,4");

    ssmall = 0.4 * (s[0]+s[1]+s[2]+s[3]);
    isur   = 0;
    for (i=0; i<4; i++)
      isur += s[i] > ssmall;

    /* types 2,3 */
    item[0] = iarmax;
    item[1] = isar[iarmax][0];
    if ( isur == 1 ) {
      surmin   = s[0];
      isurmin = 0;
      surmax   = s[0];
      isurmax = 0;
      for (i=1; i<4; i++) {
        if ( s[i] < surmin ) {
          surmin  = s[i];
    isurmin = i;
  }
        else if ( s[i] > surmax ) {
    surmax  = s[i];
    isurmax = i;
  }
      }
      dd = surmin / surmax;
      if ( dd < RAPMAX ) {
        item[1] = isar[iarmax][0];
        return(3);
      }
      else {
        item[0] = isurmax;
  item[1] = isurmin;
        return(2);
      }
    }

    /* types 1 */
    isur = 0;
    if ( s[0]+s[1] > ssmall )  isur = 1;
    if ( s[0]+s[2] > ssmall )  isur++;
    if ( s[0]+s[3] > ssmall )  isur++;

    if ( isur > 2 ) {
      dd = rapmin / rapmax;
      item[0] = iarmin;
      item[1] = idir[iarmin][0];
      if ( dd < 0.01 )  return(4);
      if ( s[0]+s[1] > ssmall ) {
        item[0] = 0;
        return(1);
      }
      if ( s[0]+s[2] > ssmall ) {
        item[0] = 1;
        return(1);
      }
      if ( s[0]+s[3] > ssmall ) {
        item[0] = 2;
        return(1);
      }
    }

//puts("default");
    item[0] = 0;
    return(1);
  }/*end chkvol*/

 dd = rapmin / rapmax;
 // printf("dd %e %e %e\n",dd,RAPMAX,0.7*RAPMAX);
  /* types 3,6,7 */
 if ( dd < RAPMAX ) { /*ie une arete 3 fois plus gde qu'une autre*/
   lmoy = 0;
   for (i=0; i<6; i++)  {h[i] = sqrt(h[i]); lmoy+=h[i];}
   lmoy *= 1./6.;

    nobtus = 0;
    naigu = 0;
    for (k=0; k<4; k++) {
      aigu = 0;
      for(i=0 ; i<3 ; i++) {
  i0 = idir[k][i];
        i1 = idir[k][inxt[i]];
        i2 = idir[k][inxt[i+1]];
  if((h[i0]>h[i1] && h[i0]>h[i2]) && (h[i0]>2.5*h[i1] || h[i0]>2.5*h[i2])) {//obtu ? > /*130*/ 150
    //be carefull isocele triangle --> opposite edge only
    if(!( fabs(1-h[i0]/h[i1]) < 0.1 || fabs(1-h[i0]/h[i2])< 0.1 ) ) {
      if(h[i0] > 0.9659258/*0.9063078*/*(h[i1]+h[i2])) break;
    }
  }
    if((h[i0]<h[i1] && h[i0]<h[i2]) && (2.5*h[i0]<h[i1] || 2.5*h[i0]<h[i2])) {//aigu ? <20
      if(h[i0] < 0.17*(h[i1]+h[i2])) aigu++;
  }
      }
      if(i<3) nobtus++;
      else if(aigu) naigu++;
    }
    //printf("%d on trouve %d aigu et %d obtu\n",iel,naigu,nobtus);
    switch(nobtus) {
    case(3):case(4):
      return(8);
    case(2):
      if(naigu==2) {
  return(6);
      } else {
  return(8);
      }
    case(1):
      return(3);
    case(0):
      if(naigu==4) {
  return(7);
      } else if(naigu==3) {
  return(5);
      } else if(naigu==2) {
  return(4);
      } else {
  //printf("%d on trouve %d aigu et %d obtu\n",iel,naigu,nobtus);
  return(9); //toutes les faces sont ok
      }
    }
 }


  item[0] = 0;
  // printf("edge %e %e -- vol %e\n",rapmin,rapmax,vol);
  //puts("default");
  return(9);
}

int badelt(pMesh mesh,pSol met) {
  pTetra   pt;
  double   kal;
  int      k,it,maxit,nd/*,item[2],typ*/;
  /*int      ntyp[10];*/
  int      list[LMAX+2],i,ilist,nconf,ns;

  it = 0;
  maxit = 3;
  do {
    nd = 0;
    ns = 0;
    //for (k=1; k<10; k++) ntyp[k]=0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      kal = /*ALPHAD **/ pt->qual;
      if ( kal > 0.0096225 /*BADKAL/ALPHAD*/ )  continue;
      //typ =  typelt(mesh,k,item);
      //ntyp[typ]++;
      nd++;
      /*treat bad elt*/
      /*1) try to swp one edge*/
      for(i=0 ; i<6 ; i++) {
        nconf = chkswpgen(mesh,k,i,&ilist,list,1.01);
        if ( nconf ) {
          ns++;
#ifdef PATTERN
          if(!swpgen(mesh,met,nconf,ilist,list)) return(-1);
#else
          if(!swpgen(mesh,met,nconf,ilist,list,NULL)) return(-1);
#endif
          break;
        }
      }
    }
    /*printf("on trouve %d bad elt\n",nd);
    for (k=0; k<=9; k++)
      if ( ntyp[k] )
        printf("  optim [%d]      = %5d  %6.2f %%\n",k,ntyp[k],100.0*ntyp[k]/nd);
    */if ( ns > 0 )
  fprintf(stdout,"     %8d edge swapped\n",ns);
   }
  while ( ++it < maxit && nd > 0 );
  return(nd);
}

/** Compute sizes of edges of the mesh, and displays histo */
int prilen(pMesh mesh, pSol met) {
  pTetra          pt;
  Hash            hash;
  double          len,avlen,dned,lmin,lmax;
  int             k,np,nq,amin,bmin,amax,bmax,ned,hl[9];
  char            ia,i0,i1,ier,i;
  static double   bd[9]= {0.0, 0.3, 0.6, 0.7071, 0.9, 1.3, 1.4142, 2.0, 5.0};
  //{0.0, 0.2, 0.5, 0.7071, 0.9, 1.111, 1.4142, 2.0, 5.0};

  memset(hl,0,9*sizeof(int));
  ned = 0;
  avlen = 0.0;
  lmax = 0.0;
  lmin = 1.e30;
  amin = amax = bmin = bmax = 0;

  /* Hash all edges in the mesh */
  if ( !hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return(0);

  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = iare[ia][0];
      i1 = iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      if(!hashEdge(mesh,&hash,np,nq,0)){
        fprintf(stdout,"%s:%d: Error: function hashEdge return 0\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
      }
    }
  }

  /* Pop edges from hash table, and analyze their length */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = iare[ia][0];
      i1 = iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      /* Remove edge from hash ; ier = 1 if edge has been found */
      ier = hashPop(&hash,np,nq);
      if( ier ) {
        ned ++;
        len = lenedg(mesh,met,np,nq);
        avlen += len;

        if( len < lmin ) {
          lmin = len;
          amin = np;
          bmin = nq;
        }

        if ( len > lmax ) {
          lmax = len;
          amax = np;
          bmax = nq;
        }

        /* Locate size of edge among given table */
        for(i=0; i<8; i++) {
          if ( bd[i] <= len && len < bd[i+1] ) {
            hl[i]++;
            break;
          }
        }
        if( i == 8 ) hl[8]++;
      }
    }
  }

  /* Display histogram */
  dned = (double)ned;
  avlen = avlen / dned;

  fprintf(stdout,"\n  -- RESULTING EDGE LENGTHS  %d\n",ned);
  fprintf(stdout,"     AVERAGE LENGTH         %12.4f\n",avlen);
  fprintf(stdout,"     SMALLEST EDGE LENGTH   %12.4f   %6d %6d\n",
          lmin,amin,bmin);
  fprintf(stdout,"     LARGEST  EDGE LENGTH   %12.4f   %6d %6d \n",
          lmax,amax,bmax);

  if ( hl[3]+hl[4]+hl[5] )
    fprintf(stdout,"   %6.2f < L <%5.2f  %8d   %5.2f %%  \n",
            bd[3],bd[6],hl[3]+hl[4]+hl[5],100.*(hl[3]+hl[4]+hl[5])/(double)ned);
  if ( hl[2]+hl[3]+hl[4] )
    fprintf(stdout,"   %6.2f < L <%5.2f  %8d   %5.2f %%  \n",
            bd[2],bd[5],hl[2]+hl[3]+hl[4],100.*(hl[2]+hl[3]+hl[4])/(double)ned);


  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"\n     HISTOGRAMM:\n");
    if ( hl[0] )
      fprintf(stdout,"     0.00 < L < 0.30  %8d   %5.2f %%  \n",
              hl[0],100.*(hl[0]/(float)ned));
    if ( lmax > 0.2 ) {
      for (k=2; k<9; k++) {
        if ( hl[k-1] > 0 )
          fprintf(stdout,"   %6.2f < L <%5.2f  %8d   %5.2f %%  \n",
                  bd[k-1],bd[k],hl[k-1],100.*(hl[k-1]/(float)ned));
      }
      if ( hl[8] )
        fprintf(stdout,"     5.   < L         %8d   %5.2f %%  \n",
                hl[8],100.*(hl[8]/(float)ned));
    }
  }

  DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
  return(1);
}

/** print mesh quality histo */
void outqua(pMesh mesh,pSol met) {
  pTetra    pt;
  double   rap,rapmin,rapmax,rapavg,med,good;
  int      i,k,iel,ok,ir,imax,nex,his[5];

  rapmin  = 2.0;
  rapmax  = 0.0;
  rapavg  = med = good = 0.0;
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
    if ( rap > 0.12 ) good++;
    if ( rap < BADKAL )  mesh->info.badkal = 1;
    rapavg += rap;
    rapmax  = MG_MAX(rapmax,rap);
    ir = MG_MIN(4,(int)(5.0*rap));
    his[ir] += 1;
  }

#ifndef DEBUG
  fprintf(stdout,"\n  -- MESH QUALITY   %d\n",mesh->ne - nex);
  fprintf(stdout,"     BEST   %8.6f  AVRG.   %8.6f  WRST.   %8.6f (%d)\n",
          rapmax,rapavg / (mesh->ne-nex),rapmin,iel);
#else
  fprintf(stdout,"     BEST   %e  AVRG.   %e  WRST.   %e (%d)\n => %d %d %d %d\n",
          rapmax,rapavg / (mesh->ne-nex),rapmin,iel,
          indPt(mesh,mesh->tetra[iel].v[0]),indPt(mesh,mesh->tetra[iel].v[1]),
          indPt(mesh,mesh->tetra[iel].v[2]),indPt(mesh,mesh->tetra[iel].v[3]));
#endif
  if ( abs(mesh->info.imprim) < 4 ){
    if (rapmin == 0){
      fprintf(stdout,"  ## WARNING: TOO BAD QUALITY FOR THE WORST ELEMENT\n");
      saveMesh(mesh);
      exit(EXIT_FAILURE);
    }
    return;
  }

  /* print histo */
  fprintf(stdout,"     HISTOGRAMM:");
  fprintf(stdout,"  %6.2f %% > 0.12\n",100.0*(good/(float)(mesh->ne-nex)));
  if ( abs(mesh->info.imprim) > 4 ) {
    fprintf(stdout,"                  %6.2f %% >  0.5\n",100.0*( med/(float)(mesh->ne-nex)));
    imax = MG_MIN(4,(int)(5.*rapmax));
    for (i=imax; i>=(int)(5*rapmin); i--) {
      fprintf(stdout,"     %5.1f < Q < %5.1f   %7d   %6.2f %%\n",
              i/5.,i/5.+0.2,his[i],100.*(his[i]/(float)(mesh->ne-nex)));
    }
  }
  if (rapmin == 0){
    fprintf(stdout,"  ## WARNING: TOO BAD QUALITY FOR THE WORST ELEMENT\n");
    saveMesh(mesh);
    exit(EXIT_FAILURE);
  }
}

/*approximation of the final number of vertex*/
int countelt(pMesh mesh,pSol sol, double *weightelt, int *npcible) {
  pTetra pt;
  double *ca,*cb,*ma,*mb,len;
  int    k,ia,ipa,ipb,iad,lon,l,nptot,iadr;
  int    *pdel,lenint,loc,nedel,longen;
  int    npbdry,isbdry;
  double   dned,dnface,dnint,dnins;
  double   dnpdel,dnadd,leninv,dnaddloc,dnpdelloc;
  int   list[LMAX];

  pdel = (int*) calloc(mesh->np,sizeof(int));
  nptot = mesh->np;
  npbdry = 0;

  //substraction of the half of the number of bdry vertex to avoid the surestimation due of the interface
  for (k=1; k<=mesh->np; k++) {
    if(mesh->point[k].tag & MG_BDY) npbdry++;
  }
  nptot -= 0.5*npbdry;

  dnadd = dnpdel = 0;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] )  continue;

    if(weightelt)
      weightelt[k] = 0;
    nedel = 0;

    for (ia=0; ia<6; ia++) {
      //lon = MMG_coquil(mesh,k,ia,&list);
      longen = coquil(mesh,k,ia,list);
      lon = longen/2;
      isbdry = longen%2;
      if(!lon) continue;
      /* if ( isbdry )  { */
      /* 	assert(longen%2); */
      /* 	//printf("coquil %d\n",longen/2); */
      /* 	continue; */
      /* } */
      //assert(!(longen%2));
      for (l=1; l<lon; l++)
	if ( list[l] < 6*k )  break;

      if ( l < lon )  {
	loc = 1;
	//continue;
      } else {
	loc = 0;
      } 

      dnaddloc = 0;
      dnpdelloc = 0;
      
      ipa = iare[ia][0];
      ipb = iare[ia][1];
      /* ca  = &mesh->point[pt->v[ipa]].c[0]; */
      /* cb  = &mesh->point[pt->v[ipb]].c[0]; */

      /* iadr = (pt->v[ipa]-1)*1 + 1; */
      /* ma   = &sol->met[iadr]; */
      /* iadr = (pt->v[ipb]-1)*1 + 1; */
      /* mb   = &sol->met[iadr]; */
      //if(sol->offset==6)
      //	len = MMG_long_ani_init(ca,cb,ma,mb);
      //else
	len = lenedg(mesh,sol,pt->v[ipa],pt->v[ipb]);

      if(len > 3) {
	lenint = ((int) len); 
	if(fabs(lenint -len) > 0.5) lenint++;
	lenint++;
	//nb de point a inserer sur l'arete
	dned = lenint - 2;
	//nb de point a l'interieur de la face si toutes les aretes sont coupees le meme nb de fois
	dnface = (lenint+1)*lenint / 2. - 3 - 3*dned;
	//nb de point a l'interieur du tetra si toutes les aretes sont coupees le meme nb de fois
	dnint = (lenint+2)*(lenint+1)*lenint / 6. - 4 - 4*dnface - 6*dned;
	//nb de point a inserer pour cette arete de ce tetra : on divise par lon
	dnins = dned*(1./lon) + (dnface/3. + dnint/6.);//(dnface/12. + dnint/6.);
	if(!isbdry) {
	  //nb points sur l'arete + lon*(1/3 nb point sur la face + 1/6 nb de point interne)
	  dnaddloc = dned + lon*(dnface/3. + dnint/6.);
	} else {
	  dnaddloc = 0.5*(dned + lon*(dnface/3. + dnint/6.));
	}
	if(!loc) {
	  dnadd += dnaddloc;
	}
      } else if(len > 2.8) {
	if(!isbdry) {
	  dnaddloc = 2.;
	} else {
	  dnaddloc = 1;
	}
	if(!loc){
	  if(!isbdry) {
	    dnadd += 2.;
	  } else {
	    dnadd++;
	  }
	}
	dnins = 2;
      } else if(len > 1.7) {
	if(!loc) {
	  if(!isbdry) dnadd += 1.;
	}
	dnins = 1;
      } else if(len < 0.6) {
	nedel = 1;
	
	leninv = 1./len;
	if(pt->v[ipa]<pt->v[ipb]) {
	  if(!pdel[pt->v[ipa]]) {
	    if(!isbdry) {
	      dnpdelloc = (leninv - 1.)/leninv;
	    } else {
	      dnpdelloc = 0.5*(leninv - 1.)/leninv;
	    }
	    if(!loc) {
	      dnpdel+=dnpdelloc;
	      pdel[pt->v[ipa]]=1;
	    }
	  } else if(!pdel[pt->v[ipb]]) {
	    if(!isbdry) {
	      dnpdelloc = (leninv - 1.)/leninv;
	    } else {
	      dnpdelloc = 0.5*(leninv - 1.)/leninv;
	    }
	    if(!loc) {
	      dnpdel +=dnpdelloc;
	      pdel[pt->v[ipb]]=1;
	    }	  
	  }
	} else {
	  if(!pdel[pt->v[ipb]]) {
	    if(!isbdry) {
	      dnpdelloc = (leninv - 1.)/leninv;
	    } else {
	      dnpdelloc = 0.5*(leninv - 1.)/leninv;
	    }
	    if(!loc) {
	      dnpdel+=dnpdelloc;
	      pdel[pt->v[ipb]]=1;
	    }
	  } else if(!pdel[pt->v[ipa]]) {
	    if(!isbdry) {
	      dnpdelloc = (leninv - 1.)/leninv;
	    } else {
	      dnpdelloc = 0.5*(leninv - 1.)/leninv;
	    }
	    if(!loc) {
	      dnpdel+=dnpdelloc;
	      pdel[pt->v[ipa]]=1;
	    }
	  }
	}
	//ndel++;
      }
      
      //pour cette arete de ce tetra :
      //PHASE 1 = dnaddloc + nedel (on compte un si arete trop petite)
      //PHASE 2 = dnaddloc
      if(weightelt)
	weightelt[k] += 1./lon*(2*dnaddloc);//1./lon*(2*dnaddloc + dnpdelloc);

    }/*for ia*/
    if(weightelt)
	weightelt[k] += nedel;
 
  } /*For k*/


  nptot += (int) dnadd - (int) dnpdel;
  *npcible = nptot;
  fprintf(stdout,"ESTIMATION OF THE FINAL NUMBER OF NODES : %8d  ADD %f  DEL %f\n",nptot,dnadd,dnpdel);

  return(1);
}
