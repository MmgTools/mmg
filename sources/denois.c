#include "mmg3d.h"

#define SHORT_MAX    0x7fff
#define  LLAMBDA     0.331
#define  LMU         0.33
#define ASIXTH        0.166667

extern Info  info;
extern char  ddb;

/** compute oriented quality of tetra defined by vertices (a,b,c,d)
   (return 0.0 when element is inverted) */
inline double orcal_poi(double a[3],double b[3],double c[3],double d[3]) {
  double     abx,aby,abz,acx,acy,acz,adx,ady,adz,bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     vol,v1,v2,v3,rap;

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

/** Attempts displacement of ALL nodes of the mesh at a position p+ t/SHORT_MAX * (optpos - p) */
int trydisp(pMesh mesh,double *optpos,short t){
  pTetra    pt;
  pPoint    p0;
  int       k,np,nd;
  double    cal,alpha;
  double    a[4][3],c[3];
  char      i;

  alpha = (double) t/SHORT_MAX;

  /* test volume degradation in each tet */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    for(i=0; i<4; i++){
      np = pt->v[i];
      p0 = &mesh->point[np];
      if(p0->flag){
        c[0] = p0->c[0] + alpha*(optpos[3*(np-1)+1] - p0->c[0]);
        c[1] = p0->c[1] + alpha*(optpos[3*(np-1)+2] - p0->c[1]);
        c[2] = p0->c[2] + alpha*(optpos[3*(np-1)+3] - p0->c[2]);
        memcpy(&(a[i][0]),c,3*sizeof(double));
      }
      else
        memcpy(&(a[i][0]),p0->c,3*sizeof(double));
    }

    cal = orcal_poi(&(a[0][0]),&(a[1][0]),&(a[2][0]),&(a[3][0]));
    if(cal < NULKAL)
      return(0);
  }

  /* update point coordinates */
  nd = 0;
  for(k=1; k<= mesh->np; k++){
    p0 = &mesh->point[k];
    if ( (!MG_VOK(p0)) || !p0->flag )  continue;
    nd++;

    c[0] = p0->c[0] + alpha*(optpos[3*(k-1)+1] - p0->c[0]);
    c[1] = p0->c[1] + alpha*(optpos[3*(k-1)+2] - p0->c[1]);
    c[2] = p0->c[2] + alpha*(optpos[3*(k-1)+3] - p0->c[2]);

    memcpy(p0->c,c,3*sizeof(double));
  }

  return(nd);
}

/** Find last valid position for the move of mesh nodes in attempt to reach optimal position
   by a dichotomic process */
int dichodisp(pMesh mesh,double *optpos){
  int   it,maxit,lastit,nd;
  short t,tm;

  it = 0;
  lastit = 0;
  maxit = 100;
  t = SHORT_MAX;
  tm = 0;

  nd = trydisp(mesh,optpos,t);
  if(nd){
    tm = t;
  }
  else{
    while(tm < SHORT_MAX && it < maxit){
      t = t >> 1;
      nd = trydisp(mesh,optpos,t);
      if(nd){
        tm += t;
        lastit = it;
      }
      else{
        if(lastit <= it-2)
          break;
      }
      it++;
    }
  }

  return(tm);
}

/** Generate a denoised optimal position by the laplacian/antilaplacian smoothing approach */
int lapantilap(pMesh mesh,double *optpos){
  pTetra      pt;
  pPoint      p0,p1;
  double      *intpos;
  int         nb,iel,np,nw,nw1,k;
  int         *adja,*w;
  char        i0,i1,ip,i,j;

  for(k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  nb = 0;

  /* Get number nb of points lying on the implicit boundary of mesh */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    adja = &mesh->adja[4*(k-1)+1];

    for(i=0; i<4; i++){
      iel = adja[i] / 4;
      if(!iel) continue;
      if(mesh->tetra[iel].ref == pt->ref) continue;

      for(j=0; j<3; j++){
        ip = idir[i][j];
        np = pt->v[ip];
        p0 = &mesh->point[np];
        if(p0->flag)
          continue;
        else{
          nb++;
          p0->flag = nb;
        }
      }
    }
  }

  w = (int*)calloc(nb+1,sizeof(int));
  assert(w);
  intpos = (double*)calloc(3*nb+1,sizeof(double));
  assert(intpos);

  /* Travel implicit boundary triangles and compute intermediate positions for Laplacian step ;
     each point is counted 4 times (twice for interior tets, twice for exterior) */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    adja = &mesh->adja[4*(k-1)+1];

    for(i=0; i<4; i++){
      iel = adja[i] / 4;
      if(!iel) continue;
      if(mesh->tetra[iel].ref == pt->ref) continue;

      for(j=0; j<3; j++){
        ip = idir[i][j];
        np = pt->v[ip];
        p0 = &mesh->point[np];

        nw = p0->flag;
        assert(nw);

        i0 = idir[i][inxt2[j]];
        i1 = idir[i][iprv2[j]];

        p1 = &mesh->point[pt->v[i0]];
        intpos[3*(nw-1)+1] += p1->c[0];
        intpos[3*(nw-1)+2] += p1->c[1];
        intpos[3*(nw-1)+3] += p1->c[2];
        w[nw]++;

        p1 = &mesh->point[pt->v[i1]];
        intpos[3*(nw-1)+1] += p1->c[0];
        intpos[3*(nw-1)+2] += p1->c[1];
        intpos[3*(nw-1)+3] += p1->c[2];
        w[nw]++;
      }
    }
  }

  for(k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    nw = p0->flag;
    if(!nw) continue;

    intpos[3*(nw-1)+1] = p0->c[0] + LLAMBDA*(intpos[3*(nw-1)+1]/w[nw] - p0->c[0]);
    intpos[3*(nw-1)+2] = p0->c[1] + LLAMBDA*(intpos[3*(nw-1)+2]/w[nw] - p0->c[1]);
    intpos[3*(nw-1)+3] = p0->c[2] + LLAMBDA*(intpos[3*(nw-1)+3]/w[nw] - p0->c[2]);
  }

  /* Final travel of implicit boundary and compute final positions for antilaplacian step */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    adja = &mesh->adja[4*(k-1)+1];

    for(i=0; i<4; i++){
      iel = adja[i] / 4;
      if(!iel) continue;
      if(mesh->tetra[iel].ref == pt->ref) continue;

      for(j=0; j<3; j++){
        ip = idir[i][j];
        np = pt->v[ip];

        i0 = idir[i][inxt2[j]];
        i1 = idir[i][iprv2[j]];

        p1 = &mesh->point[pt->v[i0]];
        nw1 = p1->flag;
        assert(nw1);

        optpos[3*(np-1)+1] += intpos[3*(nw1-1)+1];
        optpos[3*(np-1)+2] += intpos[3*(nw1-1)+2];
        optpos[3*(np-1)+3] += intpos[3*(nw1-1)+3];

        p1 = &mesh->point[pt->v[i1]];
        nw1 = p1->flag;
        assert(nw1);

        optpos[3*(np-1)+1] += intpos[3*(nw1-1)+1];
        optpos[3*(np-1)+2] += intpos[3*(nw1-1)+2];
        optpos[3*(np-1)+3] += intpos[3*(nw1-1)+3];
      }
    }
  }

  for(k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    if(!MG_VOK(p0)) continue;

    nw = p0->flag;

    if(nw){
      optpos[3*(k-1)+1] = p0->c[0] - LMU * (optpos[3*(k-1)+1]/w[nw] - p0->c[0]);
      optpos[3*(k-1)+2] = p0->c[1] - LMU * (optpos[3*(k-1)+2]/w[nw] - p0->c[1]);
      optpos[3*(k-1)+3] = p0->c[2] - LMU * (optpos[3*(k-1)+3]/w[nw] - p0->c[2]);
    }
    else{
      optpos[3*(k-1)+1] = p0->c[0];
      optpos[3*(k-1)+2] = p0->c[1];
      optpos[3*(k-1)+3] = p0->c[2];
    }
  }

  free(w);
  free(intpos);
  w=NULL;
  intpos=NULL;

  return(1);
}

/* Compute mean curvature vector at point np, whose surfacic ball is passed. Coordinates of np
   are passed, in case intermediate position of this point is simulated */
/*inline int meancur(pMesh mesh,int np,double c[3],int ilist,int *list,double h[3]){
  pTetra     pt;
  pPoint     p0,p1,p2;
  double     nt[3],n0[3],r[3][3],tAA[6],tAb[3],d[3],ux,uy,uz,kappa,dd;
  int        k,iel;
  char       iface,i0,i1,i2,j;

  ddb = np == 22316;
  memset(n0,0.0,3*sizeof(double));
  memset(h,0.0,3*sizeof(double));
*/
/* Approximation of normal at np (boulesurfvol enumerates the oriented ball) */
/* for(k=0; k<ilist; k++){
   iel = list[k] / 4;
   iface = list[k] %4;
   pt = &mesh->tetra[iel];

   for(j=0; j<3; j++){
   i0 = idir[iface][j];
   if(pt->v[i0] == np) break;
   }
   assert(j<3);

   i1 = idir[iface][inxt2[j]];
   i2 = idir[iface][iprv2[j]];

   p0 = &mesh->point[pt->v[i0]];
   p1 = &mesh->point[pt->v[i1]];
   p2 = &mesh->point[pt->v[i2]];

   nt[0] = (p1->c[1]-c[1])*(p2->c[2]-c[2]) - (p1->c[2]-c[2])*(p2->c[1]-c[1]);
   nt[1] = (p1->c[2]-c[2])*(p2->c[0]-c[0]) - (p1->c[0]-c[0])*(p2->c[2]-c[2]);
   nt[2] = (p1->c[0]-c[0])*(p2->c[1]-c[1]) - (p1->c[1]-c[1])*(p2->c[0]-c[0]);

   dd = nt[0]*nt[0] + nt[1]*nt[1] + nt[2]*nt[2];
   if(dd < EPSD) continue;
   dd = 1.0 / sqrt(dd);

   nt[0] *= dd;
   nt[1] *= dd;
   nt[2] *= dd;

   n0[0] += nt[0];
   n0[1] += nt[1];
   n0[2] += nt[2];
   }

   dd = n0[0]*n0[0] + n0[1]*n0[1] + n0[2]*n0[2];
   if(dd < EPSD)
   return(0);

   dd = 1.0 / sqrt(dd);
   n0[0] *= dd;
   n0[1] *= dd;
   n0[2] *= dd;
*/
/* rotation matrix that send n0 to third vector of canonical basis */
// rotmatrix(n0,r);

/* Fill in least-square matrix and vector */
/* memset(tAA,0.0,6*sizeof(double));
   memset(tAb,0.0,3*sizeof(double));

   for(k=0; k<ilist; k++){
   iel = list[k] / 4;
   iface = list[k] % 4;
   pt = &mesh->tetra[iel];

   for(j=0; j<3; j++){
   i0 = idir[iface][j];
   if(pt->v[i0] == np) break;
   }
   assert(j<3);

   i1 = idir[iface][inxt2[j]];

   p0 = &mesh->point[pt->v[i0]];
   p1 = &mesh->point[pt->v[i1]];

   ux = p1->c[0] - c[0];
   uy = p1->c[1] - c[1];
   uz = p1->c[2] - c[2];

   d[0] = r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
   d[1] = r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
   d[2] = r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

   if(ddb) printf("uz = %f \n",d[2]);

   tAA[0] += d[0]*d[0]*d[0]*d[0];
   tAA[1] += d[0]*d[0]*d[0]*d[1];
   tAA[2] += d[0]*d[0]*d[1]*d[1];
   tAA[3] += d[0]*d[0]*d[1]*d[1];
   tAA[4] += d[0]*d[1]*d[1]*d[1];
   tAA[5] += d[1]*d[1]*d[1]*d[1];

   tAb[0] += d[0]*d[0]*d[2];
   tAb[1] += d[0]*d[1]*d[2];
   tAb[2] += d[1]*d[1]*d[2];
   }

   tAA[1] *= 2.0;
   tAA[3] *= 4.0;
   tAA[4] *= 2.0;

   tAb[0] *= 2.0;
   tAb[1] *= 4.0;
   tAb[2] *= 2.0;

   if(!sys33sym(tAA,tAb,d))
   return(0);

   kappa = d[0] + d[2];
   h[0] = kappa * n0[0];
   h[1] = kappa * n0[1];
   h[2] = kappa * n0[2];

   if(ddb) printf("la courbure = %f %f %f alors que ilist = %D \n",h[0],h[1],h[2],ilist);

   return(1);
   }
*/
/** Compute mean curvature vector at point np, whose surfacic ball is passed. Coordinates of np
   are passed, in case intermediate position of this point is simulated */
inline int meancur(pMesh mesh,int np,double c[3],int ilist,int *list,double h[3]){
  pTetra     pt;
  pPoint     pa,pb;
  double     nt[3],ab[3],dd,area;
  int        k,iel;
  char       iface,j,i;

  memset(h,0.0,3*sizeof(double));
  area = 0.0;

  /* For each triangle of surface ball of p, compute local contribution to the mean
     curvature vector, relying on the discretization of the area's first variation formula */
  for(k=0; k<ilist; k++){
    iel = list[k] / 4;
    iface = list[k] %4;
    pt = &mesh->tetra[iel];

    for(j=0; j<3; j++){
      i = idir[iface][j];
      if(pt->v[i] == np)
        break;
    }
    assert( j < 3 );

    i = idir[iface][inxt2[j]];
    pa = &mesh->point[pt->v[i]];

    i = idir[iface][iprv2[j]];
    pb = &mesh->point[pt->v[i]];

    ab[0] = pb->c[0] - pa->c[0];
    ab[1] = pb->c[1] - pa->c[1];
    ab[2] = pb->c[2] - pa->c[2];

    /* nt = (b-a)^(p-a) / norm = normal to triangle */
    nt[0] = ab[1]*(c[2] - pa->c[2])\
      - ab[2]*(c[1] - pa->c[1]);
    nt[1] = ab[2]*(c[0] - pa->c[0])\
      - ab[0]*(c[2] - pa->c[2]);
    nt[2] = ab[0]*(c[1] - pa->c[1])\
      - ab[1]*(c[0] - pa->c[0]);

    dd = nt[0]*nt[0] + nt[1]*nt[1] + nt[2]*nt[2];

    if(dd < EPSD) continue;

    dd = sqrt(dd);
    area += (0.5*dd);

    dd = 1.0 / dd;
    nt[0] *= dd;
    nt[1] *= dd;
    nt[2] *= dd;

    h[0] += (nt[1]*ab[2] - nt[2]*ab[1]);
    h[1] += (nt[2]*ab[0] - nt[0]*ab[2]);
    h[2] += (nt[0]*ab[1] - nt[1]*ab[0]);
  }

  assert(area > EPSD);
  dd = 1.0 / area;

  h[0] *= dd;
  h[1] *= dd;
  h[2] *= dd;

  return(1);
}

/** Compute volume of the interior (ref MG_MINUS) part of a mesh */
inline double volint(pMesh mesh){
  pTetra      pt;
  pPoint      p0,p1,p2,p3;
  int         k;
  double      vol,voltet;

  vol = 0.0;

  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if(!MG_EOK(pt)) continue;
    if(pt->ref != MG_MINUS) continue;

    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    p3 = &mesh->point[pt->v[3]];

    voltet = ASIXTH*det4pt(p0->c,p1->c,p2->c,p3->c);
    vol += voltet;
  }

  return(vol);
}

/** Compute (unoriented) area of face iface in tetra iel */
inline double surftri(pMesh mesh, int iel, int iface){
  pTetra      pt;
  pPoint      pa,pb,pc;
  double      area;
  double      ab[3],ac[3],pv[3];

  pt = &mesh->tetra[iel];
  pa = &mesh->point[pt->v[idir[iface][0]]];
  pb = &mesh->point[pt->v[idir[iface][1]]];
  pc = &mesh->point[pt->v[idir[iface][2]]];

  ab[0] = pb->c[0] - pa->c[0];
  ab[1] = pb->c[1] - pa->c[1];
  ab[2] = pb->c[2] - pa->c[2];

  ac[0] = pc->c[0] - pa->c[0];
  ac[1] = pc->c[1] - pa->c[1];
  ac[2] = pc->c[2] - pa->c[2];

  pv[0] = ab[1]*ac[2] - ab[2]*ac[1];
  pv[1] = ab[2]*ac[0] - ab[0]*ac[2];
  pv[2] = ab[0]*ac[1] - ab[1]*ac[0];

  area = pv[0]*pv[0] + pv[1]*pv[1] + pv[2]*pv[2];
  area = 0.5*sqrt(area);

  return(area);
}

/** Dimension time step of mean curvature flow so that loss of perimeter between initial and
   deformed configurations is no more than defrate (in %) */
double timestepMCF(pMesh mesh,double defrate){
  pTetra      pt;
  pPoint      p0;
  double      dt,dd,locarea,area,intk2;
  double      h[3];
  int         k,l,iel,jel,np,ilistv,ilists;
  int         *adja,listv[LMAX+2],lists[LMAX+2];
  char        i,j,jface,ip;

  area = intk2 = 0.0;

  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    adja = &mesh->adja[4*(k-1)+1];

    for(i=0; i<4; i++){
      iel = adja[i] / 4;
      if(!iel) continue;
      if(mesh->tetra[iel].ref == pt->ref) continue;

      area += surftri(mesh,k,i);

      for(j=0; j<3; j++){
        ip = idir[i][j];
        np = pt->v[ip];
        p0 = &mesh->point[np];

        if(p0->flag) continue;

        p0->flag = 1;

        if(!boulesurfvolp(mesh,k,ip,i,listv,&ilistv,lists,&ilists)){
          printf("%s:%d: Error: function boulesurfvolp return 0\n",
                 __FILE__,__LINE__);
          exit(0);
        }
        if(!meancur(mesh,np,p0->c,ilists,lists,h)){
          printf("%s:%d: Error: function meancur return 0\n",__FILE__,__LINE__);
          exit(0);
        }

        locarea = 0.0;
        for(l=0; l<ilists; l++){
          jel = lists[l] / 4;
          jface = lists[l] % 4;
          locarea += surftri(mesh,jel,jface);
        }

        dd = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
        intk2 += locarea * dd;
      }
    }
  }

  area *= 0.5;
  intk2 *= ATHIRD;

  dt = defrate*area/intk2;

  return(dt);
}

/** Move each point of the implicit boundary mesh according to a mean curvature flow / anti
   mean curvature flow approach */
int bdyMCF(pMesh mesh){
  pTetra      pt,pt0,pt1;
  pPoint      p0,p1,ppt0;
  int         k,l,iel,jel,maxit,it,base,np,ilistv,ilists,nbdy,nint;
  int         *adja,listv[LMAX+2],lists[LMAX+2];
  double      defrate,dt,vol,totvol,volbeg,volend;
  double      h[3],o[3];
  char        i,j,ip,jj,isok;

  defrate = 0.002;
  maxit = 10;
  it = 0;
  base = 2;
  nint = nbdy = 0;

  pt0 = &mesh->tetra[0];
  ppt0 = &mesh->point[0];

  volbeg = volint(mesh);

  /* reset point flags */
  for(k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Initialize time step */
  dt = timestepMCF(mesh,defrate);
  printf("Le pas de temps : %f \n",dt);

  /* reset point flags */
  for(k=1; k<=mesh->ne; k++)
    mesh->tetra[k].flag = 0;

  for(k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Tag points of the external boundary to -2, and -1 for points of the implicit boundary */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    adja = &mesh->adja[4*(k-1)+1];
    for(i=0; i<4; i++){
      iel = adja[i] / 4;
      if(!iel){
        for(j=0; j<3; j++){
          ip = idir[i][j];
          p0 = &mesh->point[pt->v[ip]];
          p0->flag = -2;
        }
      }
      else if(mesh->tetra[iel].ref != pt->ref){
        for(j=0; j<3; j++){
          ip = idir[i][j];
          p0 = &mesh->point[pt->v[ip]];
          if(p0->flag != -2) p0->flag = -1;
        }
      }
    }
  }

  /* Main loop each point of the implicit boundary is moved with the MCF/antiMCF approach,
     and each internal point is moved with the barycentric approach */
  do{
    base++;

    for(k=1; k<=mesh->ne; k++){
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;

      adja = &mesh->adja[4*(k-1)+1];
      for(i=0; i<4; i++){
        iel = adja[i] / 4;
        if(!iel) continue;

        /* Case of an implicit boundary face */
        if(mesh->tetra[iel].ref != pt->ref){
          for(j=0; j<3; j++){
            ip = idir[i][j];
            np = pt->v[ip];

            p0 = &mesh->point[np];
            if(p0->flag == -2 || p0->flag == -base) continue;

            p0->flag = -base;
            if(!boulesurfvolp(mesh,k,ip,i,listv,&ilistv,lists,&ilists)){
              printf("%s:%d: Error: function boulesurfvolp return 0\n",
                     __FILE__,__LINE__);
              exit(0);
            }

            /* Forward step of mean curvature flow */
            if(!meancur(mesh,np,p0->c,ilists,lists,h)){
              printf("%s:%d: Error: function meancur return 0\n",
                     __FILE__,__LINE__);
              exit(0);
            }
            o[0] = p0->c[0] - dt*h[0];
            o[1] = p0->c[1] - dt*h[1];
            o[2] = p0->c[2] - dt*h[2];

            /* Backward step from intermediate position */
            if(!  meancur(mesh,np,o,ilists,lists,h)){
              printf("%s:%d: Error: function meancur return 0\n",
                     __FILE__,__LINE__);
              exit(0);
            }

            o[0] += dt*h[0];
            o[1] += dt*h[1];
            o[2] += dt*h[2];

            /* Check validity of resulting position */
            isok = 1;
            for(l=0; l<ilistv; l++){
              jel = listv[l] / 4;
              jj = listv[l] % 4;

              pt1 = &mesh->tetra[jel];
              memcpy(pt0,pt1,sizeof(Tetra));
              memcpy(ppt0,p0,sizeof(Point));
              memcpy(ppt0->c,o,3*sizeof(double));

              if(orcal(mesh,0)<NULKAL){
                isok = 0;
                break;
              }
            }

            if(isok){
              nbdy++;
              memcpy(p0->c,o,3*sizeof(double));
            }
          }
        }

        /* Case of an internal face */
        else{
          for(j=0; j<3; j++){
            ip = idir[i][j];
            np = pt->v[ip];

            p0 = &mesh->point[np];
            if(p0->flag < 0.0 || p0->flag == base) continue;

            p0->flag = base;

            ilistv = boulevolp(mesh,k,ip,listv);

            memset(o,0.0,3*sizeof(double));
            totvol = 0.0;
            for(l=0; l<ilistv; l++){
              jel = listv[l] / 4;
              pt1 = &mesh->tetra[jel];

              vol= det4pt(mesh->point[pt1->v[0]].c,mesh->point[pt1->v[1]].c,\
                          mesh->point[pt1->v[2]].c,mesh->point[pt1->v[3]].c);

              totvol += vol;

              for(jj=0; jj<4; jj++){
                p1 = &mesh->point[pt1->v[jj]];
                o[0] += 0.25*vol*p1->c[0];
                o[1] += 0.25*vol*p1->c[1];
                o[2] += 0.25*vol*p1->c[2];
              }
            }
            if(totvol < EPSD2) continue;
            totvol = 1.0 / totvol;

            o[0] *= totvol;
            o[1] *= totvol;
            o[2] *= totvol;

            /* Check validity of resulting position */
            isok = 1;
            for(l=0; l<ilistv; l++){
              jel = listv[l] / 4;
              jj = listv[l] % 4;

              pt1 = &mesh->tetra[jel];
              memcpy(pt0,pt1,sizeof(Tetra));
              memcpy(ppt0,p0,sizeof(Point));
              memcpy(ppt0->c,o,3*sizeof(double));

              if(orcal(mesh,0)<NULKAL){
                isok = 0;
                break;
              }
            }

            if(isok){
              nint++;
              memcpy(p0->c,o,3*sizeof(double));
            }
          }
        }
      }
    }

    //if ( abs(info.imprim) < 5 && (nint > 0 || nbdy > 0) )
    fprintf(stdout,"     %8d boundary points moved, %8d internal points moved \n",nbdy,nint);

  }
  while(++it<maxit);

  volend = volint(mesh);
  printf("Loss of volume = %E percent\n",100.0*fabs(volbeg - volend)/volbeg);

  return(1);
}

/** Extend displacement field, defined only at boundary points, to all vertices of mesh */
int ppgdisp(pMesh mesh,double *optpos){
  pTetra     pt;
  pPoint     p0,p1;
  double     l,r,norm,inorm;
  double     v[3];
  int        k,iel,np,np1,it,maxit,base,newbase;
  int        *pile,*adja,*w;
  char       i,j,ip,ip1;

  for(k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  for(k=1; k<=mesh->ne; k++){
    mesh->tetra[k].flag = 0;
  }

  pile = (int*)calloc(mesh->ne+1,sizeof(int));
  assert(pile);

  w = (int*)calloc(mesh->np+1,sizeof(int));
  assert(w);

  it = 0;
  maxit = 10;
  base = 0;
  r = 0.4;

  /* Tag points lying on implicit boundary */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    adja = &mesh->adja[4*(k-1)+1];

    for(i=0; i<4; i++){
      iel = adja[i] / 4;
      if(!iel) continue;
      if(mesh->tetra[iel].ref != pt->ref) continue;

      for(j=0; j<3; j++){
        ip = idir[i][j];
        np = pt->v[ip];
        p0 = &mesh->point[np];

        if(!p0->flag){
          p0->flag = 1;
          w[np] = 1;
        }
      }
    }
  }

  /* Main loop : diffuse given vector field from the boundary */
  do{

    base++;
    newbase = base + 1;

    for(k=1; k<=mesh->np; k++){
      p0 = &mesh->point[k];
      if ( !MG_VOK(p0) )  continue;
      if(p0->flag != base) continue;

      optpos[3*(k-1)+1] = optpos[3*(k-1)+1] / w[k];
      optpos[3*(k-1)+2] = optpos[3*(k-1)+2] / w[k];
      optpos[3*(k-1)+3] = optpos[3*(k-1)+3] / w[k];
    }

    for(k=1; k<=mesh->ne; k++){
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      if(pt->flag) continue;

      for(i=0; i<4; i++){
        np = pt->v[i];
        p0 = &mesh->point[np];
        if(p0->flag == base)
          break;
      }

      if(i == 4) continue;

      pt->flag = base;

      for(i=0; i<4; i++){
        np = pt->v[i];
        p0 = &mesh->point[np];
        if(p0->flag) continue;

        ip1 = i;
        for(j=0; j<3; j++){
          ip1 = inxt3[ip1];
          np1 = pt->v[ip1];
          p1 = &mesh->point[np1];

          if(!p1->flag || p1->flag == newbase) continue;

          w[np]++;

          v[0] = optpos[3*(np1-1)+1] - p1->c[0];
          v[1] = optpos[3*(np1-1)+2] - p1->c[1];
          v[2] = optpos[3*(np1-1)+3] - p1->c[2];

          norm = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
          norm = sqrt(norm);

          if(norm < EPSD) continue;
          inorm = 1.0 / norm;
          v[0] *= inorm;
          v[1] *= inorm;
          v[2] *= inorm;

          l = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1])\
            + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);

          l = MG_MAX(0.0,norm - r*l);

          optpos[3*(np-1)+1] += p0->c[0] + l*v[0];
          optpos[3*(np-1)+2] += p0->c[1] + l*v[1];
          optpos[3*(np-1)+3] += p0->c[2] + l*v[2];
        }

        p0->flag = newbase;
      }
    }
  }
  while( ++it < maxit );

  return(1);
}

/** Denoise boundary mesh */
int denoisbdy(pMesh mesh){
  int       ier,it,maxit;
  short     tm;
  double    alpha;
  double    *optpos;

  maxit = 2;
  optpos = (double*)calloc(3*mesh->np+1,sizeof(double));

  for(it=0; it<maxit; it++){

    ier = lapantilap(mesh,optpos);
    assert(ier);

    //ier = ppgdisp(mesh,optpos);
    assert(ier);

    tm = dichodisp(mesh,optpos);
    alpha = (double) tm/SHORT_MAX;

    memset(optpos,0,(3*mesh->np+1)*sizeof(double));

    //if ( abs(info.imprim) > 5 || info.ddebug )
    printf("     Mesh denoising : iteration %d %f from optimal position.\n",it,alpha);

    if(tm == 0)  break;
  }

  free(optpos);
  optpos=NULL;
  return(1);
}
