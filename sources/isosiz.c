#include "mmg3d.h"

extern Info   info;
extern char   ddb;

#define MAXLEN    1.0e9
#define A64TH     0.015625
#define A16TH     0.0625
#define A32TH     0.03125


/** Define isotropic size at regular point nump, whose surfacic ball is provided */
static double defsizreg(pMesh mesh,pSol met,int nump,int *lists,int ilists, double hausd) {
  pTetra       pt;
  pxTetra      pxt;
  pPoint       p0,p1;
  Tria         tt;
  Bezier       b;
  double       ux,uy,uz,det2d,h,isqhmin,isqhmax,ll,lmin,lmax,hnm,s;
  double       *n,*t,r[3][3],lispoi[3*LMAX+1],intm[3],b0[3],b1[3],c[3],tAA[6],tAb[3],d[3];
  double       kappa[2],vp[2][2];
  int          k,na,nb,ntempa,ntempb,iel,ip0;
  char         iface,i,j,i0;

  p0 = &mesh->point[nump];

  if ( !p0->xp || MG_EDG(p0->tag) || (p0->tag & MG_NOM) || (p0->tag & MG_REQ))  {
    fprintf(stdout,"    ## Func. defsizreg : wrong point qualification : xp ? %d\n",p0->xp);
    return(0);
  }
  isqhmin = 1.0 / (info.hmin*info.hmin);
  isqhmax = 1.0 / (info.hmax*info.hmax);

  n = &mesh->xpoint[p0->xp].n1[0];

  /* Step 1 : rotation matrix that sends normal n to the third coordinate vector of R^3 */
  rotmatrix(n,r);

  /* Step 2 : rotation of the oriented surfacic ball with r : lispoi[k] is the common edge
     between faces lists[k-1] and lists[k] */
  iel   = lists[0] / 4;
  iface = lists[0] % 4;
  pt    = &mesh->tetra[iel];
  lmin  = MAXLEN;
  lmax  = 0.0;

  na = nb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[idir[iface][i]] != nump ) {
      if ( !na )
        na = pt->v[idir[iface][i]];
      else
        nb = pt->v[idir[iface][i]];
    }
  }

  for (k=1; k<ilists; k++) {
    iel   = lists[k] / 4;
    iface = lists[k] % 4;
    pt    = &mesh->tetra[iel];
    ntempa = ntempb = 0;
    for (i=0; i<3; i++) {
      if ( pt->v[idir[iface][i]] != nump ) {
        if ( !ntempa )
          ntempa = pt->v[idir[iface][i]];
        else
          ntempb = pt->v[idir[iface][i]];
      }
    }
    if ( ntempa == na )
      p1 = &mesh->point[na];
    else if ( ntempa == nb )
      p1 = &mesh->point[nb];
    else if ( ntempb == na )
      p1 = &mesh->point[na];
    else {
      assert(ntempb == nb);
      p1 = &mesh->point[nb];
    }
    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

    ll = lispoi[3*k+1]*lispoi[3*k+1] + lispoi[3*k+2]*lispoi[3*k+2] + lispoi[3*k+3]*lispoi[3*k+3];
    lmin = MG_MIN(lmin,ll);
    lmax = MG_MAX(lmax,ll);

    na = ntempa;
    nb = ntempb;
  }

  /* Finish with point 0 */
  iel   = lists[0] / 4;
  iface = lists[0] % 4;
  pt    = &mesh->tetra[iel];
  ntempa = ntempb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[idir[iface][i]] != nump ) {
      if ( !ntempa )
        ntempa = pt->v[idir[iface][i]];
      else
        ntempb = pt->v[idir[iface][i]];
    }
  }
  if ( ntempa == na )
    p1 = &mesh->point[na];
  else if ( ntempa == nb )
    p1 = &mesh->point[nb];
  else if ( ntempb == na )
    p1 = &mesh->point[na];
  else {
    assert(ntempb == nb);
    p1 = &mesh->point[nb];
  }

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  lispoi[1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
  lispoi[2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
  lispoi[3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

  ll = lispoi[1]*lispoi[1] + lispoi[2]*lispoi[2] + lispoi[3]*lispoi[3];
  lmin = MG_MIN(lmin,ll);
  lmax = MG_MAX(lmax,ll);

  /* list goes modulo ilist */
  lispoi[3*ilists+1] = lispoi[1];
  lispoi[3*ilists+2] = lispoi[2];
  lispoi[3*ilists+3] = lispoi[3];

  /* At this point, lispoi contains the oriented surface ball of point p0, that has been rotated
     through r, with the convention that triangle l has edges lispoi[l]; lispoi[l+1] */
  if ( lmax/lmin > 4.0*info.hmax*info.hmax/
       (info.hmin*info.hmin) )  return(info.hmax);

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( det2d < 0.0 )  return(info.hmax);
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  if ( det2d < 0.0 )    return(info.hmax);

  /* Reconstitution of the curvature tensor at p0 in the tangent plane,
     with a quadric fitting approach */
  memset(intm,0.0,3*sizeof(double));
  memset(tAA,0.0,6*sizeof(double));
  memset(tAb,0.0,3*sizeof(double));

  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    iface = lists[k] % 4;

    tet2tri(mesh,iel,iface,&tt);

    pxt   = &mesh->xtetra[mesh->tetra[iel].xt];
    if ( !bezierCP(mesh,&tt,&b,MG_GET(pxt->ori,iface)) ) {
      fprintf(stdout,"%s:%d: Error: function bezierCP return 0\n",
              __FILE__,__LINE__);
      exit(EXIT_FAILURE);
    }

    for (i0=0; i0<3; i0++) {
      if ( tt.v[i0] == nump )  break;
    }
    assert(i0 < 3);

    for (j=0; j<10; j++) {
      c[0] = b.b[j][0] - p0->c[0];
      c[1] = b.b[j][1] - p0->c[1];
      c[2] = b.b[j][2] - p0->c[2];

      b.b[j][0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
      b.b[j][1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
      b.b[j][2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];
    }

    /* Mid-point along left edge and endpoint in the rotated frame */
    if ( i0 == 0 ) {
      memcpy(b0,&(b.b[7][0]),3*sizeof(double));
      memcpy(b1,&(b.b[8][0]),3*sizeof(double));
    }
    else if ( i0 == 1 ) {
      memcpy(b0,&(b.b[3][0]),3*sizeof(double));
      memcpy(b1,&(b.b[4][0]),3*sizeof(double));
    }
    else {
      memcpy(b0,&(b.b[5][0]),3*sizeof(double));
      memcpy(b1,&(b.b[6][0]),3*sizeof(double));
    }
    s = 0.5;

    /* At this point, the two control points are expressed in the rotated frame */
    c[0] = 3.0*s*(1.0-s)*(1.0-s)*b0[0] + 3.0*s*s*(1.0-s)*b1[0] + s*s*s*lispoi[3*k+1];
    c[1] = 3.0*s*(1.0-s)*(1.0-s)*b0[1] + 3.0*s*s*(1.0-s)*b1[1] + s*s*s*lispoi[3*k+2];
    c[2] = 3.0*s*(1.0-s)*(1.0-s)*b0[2] + 3.0*s*s*(1.0-s)*b1[2] + s*s*s*lispoi[3*k+3];

    /* Fill matric tAA and second member tAb*/
    tAA[0] += c[0]*c[0]*c[0]*c[0];
    tAA[1] += c[0]*c[0]*c[1]*c[1];
    tAA[2] += c[0]*c[0]*c[0]*c[1];
    tAA[3] += c[1]*c[1]*c[1]*c[1];
    tAA[4] += c[0]*c[1]*c[1]*c[1];
    tAA[5] += c[0]*c[0]*c[1]*c[1];

    tAb[0] += c[0]*c[0]*c[2];
    tAb[1] += c[1]*c[1]*c[2];
    tAb[2] += c[0]*c[1]*c[2];

    s = 1.0;
    /* At this point, the two control points are expressed in the rotated frame */
    c[0] = 3.0*s*(1.0-s)*(1.0-s)*b0[0] + 3.0*s*s*(1.0-s)*b1[0] + s*s*s*lispoi[3*k+1];
    c[1] = 3.0*s*(1.0-s)*(1.0-s)*b0[1] + 3.0*s*s*(1.0-s)*b1[1] + s*s*s*lispoi[3*k+2];
    c[2] = 3.0*s*(1.0-s)*(1.0-s)*b0[2] + 3.0*s*s*(1.0-s)*b1[2] + s*s*s*lispoi[3*k+3];

    /* Fill matric tAA and second member tAb*/
    tAA[0] += c[0]*c[0]*c[0]*c[0];
    tAA[1] += c[0]*c[0]*c[1]*c[1];
    tAA[2] += c[0]*c[0]*c[0]*c[1];
    tAA[3] += c[1]*c[1]*c[1]*c[1];
    tAA[4] += c[0]*c[1]*c[1]*c[1];
    tAA[5] += c[0]*c[0]*c[1]*c[1];

    tAb[0] += c[0]*c[0]*c[2];
    tAb[1] += c[1]*c[1]*c[2];
    tAb[2] += c[0]*c[1]*c[2];

    /* Mid-point along median edge and endpoint in the rotated frame */
    if ( i0 == 0 ) {
      c[0] = A64TH*(b.b[1][0] + b.b[2][0] + 3.0*(b.b[3][0] + b.b[4][0])) \
        + 3.0*A16TH*(b.b[6][0] + b.b[7][0] + b.b[9][0]) + A32TH*(b.b[5][0] + b.b[8][0]);
      c[1] = A64TH*(b.b[1][1] + b.b[2][1] + 3.0*(b.b[3][1] + b.b[4][1])) \
        + 3.0*A16TH*(b.b[6][1] + b.b[7][1] + b.b[9][1]) + A32TH*(b.b[5][1] + b.b[8][1]);
      c[2] = A64TH*(b.b[1][2] + b.b[2][2] + 3.0*(b.b[3][2] + b.b[4][2])) \
        + 3.0*A16TH*(b.b[6][2] + b.b[7][2] + b.b[9][2]) + A32TH*(b.b[5][2] + b.b[8][2]);

      d[0] = 0.125*b.b[1][0] + 0.375*(b.b[3][0] + b.b[4][0]) + 0.125*b.b[2][0];
      d[1] = 0.125*b.b[1][1] + 0.375*(b.b[3][1] + b.b[4][1]) + 0.125*b.b[2][1];
      d[2] = 0.125*b.b[1][2] + 0.375*(b.b[3][2] + b.b[4][2]) + 0.125*b.b[2][2];
    }
    else if (i0 == 1) {
      c[0] = A64TH*(b.b[0][0] + b.b[2][0] + 3.0*(b.b[5][0] + b.b[6][0])) \
        + 3.0*A16TH*(b.b[3][0] + b.b[8][0] + b.b[9][0]) + A32TH*(b.b[4][0] + b.b[7][0]);
      c[1] = A64TH*(b.b[0][1] + b.b[2][1] + 3.0*(b.b[5][1] + b.b[6][1])) \
        + 3.0*A16TH*(b.b[3][1] + b.b[8][1] + b.b[9][1]) + A32TH*(b.b[4][1] + b.b[7][1]);
      c[2] = A64TH*(b.b[0][2] + b.b[2][2] + 3.0*(b.b[5][2] + b.b[6][2])) \
        + 3.0*A16TH*(b.b[3][2] + b.b[8][2] + b.b[9][2]) + A32TH*(b.b[4][2] + b.b[7][2]);

      d[0] = 0.125*b.b[2][0] + 0.375*(b.b[5][0] + b.b[6][0]) + 0.125*b.b[0][0];
      d[1] = 0.125*b.b[2][1] + 0.375*(b.b[5][1] + b.b[6][1]) + 0.125*b.b[0][1];
      d[2] = 0.125*b.b[2][2] + 0.375*(b.b[5][2] + b.b[6][2]) + 0.125*b.b[0][2];
    }
    else {
      c[0] = A64TH*(b.b[0][0] + b.b[1][0] + 3.0*(b.b[7][0] + b.b[8][0])) \
        + 3.0*A16TH*(b.b[4][0] + b.b[5][0] + b.b[9][0]) + A32TH*(b.b[3][0] + b.b[6][0]);
      c[1] = A64TH*(b.b[0][1] + b.b[1][1] + 3.0*(b.b[7][1] + b.b[8][1])) \
        + 3.0*A16TH*(b.b[4][1] + b.b[5][1] + b.b[9][1]) + A32TH*(b.b[3][1] + b.b[6][1]);
      c[2] = A64TH*(b.b[0][2] + b.b[1][2] + 3.0*(b.b[7][2] + b.b[8][2])) \
        + 3.0*A16TH*(b.b[4][2] + b.b[5][2] + b.b[9][2]) + A32TH*(b.b[3][2] + b.b[6][2]);

      d[0] = 0.125*b.b[0][0] + 0.375*(b.b[7][0] + b.b[8][0]) + 0.125*b.b[1][0];
      d[1] = 0.125*b.b[0][1] + 0.375*(b.b[7][1] + b.b[8][1]) + 0.125*b.b[1][1];
      d[2] = 0.125*b.b[0][2] + 0.375*(b.b[7][2] + b.b[8][2]) + 0.125*b.b[1][2];
    }

    /* Fill matric tAA and second member tAb*/
    tAA[0] += c[0]*c[0]*c[0]*c[0];
    tAA[1] += c[0]*c[0]*c[1]*c[1];
    tAA[2] += c[0]*c[0]*c[0]*c[1];
    tAA[3] += c[1]*c[1]*c[1]*c[1];
    tAA[4] += c[0]*c[1]*c[1]*c[1];
    tAA[5] += c[0]*c[0]*c[1]*c[1];

    tAb[0] += c[0]*c[0]*c[2];
    tAb[1] += c[1]*c[1]*c[2];
    tAb[2] += c[0]*c[1]*c[2];

    tAA[0] += d[0]*d[0]*d[0]*d[0];
    tAA[1] += d[0]*d[0]*d[1]*d[1];
    tAA[2] += d[0]*d[0]*d[0]*d[1];
    tAA[3] += d[1]*d[1]*d[1]*d[1];
    tAA[4] += d[0]*d[1]*d[1]*d[1];
    tAA[5] += d[0]*d[0]*d[1]*d[1];

    tAb[0] += d[0]*d[0]*d[2];
    tAb[1] += d[1]*d[1]*d[2];
    tAb[2] += d[0]*d[1]*d[2];
  }

  /* solve now (a b c) = tAA^{-1} * tAb */
  if ( !sys33sym(tAA,tAb,c) )  return(info.hmax);

  intm[0] = 2.0*c[0];
  intm[1] = c[2];
  intm[2] = 2.0*c[1];

  /* At this point, intm stands for the integral matrix of Taubin's approach : vp[0] and vp[1]
     are the two pr. directions of curvature, and the two curvatures can be inferred from lambdas*/
  if( !eigensym(intm,kappa,vp) ){
    fprintf(stdout,"%s:%d: Error: function eigensym return 0\n",
            __FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }

  /* h computation : h(x) = sqrt( 9*hausd / (2 * max(kappa1(x),kappa2(x)) ) */
  kappa[0] = 2.0/9.0 * fabs(kappa[0]) / hausd;
  kappa[0] = MG_MIN(kappa[0],isqhmin);
  kappa[0] = MG_MAX(kappa[0],isqhmax);

  kappa[1] = 2.0/9.0 * fabs(kappa[1]) / hausd;
  kappa[1] = MG_MIN(kappa[1],isqhmin);
  kappa[1] = MG_MAX(kappa[1],isqhmax);

  kappa[0] = 1.0 / sqrt(kappa[0]);
  kappa[1] = 1.0 / sqrt(kappa[1]);

  h = MG_MIN(kappa[0],kappa[1]);

  /* Travel surfacic ball one last time and update non manifold point metric */
  for (k=0; k<ilists; k++) {
    iel = lists[k] / 4;
    iface = lists[k] % 4;

    for (j=0; j<3; j++) {
      i0  = idir[iface][j];
      ip0 = pt->v[i0];
      p1  = &mesh->point[ip0];
      if( !(p1->tag & MG_NOM) || MG_SIN(p1->tag) ) continue;
      assert(p1->xp);
      t = &mesh->xpoint[p1->xp].t[0];
      memcpy(c,t,3*sizeof(double));

      d[0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
      d[1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];

      hnm = intm[0]*d[0]*d[0] + 2.0*intm[1]*d[0]*d[1] + intm[2]*d[1]*d[1];
      hnm = 2.0/9.0 * fabs(hnm) / hausd;
      hnm = MG_MIN(hnm,isqhmin);
      hnm = MG_MAX(hnm,isqhmax);
      hnm = 1.0 / sqrt(hnm);
      met->m[ip0] = MG_MIN(met->m[ip0],hnm);
    }
  }
  return(h);
}

/** Define isotropic size map at all boundary vertices of the mesh,
    associated with geometric approx, and prescribe hmax at the internal vertices
    Field h of Point is used, to store the prescribed size (not inverse, squared,...) */
int defsiz_iso(pMesh mesh,pSol met) {
  pTetra    pt;
  pxTetra   pxt;
  pPoint    p0,p1;
  double    hp,v[3],b0[3],b1[3],b0p0[3],b1b0[3],p1b1[3],hausd;
  double    secder0[3],secder1[3],kappa,tau[3],gammasec[3],ntau2,intau,ps,lm,*n;
  int       lists[LMAX+2],listv[LMAX+2],ilists,ilistv,k,ip0,ip1,l;
  char      i,j,ia,ised,i0,i1;
  pPar      par;

  if ( abs(info.imprim) > 5 || info.ddebug )
    fprintf(stdout,"  ** Defining map\n");

  if ( info.hmax < 0.0 )  info.hmax = 0.5 * info.delta;

  /* alloc structure */
  if ( !met->m ) {
    met->np    = mesh->np;
    met->npmax = mesh->npmax;
    met->size  = 1;
    met->dim   = 3;
    met->m = (double*)malloc((mesh->npmax+1)*sizeof(double));
    if ( !met->m ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    /* init constant size */
    for (k=1; k<=mesh->np; k++)
      met->m[k] = info.hmax;
  }
  else {
    /* size truncation */
    for (k=1; k<=mesh->np; k++)
      met->m[k] = MG_MIN(info.hmax,MG_MAX(info.hmin,met->m[k]));
  }

  /* size at regular surface points */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
    else if ( !pt->xt )  continue;

    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<4; i++) {
      if ( !(pxt->ftag[i] & MG_BDY) ) continue;
      /* local hausdorff for triangle */
      hausd = info.hausd;
      for (l=0; l<info.npar; l++) {
        par = &info.par[l];
        if ( (par->elt == MMG5_Triangle) && (pxt->ref[i] == par->ref ) )
          hausd = par->hausd;
      }

      for (j=0; j<3; j++) {
        i0  = idir[i][j];
        ip0 = pt->v[i0];
        p0  = &mesh->point[ip0];

        if ( MG_SIN(p0->tag) || MG_EDG(p0->tag) || (p0->tag & MG_NOM) ) continue;
        if ( !boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists) )  continue;

        n   = &mesh->xpoint[p0->xp].n1[0];
        directsurfball(mesh,ip0,lists,ilists,n);
        hp  = defsizreg(mesh,met,ip0,lists,ilists,hausd);
        met->m[ip0] = MG_MIN(met->m[ip0],hp);
      }
    }
  }

  /* Travel all boundary faces to update size prescription for points on ridges/edges */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    else if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++) {
      if ( !(pxt->ftag[i] & MG_BDY) )  continue;
      else if ( !norface(mesh,k,i,v) )  continue;

      /* local hausdorff for triangle */
      hausd = info.hausd;
      for (l=0; l<info.npar; l++) {
        par = &info.par[l];
        if ( (par->elt == MMG5_Triangle) && (pxt->ref[i] == par->ref ) )
          hausd = par->hausd;
      }

      for (j=0; j<3; j++) {
        ia = iarf[i][j];
        i0 = iare[ia][0];
        i1 = iare[ia][1];
        ip0 = pt->v[i0];
        ip1 = pt->v[i1];
        p0  = &mesh->point[ip0];
        p1  = &mesh->point[ip1];
        if ( !MG_EDG(p0->tag) && !MG_EDG(p1->tag) )  continue;

        ised = MG_EDG(pxt->tag[ia]) || ( pxt->tag[ia] & MG_NOM );

        BezierEdge(mesh,ip0,ip1,b0,b1,ised,v);

        b0p0[0] = b0[0] - p0->c[0];
        b0p0[1] = b0[1] - p0->c[1];
        b0p0[2] = b0[2] - p0->c[2];

        b1b0[0] = b1[0] - b0[0];
        b1b0[1] = b1[1] - b0[1];
        b1b0[2] = b1[2] - b0[2];

        p1b1[0] = p1->c[0] - b1[0];
        p1b1[1] = p1->c[1] - b1[1];
        p1b1[2] = p1->c[2] - b1[2];

        secder0[0] = p0->c[0] + b1[0] - 2.0*b0[0];
        secder0[1] = p0->c[1] + b1[1] - 2.0*b0[1];
        secder0[2] = p0->c[2] + b1[2] - 2.0*b0[2];

        secder1[0] = p1->c[0] + b0[0] - 2.0*b1[0];
        secder1[1] = p1->c[1] + b0[1] - 2.0*b1[1];
        secder1[2] = p1->c[2] + b0[2] - 2.0*b1[2];

        kappa = 0.0;
        for (l=0; l<4; l++) {
          tau[0] = 3.0*(1.0-ATHIRD*l)*(1.0-ATHIRD*l)*b0p0[0] + 6.0*ATHIRD*l*(1.0-ATHIRD*l)*b1b0[0]\
            + 3.0*ATHIRD*l*ATHIRD*l*p1b1[0];
          tau[1] = 3.0*(1.0-ATHIRD*l)*(1.0-ATHIRD*l)*b0p0[1] + 6.0*ATHIRD*l*(1.0-ATHIRD*l)*b1b0[1]\
            + 3.0*ATHIRD*l*ATHIRD*l*p1b1[1];
          tau[2] = 3.0*(1.0-ATHIRD*l)*(1.0-ATHIRD*l)*b0p0[2] + 6.0*ATHIRD*l*(1.0-ATHIRD*l)*b1b0[2]\
            + 3.0*ATHIRD*l*ATHIRD*l*p1b1[2];

          gammasec[0] = 6.0*((1.0-ATHIRD*l)*secder0[0] + ATHIRD*l*secder1[0]);
          gammasec[1] = 6.0*((1.0-ATHIRD*l)*secder0[1] + ATHIRD*l*secder1[1]);
          gammasec[2] = 6.0*((1.0-ATHIRD*l)*secder0[2] + ATHIRD*l*secder1[2]);

          ntau2 = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];
          if ( ntau2 < EPSD )  continue;
          intau = 1.0/sqrt(ntau2);
          ntau2 = 1.0/ntau2;
          tau[0] *= intau;
          tau[1] *= intau;
          tau[2] *= intau;

          ps = gammasec[0]*tau[0] + gammasec[1]*tau[1] + gammasec[2]*tau[2];
          gammasec[0] = gammasec[0]*ntau2 - ps*ntau2*tau[0];
          gammasec[1] = gammasec[1]*ntau2 - ps*ntau2*tau[1];
          gammasec[2] = gammasec[2]*ntau2 - ps*ntau2*tau[2];
          kappa = MG_MAX(kappa,gammasec[0]*gammasec[0] + gammasec[1]*gammasec[1] + gammasec[2]*gammasec[2] );
        }
        kappa = sqrt(kappa);
        if ( kappa < EPSD )
          lm = MAXLEN;
        else
          lm = sqrt(8.0*hausd / kappa);

        if ( MG_EDG(p0->tag) && !(p0->tag & MG_NOM) && !MG_SIN(p0->tag) )
          met->m[ip0] = MG_MAX(info.hmin,MG_MIN(met->m[ip0],lm));
        if ( MG_EDG(p1->tag) && !(p1->tag & MG_NOM) && !MG_SIN(p1->tag) )
          met->m[ip1] = MG_MAX(info.hmin,MG_MIN(met->m[ip1],lm));
      }
    }
  }
  return(1);
}

/** Enforce mesh gradation by truncating size map */
int gradsiz_iso(pMesh mesh,pSol met) {
  pTetra    pt;
  pPoint    p0,p1;
  double    l,hn;
  int       ip0,ip1,it,maxit,nu,nup,k;
  char      i,j,ia,i0,i1;

  if ( abs(info.imprim) > 5 || info.ddebug )
    fprintf(stdout,"  ** Grading mesh\n");

  mesh->base = 0;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  it = nup = 0;
  maxit = 100;
  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;

      for (i=0; i<4; i++) {
        for (j=0; j<3; j++) {
          ia  = iarf[i][j];
          i0  = iare[ia][0];
          i1  = iare[ia][1];
          ip0 = pt->v[i0];
          ip1 = pt->v[i1];
          p0  = &mesh->point[ip0];
          p1  = &mesh->point[ip1];
          if ( p0->flag < mesh->base-1 && p1->flag < mesh->base-1 )  continue;

          l = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1])\
            + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
          l = sqrt(l);

          if ( met->m[ip0] < met->m[ip1] ) {
            if ( met->m[ip0] < EPSD )  continue;
            hn = met->m[ip0] + info.hgrad*l;
            if ( met->m[ip1] > hn ) {
              met->m[ip1] = hn;
              p1->flag = mesh->base;
              nu++;
            }
          }
          else {
            if ( met->m[ip1] < EPSD )  continue;
            hn = met->m[ip1] + info.hgrad*l;
            if ( met->m[ip0] > hn ) {
              met->m[ip0] = hn;
              p0->flag = mesh->base;
              nu++;
            }
          }
        }
      }
    }
    nup += nu;
  }
  while( ++it < maxit && nu > 0 );

  if ( abs(info.imprim) > 4 )
    fprintf(stdout,"     gradation: %7d updated, %d iter.\n",nup,it);
  return(1);
}
