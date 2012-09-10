#include "mmg3d.h"


/* naive (increasing) sorting algorithm, for very small tabs ; permutation is stored in perm */
inline void nsort(int n,double *val,char *perm){
  int   i,j,aux;

  for (i=0; i<n; i++)  perm[i] = i; 

  for (i=0; i<n; i++) {
    for (j=i+1; j<n; j++) {
      if ( val[perm[i]] > val[perm[j]] ) {
        aux = perm[i];
        perm[i] = perm[j];
        perm[j] = aux;
      } 
    }
  }
}

/* Compute rotation matrix that sends vector n to the third vector of canonical basis */
inline void rotmatrix(double n[3],double r[3][3]) {
  double   aa,bb,ab,ll,l,cosalpha,sinalpha;

  aa = n[0]*n[0];
  bb = n[1]*n[1];
  ab = n[0]*n[1];	  	  
  ll = aa+bb; 
  cosalpha = n[2];
  sinalpha = sqrt(1.0- MG_MIN(1.0,cosalpha*cosalpha));
  
  /* No rotation needed in this case */
  if ( ll < EPS*EPS ) {
    if ( n[2] > 0.0 ) {
      r[0][0] = 1.0 ; r[0][1] = 0.0 ; r[0][2] = 0.0;  
      r[1][0] = 0.0 ; r[1][1] = 1.0 ; r[1][2] = 0.0;  
      r[2][0] = 0.0 ; r[2][1] = 0.0 ; r[2][2] = 1.0;
    }
    else {
      r[0][0] = -1.0 ; r[0][1] = 0.0 ; r[0][2] = 0.0;  
      r[1][0] = 0.0 ; r[1][1] = 1.0 ; r[1][2] = 0.0;  
      r[2][0] = 0.0 ; r[2][1] = 0.0 ; r[2][2] = -1.0;
    }
  } 
  else {
    l = sqrt(ll);
    r[0][0] = (aa*cosalpha + bb)/ll;
    r[0][1] = ab*(cosalpha-1)/ll;
    r[0][2] = -n[0]*sinalpha/l;
    r[1][0] = r[0][1];
    r[1][1] = (bb*cosalpha + aa)/ll;
    r[1][2] = -n[1]*sinalpha/l;
    r[2][0] = n[0]*sinalpha/l;
    r[2][1] = n[1]*sinalpha/l;
    r[2][2] = cosalpha;
  }  
}

/* Compute 3 * 3 determinant : det(c1-c0,c2-c0,v) */
inline double det3pt1vec(double c0[3],double c1[3],double c2[3],double v[3]) {
  double m00,m10,m20,m01,m11,m21,det;
  
  m00 = c1[0] - c0[0] ; m01 = c2[0] - c0[0]; 
  m10 = c1[1] - c0[1] ; m11 = c2[1] - c0[1]; 
  m20 = c1[2] - c0[2] ; m21 = c2[2] - c0[2];
  det = v[0]*(m10*m21 - m20*m11) -v[1]*(m00*m21-m20*m01) + v[2]*(m00*m11-m10*m01);

  return(det);
} 

/* Compute 3 * 3 determinant : det(c1-c0,c2-c0,c3-c0) */
inline double det4pt(double c0[3],double c1[3],double c2[3],double c3[3]) {
  double m00,m10,m20,m01,m11,m21,m02,m12,m22,det;
  
  m00 = c1[0] - c0[0] ; m01 = c2[0] - c0[0]; m02 = c3[0] - c0[0];
  m10 = c1[1] - c0[1] ; m11 = c2[1] - c0[1]; m12 = c3[1] - c0[1];
  m20 = c1[2] - c0[2] ; m21 = c2[2] - c0[2]; m22 = c3[2] - c0[2];
  det = m02*(m10*m21 - m20*m11) -m12*(m00*m21-m20*m01) + m22*(m00*m11-m10*m01);
  
  return(det);
} 

/* Compute oriented volume of a tetrahedron */
inline double orvol(pPoint point,int *v) {
  pPoint  p0,p1,p2,p3;
  
  p0 = &point[v[0]];
  p1 = &point[v[1]];
  p2 = &point[v[2]];
  p3 = &point[v[3]];

  return(det4pt(p0->c,p1->c,p2->c,p3->c));
} 

/* Compute normal to face iface of tetra k, exterior to tetra k */
inline int norface(pMesh mesh,int k,int iface,double n[3]) {
  pTetra     pt;
  pPoint     p0,p1,p2;
  double     ux,uy,uz,vx,vy,vz,norm;
  
  pt = &mesh->tetra[k];
  p0 = &mesh->point[pt->v[idir[iface][0]]];
  p1 = &mesh->point[pt->v[idir[iface][1]]];
  p2 = &mesh->point[pt->v[idir[iface][2]]];
  
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  vx = p2->c[0] - p0->c[0];
  vy = p2->c[1] - p0->c[1];
  vz = p2->c[2] - p0->c[2];
  
  n[0] = uy*vz - uz*vy;
  n[1] = uz*vx - ux*vz;
  n[2] = ux*vy - uy*vx;  
  norm = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( norm < EPSD )  return(0);
  
  norm = 1.0/sqrt(norm);
  n[0] *= norm;
  n[1] *= norm;
  n[2] *= norm;
  return(1);
}

/* If need be, invert the travelling sense of surfacic ball so that it is travelled in 
the direct sense with respect to direction n anchored at point ip (ip = global num.): 
        return 2 = orientation reversed, 1 otherwise */
inline int directsurfball(pMesh mesh, int ip, int *list, int ilist, double n[3]){
  pTetra          pt;
  pPoint          p0,p1;
  int             j,aux,iel;
  unsigned char   iface;
 
  iel = list[0]/4;
  iface = list[0]%4;
  
  pt = &mesh->tetra[iel];
  
  p0 = &mesh->point[pt->v[iface]];
  p1 = &mesh->point[ip];
  
  if((p1->c[0] - p0->c[0])*n[0] + (p1->c[1] - p0->c[1])*n[1] + \
        (p1->c[2] - p0->c[2])*n[2] > 0.0) return(1);
  
  for( j=1 ; j<= (ilist-1)/2; j++ ){
    aux = list[j];
	  list[j] = list[ilist -j];
	  list[ilist -j] = aux;
  }
     
  return(2);
}

/* If need be, reorder the surfacic ball of point ip, so that its first element has
  edge (p,q) (nump,q = global num) as edge iprv2[ip] of face iface. 
     return 2 = orientation reversed, 1 otherwise */
inline int startedgsurfball(pMesh mesh, int nump, int numq, int *list, int ilist){
  pTetra          pt;
  int             iel,tmp,l;
  unsigned char   iface,ip,ipt;
  
  iel = list[0]/4;                     
  iface = list[0]%4;
  pt = &mesh->tetra[iel];
  
  for(ip=0;ip<4;ip++){
	  if(pt->v[ip] == nump) break;
  }
  assert(ip<4);
  
  pt = &mesh->tetra[iel];
  ipt = idirinv[iface][ip]; // index of ip in face iface
  ipt = inxt2[ipt];         // next index in this face
  ipt = idir[iface][ipt];  // index of this point in local num of tetra
  
  if(pt->v[ipt] == numq) return(1);
  
  else{
    ipt = idir[iface][iprv2[idirinv[iface][ip]]];
	  assert(pt->v[ipt] == numq);
	
	  tmp = list[0];
	  for(l=0;l<ilist-1;l++){
	    list[l] = list[l+1]; 
	  }
	  list[ilist-1] = tmp;
  }
  
  return(2);
}

/* Compute point located at parameter value step from point ip0, as well as interpolate
 of normals, tangent for a RIDGE edge */
inline int BezierRidge(pMesh mesh,int ip0, int ip1, double s, double *o, double *no1, double *no2, double *to){
  pPoint          p0,p1;
  double          ux,uy,uz,n01[3],n02[3],n11[3],n12[3],t0[3],t1[3];
  double          ps,ps2,b0[3],b1[3],bn[3],ll,il,ntemp[3],dd,alpha;
  
  p0 = &mesh->point[ip0];  // Ref point, from which step is counted
  p1 = &mesh->point[ip1];
  
  assert(p0->xp && p1->xp);
  
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];
  
  ll = ux*ux + uy*uy + uz*uz;
  if(ll < EPSD) return(0);
  
  il = 1.0/sqrt(ll);
  
  assert((MG_GEO & p0->tag) && (MG_GEO & p1->tag) );
  
  /* Pathological case : don't move in that case ! */
  if(MG_SIN(p0->tag) && MG_SIN(p1->tag))
    return(0);
  
  if(MG_SIN(p0->tag)){
    t0[0] = ux * il;
	  t0[1] = uy * il;
	  t0[2] = uz * il; 
  }
  else{
    memcpy(t0,&(mesh->xpoint[p0->xp].t[0]),3*sizeof(double));
    ps = t0[0]*ux + t0[1]*uy + t0[2]*uz;
    if(ps < 0.0){
      t0[0] *= -1.0;
      t0[1] *= -1.0;
      t0[2] *= -1.0;
    }
  }
  
  if(MG_SIN(p1->tag)){
    t1[0] = -ux * il;
	  t1[1] = -uy * il;
	  t1[2] = -uz * il; 
  }
  else{
    memcpy(t1,&(mesh->xpoint[p1->xp].t[0]),3*sizeof(double));
    ps = - ( t1[0]*ux + t1[1]*uy + t1[2]*uz );
    if(ps < 0.0){
      t1[0] *= -1.0;
      t1[1] *= -1.0;
      t1[2] *= -1.0;
    }
  }
  
  alpha = BezierGeod(p0->c,p1->c,t0,t1);
  
  b0[0] = p0->c[0] + alpha * t0[0];
  b0[1] = p0->c[1] + alpha * t0[1];
  b0[2] = p0->c[2] + alpha * t0[2];
  
  b1[0] = p1->c[0] + alpha * t1[0];
  b1[1] = p1->c[1] + alpha * t1[1];
  b1[2] = p1->c[2] + alpha * t1[2];
  
  o[0] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[0] + 3.0*s*(1.0-s)*(1.0-s)*b0[0] + \
              3.0*s*s*(1.0-s)*b1[0] + s*s*s*p1->c[0];
  
  o[1] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[1] + 3.0*s*(1.0-s)*(1.0-s)*b0[1] + \
              3.0*s*s*(1.0-s)*b1[1] + s*s*s*p1->c[1];	
  
  o[2] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[2] + 3.0*s*(1.0-s)*(1.0-s)*b0[2] + \
              3.0*s*s*(1.0-s)*b1[2] + s*s*s*p1->c[2];
  
  if(MG_SIN(p0->tag)){
    assert(!(MG_SIN(p1->tag)));
    memcpy(n11,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
	  memcpy(n12,&(mesh->xpoint[p1->xp].n2[0]),3*sizeof(double));
	
	  memcpy(n01,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
	  memcpy(n02,&(mesh->xpoint[p1->xp].n2[0]),3*sizeof(double));
	
  }	
  else if(MG_SIN(p1->tag)){
    assert(!(MG_SIN(p0->tag)));
  	memcpy(n01,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
	  memcpy(n02,&(mesh->xpoint[p0->xp].n2[0]),3*sizeof(double));
	
	  memcpy(n11,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
	  memcpy(n12,&(mesh->xpoint[p0->xp].n2[0]),3*sizeof(double));
  }
  else{
	  memcpy(n01,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
	  memcpy(n02,&(mesh->xpoint[p0->xp].n2[0]),3*sizeof(double));
	
	  memcpy(n11,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
	  memcpy(n12,&(mesh->xpoint[p1->xp].n2[0]),3*sizeof(double));
	
	  /* Switch normals of p1 for pairing */
	  ps = n01[0] * n11[0] + n01[1] * n11[1] + n01[2] * n11[2];
	  ps2 = n01[0] * n12[0] + n01[1] * n12[1] + n01[2] * n12[2];
	
    if(ps2 > ps){
	    memcpy(ntemp,n11,3*sizeof(double));
	    memcpy(n11,n12,3*sizeof(double)); 
	    memcpy(n12,ntemp,3*sizeof(double)); 
	  }
  }
  
  /* Normal n1 interpolation */
  ps = ux*(n01[0] + n11[0]) + uy*(n01[1] + n11[1]) + uz*(n01[2] + n11[2]);
  ps = 2.0*ps/ll;
  
  bn[0] = n01[0] + n11[0] -ps*ux;
  bn[1] = n01[1] + n11[1] -ps*uy;
  bn[2] = n01[2] + n11[2] -ps*uz;
  
  dd = bn[0]*bn[0] + bn[1]*bn[1] + bn[2]*bn[2];
  if(dd > EPSD){
    dd = 1.0/sqrt(dd);
	  bn[0] *= dd;
	  bn[1] *= dd;
	  bn[2] *= dd;
  }
  
  no1[0] = (1.0-s)*(1.0-s)*n01[0] + 2.0*s*(1.0-s)*bn[0] + s*s*n11[0];
  no1[1] = (1.0-s)*(1.0-s)*n01[1] + 2.0*s*(1.0-s)*bn[1] + s*s*n11[1];
  no1[2] = (1.0-s)*(1.0-s)*n01[2] + 2.0*s*(1.0-s)*bn[2] + s*s*n11[2];
  
  dd = no1[0]*no1[0] + no1[1]*no1[1] + no1[2]*no1[2];
  if(dd > EPSD){
    dd = 1.0/sqrt(dd);
	  no1[0] *= dd;
	  no1[1] *= dd;
	  no1[2] *= dd;
  }
  
  /* Normal n2 interpolation */
  
  ps = ux*(n02[0] + n12[0]) + uy*(n02[1] + n12[1]) + uz*(n02[2] + n12[2]);
  ps = 2.0*ps/ll;
  
  bn[0] = n02[0] + n12[0] -ps*ux;
  bn[1] = n02[1] + n12[1] -ps*uy;
  bn[2] = n02[2] + n12[2] -ps*uz;
  
  dd = bn[0]*bn[0] + bn[1]*bn[1] + bn[2]*bn[2];
  if(dd > EPSD){
    dd = 1.0/sqrt(dd);
	  bn[0] *= dd;
	  bn[1] *= dd;
	  bn[2] *= dd;
  }
  
  no2[0] = (1.0-s)*(1.0-s)*n02[0] + 2.0*s*(1.0-s)*bn[0] + s*s*n12[0];
  no2[1] = (1.0-s)*(1.0-s)*n02[1] + 2.0*s*(1.0-s)*bn[1] + s*s*n12[1];
  no2[2] = (1.0-s)*(1.0-s)*n02[2] + 2.0*s*(1.0-s)*bn[2] + s*s*n12[2];  
  
  dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
  if(dd > EPSD){
    dd = 1.0/sqrt(dd);
	  no2[0] *= dd;
	  no2[1] *= dd;
	  no2[2] *= dd;
  }
  
  to[0] = no1[1]*no2[2] - no1[2]*no2[1];
  to[1] = no1[2]*no2[0] - no1[0]*no2[2];
  to[2] = no1[0]*no2[1] - no1[1]*no2[0];
  
  dd = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];
  if(dd > EPSD){
    dd = 1.0/sqrt(dd);
	  to[0] *= dd;
	  to[1] *= dd;
	  to[2] *= dd;
  }
  
  return(1);
}

/* Compute point located at parameter value step from point ip0, as well as interpolate
   of normals, tangent for a REF edge */
inline int BezierRef(pMesh mesh,int ip0, int ip1, double s, double *o, double *no, double *to){
  pPoint          p0,p1;
  double          ux,uy,uz,n0[3],n1[3],t0[3],t1[3];
  double          ps,b0[3],b1[3],bn[3],ll,il,dd,alpha;
    
  p0 = &mesh->point[ip0];  // Ref point, from which step is counted
  p1 = &mesh->point[ip1];
  
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];
  
  ll = ux*ux + uy*uy + uz*uz;
  if(ll < EPSD) return(0);
  
  il = 1.0/sqrt(ll);
  
  assert((MG_REF & p0->tag) && (MG_REF & p1->tag) );
  
  /* Pathological case : don't move in that case ! */
  if(MG_SIN(p0->tag) && MG_SIN(p1->tag))
    return(0);
	
  /* Coordinates of the new point */
  if(MG_SIN(p0->tag)){
    t0[0] = ux * il;
	  t0[1] = uy * il;
	  t0[2] = uz * il; 
  }
  else{
    memcpy(t0,&(mesh->xpoint[p0->xp].t[0]),3*sizeof(double));
    ps = t0[0]*ux + t0[1]*uy + t0[2]*uz;
    if( ps < 0.0 ){
      t0[0] *= -1.0;
      t0[1] *= -1.0;
      t0[2] *= -1.0;
    }
  }
  
  if(MG_SIN(p1->tag)){
    t1[0] = -ux * il;
	  t1[1] = -uy * il;
	  t1[2] = -uz * il; 
  }
  else{
    memcpy(t1,&(mesh->xpoint[p1->xp].t[0]),3*sizeof(double));
    ps = - ( t1[0]*ux + t1[1]*uy + t1[2]*uz );
    if( ps < 0.0 ){
      t1[0] *= -1.0;
      t1[1] *= -1.0;
      t1[2] *= -1.0;
    }
  }
  
  alpha = BezierGeod(p0->c,p1->c,t0,t1);
  
  b0[0] = p0->c[0] + alpha * t0[0];
  b0[1] = p0->c[1] + alpha * t0[1];
  b0[2] = p0->c[2] + alpha * t0[2];
  
  b1[0] = p1->c[0] + alpha * t1[0];
  b1[1] = p1->c[1] + alpha * t1[1];
  b1[2] = p1->c[2] + alpha * t1[2];

  o[0] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[0] + 3.0*s*(1.0-s)*(1.0-s)*b0[0] + \
            3.0*s*s*(1.0-s)*b1[0] + s*s*s*p1->c[0];
	
  o[1] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[1] + 3.0*s*(1.0-s)*(1.0-s)*b0[1] + \
            3.0*s*s*(1.0-s)*b1[1] + s*s*s*p1->c[1];	
	
  o[2] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[2] + 3.0*s*(1.0-s)*(1.0-s)*b0[2] + \
            3.0*s*s*(1.0-s)*b1[2] + s*s*s*p1->c[2];
	
  /* Coordinates of the new tangent and normal */	
  if(MG_SIN(p0->tag)){
    assert(!(MG_SIN(p1->tag)));
    memcpy(n1,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));	
	  memcpy(n0,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));

  }	
  else if(MG_SIN(p1->tag)){
    assert(!(MG_SIN(p0->tag)));
  	memcpy(n0,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
	  memcpy(n1,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
  }
  else{
  	memcpy(n0,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));	
	  memcpy(n1,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
  }
  
  /* Normal interpolation */
  ps = ux*(n0[0] + n1[0]) + uy*(n0[1] + n1[1]) + uz*(n0[2] + n1[2]);
  ps = 2.0*ps/ll;
  
  bn[0] = n0[0] + n1[0] -ps*ux;
  bn[1] = n0[1] + n1[1] -ps*uy;
  bn[2] = n0[2] + n1[2] -ps*uz;
  
  dd = bn[0]*bn[0] + bn[1]*bn[1] + bn[2]*bn[2];
  if(dd > EPSD){
    dd = 1.0/sqrt(dd);
	  bn[0] *= dd;
	  bn[1] *= dd;
	  bn[2] *= dd;
  }
  
  no[0] = (1.0-s)*(1.0-s)*n0[0] + 2.0*s*(1.0-s)*bn[0] + s*s*n1[0];
  no[1] = (1.0-s)*(1.0-s)*n0[1] + 2.0*s*(1.0-s)*bn[1] + s*s*n1[1];
  no[2] = (1.0-s)*(1.0-s)*n0[2] + 2.0*s*(1.0-s)*bn[2] + s*s*n1[2];
  
  dd = no[0]*no[0] + no[1]*no[1] + no[2]*no[2];
  if(dd > EPSD){
    dd = 1.0/sqrt(dd);
	  no[0] *= dd;
	  no[1] *= dd;
	  no[2] *= dd;
  }
  
  /* Tangent interpolation : possibly flip (back) t1 */
  ps = t0[0]*t1[0] + t0[1]*t1[1] + t0[2]*t1[2];
  if(ps < 0.0){
    t1[0] *= -1.0;
    t1[1] *= -1.0;
    t1[2] *= -1.0;
  }  
  
  to[0] = (1.0-s)*t0[0] + s*t1[0]; 
  to[1] = (1.0-s)*t0[1] + s*t1[1]; 
  to[2] = (1.0-s)*t0[2] + s*t1[2]; 
     
  /* Projection of the tangent in the tangent plane defined by no */
  ps = to[0]*no[0] + to[1]*no[1] + to[2]*no[2];
  to[0] = to[0] -ps*no[0];
  to[1] = to[1] -ps*no[1];
  to[2] = to[2] -ps*no[2];
  
  dd = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];
  if(dd > EPSD){
    dd = 1.0/sqrt(dd);
	  to[0] *= dd;
	  to[1] *= dd;
	  to[2] *= dd;
  }
 
  return(1);
}


/* Compute point located at parameter value step from point ip0, as well as interpolate
 of normals, tangent for a regular edge ; v = ref vector (normal) for choice of normals if need be */
inline int BezierReg(pMesh mesh,int ip0, int ip1, double s, double v[3], double *o, double *no){
  pPoint p0,p1;
  double b0[3],b1[3],bn[3],t0[3],t1[3],np0[3],np1[3],alpha,ux,uy,uz,ps1,ps2,ll,il,dd,*n1,*n2;
  
  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];
  
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];
  
  ll = ux*ux + uy*uy + uz*uz;
  
  /* Pathological case : don't move in that case ! */
  if((MG_SIN(p0->tag) && MG_SIN(p1->tag))||(ll<EPSD)){
    o[0] = 0.5*( p0->c[0] + p1->c[0] ); 
    o[1] = 0.5*( p0->c[1] + p1->c[1] ); 
    o[2] = 0.5*( p0->c[2] + p1->c[2] ); 
    
    memcpy(no,v,3*sizeof(double));
    return(1);
  }
  
  il = 1.0 /sqrt(ll);
  
  /* Coordinates of the new tangent and normal */	
  if(MG_SIN(p0->tag)){
    if(p1->tag & MG_GEO){
      n1 = &mesh->xpoint[p1->xp].n1[0];
      n2 = &mesh->xpoint[p1->xp].n2[0];
      ps1 = n1[0]*v[0] + n1[1]*v[1] + n1[2]*v[2];
      ps2 = n2[0]*v[0] + n2[1]*v[1] + n2[2]*v[2];
      
      if(fabs(ps1) < fabs(ps2)){
        memcpy(np1,&mesh->xpoint[p1->xp].n2[0],3*sizeof(double)); 
      }
      else{
        memcpy(np1,&mesh->xpoint[p1->xp].n1[0],3*sizeof(double));
      }
    }
    else{
      memcpy(np1,&mesh->xpoint[p1->xp].n1[0],3*sizeof(double));
    }
    memcpy(np0,np1,3*sizeof(double));
  }	
  else if(MG_SIN(p1->tag)){
    if(p0->tag & MG_GEO){
      n1 = &mesh->xpoint[p0->xp].n1[0];
      n2 = &mesh->xpoint[p0->xp].n2[0];
      ps1 = n1[0]*v[0] + n1[1]*v[1] + n1[2]*v[2];
      ps2 = n2[0]*v[0] + n2[1]*v[1] + n2[2]*v[2];
      
      if(fabs(ps1) < fabs(ps2)){
        memcpy(np0,&mesh->xpoint[p0->xp].n2[0],3*sizeof(double)); 
      }
      else{
        memcpy(np0,&mesh->xpoint[p0->xp].n1[0],3*sizeof(double));
      }
    }
    else{
      memcpy(np0,&mesh->xpoint[p0->xp].n1[0],3*sizeof(double));
    }
    memcpy(np1,np0,3*sizeof(double));
  }
  else{
    if(p0->tag & MG_GEO){
      n1 = &mesh->xpoint[p0->xp].n1[0];
      n2 = &mesh->xpoint[p0->xp].n2[0];
      ps1 = n1[0]*v[0] + n1[1]*v[1] + n1[2]*v[2];
      ps2 = n2[0]*v[0] + n2[1]*v[1] + n2[2]*v[2];
      
      if(fabs(ps1) < fabs(ps2)){
        memcpy(np0,&mesh->xpoint[p0->xp].n2[0],3*sizeof(double)); 
      }
      else{
        memcpy(np0,&mesh->xpoint[p0->xp].n1[0],3*sizeof(double));
      }
    }
    else{
      memcpy(np0,&mesh->xpoint[p0->xp].n1[0],3*sizeof(double));
    }    
    
    if(p1->tag & MG_GEO){
      n1 = &mesh->xpoint[p1->xp].n1[0];
      n2 = &mesh->xpoint[p1->xp].n2[0];
      ps1 = n1[0]*v[0] + n1[1]*v[1] + n1[2]*v[2];
      ps2 = n2[0]*v[0] + n2[1]*v[1] + n2[2]*v[2];
      
      if(fabs(ps1) < fabs(ps2)){
        memcpy(np1,&mesh->xpoint[p1->xp].n2[0],3*sizeof(double)); 
      }
      else{
        memcpy(np1,&mesh->xpoint[p1->xp].n1[0],3*sizeof(double));
      }
    }
    else{
      memcpy(np1,&mesh->xpoint[p1->xp].n1[0],3*sizeof(double));
    }
  }
  
  /* Normal interpolation */
  ps1 = ux*(np0[0] + np1[0]) + uy*(np0[1] + np1[1]) + uz*(np0[2] + np1[2]);
  ps1 = 2.0*ps1/ll;
  
  bn[0] = np0[0] + np1[0] -ps1*ux;
  bn[1] = np0[1] + np1[1] -ps1*uy;
  bn[2] = np0[2] + np1[2] -ps1*uz;
  
  dd = bn[0]*bn[0] + bn[1]*bn[1] + bn[2]*bn[2];
  if(dd > EPSD){
    dd = 1.0/sqrt(dd);
	  bn[0] *= dd;
	  bn[1] *= dd;
	  bn[2] *= dd;
  }
  
  no[0] = (1.0-s)*(1.0-s)*np0[0] + 2.0*s*(1.0-s)*bn[0] + s*s*np1[0];
  no[1] = (1.0-s)*(1.0-s)*np0[1] + 2.0*s*(1.0-s)*bn[1] + s*s*np1[1];
  no[2] = (1.0-s)*(1.0-s)*np0[2] + 2.0*s*(1.0-s)*bn[2] + s*s*np1[2];
  
  dd = no[0]*no[0] + no[1]*no[1] + no[2]*no[2];
  if(dd > EPSD){
    dd = 1.0/sqrt(dd);
	  no[0] *= dd;
	  no[1] *= dd;
	  no[2] *= dd;
  }
  
  /* vertex position interpolation */
  if(!BezierTgt(p0->c,p1->c,np0,np1,t0,t1)){
    t0[0] = ux * il;
    t0[1] = uy * il;
    t0[2] = uz * il;
    
    t1[0] = - ux * il;
    t1[1] = - uy * il;
    t1[2] = - uz * il;
  }

  alpha = BezierGeod(p0->c,p1->c,t0,t1);
  
  b0[0] = p0->c[0] + alpha * t0[0];
  b0[1] = p0->c[1] + alpha * t0[1];
  b0[2] = p0->c[2] + alpha * t0[2];
  
  b1[0] = p1->c[0] + alpha * t1[0];
  b1[1] = p1->c[1] + alpha * t1[1];
  b1[2] = p1->c[2] + alpha * t1[2];
  
  o[0] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[0] + 3.0*s*(1.0-s)*(1.0-s)*b0[0] + \
              3.0*s*s*(1.0-s)*b1[0] + s*s*s*p1->c[0];
	
  o[1] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[1] + 3.0*s*(1.0-s)*(1.0-s)*b0[1] + \
              3.0*s*s*(1.0-s)*b1[1] + s*s*s*p1->c[1];	
 	
  o[2] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[2] + 3.0*s*(1.0-s)*(1.0-s)*b0[2] + \
              3.0*s*s*(1.0-s)*b1[2] + s*s*s*p1->c[2];
  
  return(1);
}

