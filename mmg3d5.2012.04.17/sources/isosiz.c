#include "mmg3d.h"

extern Info   info;

#define MAXLEN   1.0e9

/* Compute length of edge [ip1,ip2] according to the size prescription */
inline double lenedg_iso(pMesh mesh, int ip1, int ip2){
  pPoint   p1,p2;
  double   h1,h2,l,r,lnr,powr,len;
  
  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  h1 = p1->h;
  h2 = p2->h;

  l = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]) \
	  + (p2->c[2]-p1->c[2])*(p2->c[2]-p1->c[2]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;

  if ( fabs(r) < EPS )  return(l / h1);
  
  lnr   = r;
  powr  = r*r;
  lnr  -= 0.5*powr;
  powr *= r;
  lnr  += ATHIRD*powr;
  powr *= r;
  lnr  -= 0.25*powr;
  
  if ( fabs(lnr) < EPSD )  return(l / h1);
  
  len = l*lnr / (r*h1);
  if ( len < 0.0 )  return(info.hmin);
  
  return(len);
}

/* Define isotropic size map at all boundary vertices of the mesh,
     associated with geometric approx, and prescribe hmax at the internal vertices 
	 Field h of Point is used, to store the prescribed size (not inverse, squared,...) */
int defsiz_iso(pMesh mesh,pSol met) {
  pTetra        pt;
  pxTetra       pxt;
  pPoint        p0,p1;
  int           k,ip0,ip1,l,ref;
  double        sqeps,sqhmin,sqhmax,v[3],b0[3],b1[3],b0p0[3],b1b0[3],p1b1[3];
  double        secder0[3],secder1[3],kappa,tau[3],gammasec[3],ntau2,intau,ps,lm;
  char          i,j,ia,ised,i0,i1,tag;

  if ( abs(info.imprim) > 5 || info.ddebug )
    fprintf(stdout,"  ** Defining map\n");

  if ( info.hmax < 0.0 )  info.hmax = 1000 * info.hmin;  

  sqeps = info.hausd*info.hausd;
  sqhmax= info.hmax*info.hmax;
  sqhmin= info.hmin*info.hmin;
  
  /* init constant size^2 */
  for (k=1; k<=mesh->np; k++) {
    mesh->point[k].h = sqhmax;
  }
    
  /* Travel all boundary faces to update size prescription for points of the boundary */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
	  if ( !MG_EOK(pt) )  continue;
	  else if ( !pt->xt ) continue;
	  pxt = &mesh->xtetra[pt->xt];  

	  for (i=0; i<4; i++) {
	    if ( !(pxt->ftag[i] & MG_BDY) )  continue;
  	  else if ( !norface(mesh,k,i,v) )  continue;		
		
	    for (j=0; j<3; j++) {
	      ia = iarf[i][j];
		    i0 = iare[ia][0];
		    i1 = iare[ia][1];

		    hGet(&mesh->htab,pt->v[i0],pt->v[i1],&ref,&tag);
		    ised = MG_EDG(tag);

		    ip0 = pt->v[i0];
		    ip1 = pt->v[i1];
		    p0  = &mesh->point[ip0];
		    p1  = &mesh->point[ip1];

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
		      lm = 8.0*info.hausd/kappa; 
		    p0->h = MG_MAX(sqhmin,MG_MIN(p0->h,lm));
		    p1->h = MG_MAX(sqhmin,MG_MIN(p1->h,lm));
	    }	
	  }
  }
    
  /* take input metric into account */
  for (k=1; k<=mesh->np; k++) {
	  p0 = &mesh->point[k];
	  if ( !MG_VOK(p0) )  continue;
	  p0->h = sqrt(p0->h);

    // For now, this test puts everything down to hmin... ???
    /*if ( met->size == 1 && met->m && k <= mesh->npi )
	    p0->h = MG_MAX(info.hmin,MG_MIN(p0->h,met->m[k]));*/
  }

  return(1);
}

/* Enforces mesh gradation by truncating size map */
int gradsiz_iso(pMesh mesh){
  pTetra    pt;
  pPoint    p0,p1;
  double    l,hn;
  int       it,nu,nutot,k;
  char      i,j,ia,i0,i1;

  if ( abs(info.imprim) > 5 || info.ddebug )
    fprintf(stdout,"  ** Grading mesh\n"); 

  mesh->base = 0;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  it = nutot = 0;
  do {
    mesh->base++;
	  nu = 0;
    for (k=1; k<=mesh->ne; k++) {
	    pt = &mesh->tetra[k];
	    if ( !MG_EOK(pt) )  continue;
	  
	    for (i=0; i<4; i++) {
		    for (j=0; j<3; j++) {
		      ia = iarf[i][j];
          i0 = iare[ia][0];
		      i1 = iare[ia][1];
		      p0 = &mesh->point[pt->v[i0]];
		      p1 = &mesh->point[pt->v[i1]];
		  		if ( p0->flag < mesh->base-1 && p1->flag < mesh->base-1 )  continue;

		      l = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1])\
		           + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
		      l = sqrt(l);
		  
		      if ( p0->h < p1->h ) {
            hn = p0->h + info.hgrad*l;
		        if ( p1->h > hn ) {
			        p1->h = hn;
			        p1->flag = mesh->base;
			        nu++;
			      }
		      }
		      else {
			      hn = p1->h + info.hgrad*l;
		        if ( p0->h > hn ) {
              p0->h = hn;
			        p0->flag = mesh->base;
			        nu++;
			      }
		      }
		    }
	    }
	  }
	  nutot += nu;
  }
  while( nu > 0 && ++it< 100 );

  if ( abs(info.imprim) > 4 )  fprintf(stdout,"     gradation: %7d updated, %d iter.\n",nutot,it);

  return(1);
}
