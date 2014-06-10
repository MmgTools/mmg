#include <stdio.h>
#include <string.h>
#include <math.h>

/* seeking 1.e-05 accuracy */
#define  EPSD           1.e-15
#define  EPSD2          1.e-200
#define  EPS6           5.e-06
#define  EPS            1.e-06
#define  EPSX2          2.e-06
#define  MAXTOU         50

/* check if numbers are equal */ 
#define egal(x,y)   ( \
  (  ((x) == 0.0) ? (fabs(y) < EPS) : \
   ( ((y) == 0.0) ? (fabs(x) < EPS) : \
     (fabs((x)-(y)) / (fabs(x) + fabs(y)) < EPSX2) )  ) )

static double Id[3][3] = { 
  {1.0, 0.0, 0.0},
  {0.0, 1.0, 0.0},
  {0.0, 0.0, 1.0} };


/* find root(s) of polynomial:  P(x)= x^3+bx^2+cx+d 
   return 1: 3 roots,  2: 2 roots,  3: 1 root */
static int newton3(double p[4],double x[3]) {
  double     b,c,d,da,db,dc,epsd;
  double     delta,fx,dfx,dxx;
  double     fdx0,fdx1,dx0,dx1,x1,x2;
  int        it,n;

  /* coeffs polynomial, a=1 */
  b = p[2];
  c = p[1];
  d = p[0];
  n = 1;

  /* 1st derivative of f */
  da = 3.0;
  db = 2.0*b;

  /* solve 2nd order eqn */
  delta = db*db - 4.0*da*c;
  epsd  = db*db*EPSD2;

  /* inflexion (f'(x)=0, x=-b/2a) */
  x1 = -db / 6.0;

  n = 1;
  if ( delta > epsd ) {
    delta = sqrt(delta);
    dx0   = (-db + delta) / 6.0;
    dx1   = (-db - delta) / 6.0;
    /* Horner */
    fdx0 = d + dx0*(c+dx0*(b+dx0));
    fdx1 = d + dx1*(c+dx1*(b+dx1));

    if ( fabs(fdx0) < EPSD ) {
      /* dx0: double root, compute single root */
      n = 2;
      x[0] = dx0;
      x[1] = dx0;
      x[2] = -b - 2.0*dx0;
      /* check if P(x) = 0 */
      fx = d + x[2]*(c+x[2]*(b+x[2]));
      if ( fabs(fx) > EPSD ) {
#ifdef DDEBUG
        fprintf(stderr,"  ## ERR 9100, newton3: fx= %E\n",fx);
#endif
        return(0);
      }
      return(n);
    }
    else if ( fabs(fdx1) < EPSD ) {
      /* dx1: double root, compute single root */
      n = 2;
      x[0] = dx1;
      x[1] = dx1;
      x[2] = -b - 2.0*dx1;
      /* check if P(x) = 0 */
      fx = d + x[2]*(c+x[2]*(b+x[2]));
      if ( fabs(fx) > EPSD ) {
#ifdef DDEBUG
        fprintf(stderr,"  ## ERR 9100, newton3: fx= %E\n",fx);
#endif
        return(0);
      }
      return(n);
    }
  }

  else if ( fabs(delta) < epsd ) {
    /* triple root */
    n = 3;
    x[0] = x1;
    x[1] = x1;
    x[2] = x1;
    /* check if P(x) = 0 */
    fx = d + x[0]*(c+x[0]*(b+x[0]));
    if ( fabs(fx) > EPSD ) {
#ifdef DDEBUG
      fprintf(stderr,"  ## ERR 9100, newton3: fx= %E\n",fx);
#endif
      return(0);
    }
    return(n);
  }
  
  else {
#ifdef DDEBUG
    fprintf(stderr,"  ## ERR 9101, newton3: no real roots\n");
#endif
    return(0);
  }

  /* Newton method: find one root (middle)
     starting point: P"(x)=0 */
  x1  = -b / 3.0;
  dfx =  c + b*x1;
  fx  = d + x1*(c -2.0*x1*x1);
  it  = 0;
  do {
    x2 = x1 - fx / dfx;
    fx = d + x2*(c+x2*(b+x2));
    if ( fabs(fx) < EPSD ) {
      x[0] = x2;
      break;
    }
    dfx = c + x2*(db + da*x2);

    /* check for break-off condition */
    dxx = fabs((x2-x1) / x2);
    if ( dxx < 1.0e-10 ) {
      x[0] = x2;
      if ( fabs(fx) > EPSD ) {
        fprintf(stderr,"  ## ERR 9102, newton3, no root found (fx %E).\n",fx);
        return(0);
      }
      break;
    }
    else
      x1 = x2;
  }
  while ( ++it < MAXTOU );

  if ( it == MAXTOU ) {
    x[0] = x1;
    fx   = d + x1*(c+(x1*(b+x1)));
    if ( fabs(fx) > EPSD ) {
      fprintf(stderr,"  ## ERR 9102, newton3, no root found (fx %E).\n",fx);
      return(0);
    }
  }

  /* solve 2nd order equation
     P(x) = (x-sol(1))* (x^2+bb*x+cc)  */
  db    = b + x[0];
  dc    = c + x[0]*db;
  delta = db*db - 4.0*dc;
  
  if ( delta <= 0.0 ) {
    fprintf(stderr,"  ## ERR 9103, newton3, det = 0.\n");
    return(0);
  }

  delta = sqrt(delta);
  x[1] = 0.5 * (-db+delta);
  x[2] = 0.5 * (-db-delta);
  
#ifdef DDEBUG
    /* check for root accuracy */
    fx = d + x[1]*(c+x[1]*(b+x[1]));
    if ( fabs(fx) > EPSD ) {
      fprintf(stderr,"  ## ERR 9104, newton3: fx= %E  x= %E\n",fx,x[1]);
      return(0);
    }
    fx = d + x[2]*(c+x[2]*(b+x[2]));
    if ( fabs(fx) > EPSD ) {
      fprintf(stderr,"  ## ERR 9104, newton3: fx= %E  x= %E\n",fx,x[2]);
      return(0);
    }
#endif

  return(n);
}


/* find eigenvalues and vectors of a 3x3 symmetric definite
 * positive matrix 
 * return order of eigenvalues (1,2,3) or 0 if failed */
int eigenv(int symmat,double *mat,double lambda[3],double v[3][3]) {
  double    a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double    aa,bb,cc,dd,ee,ii,vx1[3],vx2[3],vx3[3],dd1,dd2,dd3;
  double    maxd,maxm,valm,p[4],w1[3],w2[3],w3[3];
  int       k,n;

  /* default */
  memcpy(v,Id,9*sizeof(double));
  if ( symmat ) {
    lambda[0] = (double)mat[0];
    lambda[1] = (double)mat[3];
    lambda[2] = (double)mat[5];

    maxm = fabs(mat[0]);
    for (k=1; k<6; k++) {
      valm = fabs(mat[k]);
      if ( valm > maxm )  maxm = valm;
    }
    /* single float accuracy */
    if ( maxm < EPS6 )  return(1);

    /* normalize matrix */
    dd  = 1.0 / maxm;
    a11 = mat[0] * dd;
    a12 = mat[1] * dd;
    a13 = mat[2] * dd;
    a22 = mat[3] * dd;
    a23 = mat[4] * dd;
    a33 = mat[5] * dd;

    /* diagonal matrix */
    maxd = fabs(a12);
    valm = fabs(a13);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a23);
    if ( valm > maxd )  maxd = valm;
    if ( maxd < EPSD )  return(1);

    a21  = a12;
    a31  = a13;
    a32  = a23;
    
    /* build characteristic polynomial
       P(X) = X^3 - trace X^2 + (somme des mineurs)X - det = 0 */
    aa = a11*a22;
    bb = a23*a32;
    cc = a12*a21;
    dd = a13*a31;
    p[0] =  a11*bb + a33*(cc-aa) + a22*dd -2.0*a12*a13*a23;
    p[1] =  a11*(a22 + a33) + a22*a33 - bb - cc - dd;
    p[2] = -a11 - a22 - a33;
    p[3] =  1.0;
  }
  else {
    lambda[0] = (double)mat[0];
    lambda[1] = (double)mat[4];
    lambda[2] = (double)mat[8];

    maxm = fabs(mat[0]);
    for (k=1; k<9; k++) {
      valm = fabs(mat[k]);
      if ( valm > maxm )  maxm = valm;
    }
    if ( maxm < EPS6 )  return(1);

    /* normalize matrix */
    dd  = 1.0 / maxm;
    a11 = mat[0] * dd;
    a12 = mat[1] * dd;
    a13 = mat[2] * dd;
    a21 = mat[3] * dd;
    a22 = mat[4] * dd;
    a23 = mat[5] * dd;
    a31 = mat[6] * dd;
    a32 = mat[7] * dd;
    a33 = mat[8] * dd;

    /* diagonal matrix */
    maxd = fabs(a12);
    valm = fabs(a13);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a23);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a21);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a31);
    if ( valm > maxd )  maxd = valm;
     valm = fabs(a32);
    if ( valm > maxd )  maxd = valm;
    if ( maxd < EPSD )  return(1);

    /* build characteristic polynomial
       P(X) = X^3 - trace X^2 + (somme des mineurs)X - det = 0 */
    aa = a22*a33 - a23*a32;
    bb = a23*a31 - a21*a33;
    cc = a21*a32 - a31*a22;
    ee = a11*a33 - a13*a31;
    ii = a11*a22 - a12*a21;
    
    p[0] =  -a11*aa - a12*bb - a13*cc;
    p[1] =  aa + ee + ii;
    p[2] = -a11 - a22 - a33;
    p[3] =  1.0;
  }

  /* solve polynomial (find roots using newton) */
  n = newton3(p,lambda);
  if ( n <= 0 )  return(0);

  /* compute eigenvectors:
     an eigenvalue belong to orthogonal of Im(A-lambda*Id) */
  v[0][0] = 1.0; v[0][1] = v[0][2] = 0.0;
  v[1][1] = 1.0; v[1][0] = v[1][2] = 0.0;
  v[2][2] = 1.0; v[2][0] = v[2][1] = 0.0;

  w1[1] = a12;  w1[2] = a13;
  w2[0] = a21;  w2[2] = a23;
  w3[0] = a31;  w3[1] = a32;

  if ( n == 1 ) {
    /* vk = crsprd(wi,wj) */
    for (k=0; k<3; k++) {
      w1[0] = a11 - lambda[k];
      w2[1] = a22 - lambda[k];
      w3[2] = a33 - lambda[k];

      /* cross product vectors in (Im(A-lambda(i) Id) ortho */
      vx1[0] = w1[1]*w3[2] - w1[2]*w3[1];
      vx1[1] = w1[2]*w3[0] - w1[0]*w3[2];
      vx1[2] = w1[0]*w3[1] - w1[1]*w3[0];
      dd1    = vx1[0]*vx1[0] + vx1[1]*vx1[1] + vx1[2]*vx1[2];

      vx2[0] = w1[1]*w2[2] - w1[2]*w2[1];
      vx2[1] = w1[2]*w2[0] - w1[0]*w2[2];
      vx2[2] = w1[0]*w2[1] - w1[1]*w2[0];
      dd2    = vx2[0]*vx2[0] + vx2[1]*vx2[1] + vx2[2]*vx2[2];

      vx3[0] = w2[1]*w3[2] - w2[2]*w3[1];
      vx3[1] = w2[2]*w3[0] - w2[0]*w3[2];
      vx3[2] = w2[0]*w3[1] - w2[1]*w3[0];
      dd3    = vx3[0]*vx3[0] + vx3[1]*vx3[1] + vx3[2]*vx3[2];

      /* find vector of max norm */
      if ( dd1 > dd2 ) {
        if ( dd1 > dd3 ) {
          dd1 = 1.0 / sqrt(dd1);
          v[k][0] = vx1[0] * dd1;
          v[k][1] = vx1[1] * dd1;
          v[k][2] = vx1[2] * dd1;
        }
	      else {
          dd3 = 1.0 / sqrt(dd3);
          v[k][0] = vx3[0] * dd3;
          v[k][1] = vx3[1] * dd3;
          v[k][2] = vx3[2] * dd3;
	      }
      }
      else {
        if ( dd2 > dd3 ) { 
          dd2 = 1.0 / sqrt(dd2);
          v[k][0] = vx2[0] * dd2;
          v[k][1] = vx2[1] * dd2;
          v[k][2] = vx2[2] * dd2;
        }
        else {
          dd3 = 1.0 / sqrt(dd3);
          v[k][0] = vx3[0] * dd3;
          v[k][1] = vx3[1] * dd3;
          v[k][2] = vx3[2] * dd3;
        }  
      }
    }
  }

  /* (vp1,vp2) double,  vp3 simple root */
  else if ( n == 2 ) {
    w1[0] = a11 - lambda[2];
    w2[1] = a22 - lambda[2];
    w3[2] = a33 - lambda[2];

    /* cross product */
    vx1[0] = w1[1]*w3[2] - w1[2]*w3[1];
    vx1[1] = w1[2]*w3[0] - w1[0]*w3[2];
    vx1[2] = w1[0]*w3[1] - w1[1]*w3[0];
    dd1 = vx1[0]*vx1[0] + vx1[1]*vx1[1] + vx1[2]*vx1[2];
 
    vx2[0] = w1[1]*w2[2] - w1[2]*w2[1];
    vx2[1] = w1[2]*w2[0] - w1[0]*w2[2];
    vx2[2] = w1[0]*w2[1] - w1[1]*w2[0];
    dd2 = vx2[0]*vx2[0] + vx2[1]*vx2[1] + vx2[2]*vx2[2];

    vx3[0] = w2[1]*w3[2] - w2[2]*w3[1];
    vx3[1] = w2[2]*w3[0] - w2[0]*w3[2];
    vx3[2] = w2[0]*w3[1] - w2[1]*w3[0];
    dd3 = vx3[0]*vx3[0] + vx3[1]*vx3[1] + vx3[2]*vx3[2];

    /* find vector of max norm */
    if ( dd1 > dd2 ) {
      if ( dd1 > dd3 ) {
        dd1 = 1.0 / sqrt(dd1);
        v[2][0] = vx1[0] * dd1;
        v[2][1] = vx1[1] * dd1;
        v[2][2] = vx1[2] * dd1;
      }
      else {
        dd3 = 1.0 / sqrt(dd3);
        v[2][0] = vx3[0] * dd3;
        v[2][1] = vx3[1] * dd3;
        v[2][2] = vx3[2] * dd3;
      }
    }
    else {
      if ( dd2 > dd3 ) {
        dd2 = 1.0 / sqrt(dd2);
        v[2][0] = vx2[0] * dd2;
        v[2][1] = vx2[1] * dd2;
        v[2][2] = vx2[2] * dd2;
      }
      else {
        dd3 = 1.0 / sqrt(dd3);
        v[2][0] = vx3[0] * dd3;
        v[2][1] = vx3[1] * dd3;
        v[2][2] = vx3[2] * dd3;
      }
    }

    /* compute v1 and v2 in Im(A-vp3*Id) */
    dd1 = w1[0]*w1[0] + w1[1]*w1[1] + w1[2]*w1[2];
    dd2 = w2[0]*w2[0] + w2[1]*w2[1] + w2[2]*w2[2];
    if ( dd1 > dd2 ) {
      dd1 = 1.0 / sqrt(dd1);
      v[0][0] = w1[0]*dd1;
      v[0][1] = w1[1]*dd1;
      v[0][2] = w1[2]*dd1;
    }
    else {
      dd2 = 1.0 / sqrt(dd2);
      v[0][0] = w2[0]*dd2;
      v[0][1] = w2[1]*dd2;
      v[0][2] = w2[2]*dd2;
    }

    /* 3rd vector orthogonal */
    v[1][0] = v[2][1]*v[0][2] - v[2][2]*v[0][1];
    v[1][1] = v[2][2]*v[0][0] - v[2][0]*v[0][2];
    v[1][2] = v[2][0]*v[0][1] - v[2][1]*v[0][0];
    dd1 = v[1][0]*v[1][0] + v[1][1]*v[1][1] + v[1][2]*v[1][2];
    dd1 = 1.0 / sqrt(dd1);
    v[1][0] *= dd1;
    v[1][1] *= dd1;
    v[1][2] *= dd1;
  }

  lambda[0] *= maxm;
  lambda[1] *= maxm;
  lambda[2] *= maxm;

  /* check accuracy */
  /*-------------------------------------------------------------------
  if ( ddebug && symmat ) {
    double  err,tmpx,tmpy,tmpz;
    float   m[6];
    int     i,j;

    k = 0;
    for (i=0; i<3; i++)
      for (j=i; j<3; j++)
        m[k++] = lambda[0]*v[i][0]*v[j][0]
               + lambda[1]*v[i][1]*v[j][1]
               + lambda[2]*v[i][2]*v[j][2];
    err = fabs(mat[0]-m[0]);
    for (i=1; i<6; i++)
      if ( fabs(m[i]-mat[i]) > err )  err = fabs(m[i]-mat[i]);

    if ( err > 1.e03*maxm ) {
      printf("\nProbleme eigenv3: err= %f\n",err*maxm);
      printf("mat depart :\n");
      printf("%13.6f  %13.6f  %13.6f\n",mat[0],mat[1],mat[2]);
      printf("%13.6f  %13.6f  %13.6f\n",mat[1],mat[3],mat[4]);
      printf("%13.6f  %13.6f  %13.6f\n",mat[2],mat[4],mat[5]);
      printf("mat finale :\n");
      printf("%13.6f  %13.6f  %13.6f\n",m[0],m[1],m[2]);
      printf("%13.6f  %13.6f  %13.6f\n",m[1],m[3],m[4]);
      printf("%13.6f  %13.6f  %13.6f\n",m[2],m[4],m[5]);
      printf("lambda : %f %f %f\n",lambda[0],lambda[1],lambda[2]);
      printf(" ordre %d\n",n);
      printf("\nOrtho:\n");
      printf("v1.v2 = %.14f\n",
             v[0][0]*v[1][0]+v[0][1]*v[1][1]+ v[0][2]*v[1][2]);
      printf("v1.v3 = %.14f\n",
             v[0][0]*v[2][0]+v[0][1]*v[2][1]+ v[0][2]*v[2][2]);
      printf("v2.v3 = %.14f\n",
             v[1][0]*v[2][0]+v[1][1]*v[2][1]+ v[1][2]*v[2][2]);
      
      printf("Consistency\n");
      for (i=0; i<3; i++) {
        tmpx = v[0][i]*m[0] + v[1][i]*m[1]
             + v[2][i]*m[2] - lambda[i]*v[0][i];
        tmpy = v[0][i]*m[1] + v[1][i]*m[3]
             + v[2][i]*m[4] - lambda[i]*v[1][i];
        tmpz = v[0][i]*m[2] + v[1][i]*m[4]
             + v[2][i]*m[5] - lambda[i]*v[2][i];
        printf(" Av %d - lambda %d *v %d = %f %f %f\n",
        i,i,i,tmpx,tmpy,tmpz);
        
        printf("w1 %f %f %f\n",w1[0],w1[1],w1[2]);
        printf("w2 %f %f %f\n",w2[0],w2[1],w2[2]);
        printf("w3 %f %f %f\n",w3[0],w3[1],w3[2]);
      }
      exit(1);
    }
  }
  -------------------------------------------------------------------*/

  return(n);
}


/* eigen value + vector extraction */
int eigen2(double *m,double *lambda,double vp[2][2]) {
  double   a,b,c,dd1,dd2,disc;

  /* init */
  lambda[0] = lambda[1] = 0.0;
  vp[0][0] = 1.0;  vp[0][1] = 0.0;
  vp[1][0] = 0.0;  vp[1][1] = 1.0;

	a  = m[0];  b = m[1];  c = m[2];

  /* eigenvalues of jacobian */
  if ( fabs(a) > EPSD2 ) {
    disc = sqrt((a-c)*(a-c) + 4*b*b) / 2.;

    lambda[0] = (a+c) / 2. + disc;
    lambda[1] = (a+c) / 2. - disc;
    if ( fabs(lambda[0]) < EPSD2 || fabs(lambda[1]) < EPSD2 )  return(1);

    vp[0][0] = b;  vp[0][1] = lambda[0]-a;
    vp[1][0] = b;  vp[1][1] = lambda[1]-a;
  
    /* normalize */
    dd1 = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);  
    dd2 = sqrt(vp[1][0]*vp[1][0] + vp[1][1]*vp[1][1]);
    if ( dd1 > EPSD2 ) {
      vp[0][0] /= dd1;  vp[0][1] /= dd1;
    }
    if ( dd2 > EPSD2 ) {
      vp[1][0] /= dd2;  vp[1][1] /= dd2;
    }
  }

  return(1);
}
