/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/
#include "mmg2d.h"

extern int ddebug;
/*simultaneous reduction*/
int simred(double *m1,double *m2,double *m) {
  double  lambda[2],hh[2],det,pp[2][2];
  double  maxd1,maxd2,ex,ey,m1i[3],n[4],pi[4];
 
  /* check diag matrices */
  if ( fabs(m1[1]) < EPSD && fabs(m2[1]) < EPSD ) {
    m[0] = M_MAX(m1[0],m2[0]);
    m[2] = M_MAX(m1[2],m2[2]);
    m[1] = 0.0;
    return(1);
  }
  if ( !MMG2_invmat(m1,m1i) )  return(0); 
  if(ddebug) printf("inv %e %e %e\n",m1[0]*m1i[0]+m1[1]*m1i[1],m1[0]*m1i[1]+m1[1]*m1i[2],
                          m1[1]*m1i[1]+m1[2]*m1i[2]);

  _MMG5_eigensym(m1,lambda,pp);
  if(ddebug) printf(" m1 lambda : %e %e -- %e %e\n",lambda[0],lambda[1],1./sqrt(lambda[0]),
                                1./sqrt(lambda[1]));

  _MMG5_eigensym(m2,lambda,pp);
  if(ddebug) printf(" m2 lambda : %e %e -- %e %e\n",lambda[0],lambda[1],1./sqrt(lambda[0]),
                                1./sqrt(lambda[1]));


  /* n = (m1)^-1*m2 : stocke en ligne*/
  n[0] = m1i[0]*m2[0] + m1i[1]*m2[1];
  n[1] = m1i[0]*m2[1] + m1i[1]*m2[2];
  n[2] = m1i[1]*m2[0] + m1i[2]*m2[1];
  n[3] = m1i[1]*m2[1] + m1i[2]*m2[2];
  
  _MMG5_eigensym(n,lambda,pp);
  if(ddebug) printf("lambda : %e %e -- %e %e\n",lambda[0],lambda[1],0.5*((n[0]+n[3])+
                                sqrt((n[0]+n[3])*(n[0]+n[3])-4*(n[0]*n[3]-n[1]*n[2])))
                                ,0.5*((n[0]+n[3])-
                                sqrt((n[0]+n[3])*(n[0]+n[3])-4*(n[0]*n[3]-n[1]*n[2]))));
  
  if(ddebug) printf("l0 %e\n",lambda[0]*lambda[0]-lambda[0]*(n[0]+n[3])+n[0]*n[3]-n[1]*n[2]);
  if(ddebug) printf("l1 %e\n",lambda[1]*lambda[1]-lambda[1]*(n[0]+n[3])+n[0]*n[3]-n[1]*n[2]);
  if ( fabs(lambda[0]-lambda[1]) < EPSD ) {
    m[0] = m[2] = lambda[0];
    m[1] = 0.0;
    return(1);
  }
  else {
    /* matrix of passage */
    if(ddebug) printf("aaa : %e == %e -- %e == %e \n",n[0]*pp[0][0]+n[1]*pp[0][1],lambda[0]*pp[0][0]
                      ,n[2]*pp[0][0]+n[3]*pp[0][1],lambda[0]*pp[0][1]  ) ; 

    if(ddebug) printf("bbb : %e == %e -- %e == %e \n",n[0]*pp[1][0]+n[1]*pp[1][1],lambda[1]*pp[1][0]
                      ,n[2]*pp[1][0]+n[3]*pp[1][1],lambda[1]*pp[1][1]  ) ; 

    det = pp[0][0]*pp[1][1]-pp[1][0]*pp[0][1];
    if(fabs(det) < EPSD) return(0);

    det = 1./det;
    pi[0] = det*pp[1][1];
    pi[1] = -det*pp[0][1];
    pi[2] = -det*pp[1][0];
    pi[3] = det*pp[0][0];

    /*eigenvalues*/
    ex = pp[0][0];
    ey = pp[0][1];
    maxd1 = ex*(m1[0]*ex+m1[1]*ey) 
          + ey*(m1[1]*ex+m1[2]*ey);
    maxd2 = ex*(m2[0]*ex+m2[1]*ey) 
          + ey*(m2[1]*ex+m2[2]*ey);
    hh[0] = M_MAX(maxd1,maxd2);
    ex = pp[1][0];
    ey = pp[1][1];
    maxd1 = ex*(m1[0]*ex+m1[1]*ey) 
          + ey*(m1[1]*ex+m1[2]*ey);
    maxd2 = ex*(m2[0]*ex+m2[1]*ey) 
          + ey*(m2[1]*ex+m2[2]*ey);
    hh[1] = M_MAX(maxd1,maxd2);
    
    if(ddebug) printf("----------------- hh : %e %e -- %e %e\n",hh[0],hh[1],1./sqrt(hh[0]),1./sqrt(hh[1]));
  
    /* compose matrix tP^-1*lambda*P^-1 */
    m[0] = pi[0]*hh[0]*pi[0] + pi[2]*hh[1]*pi[2];
    m[1] = pi[0]*hh[0]*pi[1] + pi[2]*hh[1]*pi[3];
    m[2] = pi[1]*hh[0]*pi[1] + pi[3]*hh[1]*pi[3];

    if ( ddebug ) {
      _MMG5_eigensym(m,lambda,pp);
      if ( lambda[0] < -EPSD || lambda[1] < -EPSD ) {
        fprintf(stderr,"  ## simred, not a metric !\n");
        fprintf(stderr,"  %.6f %.6f %.6f\n",
                m[0],m[1],m[2]);
        fprintf(stderr,"  Lambda %f %f \n",lambda[0],lambda[1]);
        return(0);
      }
    }

    return(1);
  }   
  
  
  
}
