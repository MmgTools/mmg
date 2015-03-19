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

/**
 * \file mmgs/roots.c
 * \brief Functions to compute roots of degree 2 and 3 polynomial.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmgs.h"
#define EPSRO     1.e-3

extern Info  info;

/* Returns the 2 complex roots of a degree 2 polynomial with real coefficients a[2]T^2 + ... + a[0]
   By convention, the real roots are stored first (same thing for multiple roots) :
   return value = number of roots, counted with multiplicity */
int rootDeg2(double complex a[3], double complex r[2]){
    double complex Delta,delta,r1,r2;

    if( cabs(a[2])<EPSD ) {
        if( cabs(a[1])<EPSD ) {
            r[0] = r[1] = 0.0;
            return(0);
        }

        else{
            r[0] = r[1] = -a[0]/a[1];
            return(1);
        }
    }

    Delta = a[1]*a[1] - 4.0*a[2]*a[0];
    /* ONE square root of Delta */
    delta = cpow(Delta,0.5);
    r1 = 0.5/a[2]*(-a[1] - delta);
    r2 = 0.5/a[2]*(-a[1] + delta);

    if( fabs(cimag(r1)) < EPSRO ){
        r[0] = creal(r1);
        r[1] = ( fabs(cimag(r2)) < EPSRO ) ? creal(r2) : r2;
    }
    else {
        r[0] = ( fabs(cimag(r2)) < EPSRO ) ? creal(r2) : r2;
        r[1] = r1;
    }

    return(2);
}

/* Returns the 3 complex roots of a degree 3 polynomial with real coefficients a[3]T^3 + ... + a[0]
   By convention, the real roots are stored first (same thing for multiple roots) :
   return value = number of roots, counted with multiplicity */
int rootDeg3(double a[4],double complex r[3]) {
    double p,q,Delta,u,v,pi,b[4];
    double complex i,j,t,tbar;

    pi = 3.14159;
    i = _Complex_I;
    j = cos(2.0*ATHIRD*pi)+I*sin(2.0*ATHIRD*pi);

    /* Case when the polynomial is actually second order */
    if( fabs(a[3]) < EPSD ) {
        /* Case it is first order : return 0.0 when root does not exist */
        if( fabs(a[2]) < EPSD ) {
            if( fabs(a[1]) < EPSD ) {
                r[0] = r[1] = r[2] = 0.0;
                return(0);
            }
            else{
                r[0] = r[1] = r[2] = -a[0];
                return(1);
            }
        }
        else{
            Delta = a[1]*a[1]-4.0*a[2]*a[0];
            if( Delta > 0.0 ) {
                r[0] = 0.5/a[2]*(-a[1] - sqrt(Delta));
                r[1] = r[2] = 0.5/a[2]*(-a[1] + sqrt(Delta));
            }
            else if( Delta == 0.0 ) {
                r[0] = r[1] = -0.5*a[1]/a[2];
            }
            else{
                r[0] = 0.5/a[2]*(-a[1] - i*sqrt(-Delta));
                r[1] = r[2] = 0.5/a[2]*(-a[1] + i*sqrt(-Delta));
            }
            return(2);
        }
    }

    /* Normalize coefficients */
    b[3] = 1.0;
    b[2] = a[2]/a[3];
    b[1] = a[1]/a[3];
    b[0] = a[0]/a[3];

    p = b[1] - ATHIRD*b[2]*b[2];
    q = b[0] - ATHIRD*b[2]*b[1] + 2.0/27.0*b[2]*b[2]*b[2];
    Delta = 4.0/27.0*p*p*p+q*q;

    if( Delta>0.0 ) {
        /* Polynomial T^2 +qT -p^3/27 admits two real roots u and v */
        u = 0.5*(-q - sqrt(Delta));
        v = 0.5*(-q + sqrt(Delta));

        u = u < 0.0 ? - pow(fabs(u),ATHIRD) :pow(fabs(u),ATHIRD);
        v = v < 0.0 ? - pow(fabs(v),ATHIRD) :pow(fabs(v),ATHIRD);

        r[0] = u+v;
        r[1] = j*j*u + j*v;
        r[2] = j*u + j*j*v;

    }
    else if (Delta<0.0){
        /* Polynomial T^2 +qT -p^3/27 admits two complex conjuguate roots u and v */
        t = -0.5*q - 0.5*i*sqrt(-Delta);
        t = cpow(t,ATHIRD);
        tbar = conj(t);

        /* Theoretically speaking, the 3 roots are real... But to make sure... */
        r[0] = creal(t+tbar);
        r[1] = creal(j*t+j*j*tbar);
        r[2] = creal(j*j*t+j*tbar);
    }
    else{
        /* Polynomial T^2 +qT -p^3/27 admits one double real root */
        u = -0.5*q;
        u = pow(u,ATHIRD);

        r[0] = -u;
        r[1] = -u;
        r[2] = 2.0*u;

    }

    r[0] -= ATHIRD*b[2];
    r[1] -= ATHIRD*b[2];
    r[2] -= ATHIRD*b[2];

    return(3);
}

