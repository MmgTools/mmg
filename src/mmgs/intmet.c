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
 * \file mmgs/intmet.c
 * \brief Metric interpolations.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

extern char ddb;

/* Compute the interpolated (2 x 2) metric from metrics m and n, at parameter s :
   mr = (1-s)*m +s*n, both metrics being expressed in the simultaneous reduction basis:
   linear interpolation of sizes */
static int intmet22(double *m,double *n,double *mr,double s) {
    double  det,imn[4],dd,sqDelta,trimn,lambda[2],vp0[2],vp1[2],dm[2],dn[2],vnorm,d0,d1,ip[4];

    /* Compute imn = M^{-1}N */
    det = m[0]*m[2] - m[1]*m[1];
    if ( fabs(det) < EPS*EPS ) {
        printf("BEWARE : function intmet : null metric det : %E \n",det);
        return(0);
    }
    det = 1.0 / det;

    imn[0] = det * ( m[2]*n[0] - m[1]*n[1]);
    imn[1] = det * ( m[2]*n[1] - m[1]*n[2]);
    imn[2] = det * (-m[1]*n[0] + m[0]*n[1]);
    imn[3] = det * (-m[1]*n[1] + m[0]*n[2]);
    dd = imn[0] - imn[3];
    sqDelta = sqrt(fabs(dd*dd + 4.0*imn[1]*imn[2]));
    trimn = imn[0] + imn[3];

    lambda[0] = 0.5 * (trimn - sqDelta);
    if ( lambda[0] < 0.0 ) {
        printf("Les valeurs propres : %f \n",lambda[0]);
        return(0);
    }

    /* First case : matrices m and n are homothetic = n = lambda0*m */
    if ( sqDelta < EPS ) {
        dd  = (1.0-s)*sqrt(lambda[0]) + s;
        dd *= dd;
        if ( dd < EPSD ) {
            mr[0] = m[0];
            mr[1] = m[1];
            mr[2] = m[2];
            return(1);
        }
        dd = lambda[0] / dd;
        mr[0] = dd * m[0];
        mr[1] = dd * m[1];
        mr[2] = dd * m[2];
        return(1);
    }

    /* Second case : both eigenvalues of imn are distinct ; theory says qf associated to m and n
       are diagonalizable in basis (vp0, vp1) - the coreduction basis */
    else {
        lambda[1] = 0.5 * (trimn + sqDelta);
        assert(lambda[1] >= 0.0);

        vp0[0] = imn[1];
        vp0[1] = (lambda[0] - imn[0]);
        vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);
        if ( vnorm < EPS ) {
            vp0[0] = (lambda[0] - imn[3]);
            vp0[1] = imn[2];
            vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);
        }
        vnorm   = 1.0 / vnorm;
        vp0[0] *= vnorm;
        vp0[1] *= vnorm;

        vp1[0] = imn[1];
        vp1[1] = (lambda[1] - imn[0]);
        vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);
        if ( vnorm < EPS ) {
            vp1[0] = (lambda[1] - imn[3]);
            vp1[1] = imn[2];
            vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);
        }
        vnorm   = 1.0 / vnorm;
        vp1[0] *= vnorm;
        vp1[1] *= vnorm;

        /* Compute diagonal values in simultaneous reduction basis */
        dm[0] = m[0]*vp0[0]*vp0[0] + 2.0*m[1]*vp0[0]*vp0[1] + m[2]*vp0[1]*vp0[1];
        dm[1] = m[0]*vp1[0]*vp1[0] + 2.0*m[1]*vp1[0]*vp1[1] + m[2]*vp1[1]*vp1[1];
        dn[0] = n[0]*vp0[0]*vp0[0] + 2.0*n[1]*vp0[0]*vp0[1] + n[2]*vp0[1]*vp0[1];
        dn[1] = n[0]*vp1[0]*vp1[0] + 2.0*n[1]*vp1[0]*vp1[1] + n[2]*vp1[1]*vp1[1];

        /* Diagonal values of the interpolated metric */
        dd  = (1.0-s)*sqrt(dn[0]) + s*sqrt(dm[0]);
        dd *= dd;
        if ( dd < EPSD ) {
            d0 = s < 0.5 ? dm[0] : dn[0];
        }
        else {
            d0 = dm[0]*dn[0] / dd;
        }
        dd = (1.0-s)*sqrt(dn[1]) + s*sqrt(dm[1]);
        dd *= dd;
        if ( dd < EPSD ) {
            d1 = s < 0.5 ? dm[1] : dn[1];
        }
        else{
            d1 = dm[1]*dn[1] / dd;
        }

        /* Intersected metric = tP^-1 diag(d0,d1)P^-1, P = (vp0, vp1) stored in columns */
        det = vp0[0]*vp1[1] - vp0[1]*vp1[0];
        if ( fabs(det) < EPS )  return(0);
        det = 1.0 / det;

        ip[0] =  vp1[1]*det;
        ip[1] = -vp1[0]*det;
        ip[2] = -vp0[1]*det;
        ip[3] =  vp0[0]*det;

        mr[0] = d0*ip[0]*ip[0] + d1*ip[2]*ip[2];
        mr[1] = d0*ip[0]*ip[1] + d1*ip[2]*ip[3];
        mr[2] = d0*ip[1]*ip[1] + d1*ip[3]*ip[3];
    }

    return(1);
}

/* Metric interpolation between two points p1 and p2 such that edge 0 = (p1p2) is ridge.
   v is a direction vector, aimed at pointing towards direction of n1 at interpolated point */
int intridmet(pMesh mesh,pSol met,int k,char i,double s,double v[3],double mr[6]) {
    pTria     pt;
    pGeom     go1,go2;
    pPoint    p1,p2;
    double   *m1,*m2,*n11,*n12,*n21,*n22,ps11,ps12,dd,hn1,hn2;
    int       ip1,ip2;
    char      i1,i2;

    pt  = &mesh->tria[k];
    i1  = inxt[i];
    i2  = iprv[i];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];
    p1  = &mesh->point[ip1];
    p2  = &mesh->point[ip2];
    m1  = &met->m[6*(ip1)+1];
    m2  = &met->m[6*(ip2)+1];

    /* Case when both endpoints are singular */
    if ( MS_SIN(p1->tag) && MS_SIN(p2->tag) ) {
        dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
        dd *= dd;
        if ( dd < EPSD ) {
            mr[0] = s < 0.5 ? m1[0] : m2[0];
            mr[1] = s < 0.5 ? m1[0] : m2[0];
            mr[2] = s < 0.5 ? m1[0] : m2[0];
        }
        else {
            mr[0] = m1[0]*m2[0] / dd;
            mr[1] = m1[0]*m2[0] / dd;
            mr[2] = m1[0]*m2[0] / dd;
        }
    }
    /* vertex p1 is singular, p2 is regular */
    else if ( MS_SIN(p1->tag) && (!MS_SIN(p2->tag)) ) {
        go2 = &mesh->geom[p2->ig];
        n21 = &go2->n1[0];
        n22 = &go2->n2[0];

        /* Interpolation of the eigenvalue associated to tangent vector */
        dd = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
        dd *= dd;
        if ( dd < EPSD ) {
            mr[0] = s < 0.5 ? m1[0] : m2[0];
        }
        else {
            mr[0] = m1[0]*m2[0] / dd;
        }
        /* Interpolation of the two other eigenvalues */
        dd  = (1-s)*sqrt(m2[1]) + s*sqrt(m1[0]);
        dd *= dd;
        if ( dd < EPSD ) {
            hn1 = s < 0.5 ? m1[0] : m2[1];
        }
        else {
            hn1 = m1[0]*m2[1] / dd;
        }
        dd = (1-s)*sqrt(m2[2]) + s*sqrt(m1[0]);
        dd *= dd;
        if ( dd < EPSD ) {
            hn2 = s < 0.5 ? m1[0] : m2[2];
        }
        else {
            hn2 = m1[0]*m2[2] / dd;
        }

        /* Decision of the ordering of hn1 and hn2 */
        ps11 = n21[0]*v[0] + n21[1]*v[1] + n21[2]*v[2];
        ps12 = n22[0]*v[0] + n22[1]*v[1] + n22[2]*v[2];
        if ( fabs(ps11) > fabs(ps12) ) {
            mr[1] = hn1;
            mr[2] = hn2;
        }
        else {
            mr[1] = hn2;
            mr[2] = hn1;
        }
    }
    /* vertex p2 is singular, p1 is regular */
    else if ( MS_SIN(p2->tag) && (!MS_SIN(p1->tag)) ) {
        go1 = &mesh->geom[p2->ig];
        n11 = &go1->n1[0];
        n12 = &go1->n2[0];

        /* Interpolation of the eigenvalue associated to tangent vector */
        dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
        dd *= dd;
        if ( dd < EPSD ) {
            mr[0] = s < 0.5 ? m1[0] : m2[0];
        }
        else {
            mr[0] = m1[0]*m2[0] / dd;
        }
        /* Interpolation of the two other eigenvalues */
        dd = (1-s)*sqrt(m2[0]) + s*sqrt(m1[1]);
        dd *= dd;
        if ( dd < EPSD ) {
            hn1 = s < 0.5 ? m1[1] : m2[0];
        }
        else {
            hn1 = m1[1]*m2[0] / dd;
        }
        dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[2]);
        dd *= dd;
        if ( dd < EPSD ) {
            hn2 = s < 0.5 ? m1[2] : m2[0];
        }
        else {
            hn2 = m1[2]*m2[0] / dd;
        }

        /* Decision of the ordering of hn1 and hn2 */
        ps11 = n11[0]*v[0] + n11[1]*v[1] + n11[2]*v[2];
        ps12 = n12[0]*v[0] + n12[1]*v[1] + n12[2]*v[2];
        if ( fabs(ps11) > fabs(ps12) ) {
            mr[1] = hn1;
            mr[2] = hn2;
        }
        else {
            mr[1] = hn2;
            mr[2] = hn1;
        }
    }
    /* p1,p2 : nonsingular vertices */
    else {
        go1 = &mesh->geom[p1->ig];
        go2 = &mesh->geom[p2->ig];

        /* Interpolation of the eigenvalue associated to tangent vector */
        dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
        dd *= dd;
        if ( dd < EPSD ) {
            mr[0] = s < 0.5 ? m1[0] : m2[0];
        }
        else {
            mr[0] = m1[0]*m2[0] / dd;
        }

        /* Pairing of normal vectors at p1 and p2 */
        n11 = &go1->n1[0];
        n12 = &go1->n2[0];
        n21 = &go2->n1[0];
        n22 = &go2->n2[0];
        ps11 = n11[0]*n21[0] + n11[1]*n21[1] + n11[2]*n21[2];
        ps12 = n11[0]*n22[0] + n11[1]*n22[1] + n11[2]*n22[2];
        if ( fabs(ps11) > fabs(ps12) ) {   //n11 and n21 go together
            dd  = (1-s)*sqrt(m2[1]) + s*sqrt(m1[1]);
            dd *= dd;
            if ( dd < EPSD ) {
                hn1 = s < 0.5 ? m1[1] : m2[1];
            }
            else {
                hn1 = m1[1]*m2[1] / dd;
            }
            dd = (1-s)*sqrt(m2[2]) + s*sqrt(m1[2]);
            dd *= dd;
            if ( dd < EPSD ) {
                hn2 = s < 0.5 ? m1[2] : m2[2];
            }
            else {
                hn2 = m1[2]*m2[2] / dd;
            }
        }
        else {
            dd  = (1-s)*sqrt(m2[2]) + s*sqrt(m1[1]);
            dd *= dd;
            if ( dd < EPSD ) {
                hn1 = s < 0.5 ? m1[1] : m2[2];
            }
            else {
                hn1 = m1[1]*m2[2] / dd;
            }
            dd  = (1-s)*sqrt(m2[1]) + s*sqrt(m1[2]);
            dd *= dd;
            if ( dd < EPSD ) {
                hn2 = s < 0.5 ? m1[2] : m2[1];
            }
            else {
                hn2 = m1[2]*m2[1] / dd;
            }
        }

        /* Now, hn1 is the eigenvalue associated to the direction at interpolated point,
           closest to n11 (hn2 -> n12) ; one may need a different orientation, and put eigenvalue of
           direction closest to v (= interpolated normal) first */
        ps11 = n11[0]*v[0] + n11[1]*v[1] + n11[2]*v[2];
        ps12 = n12[0]*v[0] + n12[1]*v[1] + n12[2]*v[2];
        if ( fabs(ps11) > fabs(ps12) ) {
            mr[1] = hn1;
            mr[2] = hn2;
        }
        else {
            mr[1] = hn2;
            mr[2] = hn1;
        }
    }
    mr[3] = 0.0;
    mr[4] = 0.0;
    mr[5] = 0.0;

    return(1);
}

/* Metric interpolation between points p1 and p2, in tria it at parameter 0 <= s0 <= 1 from p1
   result is stored in mr. edge p1p2 must not be a ridge */
int intregmet(pMesh mesh,pSol met,int k,char i,double s,double mr[6]) {
    pTria     pt;
    pPoint    p1,p2;
    Bezier    b;
    double    b1[3],b2[3],bn[3],c[3],nt[3],cold[3],n[3],nold[3],mold[6],m1[6],m2[6];
    double   *n1,*n2,step,u,r[3][3],mt1[3],mt2[3],dd;
    int       ip1,ip2,nstep,l;
    char      i1,i2;

    /* Number of steps for parallel transport */
    nstep = 4;
    pt  = &mesh->tria[k];
    nortri(mesh,pt,nt);
    i1  = inxt[i];
    i2  = iprv[i];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];
    p1  = &mesh->point[ip1];
    p2  = &mesh->point[ip2];

    if ( !bezierCP(mesh,k,&b) )  return(0);

    n1 = &b.n[i1][0];
    n2 = &b.n[i2][0];
    memcpy(bn,&b.n[i+3][0],3*sizeof(double));
    memcpy(b1,&b.b[2*i+3][0],3*sizeof(double));
    memcpy(b2,&b.b[2*i+4][0],3*sizeof(double));

    /* Parallel transport of metric at p1 to point p(s) */
    step = s / nstep;
    cold[0] = p1->c[0];
    cold[1] = p1->c[1];
    cold[2] = p1->c[2];

    nold[0] = n1[0];
    nold[1] = n1[1];
    nold[2] = n1[2];

    if ( MS_SIN(p1->tag) ) {
        memcpy(m1,&met->m[6*(ip1)+1],6*sizeof(double));
    }
    else {
        if ( MS_GEO & p1->tag ) {
            if ( !buildridmetnor(mesh,met,pt->v[i1],nt,m1) )  return(0);
        }
        else {
            memcpy(m1,&met->m[6*(ip1)+1],6*sizeof(double));
        }
        memcpy(mold,m1,6*sizeof(double));

        /* Go from point (l-1)step, to point l step */
        for (l=1; l<=nstep; l++) {
            u    = l*step;
            c[0] = p1->c[0]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[0]\
                + 3.0*u*u*(1.0-u)*b2[0] + u*u*u*p2->c[0];
            c[1] = p1->c[1]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[1]\
                + 3.0*u*u*(1.0-u)*b2[1] + u*u*u*p2->c[1];
            c[2] = p1->c[2]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[2]\
                + 3.0*u*u*(1.0-u)*b2[2] + u*u*u*p2->c[2];

            n[0] = (1.0-u)*(1.0-u)*n1[0] + 2.0*u*(1.0-u)*bn[0] + u*u*n2[0];
            n[1] = (1.0-u)*(1.0-u)*n1[1] + 2.0*u*(1.0-u)*bn[1] + u*u*n2[1];
            n[2] = (1.0-u)*(1.0-u)*n1[2] + 2.0*u*(1.0-u)*bn[2] + u*u*n2[2];
            dd   = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
            if ( dd < EPSD )  return(0);
            dd = 1.0 / sqrt(dd);
            n[0] *= dd;
            n[1] *= dd;
            n[2] *= dd;

            if ( !paratmet(cold,nold,mold,c,n,m1) )  return(0);

            memcpy(cold,c,3*sizeof(double));
            memcpy(nold,n,3*sizeof(double));
            memcpy(mold,m1,6*sizeof(double));
        }
    }

    /* Parallel transport of metric at p2 to point p(s) */
    step = (1.0-s) / nstep;
    cold[0] = p2->c[0];
    cold[1] = p2->c[1];
    cold[2] = p2->c[2];

    nold[0] = n2[0];
    nold[1] = n2[1];
    nold[2] = n2[2];

    if ( MS_SIN(p2->tag) ) {
        memcpy(m2,&met->m[6*(ip2)+1],6*sizeof(double));

        /* In this pathological case, n is empty */
        if ( MS_SIN(p1->tag) )
            memcpy(n,n2,3*sizeof(double));
    }
    else {
        if ( p2->tag & MS_GEO ) {
            if ( !buildridmetnor(mesh,met,pt->v[i2],nt,m2))  return(0);
        }
        else {
            memcpy(m2,&met->m[6*(ip2)+1],6*sizeof(double));
        }
        memcpy(mold,m2,6*sizeof(double));

        /* Go from point (l-1)step, to point l step */
        for (l=1; l<=nstep; l++) {
            u    = 1.0 - l*step;
            c[0] = p1->c[0]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[0]\
                + 3.0*u*u*(1.0-u)*b2[0] + u*u*u*p2->c[0];
            c[1] = p1->c[1]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[1]\
                + 3.0*u*u*(1.0-u)*b2[1] + u*u*u*p2->c[1];
            c[2] = p1->c[2]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[2]\
                + 3.0*u*u*(1.0-u)*b2[2] + u*u*u*p2->c[2];

            n[0] = (1.0-u)*(1.0-u)*n1[0] + 2.0*u*(1.0-u)*bn[0] + u*u*n2[0];
            n[1] = (1.0-u)*(1.0-u)*n1[1] + 2.0*u*(1.0-u)*bn[1] + u*u*n2[1];
            n[2] = (1.0-u)*(1.0-u)*n1[2] + 2.0*u*(1.0-u)*bn[2] + u*u*n2[2];
            dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
            if ( dd < EPSD )  return(0);
            dd = 1.0 / sqrt(dd);
            n[0] *= dd;
            n[1] *= dd;
            n[2] *= dd;

            if ( !paratmet(cold,nold,mold,c,n,m2) )  return(0);

            memcpy(cold,c,3*sizeof(double));
            memcpy(nold,n,3*sizeof(double));
            memcpy(mold,m2,6*sizeof(double));
        }
    }
    /* At this point, c is point p(s), n is the normal at p(s), m1 and m2 are the 3*3
       transported metric tensors from p1 and p2 to p(s) */

    /* Rotate both matrices to the tangent plane */
    if ( !rotmatrix(n,r) )  return(0);
    rmtr(r,m1,mold);
    mt1[0] = mold[0];
    mt1[1] = mold[1];
    mt1[2] = mold[3];

    rmtr(r,m2,mold);
    mt2[0] = mold[0];
    mt2[1] = mold[1];
    mt2[2] = mold[3];

    /* Interpolate both metrics expressed in the same tangent plane : cold is reused */
    if ( !intmet22(mt1,mt2,cold,s) ) {
        printf("Impossible interpolation between points : %d %d\n",pt->v[i1],pt->v[i2]);
        printf("m1 : %E %E %E %E %E %E \n",m1[0],m1[1],m1[2],m1[3],m1[4],m1[5]);
        printf("m2 : %E %E %E %E %E %E \n",m2[0],m2[1],m2[2],m2[3],m2[4],m2[5]);
        printf("mt1 : %E %E %E et det %E \n",mt1[0],mt1[1],mt1[2],mt1[0]*mt1[2]-mt1[1]*mt1[1]);
        printf("mt2 : %E %E %E et det : %E \n",mt2[0],mt2[1],mt2[2],mt2[0]*mt2[2]-mt2[1]*mt2[1]);
        hashTria(mesh);
        saveMesh(mesh);
        exit(0);
        return(0);
    }

    /* End rotating back tangent metric into canonical basis of R^3 : mr = {^t}R*cold*R
       mt1 and mt2 serve for nothing ; let them be the lines of cold * R  */
    mt1[0] = cold[0]*r[0][0] + cold[1]*r[1][0];  mt1[1] = cold[0]*r[0][1] + cold[1]*r[1][1];   mt1[2] = cold[0]*r[0][2] + cold[1]*r[1][2] ;
    mt2[0] = cold[1]*r[0][0] + cold[2]*r[1][0];  mt2[1] = cold[1]*r[0][1] + cold[2]*r[1][1];   mt2[2] = cold[1]*r[0][2] + cold[2]*r[1][2] ;

    mr[0] = r[0][0] * mt1[0] + r[1][0] * mt2[0];
    mr[1] = r[0][0] * mt1[1] + r[1][0] * mt2[1];
    mr[2] = r[0][0] * mt1[2] + r[1][0] * mt2[2];
    mr[3] = r[0][1] * mt1[1] + r[1][1] * mt2[1];
    mr[4] = r[0][1] * mt1[2] + r[1][1] * mt2[2];
    mr[5] = r[0][2] * mt1[2] + r[1][2] * mt2[2];

    return(1);
}

/* linear interpolation of sizemap along edge i of tria k */
void intmet_iso(pMesh mesh,pSol met,int k,char i,int ip,double s) {
    pTria  pt;
    int    ip1,ip2;
    char   i1,i2;

    pt  = &mesh->tria[k];
    i1  = inxt[i];
    i2  = iprv[i];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];
    met->m[ip] = s * (met->m[ip1] + met->m[ip2]);
}

void intmet_ani(pMesh mesh,pSol met,int k,char i,int ip,double s) {
    pTria    pt;
    pPoint   ppt;
    pGeom    go;
    double  *m;

    pt = &mesh->tria[k];
    m  = &met->m[6*(ip)+1];
    if ( pt->tag[i] & MS_GEO ) {
        ppt = &mesh->point[ip];
        assert(ppt->ig);
        go = &mesh->geom[ppt->ig];
        intridmet(mesh,met,k,i,s,go->n1,m);
    }
    else {
        intregmet(mesh,met,k,i,s,m);
    }
}

/* Interpolation of size features for a 3*3 metric tensor field, between points np and nq,
   at parameter s from np */
int intmet33(pMesh mesh,pSol met,int np,int nq,int ip,double s) {
    int     order;
    double  *m,*n,*mr,lambda[3],vp[3][3],mu[3],is[6],isnis[6],mt[9],P[9],dd;
    char    i;

    m  = &met->m[6*np+1];
    n  = &met->m[6*nq+1];
    mr = &met->m[6*ip+1];

    /* Compute inverse of square root of matrix M : is = P*diag(1/sqrt(lambda))*{^t}P */
    order = eigenv(1,m,lambda,vp);
    if ( !order ) return(0);

    for (i=0; i<3; i++) {
        if ( lambda[i] < EPSD ) return(0);
        lambda[i] = sqrt(lambda[i]);
        lambda[i] = 1.0 / lambda[i];
    }

    is[0] = lambda[0]*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0] + lambda[1]*vp[2][0]*vp[2][0];
    is[1] = lambda[0]*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1] + lambda[1]*vp[2][0]*vp[2][1];
    is[2] = lambda[0]*vp[0][0]*vp[0][2] + lambda[1]*vp[1][0]*vp[1][2] + lambda[1]*vp[2][0]*vp[2][2];
    is[3] = lambda[0]*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1] + lambda[1]*vp[2][1]*vp[2][1];
    is[4] = lambda[0]*vp[0][1]*vp[0][2] + lambda[1]*vp[1][1]*vp[1][2] + lambda[1]*vp[2][1]*vp[2][2];
    is[5] = lambda[0]*vp[0][2]*vp[0][2] + lambda[1]*vp[1][2]*vp[1][2] + lambda[1]*vp[2][2]*vp[2][2];

    mt[0] = n[0]*is[0] + n[1]*is[1] + n[2]*is[2];
    mt[1] = n[0]*is[1] + n[1]*is[3] + n[2]*is[4];
    mt[2] = n[0]*is[2] + n[1]*is[4] + n[2]*is[5];
    mt[3] = n[1]*is[0] + n[3]*is[1] + n[4]*is[2];
    mt[4] = n[1]*is[1] + n[3]*is[3] + n[4]*is[4];
    mt[5] = n[1]*is[2] + n[3]*is[4] + n[4]*is[5];
    mt[6] = n[2]*is[0] + n[4]*is[1] + n[5]*is[2];
    mt[7] = n[2]*is[1] + n[4]*is[3] + n[5]*is[4];
    mt[8] = n[2]*is[2] + n[4]*is[4] + n[5]*is[5];

    isnis[0] = is[0]*mt[0] + is[1]*mt[3] + is[2]*mt[6];
    isnis[1] = is[0]*mt[1] + is[1]*mt[4] + is[2]*mt[7];
    isnis[2] = is[0]*mt[2] + is[1]*mt[5] + is[2]*mt[8];
    isnis[3] = is[1]*mt[1] + is[3]*mt[4] + is[4]*mt[7];
    isnis[4] = is[1]*mt[2] + is[3]*mt[5] + is[4]*mt[8];
    isnis[5] = is[2]*mt[2] + is[4]*mt[5] + is[5]*mt[8];

    order = eigenv(1,isnis,lambda,vp);
    if ( !order ) return(0);

    /* P = is * (vp) */
    P[0] = is[0]*vp[0][0] + is[1]*vp[0][1] + is[2]*vp[0][2];
    P[1] = is[0]*vp[1][0] + is[1]*vp[1][1] + is[2]*vp[1][2];
    P[2] = is[0]*vp[2][0] + is[1]*vp[2][1] + is[2]*vp[2][2];
    P[3] = is[1]*vp[0][0] + is[3]*vp[0][1] + is[4]*vp[0][2];
    P[4] = is[1]*vp[1][0] + is[3]*vp[1][1] + is[4]*vp[1][2];
    P[5] = is[1]*vp[2][0] + is[3]*vp[2][1] + is[4]*vp[2][2];
    P[6] = is[2]*vp[0][0] + is[4]*vp[0][1] + is[5]*vp[0][2];
    P[7] = is[2]*vp[1][0] + is[4]*vp[1][1] + is[5]*vp[1][2];
    P[8] = is[2]*vp[2][0] + is[4]*vp[2][1] + is[5]*vp[2][2];

    /* At this point, theory states that ^tPMP = I, {^t}PNP=\Lambda */
    /* Linear interpolation between sizes */
    for(i=0; i<3; i++) {
        if ( lambda[i] < 0.0 ) return(0);
        dd = s*sqrt(lambda[i]) + (1.0-s);
        dd = dd*dd;
        if ( dd < EPSD )  return(0);
        mu[i] = lambda[i]/dd;
    }

    if ( !invmatg(P,mt) )  return(0);

    /* Resulting matrix = ^tP^{-1} diag(mu) P^{-1} */
    mr[0] = mu[0]*mt[0]*mt[0] + mu[1]*mt[3]*mt[3] + mu[2]*mt[6]*mt[6];
    mr[1] = mu[0]*mt[0]*mt[1] + mu[1]*mt[3]*mt[4] + mu[2]*mt[6]*mt[7];
    mr[2] = mu[0]*mt[0]*mt[2] + mu[1]*mt[3]*mt[5] + mu[2]*mt[6]*mt[8];
    mr[3] = mu[0]*mt[1]*mt[1] + mu[1]*mt[4]*mt[4] + mu[2]*mt[7]*mt[7];
    mr[4] = mu[0]*mt[1]*mt[2] + mu[1]*mt[4]*mt[5] + mu[2]*mt[7]*mt[8];
    mr[5] = mu[0]*mt[2]*mt[2] + mu[1]*mt[5]*mt[5] + mu[2]*mt[8]*mt[8];

    return(1);
}
