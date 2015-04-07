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
 * \file mmgs/mettools.c
 * \brief Algebraic tools for metric handling.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

extern Info  info;
extern char ddb;


/* Build metric tensor at a fictitious ridge point, whose normal and tangent are provided */
inline int buildridmetfic(pMesh mesh,double t[3],double n[3],double dtan,double dv,double m[6]) {
    double u[3],r[3][3];

    u[0] = n[1]*t[2] - n[2]*t[1];
    u[1] = n[2]*t[0] - n[0]*t[2];
    u[2] = n[0]*t[1] - n[1]*t[0];

    /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(p0->m[0],dv,0)*/
    r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n[0];
    r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n[1];
    r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n[2];

    m[0] = dtan*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1];
    m[1] = dtan*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1];
    m[2] = dtan*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1];
    m[3] = dtan*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1];
    m[4] = dtan*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1];
    m[5] = dtan*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1];

    return(1);
}

/* Parallel transport of a metric tensor field, attached to point c0, with normal n0,
   to point c1, with normal n1 */
int paratmet(double c0[3],double n0[3],double m[6],double c1[3],double n1[3],double mt[6]) {
    double  r[3][3],mrot[6],mtan[3],lambda[2],vp[2][2],u[3],ps,ll;
    int     ord;

    /* Take the induced metric tensor in the tangent plane by change of basis : R * M * {^t}R*/
    if ( !rotmatrix(n0,r) )  return(0);
    rmtr(r,m,mrot);
    mtan[0] = mrot[0];
    mtan[1] = mrot[1];
    mtan[2] = mrot[3];

    /* Take eigenvectors of metric tensor in tangent plane */
    ord = eigensym(mtan,lambda,vp);

    /* Eigenvector in canonical basis = {t}R*vp[0] */
    u[0] = r[0][0]*vp[0][0] + r[1][0]*vp[0][1];
    u[1] = r[0][1]*vp[0][0] + r[1][1]*vp[0][1];
    u[2] = r[0][2]*vp[0][0] + r[1][2]*vp[0][1];

    /* Projection in the tangent plane of c1 */
    ps = u[0]*n1[0] + u[1]*n1[1] + u[2]*n1[2];
    u[0] -= ps*n1[0];
    u[1] -= ps*n1[1];
    u[2] -= ps*n1[2];
    ll = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    if ( ll < EPSD )  return(0);
    ll = 1.0 / sqrt(ll);
    u[0] *= ll;
    u[1] *= ll;
    u[2] *= ll;

    /* And the transported metric is diag(lambda[0], lambda[1], 0) in basis (u,n1^u,n1) */
    r[0][0] = u[0];
    r[1][0] = u[1];
    r[2][0] = u[2];

    r[0][1] = n1[1]*u[2] - n1[2]*u[1];
    r[1][1] = n1[2]*u[0] - n1[0]*u[2];
    r[2][1] = n1[0]*u[1] - n1[1]*u[0];

    ll = r[0][1]*r[0][1] + r[1][1]*r[1][1] + r[2][1]*r[2][1];
    if ( ll < EPSD )  return(0);
    ll = 1.0 / sqrt(ll);
    r[0][1] *= ll;
    r[1][1] *= ll;
    r[2][1] *= ll;

    r[0][2] = n1[0];
    r[1][2] = n1[1];
    r[2][2] = n1[2];

    /*mt = R * diag(lambda[0], lambda[1], 0)*{^t}R */
    mt[0] = lambda[0]*r[0][0]*r[0][0] + lambda[1]*r[0][1]*r[0][1];
    mt[1] = lambda[0]*r[0][0]*r[1][0] + lambda[1]*r[0][1]*r[1][1];
    mt[2] = lambda[0]*r[0][0]*r[2][0] + lambda[1]*r[0][1]*r[2][1];
    mt[3] = lambda[0]*r[1][0]*r[1][0] + lambda[1]*r[1][1]*r[1][1];
    mt[4] = lambda[0]*r[2][0]*r[1][0] + lambda[1]*r[2][1]*r[1][1];
    mt[5] = lambda[0]*r[2][0]*r[2][0] + lambda[1]*r[2][1]*r[2][1];

    return(1);
}

/* Build metric tensor at ridge point p0, when the 'good' normal direction is given by nt */
int buildridmetnor(pMesh mesh,pSol met,int np0,double nt[3],double mr[6]) {
    pPoint p0;
    pGeom  go;
    double ps1,ps2,*n1,*n2,*t,*m,dv,u[3],r[3][3];

    p0 = &mesh->point[np0];
    if ( !(MS_GEO & p0->tag) )  return(0);
    m = &met->m[6*(np0)+1];
    t = &p0->n[0];
    go = &mesh->geom[p0->ig];

    /* Decide between the two possible configurations */
    n1 = &go->n1[0];
    n2 = &go->n2[0];

    ps1 = nt[0]*n1[0] + nt[1]*n1[1] + nt[2]*n1[2];
    ps2 = nt[0]*n2[0] + nt[1]*n2[1] + nt[2]*n2[2];

    if ( fabs(ps2) > fabs(ps1) ) {
        n1 = &go->n2[0];
        dv = m[2];
    }
    else{
        dv = m[1];
    }

    u[0] = n1[1]*t[2] - n1[2]*t[1];
    u[1] = n1[2]*t[0] - n1[0]*t[2];
    u[2] = n1[0]*t[1] - n1[1]*t[0];

    /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(m[0],dv,0)*/
    r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n1[0];
    r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n1[1];
    r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n1[2];

    mr[0] = m[0]*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1];
    mr[1] = m[0]*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1];
    mr[2] = m[0]*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1];
    mr[3] = m[0]*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1];
    mr[4] = m[0]*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1];
    mr[5] = m[0]*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1];

    return(1);
}

/* Build metric tensor at ridge point p0, when computations with respect to p1 are to be held */
int buildridmet(pMesh mesh,pSol met,int np0,double ux,double uy,double uz,double mr[6]) {
    pPoint p0;
    pGeom  go;
    double ps1,ps2,*n1,*n2,*t,*m,dv,u[3],r[3][3];

    p0 = &mesh->point[np0];
    if ( !(MS_GEO & p0->tag) )  return(0);
    m = &met->m[6*(np0)+1];
    t = &p0->n[0];
    go = &mesh->geom[p0->ig];

    /* Decide between the two possible configurations */
    n1 = &go->n1[0];
    n2 = &go->n2[0];

    ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
    ps2 = ux*n2[0] + uy*n2[1] + uz*n2[2];

    if ( fabs(ps2)<fabs(ps1) ) {
        n1 = &go->n2[0];
        dv = m[2];
    }
    else{
        dv = m[1];
    }

    u[0] = n1[1]*t[2] - n1[2]*t[1];
    u[1] = n1[2]*t[0] - n1[0]*t[2];
    u[2] = n1[0]*t[1] - n1[1]*t[0];

    /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(m[0],dv,0)*/
    r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n1[0];
    r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n1[1];
    r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n1[2];

    mr[0] = m[0]*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1];
    mr[1] = m[0]*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1];
    mr[2] = m[0]*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1];
    mr[3] = m[0]*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1];
    mr[4] = m[0]*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1];
    mr[5] = m[0]*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1];

    return(1);
}

/* Compute length of edge [i0i1] according to the prescribed metric */
double lenedg_ani(pMesh mesh,pSol met,int np0,int np1,char isedg) {
    pPoint   p0,p1;
    double   gammaprim0[3],gammaprim1[3],t[3],*n1,*n2,ux,uy,uz,ps1,ps2,l0,l1;
    double  *m0,*m1,met0[6],met1[6];

    p0 = &mesh->point[np0];
    p1 = &mesh->point[np1];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    /* computation of the two tangent vectors to the underlying curve of [i0i1] */
    if ( MS_SIN(p0->tag) ) {
        gammaprim0[0] = ux;
        gammaprim0[1] = uy;
        gammaprim0[2] = uz;
    }
    else if ( isedg ) {
        memcpy(t,p0->n,3*sizeof(double));
        ps1 = ux*t[0] + uy*t[1] + uz*t[2];
        gammaprim0[0] = ps1*t[0];
        gammaprim0[1] = ps1*t[1];
        gammaprim0[2] = ps1*t[2];
    }
    else {
        if ( MS_GEO & p0->tag ) {
            //assert(p0->ig);
            n1 = &mesh->geom[p0->ig].n1[0];
            n2 = &mesh->geom[p0->ig].n2[0];
            ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
            ps2 = ux*n2[0] + uy*n2[1] + uz*n2[2];

            if ( fabs(ps2) < fabs(ps1) ) {
                n1  = &mesh->geom[p0->ig].n2[0];
                ps1 = ps2;
            }
        }
        else if ( MS_REF & p0->tag ) {
            n1  = &mesh->geom[p0->ig].n1[0];
            ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
        }
        else {
            n1  = &(p0->n[0]);
            ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
        }
        gammaprim0[0] = ux - ps1*n1[0];
        gammaprim0[1] = uy - ps1*n1[1];
        gammaprim0[2] = uz - ps1*n1[2];
    }

    if ( MS_SIN(p1->tag) ) {
        gammaprim1[0] = -ux;
        gammaprim1[1] = -uy;
        gammaprim1[2] = -uz;
    }
    else if ( isedg ) {
        memcpy(t,p1->n,3*sizeof(double));
        ps1 = -ux*t[0] - uy*t[1] - uz*t[2];
        gammaprim1[0] = ps1*t[0];
        gammaprim1[1] = ps1*t[1];
        gammaprim1[2] = ps1*t[2];
    }
    else {
        if ( MS_GEO & p1->tag ) {
            n1 = &mesh->geom[p1->ig].n1[0];
            n2 = &mesh->geom[p1->ig].n2[0];
            ps1 = -ux*n1[0] - uy*n1[1] - uz*n1[2];
            ps2 = -ux*n2[0] - uy*n2[1] - uz*n2[2];

            if ( fabs(ps2) < fabs(ps1) ) {
                n1  = &mesh->geom[p1->ig].n2[0];
                ps1 = ps2;
            }
        }
        else if ( MS_REF & p1->tag ) {
            n1  = &mesh->geom[p1->ig].n1[0];
            ps1 = - ux*n1[0] - uy*n1[1] - uz*n1[2];
        }
        else {
            n1  = &(p1->n[0]);
            ps1 = -ux*n1[0] - uy*n1[1] - uz*n1[2];
        }
        gammaprim1[0] = - ux - ps1*n1[0];
        gammaprim1[1] = - uy - ps1*n1[1];
        gammaprim1[2] = - uz - ps1*n1[2];
    }

    /* Set metrics */
    if ( MS_SIN(p0->tag) ) {
        m0 = &met->m[6*(np0)+1];
    }
    else if ( MS_GEO & p0->tag ) {
        if ( !buildridmet(mesh,met,np0,ux,uy,uz,met0) )  return(-1.0);
        m0 = met0;
    }
    else {
        m0 = &met->m[6*(np0)+1];
    }

    if ( MS_SIN(p1->tag) ) {
        m1 = &met->m[6*(np1)+1];
    }
    else if ( MS_GEO & p1->tag ) {
        if ( !buildridmet(mesh,met,np1,ux,uy,uz,met1) )  return(-1.0);
        m1 = met1;
    }
    else {
        m1 = &met->m[6*(np1)+1];
    }

    /* computation of the length of the two tangent vectors in their respective tangent plane */
    l0 = m0[0]*gammaprim0[0]*gammaprim0[0] + m0[3]*gammaprim0[1]*gammaprim0[1] \
        + m0[5]*gammaprim0[2]*gammaprim0[2] \
        + 2.0*m0[1]*gammaprim0[0]*gammaprim0[1]  + 2.0*m0[2]*gammaprim0[0]*gammaprim0[2] \
        + 2.0*m0[4]*gammaprim0[1]*gammaprim0[2];

    l1 = m1[0]*gammaprim1[0]*gammaprim1[0] + m1[3]*gammaprim1[1]*gammaprim1[1] \
        + m1[5]*gammaprim1[2]*gammaprim1[2] \
        +2.0*m1[1]*gammaprim1[0]*gammaprim1[1]  + 2.0*m1[2]*gammaprim1[0]*gammaprim1[2] \
        + 2.0*m1[4]*gammaprim1[1]*gammaprim1[2];

    l0 = 0.5*(sqrt(l0) + sqrt(l1));
    return(l0);
}

/* Compute anisotropic volume of element iel, with respect to metric met */
double surftri_ani(pMesh mesh,pSol met,int iel) {
    pTria     pt;
    pPoint    p[3];
    Bezier    b;
    int       np[3];
    double    surf,ux,uy,uz,dens,m[3][6],J[3][2],mJ[3][2],tJmJ[2][2];
    char      i,i1,i2;

    surf = 0.0;
    pt = &mesh->tria[iel];

    for (i=0; i<3; i++) {
        np[i] = pt->v[i];
        p[i]  = &mesh->point[np[i]];
    }
    if ( !bezierCP(mesh,iel,&b) ) return(0.0);

    /* Set metric tensors at vertices of tria iel */
    for(i=0; i<3; i++) {
        i1 = inxt[i];
        i2 = iprv[i];
        ux = 0.5*(p[i1]->c[0]+p[i2]->c[0]) - p[i]->c[0];
        uy = 0.5*(p[i1]->c[1]+p[i2]->c[1]) - p[i]->c[1];
        uz = 0.5*(p[i1]->c[2]+p[i2]->c[2]) - p[i]->c[2];

        if ( MS_SIN(p[i]->tag) ) {
            memcpy(&m[i][0],&met->m[6*np[i]+1],6*sizeof(double));
        }
        else if ( p[i]->tag & MS_GEO ) {
            if ( !buildridmet(mesh,met,np[i],ux,uy,uz,&m[i][0]) )  return(0.0);
        }
        else {
            memcpy(&m[i][0],&met->m[6*np[i]+1],6*sizeof(double));
        }
    }

    /* Compute density integrand of volume at the 3 vertices of T */
    for (i=0; i<3; i++) {
        if ( i == 0 ) {
            J[0][0] = 3.0*( b.b[7][0] - b.b[0][0] ) ; J[0][1] = 3.0*( b.b[6][0] - b.b[0][0] );
            J[1][0] = 3.0*( b.b[7][1] - b.b[0][1] ) ; J[1][1] = 3.0*( b.b[6][1] - b.b[0][1] );
            J[2][0] = 3.0*( b.b[7][2] - b.b[0][2] ) ; J[2][1] = 3.0*( b.b[6][2] - b.b[0][2] );
        }
        else if ( i == 1 ) {
            J[0][0] = 3.0*( b.b[1][0] - b.b[8][0] ) ; J[0][1] = 3.0*( b.b[3][0] - b.b[8][0] );
            J[1][0] = 3.0*( b.b[1][1] - b.b[8][1] ) ; J[1][1] = 3.0*( b.b[3][1] - b.b[8][1] );
            J[2][0] = 3.0*( b.b[1][2] - b.b[8][2] ) ; J[2][1] = 3.0*( b.b[3][2] - b.b[8][2] );
        }
        else {
            J[0][0] = 3.0*( b.b[4][0] - b.b[5][0] ) ; J[0][1] = 3.0*( b.b[2][0] - b.b[5][0] );
            J[1][0] = 3.0*( b.b[4][1] - b.b[5][1] ) ; J[1][1] = 3.0*( b.b[2][1] - b.b[5][1] );
            J[2][0] = 3.0*( b.b[4][2] - b.b[5][2] ) ; J[2][1] = 3.0*( b.b[2][2] - b.b[5][2] );
        }

        mJ[0][0] = m[i][0]*J[0][0] + m[i][1]*J[1][0] + m[i][2]*J[2][0];
        mJ[1][0] = m[i][1]*J[0][0] + m[i][3]*J[1][0] + m[i][4]*J[2][0];
        mJ[2][0] = m[i][2]*J[0][0] + m[i][4]*J[1][0] + m[i][5]*J[2][0];

        mJ[0][1] = m[i][0]*J[0][1] + m[i][1]*J[1][1] + m[i][2]*J[2][1];
        mJ[1][1] = m[i][1]*J[0][1] + m[i][3]*J[1][1] + m[i][4]*J[2][1];
        mJ[2][1] = m[i][2]*J[0][1] + m[i][4]*J[1][1] + m[i][5]*J[2][1];

        /* dens = sqrt(tJacsigma * M * Jacsigma )*/
        tJmJ[0][0] = J[0][0]*mJ[0][0] + J[1][0]*mJ[1][0] + J[2][0]*mJ[2][0];
        tJmJ[0][1] = J[0][0]*mJ[0][1] + J[1][0]*mJ[1][1] + J[2][0]*mJ[2][1];
        tJmJ[1][0] = J[0][1]*mJ[0][0] + J[1][1]*mJ[1][0] + J[2][1]*mJ[2][0];
        tJmJ[1][1] = J[0][1]*mJ[0][1] + J[1][1]*mJ[1][1] + J[2][1]*mJ[2][1];

        dens = tJmJ[0][0]*tJmJ[1][1] - tJmJ[1][0]*tJmJ[0][1];
        if ( dens < 0.0 ) {
            //fprintf(stdout,"  ## Density should be positive : %E for elt %d %d %d \n",dens,pt->v[0],pt->v[1],pt->v[2]);
            //return(0.0);
        }
        surf += sqrt(fabs(dens));
    }

    surf *= ATHIRD;
    return(surf);
}

/* Compute the intersected (2 x 2) metric from metrics m and n : take simultaneous reduction,
   and proceed to truncation in sizes */
static int intersecmet22(double *m,double *n,double *mr) {
    double  det,imn[4],dd,sqDelta,trimn,lambda[2],vp0[2],vp1[2],dm[2],dn[2],vnorm,d0,d1,ip[4];
    double  isqhmin,isqhmax;

    isqhmin  = 1.0 / (info.hmin*info.hmin);
    isqhmax  = 1.0 / (info.hmax*info.hmax);

    /* Compute imn = M^{-1}N */
    det = m[0]*m[2] - m[1]*m[1];
    if ( fabs(det) < EPS*EPS ) {
        printf("  ## Function intersecmet : null metric det : %E \n",det);
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
        printf(" ## Eigenvalues : %f \n",lambda[0]);
        return(0);
    }

    /* First case : matrices m and n are homothetic : n = lambda0*m */
    if ( sqDelta < EPS ) {
        /* Diagonalize m and truncate eigenvalues : trimn, det, etc... are reused */
        dd    = m[0] - m[2];
        trimn = m[0] + m[2];
        det   = m[0]*m[2] - m[1]*m[1];

        sqDelta = sqrt(fabs(dd*dd +4*0*m[1]*m[1]));
        dm[0]   = 0.5 * (trimn + sqDelta);
        dm[1]   = 0.5 * (trimn - sqDelta);

        vp0[0] = m[1];
        vp0[1] = (dm[0]-m[0]);
        vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);

        if ( vnorm < EPS ) {
            vp0[0] = (dm[0] - m[2]);
            vp0[1] = m[1];
            vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);
            if ( vnorm < EPS ) return(0);
        }

        vnorm   = 1.0 / vnorm;
        vp0[0] *= vnorm;
        vp0[1] *= vnorm;

        vp1[0] = m[1];
        vp1[1] = (dm[1]-m[0]);
        vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);

        if ( vnorm < EPS ) {
            vp1[0] = (dm[1] - m[2]);
            vp1[1] = m[1];
            vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);
            if ( vnorm < EPS ) return(0);
        }

        vnorm   = 1.0 / vnorm;
        vp1[0] *= vnorm;
        vp1[1] *= vnorm;

        /* Eigenvalues of the resulting matrix*/
        dn[0] = MS_MAX(dm[0],lambda[0]*dm[0]);
        dn[0] = MS_MIN(isqhmin,MS_MAX(isqhmax,dn[0]));
        dn[1] = MS_MAX(dm[1],lambda[0]*dm[1]);
        dn[1] = MS_MIN(isqhmin,MS_MAX(isqhmax,dn[1]));

        /* Intersected metric = P diag(d0,d1){^t}P, P = (vp0, vp1) stored in columns */
        mr[0] = dn[0]*vp0[0]*vp0[0] + dn[1]*vp1[0]*vp1[0];
        mr[1] = dn[0]*vp0[0]*vp0[1] + dn[1]*vp1[0]*vp1[1];
        mr[2] = dn[0]*vp0[1]*vp0[1] + dn[1]*vp1[1]*vp1[1];

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

        /* Diagonal values of the intersected metric */
        d0 = MS_MAX(dm[0],dn[0]);
        d0 = MS_MIN(isqhmin,MS_MAX(d0,isqhmax));

        d1 = MS_MAX(dm[1],dn[1]);
        d1 = MS_MIN(isqhmin,MS_MAX(d1,isqhmax));

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

/* Intersect metric held in np (supported in tangent plane of np) with 3*3 metric in me */
int intextmet(pMesh mesh,pSol met,int np,double me[6]) {
    pPoint         p0;
    pGeom          go;
    double         hu,isqhmin,isqhmax,dd;
    double        *m,*n,*n1,*n2,*t,r[3][3],mrot[6],mr[3],mtan[3],metan[3],u[3],a[4];
    double complex ro[3];
    char           i;

    isqhmin = 1.0 / (info.hmin*info.hmin);
    isqhmax = 1.0 / (info.hmax*info.hmax);

    p0 = &mesh->point[np];
    m  = &met->m[6*np+1];

    /* Case of a singular point : take smallest size prescribed by met, or me in every direction */
    if ( MS_SIN(p0->tag) ) {
        /* Characteristic polynomial of me */
        a[3] = -1.0;
        a[2] = me[0]+me[3]+me[5];
        a[1] = -(me[0]*me[3]+me[0]*me[5]+me[3]*me[5]) + (me[1]*me[1]+me[2]*me[2]+me[4]*me[4]);
        a[0] = me[0]*(me[3]*me[5]-me[4]*me[4]) -me[1]*(me[1]*me[5]-me[2]*me[4]) \
            + me[2]*(me[1]*me[4]-me[2]*me[3]);

        rootDeg3(a,ro);
        hu = m[0];
        for(i=0; i<3; i++) {
            if( cimag(ro[i]) != 0.0 )
                break;
            else
                hu = MS_MAX(hu,creal(ro[i]));
        }
        hu = MS_MIN(isqhmin,hu);
        hu = MS_MAX(isqhmax,hu);
        m[0] = hu;
        m[3] = hu;
        m[5] = hu;
    }
    /* Case of a ridge point : take sizes in 3 directions t,n1,u */
    else if ( p0->tag & MS_GEO ) {
        /* Size prescribed by metric me in direction t */
        t = &p0->n[0];
        hu = me[0]*t[0]*t[0] + me[3]*t[1]*t[1] + me[5]*t[2]*t[2] \
            + 2.0*me[1]*t[0]*t[1] + 2.0*me[2]*t[0]*t[2] + 2.0*me[4]*t[1]*t[2];

        hu = MS_MIN(isqhmin,hu);
        hu = MS_MAX(isqhmax,hu);
        m[0] = MS_MAX(m[0],hu);

        /* Size prescribed by metric me in direction u1 = n1 ^ t */
        assert ( p0->ig );
        go = &mesh->geom[p0->ig];
        n1 = &go->n1[0];
        n2 = &go->n2[0];

        u[0] = n1[1]*t[2] - n1[2]*t[1];
        u[1] = n1[2]*t[0] - n1[0]*t[2];
        u[2] = n1[0]*t[1] - n1[1]*t[0];
        dd = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
        if ( dd < EPSD ) return(0);
        dd = 1.0 / sqrt(dd);

        u[0] *= dd;
        u[1] *= dd;
        u[2] *= dd;

        hu = me[0]*u[0]*u[0] + me[3]*u[1]*u[1] + me[5]*u[2]*u[2] \
            + 2.0*me[1]*u[0]*u[1] + 2.0*me[2]*u[0]*u[2] + 2.0*me[4]*u[1]*u[2];

        hu = MS_MIN(isqhmin,hu);
        hu = MS_MAX(isqhmax,hu);
        m[1] = MS_MAX(m[1],hu);

        /* Size prescribed by metric me in direction u1 = n1 ^ t */
        u[0] = n2[1]*t[2] - n2[2]*t[1];
        u[1] = n2[2]*t[0] - n2[0]*t[2];
        u[2] = n2[0]*t[1] - n2[1]*t[0];
        dd = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
        if ( dd < EPSD ) return(0);
        dd = 1.0 / sqrt(dd);

        u[0] *= dd;
        u[1] *= dd;
        u[2] *= dd;

        hu =     me[0]*u[0]*u[0] +     me[3]*u[1]*u[1] +     me[5]*u[2]*u[2] \
            + 2.0*me[1]*u[0]*u[1] + 2.0*me[2]*u[0]*u[2] + 2.0*me[4]*u[1]*u[2];

        hu = MS_MIN(isqhmin,hu);
        hu = MS_MAX(isqhmax,hu);
        m[2] = MS_MAX(m[2],hu);
    }
    /* Case of a ref, or regular point : intersect metrics in tangent plane */
    else {
        n = &p0->n[0];
        rotmatrix(n,r);

        /* Expression of both metrics in tangent plane */
        rmtr(r,m,mrot);
        mtan[0] = mrot[0];
        mtan[1] = mrot[1];
        mtan[2] = mrot[3];

        rmtr(r,me,mrot);
        metan[0] = mrot[0];
        metan[1] = mrot[1];
        metan[2] = mrot[3];

        /* Intersection of metrics in the tangent plane */
        if ( !intersecmet22(mtan,metan,mr) ) return(0);

        /* Back to the canonical basis of \mathbb{R}^3 : me = ^tR*mr*R : mtan and metan are reused */
        mtan[0]  = mr[0]*r[0][0] + mr[1]*r[1][0];  mtan[1]  = mr[0]*r[0][1] + mr[1]*r[1][1];   mtan[2]  = mr[0]*r[0][2] + mr[1]*r[1][2] ;
        metan[0] = mr[1]*r[0][0] + mr[2]*r[1][0];  metan[1] = mr[1]*r[0][1] + mr[2]*r[1][1];   metan[2] = mr[1]*r[0][2] + mr[2]*r[1][2] ;

        m[0] = r[0][0] * mtan[0] + r[1][0] * metan[0];
        m[1] = r[0][0] * mtan[1] + r[1][0] * metan[1];
        m[2] = r[0][0] * mtan[2] + r[1][0] * metan[2];
        m[3] = r[0][1] * mtan[1] + r[1][1] * metan[1];
        m[4] = r[0][1] * mtan[2] + r[1][1] * metan[2];
        m[5] = r[0][2] * mtan[2] + r[1][2] * metan[2];
    }

    return(1);
}

