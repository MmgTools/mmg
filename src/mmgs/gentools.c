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
 * \file mmgs/gentool.c
 * \brief Generic algebraic and algorithmic tools.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgs.h"

extern Info  info;

/* Compute rotation matrix that sends vector n to the third vector of canonical basis */
inline int rotmatrix(double n[3],double r[3][3]) {
    double aa,bb,ab,ll,l,cosalpha,sinalpha;

    aa = n[0]*n[0];
    bb = n[1]*n[1];
    ab = n[0]*n[1];
    ll = aa+bb;
    cosalpha = n[2];
    sinalpha = sqrt(1.0- MS_MIN(1.0,cosalpha*cosalpha));

    /* No rotation needed in this case */
    if ( ll < EPS ) {
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
    return(1);
}

double surftri_iso(pMesh mesh,pSol met,int iel) {
    pTria    pt;
    double   *a,*b,*c,abx,aby,abz,acx,acy,acz,det,n[3];

    pt = &mesh->tria[iel];

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

    return( 0.5*sqrt(det) );
}

/* Compute product R*M*tR when M is symmetric */
inline int rmtr(double r[3][3],double m[6], double mr[6]){
    double n[3][3];

    n[0][0] = m[0]*r[0][0] + m[1]*r[0][1] + m[2]*r[0][2];
    n[1][0] = m[1]*r[0][0] + m[3]*r[0][1] + m[4]*r[0][2];
    n[2][0] = m[2]*r[0][0] + m[4]*r[0][1] + m[5]*r[0][2];

    n[0][1] = m[0]*r[1][0] + m[1]*r[1][1] + m[2]*r[1][2];
    n[1][1] = m[1]*r[1][0] + m[3]*r[1][1] + m[4]*r[1][2];
    n[2][1] = m[2]*r[1][0] + m[4]*r[1][1] + m[5]*r[1][2];

    n[0][2] = m[0]*r[2][0] + m[1]*r[2][1] + m[2]*r[2][2];
    n[1][2] = m[1]*r[2][0] + m[3]*r[2][1] + m[4]*r[2][2];
    n[2][2] = m[2]*r[2][0] + m[4]*r[2][1] + m[5]*r[2][2];

    mr[0] = r[0][0]*n[0][0] + r[0][1]*n[1][0] + r[0][2]*n[2][0];
    mr[1] = r[0][0]*n[0][1] + r[0][1]*n[1][1] + r[0][2]*n[2][1];
    mr[2] = r[0][0]*n[0][2] + r[0][1]*n[1][2] + r[0][2]*n[2][2];
    mr[3] = r[1][0]*n[0][1] + r[1][1]*n[1][1] + r[1][2]*n[2][1];
    mr[4] = r[1][0]*n[0][2] + r[1][1]*n[1][2] + r[1][2]*n[2][2];
    mr[5] = r[2][0]*n[0][2] + r[2][1]*n[1][2] + r[2][2]*n[2][2];

    return(1);
}

/* Solve 3*3 symmetric system A*r = b */
inline int sys33sym(double a[6],double b[3],double r[3]) {
    double ia[6],as[6],det,m;
    int    i;

    /* Multiply matrix by a constant coefficient for stability purpose (because of the scaling) */
    m = fabs(a[0]);
    for (i=1; i<6; i++) {
        if ( fabs(a[i]) < m )  m = fabs(a[i]);
    }
    if ( m < EPSD )  return(0);
    m = 1.0 / m;

    for (i=0; i<6; i++) {
        as[i] = a[i]*m;
    }
    det = as[0]*(as[3]*as[5]-as[4]*as[4]) - as[1]*(as[1]*as[5]-as[2]*as[4]) \
        + as[2]*(as[1]*as[4]-as[2]*as[3]);
    if ( fabs(det) < EPSD )  return(0);
    det = 1.0 / det;

    ia[0] = (as[3]*as[5]-as[4]*as[4]);
    ia[1] = - (as[1]*as[5]-as[2]*as[4]);
    ia[2] = (as[1]*as[4]-as[2]*as[3]);
    ia[3] = (as[0]*as[5]-as[2]*as[2]);
    ia[4] = -(as[0]*as[4]-as[2]*as[1]);
    ia[5] = (as[0]*as[3]-as[1]*as[1]);

    r[0] = ia[0]*b[0] + ia[1]*b[1] + ia[2]*b[2];
    r[1] = ia[1]*b[0] + ia[3]*b[1] + ia[4]*b[2];
    r[2] = ia[2]*b[0] + ia[4]*b[1] + ia[5]*b[2];

    r[0] *= (det*m);
    r[1] *= (det*m);
    r[2] *= (det*m);

    return(1);
}

/* Compute eigenelements of a SYMMETRIC matrix m. Eigenvectors are orthogonal. Return order */
inline int eigensym(double m[3],double lambda[2],double vp[2][2]) {
    double   sqDelta,dd,trm,vnorm;

    dd  = m[0]-m[2];
    trm = m[0]+m[2];
    sqDelta = sqrt(dd*dd + 4.0*m[1]*m[1]);
    lambda[0] = 0.5*(trm - sqDelta);

    /* Case when m = lambda[0]*I */
    if ( sqDelta < EPS ) {
        lambda[1] = lambda[0];
        vp[0][0] = 1.0;
        vp[0][1] = 0.0;

        vp[1][0] = 0.0;
        vp[1][1] = 1.0;
        return(2);
    }
    vp[0][0] = m[1];
    vp[0][1] = (lambda[0] - m[0]);
    vnorm = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);

    if ( vnorm < EPS ) {
        vp[0][0] = (lambda[0] - m[2]);
        vp[0][1] = m[1];
        vnorm = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);
    }

    if( !(vnorm > EPSD) ) {
        printf("la matrice : %f %f %f \n",m[0],m[1],m[2]);
        printf("le vp qui deconne %f %f \n",vp[0][0],vp[0][1]);
    }
    assert(vnorm > EPSD);

    vnorm = 1.0/vnorm;
    vp[0][0] *= vnorm;
    vp[0][1] *= vnorm;

    vp[1][0] = -vp[0][1];
    vp[1][1] = vp[0][0];

    lambda[1] = m[0]*vp[1][0]*vp[1][0] + 2.0*m[1]*vp[1][0]*vp[1][1] + m[2]*vp[1][1]*vp[1][1];

    return(1);
}

/* Compute the intersected (2 x 2) metric between metrics m and n, PRESERVING the directions
   of m. Result is stored in mr*/
int intmetsavedir(double *m,double *n,double *mr) {
    int    i;
    double lambda[2],vp[2][2],siz,isqhmin;

    isqhmin = 1.0 / (info.hmin * info.hmin);
    eigensym(m,lambda,vp);

    for (i=0; i<2; i++) {
        siz = n[0]*vp[i][0]*vp[i][0] + 2.0*n[1]*vp[i][0]*vp[i][1] + n[2]*vp[i][1]*vp[i][1];
        lambda[i] = MS_MAX(lambda[i],siz);
        lambda[i] = MS_MIN(lambda[i],isqhmin);
    }
    mr[0] = lambda[0]*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
    mr[1] = lambda[0]*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
    mr[2] = lambda[0]*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];

    return(1);
}

/* Delete all triangle references in mesh */
int delref(pMesh mesh) {
    pTria    pt;
    int      k;

    for(k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        pt->ref = 0;
    }

    return(1);
}

/* Start from triangle start, and pile up triangles by adjacency, till a GEO or REF curve is met ;
   pass all references of travelled faces to ref ; putreq = 1 if boundary edges met must
   be set to MS_REQ, 0 otherwise. */
int setref(pMesh mesh,int start,int ref,int putreq) {
    pTria      pt,pt1;
    int        *list,*adja,cur,base,k,iel,jel,ilist;
    char       j,voy;

    ilist = cur = 0;
    list = (int*)calloc(mesh->nt+1,sizeof(int));
    base = ++mesh->base;

    /* Pile up triangles from start, till a GEO boundary is met */
    pt = &mesh->tria[start];
    list[ilist] = start;
    ilist++;
    assert( ilist <= mesh->nt );
    pt->flag = base;

    do {
        iel = list[cur];
        pt = &mesh->tria[iel];
        adja = &mesh->adja[3*(iel-1)+1];

        for(j=0; j<3; j++) {
            if( MS_EDG(pt->tag[j]) ) {
                if( putreq ) {
                    pt->tag[j] |= MS_REQ;
                    jel = adja[j] / 3;
                    voy = adja[j] % 3;
                    if( !jel ) continue;
                    pt1 = &mesh->tria[jel];
                    pt1->tag[voy] |= MS_REQ;
                }
                continue;
            }
            jel = adja[j] / 3;
            assert(jel);
            pt1 = &mesh->tria[jel];
            if ( pt1->flag == base )  continue;

            list[ilist] = jel;
            ilist++;
            assert( ilist <= mesh->nt );
            pt1->flag = base;
        }
        cur++;
    }
    while( cur < ilist );

    /* Set all references of triangles of list to ref */
    for (k=0; k<ilist; k++) {
        iel = list[k];
        pt  = &mesh->tria[iel];
        pt->ref = ref;
        printf("Le tria %d passe a %d \n",k,ref);
    }

    return(1);
}

/* invert 3x3 non-symmetric matrix */
int invmatg(double m[9],double mi[9]) {
    double  aa,bb,cc,det,vmin,vmax,maxx;
    int     k;

    /* check ill-conditionned matrix */
    vmin = vmax = fabs(m[0]);
    for (k=1; k<9; k++) {
        maxx = fabs(m[k]);
        if ( maxx < vmin )  vmin = maxx;
        else if ( maxx > vmax )  vmax = maxx;
    }
    if ( vmax == 0.0 )  return(0);

    /* compute sub-dets */
    aa = m[4]*m[8] - m[5]*m[7];
    bb = m[5]*m[6] - m[3]*m[8];
    cc = m[3]*m[7] - m[4]*m[6];
    det = m[0]*aa + m[1]*bb + m[2]*cc;
    if ( fabs(det) < EPSD )  return(0);
    det = 1.0 / det;

    mi[0] = aa*det;
    mi[3] = bb*det;
    mi[6] = cc*det;
    mi[1] = (m[2]*m[7] - m[1]*m[8])*det;
    mi[4] = (m[0]*m[8] - m[2]*m[6])*det;
    mi[7] = (m[1]*m[6] - m[0]*m[7])*det;
    mi[2] = (m[1]*m[5] - m[2]*m[4])*det;
    mi[5] = (m[2]*m[3] - m[0]*m[5])*det;
    mi[8] = (m[0]*m[4] - m[1]*m[3])*det;

    return(1);
}
