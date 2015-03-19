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
 * \file mmgs/anisosiz.c
 * \brief Fonctions for anisotropic size map computation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

extern Info   info;

/* Define anisotropic metric map at a SINGULARITY of the geometry, associated to the
   geometric approx of the surface. metric= alpha* Id, alpha = size */
static int defmetsin(pMesh mesh,pSol met,int it,int ip) {
    pTria         pt;
    pPoint        p0,p1;
    double       *m,n[3],ux,uy,uz,isqhmin,isqhmax,b0[3],b1[3],ps1,tau[3],ntau2,gammasec[3];
    double        c[3],kappa,maxkappa,alpha;
    int           ilist,list[LMAX+2],k,iel,idp;
    unsigned char i0,i1,i2;

    pt  = &mesh->tria[it];
    idp = pt->v[ip];
    p0  = &mesh->point[idp];

    ilist = boulet(mesh,it,ip,list);
    assert(ilist);

    isqhmin  = 1.0 / (info.hmin*info.hmin);
    isqhmax  = 1.0 / (info.hmax*info.hmax);
    maxkappa = 0.0;
    for (k=0; k<ilist; k++) {
        iel = list[k] / 3;
        i0  = list[k] % 3;
        i1  = inxt[i0];
        i2  = iprv[i0];
        pt  = &mesh->tria[iel];
        p1  = &mesh->point[pt->v[i1]];

        ux = p1->c[0] - p0->c[0];
        uy = p1->c[1] - p0->c[1];
        uz = p1->c[2] - p0->c[2];

        /* Computation of the two control points associated to edge p0p1: p0 is singular */
        nortri(mesh,pt,n);
        if ( MS_EDG(pt->tag[i2]) )
            bezierEdge(mesh,idp,pt->v[i1],b0,b1,1,n);
        else
            bezierEdge(mesh,idp,pt->v[i1],b0,b1,0,n);

        /* tangent vector */
        tau[0] = 3.0*(b0[0] - p0->c[0]);
        tau[1] = 3.0*(b0[1] - p0->c[1]);
        tau[2] = 3.0*(b0[2] - p0->c[2]);
        ntau2  = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];

        /* 2nd order derivative */
        gammasec[0] = 6.0*p0->c[0] - 12.0*b0[0] + 6.0*b1[0];
        gammasec[1] = 6.0*p0->c[1] - 12.0*b0[1] + 6.0*b1[1];
        gammasec[2] = 6.0*p0->c[2] - 12.0*b0[2] + 6.0*b1[2];
        if ( ntau2 < EPSD )  continue;
        ntau2 = 1.0 / ntau2;

        /* derivative via the normal parametrization */
        ps1  = gammasec[0]*tau[0] + gammasec[1]*tau[1] + gammasec[2]*tau[2];
        c[0] = gammasec[0] - ps1*tau[0]*ntau2;
        c[1] = gammasec[1] - ps1*tau[1]*ntau2;
        c[2] = gammasec[2] - ps1*tau[2]*ntau2;

        kappa = ntau2 * sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
        maxkappa = MS_MAX(kappa,maxkappa);
    }
    alpha = 1.0 / 8.0 * maxkappa / info.hausd;
    alpha = MS_MIN(alpha,isqhmin);
    alpha = MS_MAX(alpha,isqhmax);

    m = &met->m[6*(idp)+1];
    memset(m,0,6*sizeof(double));
    m[0] = m[3] = m[5] = alpha;

    return(1);
}

/* Compute metric tensor associated to a ridge point : convention is a bit weird here :
   p->m[0] is the specific size in direction t,
   p->m[1] is the specific size in direction u1 = n1 ^ t
   p->m[2] is the specific size in direction u2 = n2 ^ t,
   and at each time, metric tensor has to be recomputed, depending on the side */
static int defmetrid(pMesh mesh,pSol met,int it,int ip) {
    pTria          pt;
    pPoint         p0,p1,p2;
    Bezier         b;
    int            k,iel,idp,ilist1,ilist2,ilist,*list,list1[LMAX+2],list2[LMAX+2],iprid[2],ier;
    double        *m,isqhmin,isqhmax,*n1,*n2,*n,*t,kappacur,b0[3],b1[3],n0[3],tau[3],trot[2],u[2];
    double         l,ll,ps,gammasec[3],c[3],r[3][3],lispoi[3*LMAX+1],ux,uy,uz,det,bcu[3];
    double         detg,detd,Jacb[3][2],Hb[3][3],lambda[2];
    unsigned char  i,i0,i1,i2;

    pt  = &mesh->tria[it];
    idp = pt->v[ip];
    p0  = &mesh->point[idp];

    isqhmin = 1.0 / (info.hmin*info.hmin);
    isqhmax = 1.0 / (info.hmax*info.hmax);

    n1 = &mesh->geom[p0->ig].n1[0];
    n2 = &mesh->geom[p0->ig].n2[0];
    t  = p0->n;

    m = &met->m[6*(idp)+1];
    memset(m,0,6*sizeof(double));
    m[0] = isqhmax;
    m[1] = isqhmax;
    m[2] = isqhmax;

    ier = bouletrid(mesh,it,ip,&ilist1,list1,&ilist2,list2,&iprid[0],&iprid[1]);
    assert(ier);

    /* Specific size in direction of t */
    for (i=0; i<2; i++) {
        kappacur = 0.0;
        bezierEdge(mesh,idp,iprid[i],b0,b1,1,n0);

        tau[0] = 3.0*(b0[0] - p0->c[0]);
        tau[1] = 3.0*(b0[1] - p0->c[1]);
        tau[2] = 3.0*(b0[2] - p0->c[2]);
        ll = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];
        if ( ll < EPSD )  continue;
        l = 1.0 / sqrt(ll);
        tau[0] *= l;
        tau[1] *= l;
        tau[2] *= l;

        gammasec[0] = 6.0*p0->c[0] -12.0*b0[0] + 6.0*b1[0];
        gammasec[1] = 6.0*p0->c[1] -12.0*b0[1] + 6.0*b1[1];
        gammasec[2] = 6.0*p0->c[2] -12.0*b0[2] + 6.0*b1[2];

        ps = tau[0]*gammasec[0] + tau[1]*gammasec[1] + tau[2]*gammasec[2];
        c[0] = gammasec[0] - ps*tau[0];
        c[1] = gammasec[1] - ps*tau[1];
        c[2] = gammasec[2] - ps*tau[2];

        kappacur = MS_MAX(0.0,1.0/ll*sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]));
        kappacur = 1.0/8.0*kappacur/info.hausd;
        kappacur = MS_MIN(kappacur,isqhmin);
        kappacur = MS_MAX(kappacur,isqhmax);
        m[0] = MS_MAX(m[0],kappacur);
    }

    /* Characteristic sizes in directions u1 and u2 */
    for (i=0; i<2; i++) {
        if ( i==0 ) {
            n = n1;
            ilist = ilist1;
            list  = &list1[0];
        }
        else {
            n = n2;
            ilist = ilist2;
            list  = &(list2[0]);
        }
        assert(rotmatrix(n,r));

        /* Apply rotation to the half-ball under consideration */
        i1 = 0;
        for (k=0; k<ilist; k++) {
            iel = list[k] / 3;
            i0  = list[k] % 3;
            i1  = inxt[i0];
            pt = &mesh->tria[iel];
            p1 = &mesh->point[pt->v[i1]];

            ux = p1->c[0] - p0->c[0];
            uy = p1->c[1] - p0->c[1];
            uz = p1->c[2] - p0->c[2];

            lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
            lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
            lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
        }

        /* last point : the half-ball is open : ilist tria, and ilist +1 points ;
           lists are enumerated in direct order */
        i2 = inxt[i1];
        p2 = &mesh->point[pt->v[i2]];

        ux = p2->c[0] - p0->c[0];
        uy = p2->c[1] - p0->c[1];
        uz = p2->c[2] - p0->c[2];

        lispoi[3*ilist+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
        lispoi[3*ilist+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
        lispoi[3*ilist+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

        /* At this point, lispoi contains all the points of the half-ball of p0, rotated
           so that t_{p_0}S = [z = 0] */

        /* Rotated tangent vector (trot[2] = 0), and third direction */
        trot[0] = r[0][0]*t[0] + r[0][1]*t[1] + r[0][2]*t[2];
        trot[1] = r[1][0]*t[0] + r[1][1]*t[1] + r[1][2]*t[2];

        u[0] = -trot[1];
        u[1] =  trot[0];

        /* Travel half-ball at p0 and stop at first triangle containing u */
        for (k=0; k<ilist; k++) {
            detg = lispoi[3*k+1]*u[1] - lispoi[3*k+2]*u[0];
            detd = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
            if ( detg > 0.0 && detd > 0.0 )  break;
        }

        /* If triangle not found, try with -u */
        if ( k == ilist ) {
            u[0] *= -1.0;
            u[1] *= -1.0;

            for (k=0; k<ilist; k++) {
                detg = lispoi[3*k+1]*u[1] - lispoi[3*k+2]*u[0];
                detd = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
                if ( detg > 0.0 && detd > 0.0 )  break;
            }
        }
        if ( k == ilist )  continue;

        iel = list[k] / 3;
        i0  = list[k] % 3;
        i1  = inxt[i0];
        i2  = iprv[i0];
        pt = &mesh->tria[iel];
        if ( !bezierCP(mesh,iel,&b) )  continue;

        /* Barycentric coordinates of vector u in tria iel */
        detg = lispoi[3*k+1]*u[1] - lispoi[3*k+2]*u[0];
        detd = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
        det = detg + detd;
        if ( det < EPSD )  continue;

        det = 1.0 / det;
        bcu[0] = 0.0;
        bcu[1] = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
        bcu[1] *= det;
        assert(bcu[1] <= 1.0);
        bcu[2] = 1.0 - bcu[1];

        /* Computation of tangent vector and second derivative of curve t \mapsto b(tbcu)
           (not in rotated frame) */
        if ( i0 == 0 ) { // w = 1, u,v = 0
            lambda[0] = bcu[1];
            lambda[1] = bcu[2];

            Jacb[0][0] = 3.0*(b.b[7][0]-b.b[0][0]);
            Jacb[1][0] = 3.0*(b.b[7][1]-b.b[0][1]);
            Jacb[2][0] = 3.0*(b.b[7][2]-b.b[0][2]);

            Jacb[0][1] = 3.0*(b.b[6][0]-b.b[0][0]);
            Jacb[1][1] = 3.0*(b.b[6][1]-b.b[0][1]);
            Jacb[2][1] = 3.0*(b.b[6][2]-b.b[0][2]);

            /* Hb[i] = hessian matrix of i-th component of b at point p0 */
            Hb[0][0] = 6.0*(b.b[0][0] - 2.0*b.b[7][0] + b.b[8][0]);
            Hb[1][0] = 6.0*(b.b[0][1] - 2.0*b.b[7][1] + b.b[8][1]);
            Hb[2][0] = 6.0*(b.b[0][2] - 2.0*b.b[7][2] + b.b[8][2]);

            Hb[0][1] = 6.0*(b.b[0][0] - b.b[7][0] - b.b[6][0] + b.b[9][0]);
            Hb[1][1] = 6.0*(b.b[0][1] - b.b[7][1] - b.b[6][1] + b.b[9][1]);
            Hb[2][1] = 6.0*(b.b[0][2] - b.b[7][2] - b.b[6][2] + b.b[9][2]);

            Hb[0][2] = 6.0*(b.b[0][0] - 2.0*b.b[6][0] + b.b[5][0]);
            Hb[1][2] = 6.0*(b.b[0][1] - 2.0*b.b[6][1] + b.b[5][1]);
            Hb[2][2] = 6.0*(b.b[0][2] - 2.0*b.b[6][2] + b.b[5][2]);
        }
        else if ( i0 == 1 ) {  //w = v = 0, u = 1;
            lambda[0] = bcu[0];
            lambda[1] = bcu[1];

            Jacb[0][0] = 3.0*(b.b[1][0]-b.b[8][0]);
            Jacb[1][0] = 3.0*(b.b[1][1]-b.b[8][1]);
            Jacb[2][0] = 3.0*(b.b[1][2]-b.b[8][2]);

            Jacb[0][1] = 3.0*(b.b[3][0]-b.b[8][0]);
            Jacb[1][1] = 3.0*(b.b[3][1]-b.b[8][1]);
            Jacb[2][1] = 3.0*(b.b[3][2]-b.b[8][2]);

            Hb[0][0] = 6.0*(b.b[1][0] - 2.0*b.b[8][0] + b.b[7][0]);
            Hb[1][0] = 6.0*(b.b[1][1] - 2.0*b.b[8][1] + b.b[7][1]);
            Hb[2][0] = 6.0*(b.b[1][2] - 2.0*b.b[8][2] + b.b[7][2]);

            Hb[0][1] = 6.0*(b.b[7][0] - b.b[8][0] - b.b[9][0] + b.b[3][0]);
            Hb[1][1] = 6.0*(b.b[7][1] - b.b[8][1] - b.b[9][1] + b.b[3][1]);
            Hb[2][1] = 6.0*(b.b[7][2] - b.b[8][2] - b.b[9][2] + b.b[3][2]);

            Hb[0][2] = 6.0*(b.b[4][0] - 2.0*b.b[9][0] + b.b[7][0]);
            Hb[1][2] = 6.0*(b.b[4][1] - 2.0*b.b[9][1] + b.b[7][1]);
            Hb[2][2] = 6.0*(b.b[4][2] - 2.0*b.b[9][2] + b.b[7][2]);
        }
        else {   //w =u = 0, v =1
            lambda[0] = bcu[2];
            lambda[1] = bcu[0];

            Jacb[0][0] = 3.0*(b.b[4][0]-b.b[5][0]);
            Jacb[1][0] = 3.0*(b.b[4][1]-b.b[5][1]);
            Jacb[2][0] = 3.0*(b.b[4][2]-b.b[5][2]);

            Jacb[0][1] = 3.0*(b.b[2][0]-b.b[5][0]);
            Jacb[1][1] = 3.0*(b.b[2][1]-b.b[5][1]);
            Jacb[2][1] = 3.0*(b.b[2][2]-b.b[5][2]);

            Hb[0][0] = 6.0*(b.b[3][0] - 2.0*b.b[9][0] + b.b[6][0]);
            Hb[1][0] = 6.0*(b.b[3][1] - 2.0*b.b[9][1] + b.b[6][1]);
            Hb[2][0] = 6.0*(b.b[3][2] - 2.0*b.b[9][2] + b.b[6][2]);

            Hb[0][1] = 6.0*(b.b[4][0] - b.b[5][0] - b.b[9][0] + b.b[6][0]);
            Hb[1][1] = 6.0*(b.b[4][1] - b.b[5][1] - b.b[9][1] + b.b[6][1]);
            Hb[2][1] = 6.0*(b.b[4][2] - b.b[5][2] - b.b[9][2] + b.b[6][2]);

            Hb[0][2] = 6.0*(b.b[2][0] - 2.0*b.b[5][0] + b.b[6][0]);
            Hb[1][2] = 6.0*(b.b[2][1] - 2.0*b.b[5][1] + b.b[6][1]);
            Hb[2][2] = 6.0*(b.b[2][2] - 2.0*b.b[5][2] + b.b[6][2]);
        }

        /* tau = jacb *(lambda0,lambda1)*/
        tau[0] = Jacb[0][0]*lambda[0] + Jacb[0][1]*lambda[1];
        tau[1] = Jacb[1][0]*lambda[0] + Jacb[1][1]*lambda[1];
        tau[2] = Jacb[2][0]*lambda[0] + Jacb[2][1]*lambda[1];
        ll = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];
        if ( ll < EPSD )  continue;

        l = 1.0 / sqrt(ll);
        tau[0] *= l;
        tau[1] *= l;
        tau[2] *= l;

        gammasec[0] = Hb[0][0]*lambda[0]*lambda[0] + 2.0*Hb[0][1]*lambda[0]*lambda[1] + Hb[0][2]*lambda[1]*lambda[1];
        gammasec[1] = Hb[1][0]*lambda[0]*lambda[0] + 2.0*Hb[1][1]*lambda[0]*lambda[1] + Hb[1][2]*lambda[1]*lambda[1];
        gammasec[2] = Hb[2][0]*lambda[0]*lambda[0] + 2.0*Hb[2][1]*lambda[0]*lambda[1] + Hb[2][2]*lambda[1]*lambda[1];

        ps = tau[0]*gammasec[0] + tau[1]*gammasec[1] + tau[2]*gammasec[2];
        c[0] = gammasec[0] - ps*tau[0];
        c[1] = gammasec[1] - ps*tau[1];
        c[2] = gammasec[2] - ps*tau[2];

        kappacur = MS_MAX(0.0,1.0/ll*sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]));
        kappacur = 1.0/8.0 * kappacur/info.hausd;
        kappacur = MS_MIN(kappacur,isqhmin);
        kappacur = MS_MAX(kappacur,isqhmax);

        m[i+1] = MS_MAX(m[i+1],kappacur);
    }

    return(1);
}

/* Define anisotropic metric map at a REF point of the geometry, associated to the
   geometric approx of the surface.*/
static int defmetref(pMesh mesh,pSol met,int it,int ip) {
    pTria         pt;
    pPoint        p0,p1;
    Bezier        b;
    int           ilist,list[LMAX+2],k,iel,ipref[2],idp;
    double        *m,isqhmin,isqhmax,*n,*t,l,ll,r[3][3],lispoi[3*LMAX+1];
    double        ux,uy,uz,det2d,intm[3],tau[2],b0[3],b1[3],kappa[2],vp[2][2],c[3];
    double        ps1,gammasec[3],kappacur,*t1,tAA[6],tAb[3],d[3];
    unsigned char i0,i1,i2,i,j;

    ipref[0] = ipref[1] = 0;
    pt  = &mesh->tria[it];
    idp = pt->v[ip];
    p0  = &mesh->point[idp];
    ilist = boulet(mesh,it,ip,list);
    assert(ilist);

    isqhmin = 1.0 / (info.hmin*info.hmin);
    isqhmax = 1.0 / (info.hmax*info.hmax);

    /* Computation of the rotation matrix T_p0 S -> [z = 0] */
    n  = &mesh->geom[p0->ig].n1[0];
    assert(rotmatrix(n,r));
    m = &met->m[6*(idp)+1];

    /* Apply rotation \circ translation to the whole ball */
    for (k=0; k<ilist; k++) {
        iel = list[k] / 3;
        i0  = list[k] % 3;
        i1  = inxt[i0];
        i2  = iprv[i0];
        pt = &mesh->tria[iel];
        p1 = &mesh->point[pt->v[i1]];

        /* Store the two ending points of ref curves */
        if ( MS_REF & pt->tag[i1] ) {
            if ( !ipref[0] ) {
                ipref[0] = pt->v[i2];
            }
            else if ( !ipref[1] && (pt->v[i2] != ipref[0]) ) {
                ipref[1] = pt->v[i2];
            }
            else if ( (pt->v[i2] != ipref[0]) && (pt->v[i2] != ipref[1]) ) {
                printf("Problem (func defmetref) : three adjacent ref at a non singular point\n");
                exit(0);
            }
        }

        if ( MS_REF & pt->tag[i2] ) {
            if ( !ipref[0] ) {
                ipref[0] = pt->v[i1];
            }
            else if ( !ipref[1] && (pt->v[i1] != ipref[0]) ) {
                ipref[1] = pt->v[i1];
            }
            else if ( (pt->v[i1] != ipref[0]) && (pt->v[i1] != ipref[1]) ) {
                printf("Problem (func defmetref) : three adjacent ref at a non singular point\n");
                exit(0);
            }
        }

        ux = p1->c[0] - p0->c[0];
        uy = p1->c[1] - p0->c[1];
        uz = p1->c[2] - p0->c[2];

        lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
        lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
        lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
    }

    /* list goes modulo ilist */
    lispoi[3*ilist+1] =  lispoi[1];
    lispoi[3*ilist+2] =  lispoi[2];
    lispoi[3*ilist+3] =  lispoi[3];

    /* Check all projections over tangent plane. */
    for (k=0; k<ilist-1; k++) {
        det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
        if ( det2d < 0.0 ) {
            printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
            return(0);
        }
    }
    det2d = lispoi[3*(ilist-1)+1]*lispoi[3*0+2] - lispoi[3*(ilist-1)+2]*lispoi[3*0+1];
    if ( det2d < 0.0 ) {
        printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
        return(0);
    }
    assert(ipref[0] && ipref[1]);

    /* At this point, lispoi contains all the points of the ball of p0, rotated
       so that t_{p_0}S = [z = 0], ipref1 and ipref2 are the indices of other ref points. */

    /* Second step : reconstitution of the curvature tensor at p0 in the tangent plane,
       with a quadric fitting approach */
    memset(intm,0.0,3*sizeof(double));
    memset(tAA,0.0,6*sizeof(double));
    memset(tAb,0.0,3*sizeof(double));

    for (k=0; k<ilist; k++) {
        /* Approximation of the curvature in the normal section associated to tau : by assumption,
           p1 is either regular, either on a ridge (or a singularity), but p0p1 is not ridge*/
        iel = list[k] / 3;
        i0  = list[k] % 3;
        i1  = inxt[i0];
        pt = &mesh->tria[iel];
        bezierCP(mesh,iel,&b);

        for(j=0; j<10; j++){
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
        else if(i0 == 1){
            memcpy(b0,&(b.b[3][0]),3*sizeof(double));
            memcpy(b1,&(b.b[4][0]),3*sizeof(double));
        }
        else{
            memcpy(b0,&(b.b[5][0]),3*sizeof(double));
            memcpy(b1,&(b.b[6][0]),3*sizeof(double));
        }

        /* At this point, the two control points are expressed in the rotated frame */
        c[0] = 3.0/8.0*b0[0] + 3.0/8.0*b1[0] + 1.0/8.0*lispoi[3*k+1];
        c[1] = 3.0/8.0*b0[1] + 3.0/8.0*b1[1] + 1.0/8.0*lispoi[3*k+2];
        c[2] = 3.0/8.0*b0[2] + 3.0/8.0*b1[2] + 1.0/8.0*lispoi[3*k+3];

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

        tAA[0] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1];
        tAA[1] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2];
        tAA[2] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2];
        tAA[3] += lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2];
        tAA[4] += lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2];
        tAA[5] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2];

        tAb[0] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+3];
        tAb[1] += lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+3];
        tAb[2] += lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+3];

        /* Mid-point along median edge and endpoint in the rotated frame */
        if(i0 == 0){
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
        else if(i0 == 1){
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
        else{
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
    if ( !sys33sym(tAA,tAb,c) )  return(0);

    intm[0] = 2.0*c[0];
    intm[1] = c[2];
    intm[2] = 2.0*c[1];

    /* At this point, intm stands for the integral matrix of Taubin's approach : vp[0] and vp[1]
       are the two pr. directions of curvature, and the two curvatures can be inferred from lambdas*/
    assert(eigensym(intm,kappa,vp));

    /* Truncation of eigenvalues */
    kappa[0] = 2.0/9.0 * fabs(kappa[0])/info.hausd;
    kappa[0] = MS_MIN(kappa[0],isqhmin);
    kappa[0] = MS_MAX(kappa[0],isqhmax);

    kappa[1] = 2.0/9.0 * fabs(kappa[1])/info.hausd;
    kappa[1] = MS_MIN(kappa[1],isqhmin);
    kappa[1] = MS_MAX(kappa[1],isqhmax);

    /* Send back the metric to the canonical basis of tangent plane :
       diag(lambda) = {^t}vp * M * vp, M = vp * diag(lambda) * {^t}vp */
    intm[0] = kappa[0]*vp[0][0]*vp[0][0] + kappa[1]*vp[1][0]*vp[1][0];
    intm[1] = kappa[0]*vp[0][0]*vp[0][1] + kappa[1]*vp[1][0]*vp[1][1];
    intm[2] = kappa[0]*vp[0][1]*vp[0][1] + kappa[1]*vp[1][1]*vp[1][1];

    /* Now express metric with respect to the approx of the underlying ref curve */
    t = &p0->n[0];
    kappacur = 0.0;

    for (i=0; i<2; i++) {
        p1 = &mesh->point[ipref[i]];
        ux = p1->c[0] - p0->c[0];
        uy = p1->c[1] - p0->c[1];
        uz = p1->c[2] - p0->c[2];

        ps1 =  ux*t[0] + uy*t[1] + uz*t[2];
        c[0] = ATHIRD*ps1*t[0];
        c[1] = ATHIRD*ps1*t[1];
        c[2] = ATHIRD*ps1*t[2];

        b0[0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
        b0[1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
        b0[2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];

        if ( (MS_CRN & p1->tag) || (MS_NOM & p1->tag) ) {
            c[0] = p1->c[0] - ATHIRD*ux;
            c[1] = p1->c[1] - ATHIRD*uy;
            c[2] = p1->c[2] - ATHIRD*uz;
        }
        else {
            assert(MS_REF & p1->tag);
            t1 = &(p1->n[0]);
            ps1 =  -(ux*t1[0] + uy*t1[1] + uz*t1[2]);
            c[0] = p1->c[0] + ATHIRD*ps1*t1[0];
            c[1] = p1->c[1] + ATHIRD*ps1*t1[1];
            c[2] = p1->c[2] + ATHIRD*ps1*t1[2];
        }
        c[0] -= p0->c[0];
        c[1] -= p0->c[1];
        c[2] -= p0->c[2];

        b1[0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
        b1[1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
        b1[2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];

        /* Everything is expressed in the rotated frame */
        tau[0] = 3.0*b0[0];
        tau[1] = 3.0*b0[1];
        ll = tau[0]*tau[0] + tau[1]*tau[1];
        if ( ll < EPSD ) {
            kappacur = isqhmax;
            continue;
        }
        l = 1.0 / sqrt(ll);
        tau[0] *= l;
        tau[1] *= l;

        gammasec[0] = -12.0*b0[0] + 6.0*b1[0];
        gammasec[1] = -12.0*b0[1] + 6.0*b1[1];
        gammasec[2] = -12.0*b0[2] + 6.0*b1[2];

        ps1 = tau[0]*gammasec[0] + tau[1]*gammasec[1];
        c[0] = gammasec[0] - ps1*tau[0];
        c[1] = gammasec[1] - ps1*tau[1];
        c[2] = gammasec[2];

        kappacur = MS_MAX(kappacur,MS_MAX(0.0,1.0/ll*fabs(c[2]))); // p.s. with normal at p0
    }

    /* Rotation of tangent vector : tau is reused */
    c[0] =  r[0][0]*t[0] + r[0][1]*t[1] + r[0][2]*t[2];
    c[1] =  r[1][0]*t[0] + r[1][1]*t[1] + r[1][2]*t[2];
    c[2] =  r[2][0]*t[0] + r[2][1]*t[1] + r[2][2]*t[2];
    memcpy(tau,&c[0],2*sizeof(double));

    /* Truncation of curvature */
    kappacur = 1.0/8.0 * kappacur/info.hausd;
    kappacur = MS_MIN(kappacur,isqhmin);
    kappacur = MS_MAX(kappacur,isqhmax);

    /* The associated matrix in basis (rt, orth rt) */
    c[0] = kappacur*tau[0]*tau[0] + isqhmax*tau[1]*tau[1];
    c[1] = (kappacur - isqhmax)*tau[0]*tau[1];
    c[2] = kappacur*tau[1]*tau[1] + isqhmax*tau[0]*tau[0];

    /* Reuse b0 for commodity */
    assert(intmetsavedir(c,intm,b0));
    memcpy(intm,b0,3*sizeof(double));

    /* At this point, intm (with 0 in the z direction)  is the desired metric, except
       it is expressed in the rotated bc, that is intm = R * metric in bc * ^t R,
       so metric in bc = ^tR*intm*R */

    /* b0 and b1 serve now for nothing : let them be the lines of matrix intm*R  */
    b0[0] = intm[0]*r[0][0] + intm[1]*r[1][0];  b0[1] = intm[0]*r[0][1] + intm[1]*r[1][1];  b0[2] = intm[0]*r[0][2] + intm[1]*r[1][2];
    b1[0] = intm[1]*r[0][0] + intm[2]*r[1][0];  b1[1] = intm[1]*r[0][1] + intm[2]*r[1][1];  b1[2] = intm[1]*r[0][2] + intm[2]*r[1][2];

    m[0] = r[0][0] * b0[0] + r[1][0] * b1[0];
    m[1] = r[0][0] * b0[1] + r[1][0] * b1[1];
    m[2] = r[0][0] * b0[2] + r[1][0] * b1[2];

    m[3] = r[0][1] * b0[1] + r[1][1] * b1[1];
    m[4] = r[0][1] * b0[2] + r[1][1] * b1[2];

    m[5] = r[0][2] * b0[2] + r[1][2] * b1[2];

    return(1);
}

/* Define anisotropic metric map at a REGULAR vertex of the mesh, associated to the
   geometric approx of the surface.*/
static int defmetreg(pMesh mesh,pSol met,int it,int ip) {
    pTria          pt;
    pPoint         p0,p1;
    Bezier         b;
    int            ilist,list[LMAX+2],k,iel,idp;
    double        *n,*m,r[3][3],ux,uy,uz,lispoi[3*LMAX+1];
    double         det2d,intm[3],b0[3],b1[3],c[3],isqhmin,isqhmax;
    double         kappa[2],vp[2][2],tAA[6],tAb[3],d[3];
    unsigned char  i0,i1,j;

    pt  = &mesh->tria[it];
    idp = pt->v[ip];
    p0  = &mesh->point[idp];

    ilist = boulet(mesh,it,ip,list);
    assert(ilist);

    isqhmin = 1.0 / (info.hmin*info.hmin);
    isqhmax = 1.0 / (info.hmax*info.hmax);

    /* Computation of the rotation matrix T_p0 S -> [z = 0] */
    n  = &p0->n[0];
    assert(rotmatrix(n,r));
    m = &met->m[6*(idp)+1];

    /* Apply rotation \circ translation to the whole ball */
    for (k=0; k<ilist; k++) {
        iel = list[k] / 3;
        i0  = list[k] % 3;
        i1  = inxt[i0];
        pt = &mesh->tria[iel];
        p1 = &mesh->point[pt->v[i1]];

        ux = p1->c[0] - p0->c[0];
        uy = p1->c[1] - p0->c[1];
        uz = p1->c[2] - p0->c[2];

        lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
        lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
        lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
    }

    /* list goes modulo ilist */
    lispoi[3*ilist+1] =  lispoi[1];
    lispoi[3*ilist+2] =  lispoi[2];
    lispoi[3*ilist+3] =  lispoi[3];

    /* Check all projections over tangent plane. */
    for (k=0; k<ilist-1; k++) {
        det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
        if ( det2d < 0.0 ) {
            printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
            return(0);
        }
    }
    det2d = lispoi[3*(ilist-1)+1]*lispoi[3*0+2] - lispoi[3*(ilist-1)+2]*lispoi[3*0+1];
    if ( det2d < 0.0 ) {
        printf("PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
        return(0);
    }

    /* At this point, lispoi contains all the points of the ball of p0, rotated
       so that t_{p_0}S = [z = 0] */

    /* Second step : reconstitution of the curvature tensor at p0 in the tangent plane,
       with a quadric fitting approach */
    memset(intm,0.0,3*sizeof(double));
    memset(tAA,0.0,6*sizeof(double));
    memset(tAb,0.0,3*sizeof(double));

    for (k=0; k<ilist; k++) {
        /* Approximation of the curvature in the normal section associated to tau : by assumption,
           p1 is either regular, either on a ridge (or a singularity), but p0p1 is not ridge*/
        iel = list[k] / 3;
        i0  = list[k] % 3;
        i1  = inxt[i0];
        pt = &mesh->tria[iel];
        bezierCP(mesh,iel,&b);

        for(j=0; j<10; j++){
            c[0] = b.b[j][0] - p0->c[0];
            c[1] = b.b[j][1] - p0->c[1];
            c[2] = b.b[j][2] - p0->c[2];

            b.b[j][0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
            b.b[j][1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
            b.b[j][2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];
        }

        /* Mid-point along left edge and endpoint in the rotated frame */
        if(i0 == 0){
            memcpy(b0,&(b.b[7][0]),3*sizeof(double));
            memcpy(b1,&(b.b[8][0]),3*sizeof(double));
        }
        else if(i0 == 1){
            memcpy(b0,&(b.b[3][0]),3*sizeof(double));
            memcpy(b1,&(b.b[4][0]),3*sizeof(double));
        }
        else{
            memcpy(b0,&(b.b[5][0]),3*sizeof(double));
            memcpy(b1,&(b.b[6][0]),3*sizeof(double));
        }

        /* At this point, the two control points are expressed in the rotated frame */
        c[0] = 3.0/8.0*b0[0] + 3.0/8.0*b1[0] + 1.0/8.0*lispoi[3*k+1];
        c[1] = 3.0/8.0*b0[1] + 3.0/8.0*b1[1] + 1.0/8.0*lispoi[3*k+2];
        c[2] = 3.0/8.0*b0[2] + 3.0/8.0*b1[2] + 1.0/8.0*lispoi[3*k+3];

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

        tAA[0] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1];
        tAA[1] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2];
        tAA[2] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2];
        tAA[3] += lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2];
        tAA[4] += lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2];
        tAA[5] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2];

        tAb[0] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+3];
        tAb[1] += lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+3];
        tAb[2] += lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+3];

        /* Mid-point along median edge and endpoint in the rotated frame */
        if(i0 == 0){
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
        else if(i0 == 1){
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
        else{
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
    if ( !sys33sym(tAA,tAb,c) ) {
        //printf(" La matrice %f %f %f %f %f %f \n",tAA[0],tAA[1],tAA[2],tAA[3],tAA[4],tAA[5]);
        return(0);
    }

    intm[0] = 2.0*c[0];
    intm[1] = c[2];
    intm[2] = 2.0*c[1];

    /* At this point, intm stands for the integral matrix of Taubin's approach : vp[0] and vp[1]
       are the two pr. directions of curvature, and the two curvatures can be inferred from lambdas*/
    assert(eigensym(intm,kappa,vp));

    /* Truncation of eigenvalues */
    kappa[0] = 2.0/9.0 * fabs(kappa[0])/info.hausd;
    kappa[0] = MS_MIN(kappa[0],isqhmin);
    kappa[0] = MS_MAX(kappa[0],isqhmax);

    kappa[1] = 2.0/9.0 * fabs(kappa[1])/info.hausd;
    kappa[1] = MS_MIN(kappa[1],isqhmin);
    kappa[1] = MS_MAX(kappa[1],isqhmax);

    /* Send back the metric to the canonical basis of tangent plane :
       diag(lambda) = {^t}vp * M * vp, M = vp * diag(lambda) * {^t}vp */
    intm[0] = kappa[0]*vp[0][0]*vp[0][0] + kappa[1]*vp[1][0]*vp[1][0];
    intm[1] = kappa[0]*vp[0][0]*vp[0][1] + kappa[1]*vp[1][0]*vp[1][1];
    intm[2] = kappa[0]*vp[0][1]*vp[0][1] + kappa[1]*vp[1][1]*vp[1][1];

    /* At this point, intm (with 0 in the z direction)  is the desired metric, except
       it is expressed in the rotated bc, that is intm = R * metric in bc * ^t R,
       so metric in bc = ^tR*intm*R */

    /* b0 and b1 serve now for nothing : let them be the lines of matrix intm*R  */
    b0[0] = intm[0]*r[0][0] + intm[1]*r[1][0] ;   b0[1] = intm[0]*r[0][1] + intm[1]*r[1][1] ;    b0[2] = intm[0]*r[0][2] + intm[1]*r[1][2] ;
    b1[0] = intm[1]*r[0][0] + intm[2]*r[1][0] ;   b1[1] = intm[1]*r[0][1] + intm[2]*r[1][1] ;    b1[2] = intm[1]*r[0][2] + intm[2]*r[1][2] ;
    //last line = 0.0;

    m[0] = r[0][0] * b0[0] + r[1][0] * b1[0];
    m[1] = r[0][0] * b0[1] + r[1][0] * b1[1];
    m[2] = r[0][0] * b0[2] + r[1][0] * b1[2];

    m[3] = r[0][1] * b0[1] + r[1][1] * b1[1];
    m[4] = r[0][1] * b0[2] + r[1][1] * b1[2];

    m[5] = r[0][2] * b0[2] + r[1][2] * b1[2];

    /* Security check : normal in the kernel */
    /*if((fabs(p0->m[0]*n[0] + p0->m[1]*n[1] + p0->m[2]*n[2] ) > 1.0e-4)){
      printf("VALEUR ETRANGE... %f \n",fabs(p0->m[0]*n[0] + p0->m[1]*n[1] + p0->m[2]*n[2] ));
      }
      if((fabs(p0->m[1]*n[0] + p0->m[3]*n[1] + p0->m[4]*n[2] ) > 1.0e-4)){
      printf("VALEUR ETRANGE... %f \n",fabs(p0->m[1]*n[0] + p0->m[3]*n[1] + p0->m[4]*n[2] ));
      }

      if((fabs(p0->m[2]*n[0] + p0->m[4]*n[1] + p0->m[5]*n[2] ) > 1.0e-4)){
      printf("VALEUR ETRANGE... %f \n",fabs(p0->m[2]*n[0] + p0->m[4]*n[1] + p0->m[5]*n[2]));
      } */

    /*if(ddb) {
      printf("La matrice %f %f %f\n",p0->m[0],p0->m[1],p0->m[2]);
      printf("            %f %f %f\n",p0->m[1],p0->m[3],p0->m[4]);
      printf("            %f %f %f\n",p0->m[2],p0->m[4],p0->m[5]);

      }*/
    return(1);
}

int defsiz_ani(pMesh mesh,pSol met) {
    pTria    pt;
    pPoint   ppt;
    double  *m,*n,mm[6],r[3][3],isqhmax;
    int      k;
    char     i,ismet;

    if ( abs(info.imprim) > 5 || info.ddebug )
        fprintf(stdout,"  ** Defining map\n");

    ismet = (met->m > 0);
    if ( !met->m ) {
        met->np    = mesh->np;
        met->npmax = mesh->npmax;
        met->m = calloc(6*(mesh->npmax+1)+1,sizeof(double));
        assert(met->m);
    }
    if ( info.hmax < 0.0 )  info.hmax = 0.5 * info.delta;

    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        ppt->flag = 0;
    }

    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) || pt->ref < 0 )  continue;

        for (i=0; i<3; i++) {
            ppt = &mesh->point[pt->v[i]];
            if ( ppt->flag || !MS_VOK(ppt) )  continue;
            if ( ismet )  memcpy(mm,&met->m[6*(pt->v[i])+1],6*sizeof(double));

            if ( MS_SIN(ppt->tag) ) {
                if ( !defmetsin(mesh,met,k,i) )  continue;
            }
            else if ( ppt->tag & MS_GEO ) {
                if ( !defmetrid(mesh,met,k,i))  continue;
            }
            else if ( (ppt->tag & MS_REF) && (!(ppt->tag & MS_GEO)) ) {
                if ( !defmetref(mesh,met,k,i) )  continue;
            }
            else if ( ppt->tag )  continue;
            else {
                if ( !defmetreg(mesh,met,k,i) )  continue;
            }
/* A FAIRE */
            if ( ismet )  intextmet(mesh,met,pt->v[i],mm);
            ppt->flag = 1;
        }
    }

    /* search for unintialized metric */
    isqhmax = 1.0 / (info.hmax*info.hmax);
    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MS_VOK(ppt) || ppt->flag == 1 )  continue;

        m = &met->m[6*(k)+1];
        memset(m,0,6*sizeof(double));
        if ( MS_SIN(ppt->tag) ) {
            m[0] = m[3] = m[5] = isqhmax;
        }
        else if ( ppt->tag & MS_GEO ) {
            m[0] = m[1] = m[2] = isqhmax;
        }
        else {
            n = ppt->tag & MS_REF ? &mesh->geom[ppt->ig].n1[0] : ppt->n;
            rotmatrix(n,r);
            m[0] = isqhmax*(r[0][0]*r[0][0]+r[1][0]*r[1][0]);
            m[1] = isqhmax*(r[0][0]*r[0][1]+r[1][0]*r[1][1]);
            m[2] = isqhmax*(r[0][0]*r[0][2]+r[1][0]*r[1][2]);
            m[3] = isqhmax*(r[0][1]*r[0][1]+r[1][1]*r[1][1]);
            m[4] = isqhmax*(r[0][1]*r[0][2]+r[1][1]*r[1][2]);
            m[5] = isqhmax*(r[0][2]*r[0][2]+r[1][2]*r[1][2]);
        }
        ppt->flag = 1;
    }

    return(1);
}

/* Enforces gradation of metric in one extremity of edge i in tria k with respect to the other,
   along the direction of the associated support curve
   Return -1 if no gradation is needed, else index of graded point */
static int grad2met(pMesh mesh, pSol met, int iel, int i){
    pTria    pt;
    pPoint   p1,p2;
    double   *mm1,*mm2,*nn1,*nn2,ps1,ps2,ux,uy,uz,m1[6],m2[6],n1[3],n2[3],nt[3];
    double   r1[3][3],r2[3][3],t1[3],t2[3],c[3],mtan1[3],mtan2[3],mr[6],l1,l2,l,dd;
    double   lambda[2],vp[2][2],alpha,beta,mu;
    int      np1,np2;
    char     i1,i2,ichg;

    pt = &mesh->tria[iel];

    i1 = inxt[i];
    i2 = iprv[i];
    np1 = pt->v[i1];
    np2 = pt->v[i2];

    p1 = &mesh->point[np1];
    p2 = &mesh->point[np2];

    ux = p2->c[0] - p1->c[0];
    uy = p2->c[1] - p1->c[1];
    uz = p2->c[2] - p1->c[2];

    mm1 = &met->m[6*(np1)+1];
    mm2 = &met->m[6*(np2)+1];

    if( !nortri(mesh,pt,nt) )
        return(-1);

    /* Recover normal and metric associated to p1 */
    if( MS_SIN(p1->tag) ){
        memcpy(n1,nt,3*sizeof(double));
        memcpy(m1,mm1,6*sizeof(double));
    }
    else if( MS_GEO & p1->tag ){
        nn1 = &mesh->geom[p1->ig].n1[0];
        nn2 = &mesh->geom[p1->ig].n2[0];
        ps1 = nt[0]*nn1[0] + nt[1]*nn1[1] + nt[2]*nn1[2];
        ps2 = nt[0]*nn2[0] + nt[1]*nn2[1] + nt[2]*nn2[2];

        if( fabs(ps1) < fabs(ps2))
            memcpy(n1,nn2,3*sizeof(double));
        else
            memcpy(n1,nn1,3*sizeof(double));

        if( !buildridmet(mesh,met,np1,ux,uy,uz,m1) )
            return(-1);
    }
    else if( MS_REF & p1->tag ){
        memcpy(n1,&(mesh->geom[p1->ig].n1[0]),3*sizeof(double));
        memcpy(m1,mm1,6*sizeof(double));
    }
    else{
        memcpy(n1,p1->n,3*sizeof(double));
        memcpy(m1,mm1,6*sizeof(double));
    }

    /* Recover normal and metric associated to p2 */
    if ( MS_SIN(p2->tag) ) {
        memcpy(n2,nt,3*sizeof(double));
        memcpy(m2,mm2,6*sizeof(double));
    }
    else if ( MS_GEO & p2->tag ) {
        nn1 = &mesh->geom[p2->ig].n1[0];
        nn2 = &mesh->geom[p2->ig].n2[0];
        ps1 = nt[0]*nn1[0] + nt[1]*nn1[1] + nt[2]*nn1[2];
        ps2 = nt[0]*nn2[0] + nt[1]*nn2[1] + nt[2]*nn2[2];

        if( fabs(ps1) < fabs(ps2))
            memcpy(n2,nn2,3*sizeof(double));
        else
            memcpy(n2,nn1,3*sizeof(double));

        if( !buildridmet(mesh,met,np2,ux,uy,uz,m2) )
            return(-1);
    }
    else if( MS_REF & p2->tag ){
        memcpy(n2,&(mesh->geom[p2->ig].n1[0]),3*sizeof(double));
        memcpy(m2,mm2,6*sizeof(double));
    }
    else{
        memcpy(n2,p2->n,3*sizeof(double));
        memcpy(m2,mm2,6*sizeof(double));
    }

    /* Rotation matrices mapping n1/n2 to e_3 */
    rotmatrix(n1,r1);
    rotmatrix(n2,r2);

    /* Geodesic length of support curve to edge i */
    ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
    t1[0] = ux - ps1*n1[0];
    t1[1] = uy - ps1*n1[1];
    t1[2] = uz - ps1*n1[2];

    ps2 = - (ux*n2[0] + uy*n2[1] + uz*n2[2]);
    t2[0] = -ux - ps1*n2[0];
    t2[1] = -uy - ps1*n2[1];
    t2[2] = -uz - ps1*n2[2];

    l1 = m1[0]*t1[0]*t1[0] + m1[3]*t1[1]*t1[1] + m1[5]*t1[2]*t1[2] \
        + 2.0 * ( m1[1]*t1[0]*t1[1] + m1[2]*t1[0]*t1[2] + m1[4]*t1[1]*t1[2] ) ;
    l2 = m2[0]*t2[0]*t2[0] + m2[3]*t2[1]*t2[1] + m2[5]*t2[2]*t2[2] \
        + 2.0 * ( m2[1]*t2[0]*t2[1] + m2[2]*t2[0]*t2[2] + m2[4]*t2[1]*t2[2] ) ;
    l = 0.5* ( sqrt(l1) + sqrt(l2) ) ;

    l = sqrt(ux*ux+uy*uy+uz*uz);

    /* Characteristic sizes in direction of support curve */
    rmtr(r1,m1,mr);
    mtan1[0] = mr[0];
    mtan1[1] = mr[1];
    mtan1[2] = mr[3];
    c[0] = r1[0][0]*ux + r1[0][1]*uy + r1[0][2]*uz;
    c[1] = r1[1][0]*ux + r1[1][1]*uy + r1[1][2]*uz;
    c[2] = r1[2][0]*ux + r1[2][1]*uy + r1[2][2]*uz;
    memcpy(t1,c,3*sizeof(double));
    dd = t1[0]*t1[0] + t1[1]*t1[1];
    if(dd < EPSD2)
        return(-1);

    dd = 1.0/sqrt(dd);
    t1[0] *= dd;
    t1[1] *= dd;
    ps1 = mtan1[0]*t1[0]*t1[0] + 2.0*mtan1[1]*t1[0]*t1[1] + mtan1[2]*t1[1]*t1[1];
    ps1 = sqrt(ps1);

    rmtr(r2,m2,mr);
    mtan2[0] = mr[0];
    mtan2[1] = mr[1];
    mtan2[2] = mr[3];
    c[0] = - ( r2[0][0]*ux + r2[0][1]*uy + r2[0][2]*uz );
    c[1] = - ( r2[1][0]*ux + r2[1][1]*uy + r2[1][2]*uz );
    c[2] = - ( r2[2][0]*ux + r2[2][1]*uy + r2[2][2]*uz );
    memcpy(t2,c,3*sizeof(double));

    dd = t2[0]*t2[0] + t2[1]*t2[1];
    if(dd < EPSD2)
        return(-1);

    dd = 1.0/sqrt(dd);
    t2[0] *= dd;
    t2[1] *= dd;
    ps2 = mtan2[0]*t2[0]*t2[0] + 2.0*mtan2[1]*t2[0]*t2[1] + mtan2[2]*t2[1]*t2[1];
    ps2 = sqrt(ps2);

    /* Metric in p1 has to be changed */
    if( ps2 > ps1 ){
        alpha = ps2 /(1.0+info.hgrad*l*ps2);
        if( ps1 >= alpha -EPS )
            return(-1);

        eigensym(mtan1,lambda,vp);
        c[0] = t1[0]*vp[0][0] + t1[1]*vp[0][1];
        c[1] = t1[0]*vp[1][0] + t1[1]*vp[1][1];

        if( fabs(c[0]) > fabs(c[1]) ){
            ichg  = 0;
            beta = (alpha*alpha - ps1*ps1)/(c[0]*c[0]);
            mu = lambda[0] + beta ;
            mtan1[0] = mu*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
            mtan1[1] = mu*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
            mtan1[2] = mu*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];
        }
        else{
            ichg = 1;
            beta = (alpha*alpha - ps1*ps1)/(c[1]*c[1]);
            mu = lambda[1] + beta;
            mtan1[0] = lambda[0]*vp[0][0]*vp[0][0] + mu*vp[1][0]*vp[1][0];
            mtan1[1] = lambda[0]*vp[0][0]*vp[0][1] + mu*vp[1][0]*vp[1][1];
            mtan1[2] = lambda[0]*vp[0][1]*vp[0][1] + mu*vp[1][1]*vp[1][1];
        }

        /* Metric update */
        if( MS_SIN(p1->tag) ){
            mm1[0] += 0.5*beta;
            mm1[3] += 0.5*beta;
            mm1[5] += 0.5*beta;
        }
        else if( p1->tag & MS_GEO ){
            c[0] = fabs(mm1[0]-lambda[ichg]);
            c[1] = fabs(mm1[1]-lambda[ichg]);
            c[2] = fabs(mm1[2]-lambda[ichg]);
            if( c[0] < c[1] ){
                if( c[0] < c[2] ){
                    mm1[0] += beta;
                }
                else{
                    mm1[2] += beta;
                }
            }
            else{
                if( c[1] < c[2] ){
                    mm1[1] += beta;
                }
                else{
                    mm1[2] += beta;
                }
            }
        }
        else{
            /* Reuse t1 and t2 */
            t1[0] = mtan1[0]*r1[0][0] + mtan1[1]*r1[1][0];  t1[1] = mtan1[0]*r1[0][1] + mtan1[1]*r1[1][1];  t1[2] = mtan1[0]*r1[0][2] + mtan1[1]*r1[1][2];
            t2[0] = mtan1[1]*r1[0][0] + mtan1[2]*r1[1][0];  t2[1] = mtan1[1]*r1[0][1] + mtan1[2]*r1[1][1];  t2[2] = mtan1[1]*r1[0][2] + mtan1[2]*r1[1][2];

            m1[0] = r1[0][0]*t1[0] + r1[1][0]*t2[0];
            m1[1] = r1[0][0]*t1[1] + r1[1][0]*t2[1];
            m1[2] = r1[0][0]*t1[2] + r1[1][0]*t2[2];
            m1[3] = r1[0][1]*t1[1] + r1[1][1]*t2[1];
            m1[4] = r1[0][1]*t1[2] + r1[1][1]*t2[2];
            m1[5] = r1[0][2]*t1[2] + r1[1][2]*t2[2];

            memcpy(mm1,m1,6*sizeof(double));
        }
        return(i1);
    }
    /* Metric in p2 has to be changed */
    else{
        alpha = ps1 /(1.0+info.hgrad*l*ps1);
        if( ps2 >= alpha - EPS)
            return(-1);

        eigensym(mtan2,lambda,vp);
        c[0] = t2[0]*vp[0][0] + t2[1]*vp[0][1];
        c[1] = t2[0]*vp[1][0] + t2[1]*vp[1][1];

        if( fabs(c[0]) > fabs(c[1]) ){
            ichg = 0;
            beta = (alpha*alpha - ps2*ps2)/(c[0]*c[0]);
            mu = lambda[0] + beta;
            mtan2[0] = mu*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
            mtan2[1] = mu*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
            mtan2[2] = mu*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];
        }
        else{
            ichg = 1;
            beta = (alpha*alpha - ps2*ps2)/(c[1]*c[1]);
            mu = lambda[1] + beta;
            mtan2[0] = lambda[0]*vp[0][0]*vp[0][0] + mu*vp[1][0]*vp[1][0];
            mtan2[1] = lambda[0]*vp[0][0]*vp[0][1] + mu*vp[1][0]*vp[1][1];
            mtan2[2] = lambda[0]*vp[0][1]*vp[0][1] + mu*vp[1][1]*vp[1][1];
        }

        /* Metric update */
        if( MS_SIN(p2->tag) ){
            mm2[0] += 0.5*beta;
            mm2[3] += 0.5*beta;
            mm2[5] += 0.5*beta;
        }
        else if( p2->tag & MS_GEO ){
            c[0] = fabs(mm2[0]-lambda[ichg]);
            c[1] = fabs(mm2[1]-lambda[ichg]);
            c[2] = fabs(mm2[2]-lambda[ichg]);
            if( c[0] < c[1] ){
                if( c[0] < c[2] ){
                    mm2[0] += beta;
                }
                else{
                    mm2[2] += beta;
                }
            }
            else{
                if( c[1] < c[2] ){
                    mm2[1] += beta;
                }
                else{
                    mm2[2] += beta;
                }
            }
        }
        else{
            /* Reuse t1 and t2 */
            t1[0] = mtan2[0]*r2[0][0] + mtan2[1]*r2[1][0];  t1[1] = mtan2[0]*r2[0][1] + mtan2[1]*r2[1][1];  t1[2] = mtan2[0]*r2[0][2] + mtan2[1]*r2[1][2];
            t2[0] = mtan2[1]*r2[0][0] + mtan2[2]*r2[1][0];  t2[1] = mtan2[1]*r2[0][1] + mtan2[2]*r2[1][1];  t2[2] = mtan2[1]*r2[0][2] + mtan2[2]*r2[1][2];

            m2[0] = r2[0][0]*t1[0] + r2[1][0]*t2[0];
            m2[1] = r2[0][0]*t1[1] + r2[1][0]*t2[1];
            m2[2] = r2[0][0]*t1[2] + r2[1][0]*t2[2];
            m2[3] = r2[0][1]*t1[1] + r2[1][1]*t2[1];
            m2[4] = r2[0][1]*t1[2] + r2[1][1]*t2[2];
            m2[5] = r2[0][2]*t1[2] + r2[1][2]*t2[2];

            memcpy(mm2,m2,6*sizeof(double));
        }

        return(i2);
    }
}

/* Enforces mesh gradation by truncating metric field */
int gradsiz_ani(pMesh mesh,pSol met) {
    pTria   pt;
    pPoint  p1,p2;
    double  *m,mv;
    int     k,it,nup,nu,maxit;
    char    i,ier,i1,i2;

    if ( abs(info.imprim) > 5 || info.ddebug )
        fprintf(stdout,"  ** Anisotropic mesh gradation\n");

    mesh->base = 0;
    for (k=1; k<=mesh->np; k++)
        mesh->point[k].flag = mesh->base;

    /* First step : make ridges iso */
    for (k=1; k<= mesh->np; k++) {
        p1 = &mesh->point[k];
        if ( !MS_VOK(p1) ) continue;
        if ( MS_SIN(p1->tag) ) continue;
        if ( !(p1->tag & MS_GEO) ) continue;

        m = &met->m[6*k+1];
        mv = MS_MAX(m[0],MS_MAX(m[1],m[2]));
        m[0] = mv;
        m[1] = mv;
        m[2] = mv;
    }

    /* Second step : standard gradation procedure */
    it = nup = 0;
    maxit = 100;
    do {
        mesh->base++;
        nu = 0;
        for (k=1; k<=mesh->nt; k++) {
            pt = &mesh->tria[k];
            if ( !MS_EOK(pt) )  continue;

            for (i=0; i<3; i++) {
                i1 = inxt[i];
                i2 = iprv[i];
                p1 = &mesh->point[pt->v[i1]];
                p2 = &mesh->point[pt->v[i2]];

                if ( p1->flag < mesh->base-1 && p2->flag < mesh->base-1 )  continue;
                ier = grad2met(mesh,met,k,i);
                if ( ier == i1 ) {
                    p1->flag = mesh->base;
                    nu++;
                }
                else if ( ier == i2 ) {
                    p2->flag = mesh->base;
                    nu++;
                }
            }
        }
        nup += nu;
    }
    while( ++it < maxit && nu > 0 );

    if ( abs(info.imprim) > 4 )  fprintf(stdout,"     gradation: %7d updated, %d iter.\n",nup,it);
    return(1);
}

