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
 * \file mmgs/swap.c
 * \brief Functions for swapping process.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmgs.h"

extern Info  info;

/* Check whether edge i of triangle k should be swapped for geometric approximation purposes */
int chkswp(pMesh mesh,pSol met,int k,int i,char typchk) {
    pTria    pt,pt0,pt1;
    pPoint   p[3],q;
    double   np[3][3],nq[3],*nr1,*nr2,nt[3],ps,ps2,*n1,*n2,dd,c1[3],c2[3];
    double   cosn1,cosn2,calnat,calchg,cal1,cal2,cosnat,coschg,ux,uy,uz,ll,loni,lona;
    int     *adja,j,kk,ip0,ip1,ip2,iq;
    char     ii,i1,i2,jj;

    pt0 = &mesh->tria[0];
    pt  = &mesh->tria[k];
    i1 = inxt[i];
    i2 = iprv[i];
    if ( MS_EDG(pt->tag[i]) || MS_SIN(pt->tag[i]) )  return(0);
    else if ( MS_SIN(pt->tag[i1]) )  return(0);

    ip0  = pt->v[i];
    ip1  = pt->v[i1];
    ip2  = pt->v[i2];
    p[0] = &mesh->point[ip0];
    p[1] = &mesh->point[ip1];
    p[2] = &mesh->point[ip2];

    adja = &mesh->adja[3*(k-1)+1];
    if ( !adja[i] )  return(0);

    kk = adja[i] / 3;
    ii = adja[i] % 3;
    jj = inxt[ii];
    pt1 = &mesh->tria[kk];
    if ( MS_SIN(pt1->tag[jj]) )  return(0);

    iq = pt1->v[ii];
    q  = &mesh->point[iq];

    /* check length */
    if ( typchk == 2 && met->m ) {
        loni = lenedg(mesh,met,ip1,ip2,0);
        lona = lenedg(mesh,met,ip0,iq,0);
        if ( loni > 1.0 )  loni = MS_MIN(1.0 / loni,LSHRT);
        if ( lona > 1.0 )  lona = 1.0 / lona;
        if ( lona < loni )  return(0);
    }

    /* check non convexity */
    norpts(p[0],p[1],q,c1);
    norpts(p[0],q,p[2],c2);
    ps = c1[0]*c2[0] + c1[1]*c2[1] + c1[2]*c2[2];
    if ( ps < ANGEDG )   return(0);

    /* normal recovery at points p[0],p[1],p[2],q */
    for (j=0; j<3; j++) {
        if ( MS_SIN(p[j]->tag) ) {
            nortri(mesh,pt,np[j]);
        }
        else if ( MS_EDG(p[j]->tag) ) {
            nortri(mesh,pt,nt);
            nr1  = &mesh->geom[p[j]->ig].n1[0];
            nr2  = &mesh->geom[p[j]->ig].n2[0];
            ps  = nr1[0]*nt[0] + nr1[1]*nt[1] + nr1[2]*nt[2];
            ps2 = nr2[0]*nt[0] + nr2[1]*nt[1] + nr2[2]*nt[2];
            if ( fabs(ps) > fabs(ps2) )
                memcpy(&np[j],nr1,3*sizeof(double));
            else
                memcpy(&np[j],nr2,3*sizeof(double));
        }
        else
            memcpy(&np[j],p[j]->n,3*sizeof(double));
    }

    if ( MS_SIN(q->tag) ) {
        nortri(mesh,pt,nq);
    }
    else if ( MS_EDG(q->tag) ) {
        nortri(mesh,pt,nt);
        nr1  = &mesh->geom[q->ig].n1[0];
        nr2  = &mesh->geom[q->ig].n2[0];
        ps  = nr1[0]*nt[0] + nr1[1]*nt[1] + nr1[2]*nt[2];
        ps2 = nr2[0]*nt[0] + nr2[1]*nt[1] + nr2[2]*nt[2];
        if ( fabs(ps) > fabs(ps2) )
            memcpy(&nq,nr1,3*sizeof(double));
        else
            memcpy(&nq,nr2,3*sizeof(double));
    }
    else
        memcpy(&nq,q->n,3*sizeof(double));

    /* Estimate of the Hausdorff distance between approximation and underlying surface
       when using the 'natural' edge [i1,i2] */
    ux = p[2]->c[0] - p[1]->c[0];
    uy = p[2]->c[1] - p[1]->c[1];
    uz = p[2]->c[2] - p[1]->c[2];

    ll = ux*ux + uy*uy + uz*uz;
    if ( ll < EPS )  return(0); /* no change for short edge */

    n1 = np[1];
    n2 = np[2];

    ps = ux*n1[0] + uy*n1[1] + uz*n1[2];
    c1[0] = (2.0*p[1]->c[0] + p[2]->c[0] - ps*n1[0]) / 3.0 - p[1]->c[0];
    c1[1] = (2.0*p[1]->c[1] + p[2]->c[1] - ps*n1[1]) / 3.0 - p[1]->c[1];
    c1[2] = (2.0*p[1]->c[2] + p[2]->c[2] - ps*n1[2]) / 3.0 - p[1]->c[2];

    ps = -(ux*n2[0] + uy*n2[1] + uz*n2[2]);
    c2[0] = (2.0*p[2]->c[0] + p[1]->c[0] - ps*n2[0]) / 3.0 - p[2]->c[0];
    c2[1] = (2.0*p[2]->c[1] + p[1]->c[1] - ps*n2[1]) / 3.0 - p[2]->c[1];
    c2[2] = (2.0*p[2]->c[2] + p[1]->c[2] - ps*n2[2]) / 3.0 - p[2]->c[2];

    /* squared cosines */
    ps = c1[0]*ux + c1[1]*uy + c1[2]*uz;
    ps *= ps;
    dd = c1[0]*c1[0] + c1[1]*c1[1] + c1[2]*c1[2];
    cosn1  =  ps / (dd*ll);
    cosn1 *= (1.0-cosn1);
    cosn1 *= (0.25*ll);

    ps = -c2[0]*ux - c2[1]*uy - c2[2]*uz;
    ps *= ps;
    dd = c2[0]*c2[0]+c2[1]*c2[1]+c2[2]*c2[2];
    cosn2  =  ps / (dd*ll);
    cosn2 *= (1.0-cosn2);
    cosn2 *= (0.25*ll);

    cosnat = MS_MAX(fabs(cosn1),fabs(cosn2));
    cosnat = cosnat < EPS ? 0.0 : cosnat;

    /* Estimate of the Hausdorff distance between approximation and underlying surface
       when using the 'swapped' edge [i0,q] */
    ux = q->c[0] - p[0]->c[0];
    uy = q->c[1] - p[0]->c[1];
    uz = q->c[2] - p[0]->c[2];

    ll = ux*ux + uy*uy + uz*uz;
    if ( ll < EPS )  return(0);

    n1 = np[0];
    n2 = nq;

    ps = ux*n1[0] + uy*n1[1] + uz*n1[2];
    c1[0] = (2.0*p[0]->c[0] + q->c[0] - ps*n1[0]) / 3.0 - p[0]->c[0];
    c1[1] = (2.0*p[0]->c[1] + q->c[1] - ps*n1[1]) / 3.0 - p[0]->c[1];
    c1[2] = (2.0*p[0]->c[2] + q->c[2] - ps*n1[2]) / 3.0 - p[0]->c[2];

    ps = -(ux*n2[0] + uy*n2[1] + uz*n2[2]);
    c2[0] = (2.0*q->c[0] + p[0]->c[0] - ps*n2[0]) / 3.0 - q->c[0];
    c2[1] = (2.0*q->c[1] + p[0]->c[1] - ps*n2[1]) / 3.0 - q->c[1];
    c2[2] = (2.0*q->c[2] + p[0]->c[2] - ps*n2[2]) / 3.0 - q->c[2];

    /* squared cosines */
    ps = c1[0]*ux + c1[1]*uy + c1[2]*uz;
    ps *= ps;
    dd = c1[0]*c1[0] + c1[1]*c1[1] + c1[2]*c1[2];
    cosn1  =  ps / (dd*ll);
    cosn1 *= (1.0-cosn1);
    cosn1 *= (0.25*ll);

    ps = -c2[0]*ux - c2[1]*uy - c2[2]*uz;
    ps *= ps;
    dd = c2[0]*c2[0]+c2[1]*c2[1]+c2[2]*c2[2];
    cosn2  =  ps / (dd*ll);
    cosn2 *= (1.0-cosn2);
    cosn2 *= (0.25*ll);

    coschg = MS_MAX(fabs(cosn1),fabs(cosn2));
    coschg = coschg < EPS ? 0.0 : coschg;

    /* swap if Hausdorff contribution of the swapped edge is less than existing one */
    if ( coschg > info.hausd*info.hausd )  return(0);
    else if ( coschg < info.hausd*info.hausd && cosnat > info.hausd*info.hausd )  return(1);

    if ( typchk == 2 && met->m ) {
        pt0->v[0]= ip0;  pt0->v[1]= ip1;  pt0->v[2]= ip2;
        cal1 = calelt(mesh,met,0);
        pt0->v[0]= ip1;  pt0->v[1]= iq;   pt0->v[2]= ip2;
        cal2 = calelt(mesh,met,0);
        calnat = MS_MIN(cal1,cal2);
        pt0->v[0]= ip0;  pt0->v[1]= ip1;  pt0->v[2]= iq;
        cal1 = calelt(mesh,met,0);
        pt0->v[0]= ip0;  pt0->v[1]= iq;   pt0->v[2]= ip2;
        cal2 = calelt(mesh,met,0);
        calchg = MS_MIN(cal1,cal2);
    }
    else {
        pt0->v[0]= ip0;  pt0->v[1]= ip1;  pt0->v[2]= ip2;
        cal1 = calelt_iso(mesh,met,0);
        pt0->v[0]= ip1;  pt0->v[1]= iq;   pt0->v[2]= ip2;
        cal2 = calelt_iso(mesh,met,0);
        calnat = MS_MIN(cal1,cal2);
        pt0->v[0]= ip0;  pt0->v[1]= ip1;  pt0->v[2]= iq;
        cal1 = calelt_iso(mesh,met,0);
        pt0->v[0]= ip0;  pt0->v[1]= iq;   pt0->v[2]= ip2;
        cal2 = calelt_iso(mesh,met,0);
        calchg = MS_MIN(cal1,cal2);
    }
    return(calchg > 1.01 * calnat);
}

int swapar(pMesh mesh,int k,int i) {
    pTria    pt,pt1;
    int     *adja,adj,k11,k21;
    char     i1,i2,j,jj,j2,v11,v21;

    pt   = &mesh->tria[k];
    if ( MS_EDG(pt->tag[i]) || MS_SIN(pt->tag[i]) )  return(0);

    adja = &mesh->adja[3*(k-1)+1];
    assert(adja[i]);

    adj = adja[i] / 3;
    j   = adja[i] % 3;
    pt1 = &mesh->tria[adj];

    /* simulation */
    i1 = inxt[i];
    i2 = iprv[i];

    /* update structure */
    k11 = adja[i1] / 3;
    v11 = adja[i1] % 3;
    if ( k11 < 1 )  return(0);
    adja = &mesh->adja[3*(adj-1)+1];
    jj  = inxt[j];
    j2  = iprv[j];
    k21 = adja[jj] / 3;
    v21 = adja[jj] % 3;
    if ( k21 < 1 )  return(0);

    pt->v[i2]  = pt1->v[j];
    pt1->v[j2] = pt->v[i];

    /* update info */
    pt->tag[i] = pt1->tag[jj];
    pt->edg[i] = pt1->edg[jj];
    pt->base   = mesh->base;
    pt1->tag[j] = pt->tag[i1];
    pt1->edg[j] = pt->edg[i1];
    pt->tag[i1] = 0;
    pt->edg[i1] = 0;
    pt1->tag[jj] = 0;
    pt1->edg[jj] = 0;
    pt1->base    = mesh->base;

    /* update adjacent */
    mesh->adja[3*(k-1)+1+i]     = 3*k21+v21;
    mesh->adja[3*(k21-1)+1+v21] = 3*k+i;
    mesh->adja[3*(k-1)+1+i1]    = 3*adj+jj;
    mesh->adja[3*(adj-1)+1+jj]  = 3*k+i1;
    mesh->adja[3*(k11-1)+1+v11] = 3*adj+j;
    mesh->adja[3*(adj-1)+1+j]   = 3*k11+v11;

    return(1);
}


/* flip edge i of tria k */
int litswp(pMesh mesh,int k,char i,double kali) {
    pTria    pt,pt0,pt1;
    pPoint   a,b,c,d;
    double   kalf,kalt,ps,n1[3],n2[3];
    int     *adja,ia,ib,ic,id,kk;
    char     ii,i1,i2;

    pt0 = &mesh->tria[0];
    pt  = &mesh->tria[k];
    if ( !MS_EOK(pt) || MS_EDG(pt->tag[i]) )  return(0);

    i1 = inxt[i];
    i2 = iprv[i];
    ia = pt->v[i];
    ib = pt->v[i1];
    ic = pt->v[i2];
    a  = &mesh->point[ia];
    b  = &mesh->point[ib];
    c  = &mesh->point[ic];

    adja = &mesh->adja[3*(k-1)+1];
    kk  = adja[i] / 3;
    ii  = adja[i] % 3;
    pt1 = &mesh->tria[kk];
    if ( MS_SIN(pt1->tag[ii]) )  return(0);
    id = pt1->v[ii];
    d  = &mesh->point[id];

    /* check non convexity */
    norpts(a,b,d,n1);
    norpts(a,d,c,n2);
    ps = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
    if ( ps < ANGEDG )  return(0);

    /* check quality */
    pt0->v[0] = id;  pt0->v[1] = ic;  pt0->v[2] = ib;
    kalt = calelt(mesh,NULL,0);
    kali = MS_MIN(kali,kalt);
    pt0->v[0] = ia;  pt0->v[1] = id;  pt0->v[2] = ic;
    kalt = calelt(mesh,NULL,0);
    pt0->v[0] = ia;  pt0->v[1] = ib;  pt0->v[2] = id;
    kalf = calelt(mesh,NULL,0);
    kalf = MS_MIN(kalf,kalt);
    if ( kalf > 1.02 * kali ) {
        swapar(mesh,k,i);
        return(1);
    }
    return(0);
}


/* attempt to swap any edge below quality value
   list goes from 0 to ilist-1 */
int swpedg(pMesh mesh,pSol met,int *list,int ilist,char typchk) {
    int      k,ns,iel;
    char     i,i1;

    k  = 0;
    ns = 0;
    do {
        iel = list[k] / 3;
        i   = list[k] % 3;
        i1  = inxt[i];
        if ( chkswp(mesh,met,iel,i1,typchk) ) {
            ns += swapar(mesh,iel,i1);
            k++;
        }
        k++;
    }
    while ( k < ilist );

    return(ns);
}
