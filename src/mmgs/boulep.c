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
 * \file mmgs/boulep.c
 * \brief Functions for ball of points computation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

/* find all triangles sharing P, list[0] = start
   do not stop when crossing ridge */
int boulet(pMesh mesh,int start,int ip,int *list) {
    pTria    pt;
    pPoint   ppt;
    int     *adja,k,ilist,idp;
    char     i,i1,i2;

    pt = &mesh->tria[start];
    if ( !MS_EOK(pt) )  return(0);
    idp = pt->v[ip];
    ppt = &mesh->point[pt->v[ip]];
    if ( ppt->tag & MS_NOM )  return(0);
    ilist = 0;

    /* store neighbors */
    k = start;
    i = ip;
    do {
        if ( ilist > LMAX-2 )  return(-ilist);
        list[ilist] = 3*k + i;
        ++ilist;

        adja = &mesh->adja[3*(k-1)+1];
        i1 = inxt[i];
        k  = adja[i1] / 3;
        i  = adja[i1] % 3;
        i  = inxt[i];
    }
    while ( k && k != start );
    if ( k > 0 )  return(ilist);

    /* check if boundary hit */
    k = start;
    i = ip;
    do {
        adja = &mesh->adja[3*(k-1)+1];
        i2 = iprv[i];
        k  = adja[i2] / 3;
        if ( k == 0 )  break;
        i  = adja[i2] % 3;
        i  = iprv[i];

        if ( ilist > LMAX-2 )  return(-ilist);
        list[ilist] = 3*k + i;
        ilist++;
    }
    while ( k );

    return(ilist);
}

/* find all triangles sharing P, list[0] = start ; do not stop when crossing ridge ;
   check whether resulting configuration is manifold */
int boulechknm(pMesh mesh,int start,int ip,int *list) {
    pTria    pt;
    pPoint   ppt;
    int     *adja,k,ilist,idp,base,iel;
    char     i,i1,i2,ia,iq,voy;

    base = ++mesh->base;

    pt = &mesh->tria[start];
    ia = iprv[ip];
    iq = inxt[ip];
    if ( !MS_EOK(pt) )  return(0);
    idp = pt->v[ip];
    ppt = &mesh->point[pt->v[ip]];
    if ( ppt->tag & MS_NOM )  return(0);
    ilist = 0;

    /* store neighbors */
    k = start;
    i = ip;
    do {
        if ( ilist > LMAX-2 )  return(-ilist);
        list[ilist] = 3*k + i;
        ++ilist;

        pt = &mesh->tria[k];

        i1 = inxt[i];
        i2 = iprv[i];
        ppt = &mesh->point[pt->v[i1]];
        ppt->s = base;
        ppt = &mesh->point[pt->v[i2]];
        ppt->s = base;

        adja = &mesh->adja[3*(k-1)+1];
        k  = adja[i1] / 3;
        i  = adja[i1] % 3;
        i  = inxt[i];
    }
    while ( k && k != start );

    /* check if boundary hit */
    if ( k <= 0 ) {
        k = start;
        i = ip;
        do {
            adja = &mesh->adja[3*(k-1)+1];
            i1 = inxt[i];
            i2 = iprv[i];

            pt = &mesh->tria[k];
            ppt = &mesh->point[pt->v[i1]];
            ppt->s = base;
            ppt = &mesh->point[pt->v[i2]];
            ppt->s = base;

            k  = adja[i2] / 3;
            if ( k == 0 )  break;
            i  = adja[i2] % 3;
            i  = iprv[i];

            if ( ilist > LMAX-2 )  return(-ilist);
            list[ilist] = 3*k + i;
            ilist++;
        }
        while ( k );
    }

    pt = &mesh->tria[start];
    i1 = inxt[ip];
    i2 = iprv[ip];
    ppt = &mesh->point[pt->v[i1]];
    ppt->s = 0;
    ppt = &mesh->point[pt->v[i2]];
    ppt->s = 0;

    adja = &mesh->adja[3*(start-1)+1];
    iel = adja[ia] / 3;
    voy = adja[ia] % 3;

    if( iel ) {
        pt = &mesh->tria[iel];
        ppt = &mesh->point[pt->v[voy]];
        ppt->s = 0;
    }

    /* store neighbors */
    k = start;
    i = iq;
    do {
        pt = &mesh->tria[k];

        i1 = inxt[i];
        i2 = iprv[i];
        ppt = &mesh->point[pt->v[i1]];
        if ( ppt->s == base ) return(0);
        ppt = &mesh->point[pt->v[i2]];
        if ( ppt->s == base ) return(0);

        adja = &mesh->adja[3*(k-1)+1];
        k  = adja[i1] / 3;
        i  = adja[i1] % 3;
        i  = inxt[i];
    }
    while ( k && k != start );
    if( k > 0 ) return(ilist);

    /* check if boundary hit */
    k = start;
    i = iq;
    do {
        adja = &mesh->adja[3*(k-1)+1];
        i1 = inxt[i];
        i2 = iprv[i];

        pt = &mesh->tria[k];
        ppt = &mesh->point[pt->v[i1]];
        if ( ppt->s == base ) return(0);
        ppt = &mesh->point[pt->v[i2]];
        if ( ppt->s == base ) return(0);

        k  = adja[i2] / 3;
        if ( k == 0 )  break;
        i  = adja[i2] % 3;
        i  = iprv[i];
    }
    while ( k );

    return(ilist);
}

/* return average normal of triangles sharing P */
int boulen(pMesh mesh,int start,int ip,double *nn) {
    pTria    pt;
    double   n[3],dd;
    int     *adja,k,ier;
    char     i,i1,i2;

    pt = &mesh->tria[start];
    if ( !MS_EOK(pt) )  return(0);
    nn[0] = nn[1] = nn[2] = 0.0;

    /* store neighbors */
    k  = start;
    i  = ip;
    i1 = inxt[i];
    do {
        pt = &mesh->tria[k];
        nortri(mesh,pt,n);
        nn[0] += n[0];  nn[1] += n[1];  nn[2] += n[2];

        if ( pt->tag[i1] & MS_GEO ) {
            k = 0;
            break;
        }
        adja = &mesh->adja[3*(k-1)+1];
        k  = adja[i1] / 3;
        i2 = adja[i1] % 3;
        i1 = iprv[i2];
    }
    while ( k && k != start );

    if ( k == 0 ) {
        k  = start;
        i  = ip;
        i2 = iprv[i];
        pt = &mesh->tria[k];
        do {
            if ( pt->tag[i2] & MS_GEO )  break;

            adja = &mesh->adja[3*(k-1)+1];
            k  = adja[i2] / 3;
            if ( k == 0 )  break;
            i1 = adja[i2] % 3;
            i2 = inxt[i1];
            pt = &mesh->tria[k];

            ier = nortri(mesh,pt,n);
            nn[0] += n[0];  nn[1] += n[1];  nn[2] += n[2];
        }
        while ( k && k != start );
    }

    /* normalize */
    dd = nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2];
    if ( dd > EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        nn[0] *= dd;
        nn[1] *= dd;
        nn[2] *= dd;
        return(1);
    }

    return(0);
}


/* return all vertices connected to ip, list[0] = ip */
int boulep(pMesh mesh,int start,int ip,int *list) {
    pTria    pt;
    int     *adja,k,ilist;
    char     i,i1,i2;

    pt = &mesh->tria[start];
    if ( !MS_EOK(pt) )  return(0);
    list[0] = pt->v[ip];
    ilist   = 0;

    /* store neighbors */
    k  = start;
    i  = ip;
    i1 = inxt[i];
    i2 = iprv[i];
    do {
        if ( ilist > LMAX-2 )  return(-ilist);
        ilist++;
        list[ilist] = pt->v[i2];

        adja = &mesh->adja[3*(k-1)+1];
        k  = adja[i1] / 3;
        i2 = adja[i1] % 3;
        i1 = iprv[i2];
        pt = &mesh->tria[k];
    }
    while ( k && k != start );
    if ( k > 0 )  return(ilist);

    /* reverse loop */
    k  = start;
    i  = ip;
    pt = &mesh->tria[k];
    i1 = inxt[i];
    i2 = inxt[i1];
    do {
        if ( ilist > LMAX-2 )  return(-ilist);
        ilist++;
        list[ilist] = pt->v[i1];

        adja = &mesh->adja[3*(k-1)+1];
        k  = adja[i2] / 3;
        i1 = adja[i2] % 3;
        i2 = iprv[i1];
        pt = &mesh->tria[k];
    }
    while ( k > 0 );

    return(ilist);
}


/* return tangent to curve at ip */
int boulec(pMesh mesh,int start,int ip,double *tt) {
    pTria    pt;
    pPoint   p0,p1,p2;
    double   dd;
    int     *adja,k;
    char     i,i1,i2;

    pt = &mesh->tria[start];
    if ( !MS_EOK(pt) )  return(0);
    p0 = &mesh->point[pt->v[ip]];
    if ( !MS_EDG(p0->tag) )  return(0);

    /* check other triangle vertices */
    k  = start;
    i  = ip;
    i1 = inxt[i];
    i2 = iprv[i];
    p1 = p2 = 0;
    do {
        pt = &mesh->tria[k];
        if ( MS_EDG(pt->tag[i1]) ) {
            p1 = &mesh->point[pt->v[i2]];
            k  = 0;
            break;
        }
        adja = &mesh->adja[3*(k-1)+1];
        k  = adja[i1] / 3;
        i2 = adja[i1] % 3;
        i1 = iprv[i2];
    }
    while ( k && k != start );

    /* check if open boundary hit */
    if ( k == 0 ) {
        k  = start;
        i  = ip;
        i1 = inxt[i];
        i2 = iprv[i];
        do {
            pt = &mesh->tria[k];
            if ( MS_EDG(pt->tag[i2]) ) {
                p2 = &mesh->point[pt->v[i1]];
                break;
            }
            adja = &mesh->adja[3*(k-1)+1];
            k  = adja[i2] / 3;
            i1 = adja[i2] % 3;
            i2 = inxt[i1];
        }
        while ( k );
    }

    if ( !p1 || !p2 )
        return(0);
    else if ( p1 == p2 )
        p2 = p0;

    /* tangent approx */
    tt[0] = p2->c[0] - p1->c[0];
    tt[1] = p2->c[1] - p1->c[1];
    tt[2] = p2->c[2] - p1->c[2];
    dd = tt[0]*tt[0] + tt[1]*tt[1] + tt[2]*tt[2];
    if ( dd > EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        tt[0] *= dd;
        tt[1] *= dd;
        tt[2] *= dd;
    }

    return(1);
}


/* store edges and return number (nref+ngeo) incident to ip */
int bouler(pMesh mesh,int start,int ip,int *list,int *ng,int *nr) {
    pTria    pt;
    int     *adja,k,ns;
    char     i,i1,i2;

    pt  = &mesh->tria[start];
    if ( !MS_EOK(pt) )  return(0);

    /* check other triangle vertices */
    k  = start;
    i  = ip;
    *ng = *nr = ns = 0;
    do {
        i1 = inxt[i];
        if ( MS_EDG(pt->tag[i1]) ) {
            i2 = iprv[i];
            ns++;
            if ( pt->tag[i1] & MS_GEO )
                *ng = *ng + 1;
            else
                *nr = *nr + 1;
            list[ns] = pt->v[i2];
            if ( ns > LMAX-2 )  return(-ns);
        }
        adja = &mesh->adja[3*(k-1)+1];
        k  = adja[i1] / 3;
        i  = adja[i1] % 3;
        i  = inxt[i];
        pt = &mesh->tria[k];
    }
    while ( k && k != start );

    /* reverse loop */
    if ( k != start ) {
        k = start;
        i = ip;
        do {
            pt = &mesh->tria[k];
            i2 = iprv[i];
            if ( MS_EDG(pt->tag[i2]) ) {
                i1 = inxt[i];
                ns++;
                if ( pt->tag[i2] & MS_GEO )
                    *ng = *ng + 1;
                else
                    *nr = *nr + 1;
                list[ns] = pt->v[i1];
                if ( ns > LMAX-2 )  return(-ns);
            }
            adja = &mesh->adja[3*(k-1)+1];
            k  = adja[i2] / 3;
            i  = adja[i2] % 3;
            i  = iprv[i];
        }
        while ( k && k != start );
    }

    return(ns);
}

/* Computation of the two balls of a ridge point : list1 is associated to normal n1's side
   ip0, ip1 = indices of the 2 ending point of the ridge
   Both lists are returned enumerated in direct order  */
int bouletrid(pMesh mesh,int start,int ip,int *il1,int *l1,int *il2,int *l2,int *ip0,int *ip1) {
    pTria           pt;
    pPoint          ppt;
    int             idp,k,kold,*adja,iel,*ilist1,*ilist2,*list1,*list2,aux;
    unsigned char   i,iold,i1,i2,ipn;
    double          *n1,*n2,nt[3],ps1,ps2;

    pt = &mesh->tria[start];
    if ( !MS_EOK(pt) )  return(0);

    idp = pt->v[ip];
    ppt = &mesh->point[idp];
    assert( ppt->tag & MS_GEO );

    /* set pointers: first manifold is on side of triangle */
    if ( !nortri(mesh,pt,nt) )  return(0);

    n1 = &(mesh->geom[ppt->ig].n1[0]);
    n2 = &(mesh->geom[ppt->ig].n2[0]);
    ps1 = n1[0]*nt[0] + n1[1]*nt[1] + n1[2]*nt[2];
    ps2 = n2[0]*nt[0] + n2[1]*nt[1] + n2[2]*nt[2];

    if ( fabs(ps1) < fabs(ps2) ) {
        list1  = l2;
        list2  = l1;
        ilist1 = il2;
        ilist2 = il1;
    }
    else {
        list1  = l1;
        list2  = l2;
        ilist1 = il1;
        ilist2 = il2;
    }
    *ilist1 = 0;

    /* First ball, first side (via i1)*/
    k = start;
    i = ip;
    do {
        pt   = &mesh->tria[k];
        adja = &mesh->adja[3*(k-1)+1];
        i1 = inxt[i];
        i2 = iprv[i];
        kold = k;
        iold = i;
        k = adja[i1] / 3;
        i = adja[i1] % 3;
        i = inxt[i];
    }
    while ( k && k != start && !(pt->tag[i1] & MS_GEO) );
    *ip0 = pt->v[i2];

    /* Store the needed elements to start back in the new area,
       and complete first ball, second side (via i2) */
    iel = k;
    ipn = i;
    k  = kold;
    i  = iold;
    do {
        pt   = &mesh->tria[k];
        adja = &mesh->adja[3*(k-1)+1];
        if ( (*ilist1) > LMAX-2 )  return(-(*ilist1));
        list1[(*ilist1)] = 3*k+i;
        (*ilist1)++;
        i1 = inxt[i];
        i2 = iprv[i];
        k = adja[i2] / 3;
        i = adja[i2] % 3;
        i = iprv[i];
    }
    while ( k && !(MS_GEO & pt->tag[i2]) );
    *ip1 = pt->v[i1];

    /* Invert the order in list1, so that it is enumerated in DIRECT order */
    for (k=0; k<(*ilist1) / 2;k++) {
        aux = list1[k];
        list1[k] = list1[*ilist1-1-k];
        list1[*ilist1-1-k] = aux;
    }
    /* At this point, either something has been stored in iel or an open boundary has been hit */
    *ilist2 = 0;
    if ( !iel )  return(1);

    /* Else, start back from the hit boundary, until another boundary is hit */
    k  = iel;
    i  = ipn;
    do {
        pt   = &mesh->tria[k];
        adja = &mesh->adja[3*(k-1)+1];
        if ( *ilist2 > LMAX-2 )  return(-(*ilist2));
        list2[*ilist2] = 3*k+i;
        (*ilist2)++;
        i1 = inxt[i];
        k = adja[i1] / 3;
        i = adja[i1] % 3;
        i = inxt[i];
    }
    while ( k && !(MS_GEO & pt->tag[i1]) );

    if ( !(MS_GEO & pt->tag[i1]) )
        return(0);
    else
        return(1);
}
