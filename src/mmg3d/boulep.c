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
 * \file mmg3d/boulep.c
 * \brief Functions for ball of points computation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmg3d.h"

extern MMG5_Info  info;

/** return average normal of triangles sharing P without crossing ridge */
int _MMG5_boulen(MMG5_pMesh mesh,int start,int ip,double *nn) {
    MMG5_pTria    pt;
    double   n[3],dd;
    int     *adja,k;
    char     i,i1,i2;

    pt = &mesh->tria[start];
    if ( !MG_EOK(pt) )  return(0);
    nn[0] = nn[1] = nn[2] = 0.0;

    /* store neighbors */
    k  = start;
    i  = ip;
    i1 = _MMG5_inxt2[i];
    do {
        pt = &mesh->tria[k];
        _MMG5_nortri(mesh,pt,n);
        nn[0] += n[0];  nn[1] += n[1];  nn[2] += n[2];

        if ( pt->tag[i1] & MG_GEO ) {
            k = 0;
            break;
        }
        adja = &mesh->adjt[3*(k-1)+1];
        k  = adja[i1] / 3;
        i2 = adja[i1] % 3;
        i1 = _MMG5_iprv2[i2];
    }
    while ( k && k != start );

    if ( k == 0 ) {
        k  = start;
        i  = ip;
        i2 = _MMG5_iprv2[i];
        pt = &mesh->tria[k];
        do {
            if ( pt->tag[i2] & MG_GEO )  break;

            adja = &mesh->adjt[3*(k-1)+1];
            k  = adja[i2] / 3;
            if ( k == 0 )  break;
            i1 = adja[i2] % 3;
            i2 = _MMG5_inxt2[i1];
            pt = &mesh->tria[k];

            _MMG5_nortri(mesh,pt,n);

            nn[0] += n[0];  nn[1] += n[1];  nn[2] += n[2];
        }
        while ( k && k != start );
    }

    /* normalize */
    dd = nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2];
    if ( dd > _MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        nn[0] *= dd;
        nn[1] *= dd;
        nn[2] *= dd;
        return(1);
    }

    return(0);
}


/** return tangent to curve at ip */
int _MMG5_boulec(MMG5_pMesh mesh,int start,int ip,double *tt) {
    MMG5_pTria    pt;
    MMG5_pPoint   p0,p1,p2;
    double   dd;
    int     *adja,k;
    char     i,i1,i2;

    pt = &mesh->tria[start];
    if ( !MG_EOK(pt) )       return(0);
    p0 = &mesh->point[pt->v[ip]];
    if ( !MG_EDG(p0->tag) )  return(0);

    /* check other triangle vertices */
    k  = start;
    i  = ip;
    i1 = _MMG5_inxt2[i];
    i2 = _MMG5_iprv2[i];
    p1 = p2 = 0;
    do {
        pt = &mesh->tria[k];
        if ( MG_EDG(pt->tag[i1]) ) {
            p1 = &mesh->point[pt->v[i2]];
            k  = 0;
            break;
        }
        adja = &mesh->adjt[3*(k-1)+1];
        k  = adja[i1] / 3;
        i2 = adja[i1] % 3;
        i1 = _MMG5_iprv2[i2];
    }
    while ( k && k != start );

    /* check if open boundary hit */
    if ( k == 0 ) {
        k  = start;
        i  = ip;
        i1 = _MMG5_inxt2[i];
        i2 = _MMG5_iprv2[i];
        do {
            pt = &mesh->tria[k];
            if ( MG_EDG(pt->tag[i2]) ) {
                p2 = &mesh->point[pt->v[i1]];
                break;
            }
            adja = &mesh->adjt[3*(k-1)+1];
            k  = adja[i2] / 3;
            i1 = adja[i2] % 3;
            i2 = _MMG5_inxt2[i1];
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
    if ( dd > _MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        tt[0] *= dd;
        tt[1] *= dd;
        tt[2] *= dd;
    }

    return(1);
}


/** store edges and return number (ref+geo) incident to ip */
int _MMG5_bouler(MMG5_pMesh mesh,int start,int ip,int *list,int *ng,int *nr) {
    MMG5_pTria    pt;
    int     *adja,k,ns;
    char     i,i1,i2;

    pt  = &mesh->tria[start];
    if ( !MG_EOK(pt) )  return(0);
    /* check other triangle vertices */
    k  = start;
    i  = ip;
    *ng = *nr = ns = 0;
    do {
        i1 = _MMG5_inxt2[i];
        if ( MG_EDG(pt->tag[i1])) {
            i2 = _MMG5_iprv2[i];
            if ( pt->tag[i1] & MG_GEO )
                *ng = *ng + 1;
            else if ( pt->tag[i1] & MG_REF )
                *nr = *nr + 1;
            ns++;
            list[ns] = pt->v[i2];
            if ( ns > _MMG5_LMAX-2 )  return(-ns);
        }
        adja = &mesh->adjt[3*(k-1)+1];
        k  = adja[i1] / 3;
        i  = adja[i1] % 3;
        i  = _MMG5_inxt2[i];
        pt = &mesh->tria[k];
    }
    while ( k && k != start );

    /* reverse loop */
    if ( k != start ) {
        k = start;
        i = ip;
        do {
            pt = &mesh->tria[k];
            i2 = _MMG5_iprv2[i];
            if ( MG_EDG(pt->tag[i2]) ) {
                i1 = _MMG5_inxt2[i];
                if ( pt->tag[i2] & MG_GEO )
                    *ng = *ng + 1;
                else if ( pt->tag[i1] & MG_REF )
                    *nr = *nr + 1;
                ns++;
                list[ns] = pt->v[i1];
                if ( ns > _MMG5_LMAX-2 )  return(-ns);
            }
            adja = &mesh->adjt[3*(k-1)+1];
            k = adja[i2] / 3;
            i = adja[i2] % 3;
            i = _MMG5_iprv2[i];
        }
        while ( k && k != start );
    }
    return(ns);
}

/** Return volumic ball (i.e. filled with tetrahedra) of point ip in tetra start.
    Results are stored under the form 4*kel + jel , kel = number of the tetra, jel = local
    index of p within kel */
int _MMG5_boulevolp (MMG5_pMesh mesh, int start, int ip, int * list){
    MMG5_pTetra  pt,pt1;
    int    *adja,nump,ilist,base,cur,k,k1;
    char    j,l,i;

    base = ++mesh->base;
    pt   = &mesh->tetra[start];
    nump = pt->v[ip];
    ilist = 0;

    /* Store initial tetrahedron */
    pt->flag = base;
    list[ilist] = 4*start + ip;
    ilist++;

    /* Explore list and travel by adjacency through elements sharing p */
    cur = 0;
    while ( cur < ilist ) {
        k = list[cur] / 4;
        i = list[cur] % 4; // index of point p in tetra k
        adja = &mesh->adja[4*(k-1)+1];

        for (l=0; l<3; l++) {
            i  = _MMG5_inxt3[i];
            k1 = adja[i] / 4;
            if ( !k1 )  continue;
            pt1 = &mesh->tetra[k1];
            if ( pt1->flag == base )  continue;
            pt1->flag = base;
            for (j=0; j<4; j++)
                if ( pt1->v[j] == nump )  break;
            assert(j<4);
            /* overflow */
            if ( ilist > _MMG5_LMAX-3 )  return(0);
            list[ilist] = 4*k1+j;
            ilist++;
        }
        cur++;
    }
    return(ilist);
}

/** Define normal and tangent vectors at a non manifold point (ip in start, supported by
    face iface), enumerating its (outer)surfacic ball ; return sng = whether point is singular
    or not */
int _MMG5_boulenm(MMG5_pMesh mesh,int start,int ip,int iface,
                  double n[3],double t[3]) {
    MMG5_pTetra   pt;
    MMG5_pPoint   p0,p1,ppt;
    double   dd,nt[3],l0,l1;
    int      base,nump,nr,nnm,k,piv,na,nb,adj,nvstart,fstart,aux,ip0,ip1;
    int     *adja;
    char     iopp,ipiv,indb,inda,i,ipa,ipb,isface,tag;
    char     indedg[4][4] = { {-1,0,1,2}, {0,-1,3,4}, {1,3,-1,5}, {2,4,5,-1} };

    base = ++mesh->base;
    nr  = nnm = 0;
    ip0 = ip1 = 0;

    memset(n,0.0,3*sizeof(double));
    memset(t,0.0,3*sizeof(double));

    pt   = &mesh->tetra[start];
    nump = pt->v[ip];
    k    = start;

    na   = pt->v[ip];
    nb   = pt->v[_MMG5_idir[iface][_MMG5_inxt2[_MMG5_idirinv[iface][ip]]]];
    piv  = pt->v[_MMG5_idir[iface][_MMG5_iprv2[_MMG5_idirinv[iface][ip]]]];

    iopp   = iface;
    fstart = 4*k+iopp;
    do {
        /* computation of normal and tangent at nump */
        if ( _MMG5_norface(mesh,k,iopp,nt) ) {
            n[0] += nt[0];
            n[1] += nt[1];
            n[2] += nt[2];
        }

        if ( pt->xt ) {
            for ( inda=0; inda<4; inda++ ){
                if ( pt->v[inda]==na ) break;
            }
            for ( indb=0; indb<4; indb++ ){
                if ( pt->v[indb]==nb ) break;
            }
            assert( (inda < 4) && (indb < 4));
            tag = mesh->xtetra[pt->xt].tag[indedg[inda][indb]];
        }

        else  tag = 0;

        if ( MG_EDG(tag) && !(tag & MG_NOM) )
            nr++;
        else if ( tag & MG_NOM ) {
            nnm++;
            if ( !ip0 )
                ip0 = nb;
            else
                ip1 = nb;
        }

        /* A boundary face has been hit : change travel edge */
        aux     = nb;
        nb      = piv;
        piv     = aux;
        nvstart = k;
        adj     = k;

        /* Now unfold shell of edge (na,nb) starting from k (included) */
        do {
            k = adj;
            pt = &mesh->tetra[k];
            adja = &mesh->adja[4*(k-1)+1];
            if ( pt->flag != base ) {
                for (i=0; i<4; i++)
                    if ( pt->v[i] == nump )  break;
                assert(i<4);
                pt->flag = base;
            }

            /* identification of edge number in tetra k */
            for (i=0; i<6; i++) {
                ipa = _MMG5_iare[i][0];
                ipb = _MMG5_iare[i][1];
                if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                     (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
            }
            assert(i<6);

            /* set sense of travel */
            if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
                adj = adja[ _MMG5_ifar[i][0] ] / 4;
                ipiv = _MMG5_ifar[i][1];
                iopp = _MMG5_ifar[i][0];
                piv = pt->v[ipiv];
            }
            else {
                adj = adja[ _MMG5_ifar[i][1] ] / 4;
                ipiv = _MMG5_ifar[i][0];
                iopp = _MMG5_ifar[i][1];
                piv = pt->v[ipiv];
            }
            isface = (adja[iopp] == 0);
        }
        while ( adj && (adj != nvstart) && !isface );
    }
    while ( 4*k+iopp != fstart );

    if ( (nr > 0 && nnm > 0) || nnm != 2 )  return(0);

    dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
    if ( dd > _MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        n[0] *= dd;
        n[1] *= dd;
        n[2] *= dd;
    }
    assert( ip0 && ip1 );
    if ( ip0 == ip1 )  return(0);

    p0 = &mesh->point[ip0];
    p1 = &mesh->point[ip1];
    ppt = &mesh->point[nump];

    l0 = (ppt->c[0] - p0->c[0])*(ppt->c[0] - p0->c[0]) \
        + (ppt->c[1] - p0->c[1])*(ppt->c[1] - p0->c[1]) + (ppt->c[2] - p0->c[2])*(ppt->c[2] - p0->c[2]);
    l1 = (ppt->c[0] - p1->c[0])*(ppt->c[0] - p1->c[0]) \
        + (ppt->c[1] - p1->c[1])*(ppt->c[1] - p1->c[1]) + (ppt->c[2] - p1->c[2])*(ppt->c[2] - p1->c[2]);
    l0 = sqrt(l0);
    l1 = sqrt(l1);

    if ( (l0 < _MMG5_EPSD2) || (l1 < _MMG5_EPSD2) ) {
        t[0] = p1->c[0] - p0->c[0];
        t[1] = p1->c[1] - p0->c[1];
        t[2] = p1->c[2] - p0->c[2];
    }
    else if ( l0 < l1 ) {
        dd = l0 / l1;
        t[0] = dd*(p1->c[0] - ppt->c[0]) + ppt->c[0] - p0->c[0];
        t[1] = dd*(p1->c[1] - ppt->c[1]) + ppt->c[1] - p0->c[1];
        t[2] = dd*(p1->c[2] - ppt->c[2]) + ppt->c[2] - p0->c[2];
    }
    else {
        dd = l1 / l0;
        t[0] = dd*(p0->c[0] - ppt->c[0]) + ppt->c[0] - p1->c[0];
        t[1] = dd*(p0->c[1] - ppt->c[1]) + ppt->c[1] - p1->c[1];
        t[2] = dd*(p0->c[2] - ppt->c[2]) + ppt->c[2] - p1->c[2];
    }
    dd = t[0]*n[0] + t[1]*n[1] + t[2]*n[2];
    t[0] -= dd*n[0];
    t[1] -= dd*n[1];
    t[2] -= dd*n[2];

    dd = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
    if ( dd > _MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        t[0] *= dd;
        t[1] *= dd;
        t[2] *= dd;
    }

    return(1);
}

/** Return volumic ball of a surfacic point p, as well as the part of its surfacic ball
    supported in the outer boundary starting from tet start, with point ip, and face if in tetra
    volumic ball ; list[k] = 4*number of tet + index of point
    surfacic ball : list[k] = 4*number of tet + index of FACE */
int _MMG5_bouleext(MMG5_pMesh mesh, int start, int ip, int iface, int *listv, int *ilistv, int *lists, int*ilists){
    MMG5_pTetra pt,pt1;
    int base,nump,k,k1,*adja,piv,na,nb,adj,cur,nvstart,fstart,aux;
    char iopp,ipiv,i,j,l,ipa,ipb,isface;

    base = ++mesh->base;
    *ilists = 0;
    *ilistv = 0;

    pt = &mesh->tetra[start];
    nump = pt->v[ip];
    k = start;

    na = pt->v[ip];
    nb = pt->v[_MMG5_idir[iface][_MMG5_inxt2[_MMG5_idirinv[iface][ip]]]];
    piv = pt->v[_MMG5_idir[iface][_MMG5_iprv2[_MMG5_idirinv[iface][ip]]]];

    iopp = iface;
    fstart = 4*k+iopp;

    do {
        /* A boundary face has been hit : change travel edge */
        lists[(*ilists)] = 4*k+iopp;
        (*ilists)++;
        if ( *ilists >= _MMG5_LMAX ) {
            fprintf(stdout,"  ## Warning: problem in surface remesh process.");
            fprintf(stdout," Surface ball of point %d contains too many elts.\n",
                    _MMG5_indPt(mesh,nump));
            fprintf(stdout,"  ##          Try to modify the hausdorff number,");
            fprintf(stdout," or/and the maximum mesh.\n");
            return(-1);
        }

        aux = nb;
        nb = piv;
        piv = aux;
        nvstart = k;
        adj = k;

        /* Now unfold shell of edge (na,nb) starting from k (included)*/
        do {
            k = adj;
            pt = &mesh->tetra[k];
            adja = &mesh->adja[4*(k-1)+1];
            if ( pt->flag != base ) {
                for (i=0; i<4; i++)
                    if ( pt->v[i] == nump )  break;
                assert(i<4);
                listv[(*ilistv)] = 4*k+i;
                (*ilistv)++;
                pt->flag = base;
            }

            /* identification of edge number in tetra k */
            for (i=0; i<6; i++) {
                ipa = _MMG5_iare[i][0];
                ipb = _MMG5_iare[i][1];
                if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                     (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
            }
            assert(i<6);

            /* set sense of travel */
            if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
                adj = adja[ _MMG5_ifar[i][0] ] / 4;
                ipiv = _MMG5_ifar[i][1];
                iopp = _MMG5_ifar[i][0];
                piv = pt->v[ipiv];
            }
            else {
                adj = adja[ _MMG5_ifar[i][1] ] / 4;
                ipiv = _MMG5_ifar[i][0];
                iopp = _MMG5_ifar[i][1];
                piv = pt->v[ipiv];
            }
            isface = (adja[iopp] == 0);
        }
        while ( adj && (adj != nvstart) && !isface );
    }
    while ( 4*k+iopp != fstart );

    /** Now, surfacic ball is complete ; finish travel of volumic ball */
    cur = 0;
    while ( cur < (*ilistv) ) {
        k = listv[cur]/4;
        i = listv[cur]%4;
        adja = &mesh->adja[4*(k-1)+1];

        for (l=0; l<3; l++) {
            i  = _MMG5_inxt3[i];
            k1 = adja[i]/4;
            if ( !k1 )  continue;
            pt1 = &mesh->tetra[k1];
            if ( pt1->flag == base )  continue;
            pt1->flag = base;

            for (j=0; j<4; j++)
                if ( pt1->v[j] == nump )  break;
            assert(j<4);
            /* overflow */
            if ( *ilistv > _MMG5_LMAX-3 ) {
                fprintf(stdout,"  ## Warning: problem in remesh process.");
                fprintf(stdout," Volumic ball of point %d contains too many elts.\n",
                        _MMG5_indPt(mesh,nump));
                fprintf(stdout,"  ##          Try to modify the hausdorff number,");
                fprintf(stdout," or/and the maximum mesh.\n");
                return(-1);
            }
            listv[(*ilistv)] = 4*k1+j;
            (*ilistv)++;
        }
        cur++;
    }

    return(1);
}

/** Return volumic ball of a SURFACE point p, as well as its surfacic ball, starting from tetra
    start, with point ip, and face if in tetra
    volumic ball ; list[k] = 4*number of tet + index of point
    surfacic ball : list[k] = 4*number of tet + index of FACE */
int _MMG5_boulesurfvolp(MMG5_pMesh mesh,int start,int ip,int iface,
                        int *listv,int *ilistv,int *lists,int*ilists) {
    MMG5_pTetra  pt,pt1;
    MMG5_pxTetra pxt;
    int base,nump,k,k1,*adja,piv,na,nb,adj,cur,nvstart,fstart,aux;
    char iopp,ipiv,i,j,l,ipa,ipb,isface;

    base = ++mesh->base;
    *ilists = 0;
    *ilistv = 0;

    pt = &mesh->tetra[start];
    nump = pt->v[ip];
    k = start;

    na = pt->v[ip];
    nb = pt->v[_MMG5_idir[iface][_MMG5_inxt2[_MMG5_idirinv[iface][ip]]]];
    piv = pt->v[_MMG5_idir[iface][_MMG5_iprv2[_MMG5_idirinv[iface][ip]]]];

    iopp = iface;
    fstart = 4*k+iopp;

    do {
        /* A boundary face has been hit : change travel edge */
        lists[(*ilists)] = 4*k+iopp;
        (*ilists)++;
        if ( *ilists >= _MMG5_LMAX ) {
            fprintf(stdout,"  ## Warning: problem in surface remesh process.");
            fprintf(stdout," Surface ball of point %d contains too many elts.\n",
                    _MMG5_indPt(mesh,nump));
            fprintf(stdout,"  ##          Try to modify the hausdorff number,");
            fprintf(stdout," or/and the maximum mesh.\n");
            return(-1);
        }

        aux = nb;
        nb = piv;
        piv = aux;
        nvstart = k;
        adj = k;

        /* Now unfold shell of edge (na,nb) starting from k (included)*/
        do {
            k = adj;
            pt = &mesh->tetra[k];
            adja = &mesh->adja[4*(k-1)+1];
            if ( pt->flag != base ) {
                for (i=0; i<4; i++)
                    if ( pt->v[i] == nump )  break;
                assert(i<4);
                listv[(*ilistv)] = 4*k+i;
                (*ilistv)++;
                pt->flag = base;
            }

            /* identification of edge number in tetra k */
            for (i=0; i<6; i++) {
                ipa = _MMG5_iare[i][0];
                ipb = _MMG5_iare[i][1];
                if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                     (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
            }
            assert(i<6);

            /* set sense of travel */
            if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
                adj = adja[ _MMG5_ifar[i][0] ] / 4;
                ipiv = _MMG5_ifar[i][1];
                iopp = _MMG5_ifar[i][0];
                piv = pt->v[ipiv];
            }
            else {
                adj = adja[ _MMG5_ifar[i][1] ] / 4;
                ipiv = _MMG5_ifar[i][0];
                iopp = _MMG5_ifar[i][1];
                piv = pt->v[ipiv];
            }
            isface = 0;
            if(pt->xt){
                pxt = &mesh->xtetra[pt->xt];
                isface = (MG_BDY & pxt->ftag[iopp]);
            }
        }
        while ( adj && (adj != nvstart) && !isface );
    }
    while ( 4*k+iopp != fstart );

    /* Now, surfacic ball is complete ; finish travel of volumic ball */
    cur = 0;  // Check numerotation
    while ( cur < (*ilistv) ) {
        k = listv[cur]/4;
        i = listv[cur]%4; // index of point p in tetra k
        adja = &mesh->adja[4*(k-1)+1];

        for (l=0; l<3; l++) {
            i  = _MMG5_inxt3[i];
            k1 = adja[i]/4;
            if ( !k1 )  continue;
            pt1 = &mesh->tetra[k1];
            if ( pt1->flag == base )  continue;
            pt1->flag = base;

            for (j=0; j<4; j++)
                if ( pt1->v[j] == nump )  break;
            assert(j<4);

            /* overflow */
            if ( *ilistv > _MMG5_LMAX-3 ) {
                fprintf(stdout,"  ## Warning: problem in remesh process.");
                fprintf(stdout," Volumic ball of point %d contains too many elts.\n",
                        _MMG5_indPt(mesh,nump));
                fprintf(stdout,"  ##          Try to modify the hausdorff number,");
                fprintf(stdout," or/and the maximum mesh.\n");
                return(-1);
            }
            listv[(*ilistv)] = 4*k1+j;
            (*ilistv)++;
        }
        cur++;
    }

    return(1);
}

/** Get tag of edge ia in tetra start by travelling its shell until meeting a boundary face */
static inline int
_MMG5_gettag(MMG5_pMesh mesh,int start,int ia,int *tag,int *edg) {
    MMG5_pTetra        pt;
    MMG5_pxTetra       pxt;
    int           na,nb,*adja,adj,piv;
    unsigned char i,ipa,ipb;

    if ( start < 1 )  return(0);
    pt = &mesh->tetra[start];
    if ( !MG_EOK(pt) ) return(0);

    na   = pt->v[ _MMG5_iare[ia][0] ];
    nb   = pt->v[ _MMG5_iare[ia][1] ];

    adja = &mesh->adja[4*(start-1)+1];
    adj = adja[_MMG5_ifar[ia][0]] / 4;
    piv = pt->v[_MMG5_ifar[ia][1]];

    if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        if ( (pxt->ftag[_MMG5_ifar[ia][0]] & MG_BDY) || (pxt->ftag[_MMG5_ifar[ia][1]] & MG_BDY) ) {
            *tag = pxt->tag[ia];
            *edg = pxt->edg[ia];
            return(1);
        }
    }

    while ( adj && (adj != start) ) {
        pt = &mesh->tetra[adj];
        /* identification of edge number in tetra adj */
        for (i=0; i<6; i++) {
            ipa = _MMG5_iare[i][0];
            ipb = _MMG5_iare[i][1];
            if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                 (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
        }
        assert(i<6);
        if ( pt->xt ) {
            pxt = &mesh->xtetra[pt->xt];
            if ( (pxt->ftag[_MMG5_ifar[i][0]] & MG_BDY) || (pxt->ftag[_MMG5_ifar[i][1]] & MG_BDY) ) {
                *tag = pxt->tag[i];
                *edg = pxt->edg[i];
                return(1);
            }
        }

        /* set new triangle for travel */
        adja = &mesh->adja[4*(adj-1)+1];
        if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
            adj = adja[ _MMG5_ifar[i][0] ] / 4;
            piv = pt->v[ _MMG5_ifar[i][1] ];
        }
        else {
            adj = adja[ _MMG5_ifar[i][1] ] /4;
            piv = pt->v[ _MMG5_ifar[i][0] ];
        }
    }
    return(0);
}

/** Set tag and edg of edge ia (if need be) in tetra start by travelling its shell */
inline int
_MMG5_settag(MMG5_pMesh mesh,int start,int ia,int tag,int edg) {
    MMG5_pTetra        pt;
    MMG5_pxTetra       pxt;
    int           na,nb,*adja,adj,piv;
    unsigned char i,ipa,ipb;

    assert( start >= 1 );
    pt = &mesh->tetra[start];
    assert ( MG_EOK(pt) );

    na   = pt->v[ _MMG5_iare[ia][0] ];
    nb   = pt->v[ _MMG5_iare[ia][1] ];

    adja = &mesh->adja[4*(start-1)+1];
    adj = adja[_MMG5_ifar[ia][0]] / 4;
    piv = pt->v[_MMG5_ifar[ia][1]];

    if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        if ( (pxt->ftag[_MMG5_ifar[ia][0]] & MG_BDY) ||
             (pxt->ftag[_MMG5_ifar[ia][1]] & MG_BDY) ) {
            pxt->tag[ia] |= tag;
            pxt->edg[ia]  = MG_MAX(pxt->edg[ia],edg);
        }
    }
    while ( adj && (adj != start) ) {
        pt = &mesh->tetra[adj];
        /* identification of edge number in tetra adj */
        for (i=0; i<6; i++) {
            ipa = _MMG5_iare[i][0];
            ipb = _MMG5_iare[i][1];
            if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                 (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
        }
        assert(i<6);
        if ( pt->xt ) {
            pxt = &mesh->xtetra[pt->xt];
            if ( (pxt->ftag[_MMG5_ifar[i][0]] & MG_BDY) ||
                 (pxt->ftag[_MMG5_ifar[i][1]] & MG_BDY) ) {
                pxt->tag[i] |= tag;
                pxt->edg[i]  = MG_MAX(pxt->edg[i],edg);
            }
        }
        /* set new triangle for travel */
        adja = &mesh->adja[4*(adj-1)+1];
        if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
            adj = adja[ _MMG5_ifar[i][0] ] / 4;
            piv = pt->v[ _MMG5_ifar[i][1] ];
        }
        else {
            adj = adja[ _MMG5_ifar[i][1] ] /4;
            piv = pt->v[ _MMG5_ifar[i][0] ];
        }
    }

    /* If all shell has been travelled, stop, else, travel it the other sense */
    if ( adj == start )  return(1);
    assert(!adj);

    pt = &mesh->tetra[start];
    adja = &mesh->adja[4*(start-1)+1];
    adj = adja[_MMG5_ifar[ia][1]] / 4;
    piv = pt->v[_MMG5_ifar[ia][0]];

    while ( adj && (adj != start) ) {
        pt = &mesh->tetra[adj];
        /* identification of edge number in tetra adj */
        for (i=0; i<6; i++) {
            ipa = _MMG5_iare[i][0];
            ipb = _MMG5_iare[i][1];
            if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                 (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
        }
        assert(i<6);
        if ( pt->xt ) {
            pxt = &mesh->xtetra[pt->xt];
            if ( (pxt->ftag[_MMG5_ifar[i][0]] & MG_BDY) ||
                 (pxt->ftag[_MMG5_ifar[i][1]] & MG_BDY) ) {
                pxt->tag[i] |= tag;
                pxt->edg[i]  = MG_MAX(pxt->edg[i],edg);
            }
        }
        /* set new triangle for travel */
        adja = &mesh->adja[4*(adj-1)+1];
        if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
            adj = adja[ _MMG5_ifar[i][0] ] / 4;
            piv = pt->v[ _MMG5_ifar[i][1] ];
        }
        else {
            adj = adja[ _MMG5_ifar[i][1] ] /4;
            piv = pt->v[ _MMG5_ifar[i][0] ];
        }
    }
    return(1);
}

/** Find all tets sharing edge ia of tetra start
    return 2*ilist if shell is closed, 2*ilist +1 otherwise
    return 0 if one of the tet of the shell is required */
int _MMG5_coquil(MMG5_pMesh mesh,int start,int ia,int * list) {
    MMG5_pTetra  pt;
    int     ilist,*adja,piv,adj,na,nb,ipa,ipb;
    char    i;

    assert ( start >= 1 );
    pt = &mesh->tetra[start];
    assert ( MG_EOK(pt) );

    na   = pt->v[ _MMG5_iare[ia][0] ];
    nb   = pt->v[ _MMG5_iare[ia][1] ];
    ilist = 0;
    list[ilist] = 6*start+ia;
    ilist++;

    adja = &mesh->adja[4*(start-1)+1];
    adj = adja[_MMG5_ifar[ia][0]] / 4; // start travelling by face (ia,0)
    piv = pt->v[_MMG5_ifar[ia][1]];

    while ( adj && (adj != start) ) {
        pt = &mesh->tetra[adj];
        if ( pt->tag & MG_REQ )  return(0);
        /* identification of edge number in tetra adj */
        for (i=0; i<6; i++) {
            ipa = _MMG5_iare[i][0];
            ipb = _MMG5_iare[i][1];
            if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                 (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
        }
        assert(i<6);
        list[ilist] = 6*adj +i;
        ilist++;
        /* overflow */
        if ( ilist > _MMG5_LMAX-3 ) {
            fprintf(stdout,"  ## Warning: problem in remesh process.");
            fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
                    _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
            fprintf(stdout,"  ##          Try to modify the hausdorff number,");
            fprintf(stdout," or/and the maximum mesh.\n");
            return(-1);
        }

        /* set new triangle for travel */
        adja = &mesh->adja[4*(adj-1)+1];
        if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
            adj = adja[ _MMG5_ifar[i][0] ] / 4;
            piv = pt->v[ _MMG5_ifar[i][1] ];
        }
        else {
            assert(pt->v[ _MMG5_ifar[i][1] ] == piv );
            adj = adja[ _MMG5_ifar[i][1] ] /4;
            piv = pt->v[ _MMG5_ifar[i][0] ];
        }
    }

    /* At this point, the first travel, in one direction, of the shell is
       complete. Now, analyze why the travel ended. */
    if ( adj == start )  return(2*ilist);
    assert(!adj); // a boundary has been detected

    adj = list[ilist-1] / 6;
    i   = list[ilist-1] % 6;
    ilist = 0;

    /* Start back everything from this tetra adj */
    list[ilist] = 6*adj + i;
    ilist++;
    /* overflow */
    if ( ilist > _MMG5_LMAX-3 ) {
        fprintf(stdout,"  ## Warning: problem in remesh process.");
        fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
                _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
        fprintf(stdout,"  ##          Try to modify the hausdorff number,");
        fprintf(stdout," or/and the maximum mesh.\n");
        return(-1);
    }

    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
        adj = adja[ _MMG5_ifar[i][0] ] / 4;
        piv = pt->v[ _MMG5_ifar[i][1] ];
    }
    else {
        adj = adja[ _MMG5_ifar[i][1] ] /4;
        piv = pt->v[ _MMG5_ifar[i][0] ];
    }

    while ( adj ) {
        pt = &mesh->tetra[adj];
        if ( pt->tag & MG_REQ )  return(0);
        /* identification of edge number in tetra adj */
        for (i=0; i<6; i++) {
            ipa = _MMG5_iare[i][0];
            ipb = _MMG5_iare[i][1];
            if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                 (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
        }
        assert(i<6);
        list[ilist] = 6*adj +i;
        ilist++;
        /* overflow */
        if ( ilist > _MMG5_LMAX-2 ) {
            fprintf(stdout,"  ## Warning: problem in surface remesh process.");
            fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
                    _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
            fprintf(stdout,"  ##          Try to modify the hausdorff number,");
            fprintf(stdout," or/and the maximum mesh.\n");
            return(-1);
        }

        /* set new triangle for travel */
        adja = &mesh->adja[4*(adj-1)+1];
        if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
            adj = adja[ _MMG5_ifar[i][0] ] / 4;
            piv = pt->v[ _MMG5_ifar[i][1] ];
        }
        else {
            adj = adja[ _MMG5_ifar[i][1] ] /4;
            piv = pt->v[ _MMG5_ifar[i][0] ];
        }
    }
    assert(!adj);
    return( 2*ilist+1 );
}

/** Identify whether edge ia in start is a boundary edge by unfolding its shell */
int _MMG5_srcbdy(MMG5_pMesh mesh,int start,int ia) {
    MMG5_pTetra      pt;
    MMG5_pxTetra     pxt;
    int         na,nb,adj,piv,*adja;
    char        ipa,ipb,iadj,i;

    pt = &mesh->tetra[start];
    na = pt->v[_MMG5_iare[ia][0]];
    nb = pt->v[_MMG5_iare[ia][1]];

    adja = &mesh->adja[4*(start-1)+1];
    iadj = _MMG5_ifar[ia][0];

    if(pt->xt){
        pxt = &mesh->xtetra[pt->xt];
        if( pxt->ftag[iadj] & MG_BDY )
            return(1);
    }

    adj = adja[iadj] / 4;
    piv = pt->v[_MMG5_ifar[ia][1]];

    while( adj && ( adj != start ) ) {
        pt = &mesh->tetra[adj];

        /* identification of edge number in tetra adj */
        for(i=0; i<6; i++) {
            ipa = _MMG5_iare[i][0];
            ipb = _MMG5_iare[i][1];
            if( ( pt->v[ipa] == na && pt->v[ipb] == nb ) ||
                ( pt->v[ipa] == nb && pt->v[ipb] == na ))
                break;

        }
        assert(i<6);

        /* set new triangle for travel */
        adja = &mesh->adja[4*(adj-1)+1];
        if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
            iadj = _MMG5_ifar[i][0];
            adj = adja[ iadj ] / 4;
            piv = pt->v[ _MMG5_ifar[i][1] ];
        }
        else {
            iadj = _MMG5_ifar[i][1];
            adj = adja[ iadj ] /4;
            piv = pt->v[ _MMG5_ifar[i][0] ];
        }

        if(pt->xt){
            pxt = &mesh->xtetra[pt->xt];
            if( pxt->ftag[iadj] & MG_BDY )
                return(1);
        }
    }

    return(0);
}

/** print an error message if _MMG5_coquilFace detect a boundary topology problem */
static inline void
_MMG5_errorMessage(MMG5_pMesh mesh, int k1, int k2) {
    MMG5_pPoint ppt;
    MMG5_pTetra pt;
    int    np, ne, k, kel1, kel2;

    np = ne = kel1 = kel2 = 0;
    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( MG_VOK(ppt) )  ppt->tmp = ++np;
    }
    for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( MG_EOK(pt) ) {
            ne++;
            if ( k == k1 ) kel1 = ne;
            if ( k == k2 ) kel2 = ne;
        }
    }

    fprintf(stdout,"  ## Error: problem in surface remesh process");
    fprintf(stdout," (potential creation of a lonely boundary face):\n");

    if ( kel1 != 0 ) {
        pt = &mesh->tetra[k1];
        fprintf(stdout,"            look at elt %d:",kel1);
        fprintf(stdout," %d %d %d %d.\n", mesh->point[pt->v[0]].tmp,
                mesh->point[pt->v[1]].tmp,mesh->point[pt->v[2]].tmp,
                mesh->point[pt->v[3]].tmp);
        fprintf(stdout,"adj %d %d %d %d\n",(&mesh->adja[3*(kel1-1)+1])[0],
                (&mesh->adja[3*(kel1-1)+1])[1],(&mesh->adja[3*(kel1-1)+1])[2],
                (&mesh->adja[3*(kel1-1)+1])[3]);
        fprintf(stdout,"req %d %d %d %d\n",mesh->point[pt->v[0]].tag & MG_REQ,
                mesh->point[pt->v[1]].tag & MG_REQ,
                mesh->point[pt->v[2]].tag & MG_REQ,mesh->point[pt->v[3]].tag & MG_REQ);
    } else if ( kel2 != 0 ) {
        fprintf(stdout,"            look at elt %d:",kel2);
        mesh->tetra[kel2].ref=5;
        fprintf(stdout," %d %d %d %d.\n", mesh->point[pt->v[0]].tmp,
                mesh->point[pt->v[1]].tmp,mesh->point[pt->v[2]].tmp,
                mesh->point[pt->v[3]].tmp);
    }
    fprintf(stdout,"  ##        Try to modify the hausdorff number,");
    fprintf(stdout," the maximum mesh size or/and the value of angle detection.\n");
    fprintf(stdout," You can also try to run with -noswap option but probably");
    fprintf(stdout," the final mesh will have poor quality.\n");
}

/** Find all tets sharing edge ia of tetra start, and stores boundary faces when met
    it1 & it2 = 6*iel + iface, iel = index of tetra, iface = index of face in tetra
    return 2*ilist if shell is closed, 2*ilist +1 otherwise */
int _MMG5_coquilface(MMG5_pMesh mesh,int start,int ia,int *list,int *it1,int *it2) {
    MMG5_pTetra   pt;
    MMG5_pxTetra  pxt;
    int     *adja,piv,adj,na,nb,ipa,ipb,ilist,pradj;
    char     i,iface,isbdy;

    pt = &mesh->tetra[start];

    na   = pt->v[ _MMG5_iare[ia][0] ];
    nb   = pt->v[ _MMG5_iare[ia][1] ];

    ilist = 0;
    list[ilist] = 6*start+ia;
    ilist++;

    *it1 = 0;
    *it2 = 0;

    adja = &mesh->adja[4*(start-1)+1];
    adj = adja[_MMG5_ifar[ia][0]] / 4;
    piv = pt->v[_MMG5_ifar[ia][1]];

    pxt = &mesh->xtetra[pt->xt];

    iface = _MMG5_ifar[ia][1];
    isbdy = pxt->ftag[iface];
    if ( isbdy )
        *it1 = 4*start + iface;

    while ( adj && (adj != start) ) {
        pt = &mesh->tetra[adj];
        pxt = 0;
        if ( pt->xt )
            pxt = &mesh->xtetra[pt->xt];

        /* identification of edge number in tetra adj */
        for (i=0; i<6; i++) {
            ipa = _MMG5_iare[i][0];
            ipb = _MMG5_iare[i][1];
            if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                 (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
        }
        assert(i<6);
        list[ilist] = 6*adj +i;
        ilist++;
        /* overflow */
        if ( ilist > _MMG5_LMAX-2 ) {
            fprintf(stdout,"  ## Warning: problem in surface remesh process.");
            fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
                    _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
            fprintf(stdout,"  ##          Try to modify the hausdorff number,");
            fprintf(stdout," or/and the maximum mesh.\n");
            return(-1);
        }

        /* set new tetra for travel */
        pradj = adj;
        adja = &mesh->adja[4*(adj-1)+1];
        if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
            adj = adja[ _MMG5_ifar[i][0] ] / 4;
            piv = pt->v[ _MMG5_ifar[i][1] ];
            iface = _MMG5_ifar[i][1];
        }
        else {
            assert(pt->v[ _MMG5_ifar[i][1] ] == piv );
            adj = adja[ _MMG5_ifar[i][1] ] /4;
            piv = pt->v[ _MMG5_ifar[i][0] ];
            iface = _MMG5_ifar[i][0];
        }
        isbdy = pt->xt ? pxt->ftag[iface] : 0;

        if ( isbdy ) {
            if ( *it1 == 0 )
                *it1 = 4*pradj+iface;
            else {
                if ( *it2 ) {
                    // Algiane: (assert commente) si on a 3 tetras de refs differentes dans la _MMG5_coquille??
                    printf("  ## Warning: you have more than 2 boundaries in the shell of your edge.\n");
                    printf("  Problem may occur during remesh process.\n");
                }
                //assert( *it2 == 0 );
                *it2 = 4*pradj+iface;
            }
        }
    }

    /* At this point, the first travel, in one direction, of the shell is complete. Now, analyze why
       the travel ended. */
    if ( adj == start ) {
        if ( (!(*it1) || !(*it2)) || ((*it1) == (*it2)) ) {
            _MMG5_errorMessage(mesh, (*it1)/4, (*it2)/4);
            return(-1);
        }
        return(2*ilist);
    }

    /* A boundary has been detected : slightly different configuration */
    assert(!adj);
    adj = list[ilist-1] / 6;
    i   = list[ilist-1] % 6;
    ilist = 0;

    /* Start back everything from this tetra adj */
    pradj = adj;
    /* overflow */
    if ( ilist > _MMG5_LMAX-2 ) {
        fprintf(stdout,"  ## Warning: problem in surface remesh process.");
        fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
                _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
        fprintf(stdout,"  ##          Try to modify the hausdorff number,");
        fprintf(stdout," or/and the maximum mesh.\n");
        return(-1);
    }

    pxt = 0;
    if ( pt->xt )
        pxt = &mesh->xtetra[pt->xt];

    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
        iface = _MMG5_ifar[i][1];
    }
    else {
        iface = _MMG5_ifar[i][0];
    }
    isbdy = pt->xt ? pxt->ftag[iface] : 0;
    if ( isbdy ) {
        *it1 = 4*pradj + iface;
    }

    while ( adj ) {
        pt = &mesh->tetra[adj];

        /* identification of edge number in tetra adj */
        for (i=0; i<6; i++) {
            ipa = _MMG5_iare[i][0];
            ipb = _MMG5_iare[i][1];
            if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                 (pt->v[ipa] == nb && pt->v[ipb] == na) )  break;
        }
        assert(i<6);
        list[ilist] = 6*adj +i;
        ilist++;
        /* overflow */
        if ( ilist > _MMG5_LMAX-2 ) {
            fprintf(stdout,"  ## Warning: problem in surface remesh process.");
            fprintf(stdout," Coquil of edge %d-%d contains too many elts.\n",
                    _MMG5_indPt(mesh,na),_MMG5_indPt(mesh,nb));
            fprintf(stdout,"  ##          Try to modify the hausdorff number,");
            fprintf(stdout," or/and the maximum mesh.\n");
            return(-1);
        }

        /* set new tetra for travel */
        pradj = adj;
        adja = &mesh->adja[4*(adj-1)+1];
        if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
            adj = adja[ _MMG5_ifar[i][0] ] / 4;
            piv = pt->v[ _MMG5_ifar[i][1] ];
            iface = _MMG5_ifar[i][0];
        }
        else {
            adj = adja[ _MMG5_ifar[i][1] ] /4;
            piv = pt->v[ _MMG5_ifar[i][0] ];
            iface = _MMG5_ifar[i][1];
        }
    }

    assert(!adj);
    *it2 = 4*pradj + iface;

    if ( (!(*it1) || !(*it2)) || ((*it1) == (*it2)) ) {
        _MMG5_errorMessage(mesh, (*it1)/4, (*it2)/4);
        return(-1);
    }
    return ( 2*ilist+1 );
}

