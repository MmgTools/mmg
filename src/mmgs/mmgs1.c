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
 * \file mmgs/mmgs1.c
 * \brief Perform surface mesh adaptation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgs.h"

extern Info info;
char ddb;


/* check if edge need to be split and return a binary coding the numbers of the edges of tria iel
   that should be split according to a hausdorff distance criterion */
int chkedg(pMesh mesh,int iel) {
    pTria    pt;
    pPoint   p[3];
    double   n[3][3],t[3][3],nt[3],c1[3],c2[3],*n1,*n2,t1[3],t2[3];
    double   ps,ps2,cosn,ux,uy,uz,ll,li,dd;
    char     i,i1,i2;

    pt   = &mesh->tria[iel];
    p[0] = &mesh->point[pt->v[0]];
    p[1] = &mesh->point[pt->v[1]];
    p[2] = &mesh->point[pt->v[2]];

    /* normal recovery */
    for (i=0; i<3; i++) {
        if ( MS_SIN(p[i]->tag) ) {
            nortri(mesh,pt,n[i]);
        }
        else if ( MS_EDG(p[i]->tag) ) {
            nortri(mesh,pt,nt);
            n1  = &mesh->geom[p[i]->ig].n1[0];
            n2  = &mesh->geom[p[i]->ig].n2[0];
            ps  = n1[0]*nt[0] + n1[1]*nt[1] + n1[2]*nt[2];
            ps2 = n2[0]*nt[0] + n2[1]*nt[1] + n2[2]*nt[2];
            if ( fabs(ps) > fabs(ps2) )
                memcpy(&n[i],n1,3*sizeof(double));
            else
                memcpy(&n[i],n2,3*sizeof(double));
            memcpy(&t[i],p[i]->n,3*sizeof(double));
        }
        else
            memcpy(&n[i],p[i]->n,3*sizeof(double));
    }

    /* analyze edges */
    for (i=0; i<3; i++) {
        i1 = inxt[i];
        i2 = iprv[i];

        /* check length */
        ux = p[i2]->c[0] - p[i1]->c[0];
        uy = p[i2]->c[1] - p[i1]->c[1];
        uz = p[i2]->c[2] - p[i1]->c[2];
        ll = ux*ux + uy*uy + uz*uz;
        if ( ll > info.hmax*info.hmax ) {
            MS_SET(pt->flag,i);
            continue;
        }
        else if ( !MS_EDG(pt->tag[i]) && p[i1]->tag > MS_NOTAG && p[i2]->tag > MS_NOTAG ) {
            MS_SET(pt->flag,i);
            continue;
        }

        /* Hausdorff w/r tangent direction */
        if ( MS_EDG(pt->tag[i]) ) {
            if ( MS_SIN(p[i1]->tag) ) {
                li = 1.0 / sqrt(ll);
                t1[0] = li*ux;
                t1[1] = li*uy;
                t1[2] = li*uz;
            }
            else{
                memcpy(t1,t[i1],3*sizeof(double));
            }

            if ( MS_SIN(p[i2]->tag) ) {
                li = 1.0 / sqrt(ll);
                t2[0] = li*ux;
                t2[1] = li*uy;
                t2[2] = li*uz;
            }
            else{
                memcpy(t2,t[i2],3*sizeof(double));
            }

            ps = t1[0]*ux + t1[1]*uy + t1[2]*uz;
            ps *= ps;
            cosn = ps/ll ;
            cosn *= (1.0-cosn);
            cosn *= (0.25*ll);
            if ( cosn > info.hausd*info.hausd ) {
                MS_SET(pt->flag,i);
                continue;
            }

            ps = -(t2[0]*ux + t2[1]*uy + t2[2]*uz);
            ps *= ps;
            cosn = ps/ll ;
            cosn *= (1.0-cosn);
            cosn *= (0.25*ll);
            if ( cosn > info.hausd*info.hausd ) {
                MS_SET(pt->flag,i);
                continue;
            }
        }
        else {
            n1 = n[i1];
            n2 = n[i2];

            ps = ux*n1[0] + uy*n1[1] + uz*n1[2];
            c1[0] = (2.0*p[i1]->c[0] + p[i2]->c[0] - ps*n1[0]) / 3.0 - p[i1]->c[0];
            c1[1] = (2.0*p[i1]->c[1] + p[i2]->c[1] - ps*n1[1]) / 3.0 - p[i1]->c[1];
            c1[2] = (2.0*p[i1]->c[2] + p[i2]->c[2] - ps*n1[2]) / 3.0 - p[i1]->c[2];

            ps = -(ux*n2[0] + uy*n2[1] + uz*n2[2]);
            c2[0] = (2.0*p[i2]->c[0] + p[i1]->c[0] - ps*n2[0]) / 3.0 - p[i2]->c[0];
            c2[1] = (2.0*p[i2]->c[1] + p[i1]->c[1] - ps*n2[1]) / 3.0 - p[i2]->c[1];
            c2[2] = (2.0*p[i2]->c[2] + p[i1]->c[2] - ps*n2[2]) / 3.0 - p[i2]->c[2];

            /* squared cosines */
            ps = c1[0]*ux + c1[1]*uy + c1[2]*uz;
            ps *= ps;
            dd = c1[0]*c1[0] + c1[1]*c1[1] + c1[2]*c1[2];
            cosn  =  ps / (dd*ll);
            cosn *= (1.0-cosn);
            cosn *= (0.25*ll);
            if ( cosn > info.hausd*info.hausd ) {
                MS_SET(pt->flag,i);
                continue;
            }

            ps = -c2[0]*ux - c2[1]*uy - c2[2]*uz;
            ps *= ps;
            dd = c2[0]*c2[0]+c2[1]*c2[1]+c2[2]*c2[2];
            cosn  =  ps / (dd*ll);
            cosn *= (1.0-cosn);
            cosn *= (0.25*ll);
            if ( cosn > info.hausd*info.hausd ) {
                MS_SET(pt->flag,i);
                continue;
            }
        }
    }

    return(pt->flag);
}

static int swpmsh(pMesh mesh,pSol met,char typchk) {
    pTria    pt;
    int      k,it,ns,nns,maxit;
    char     i;

    it = nns = 0;
    maxit = 2;
    mesh->base++;
    do {
        ns = 0;
        for (k=1; k<=mesh->nt; k++) {
            pt = &mesh->tria[k];
            if ( !MS_EOK(pt) || pt->ref < 0 )   continue;

            for (i=0; i<3; i++) {
                if ( MS_SIN(pt->tag[i]) || MS_EDG(pt->tag[i]) )  continue;
                else if ( chkswp(mesh,met,k,i,typchk) ) {
                    ns += swapar(mesh,k,i);
                    break;
                }
            }
        }
        nns += ns;
    }
    while ( ns > 0 && ++it < maxit );
    if ( (abs(info.imprim) > 5 || info.ddebug) && nns > 0 )
        fprintf(stdout,"     %8d edge swapped\n",nns);

    return(nns);
}

/* Analyze triangles and move points to make mesh more uniform */
static int movtri(pMesh mesh,pSol met,int maxit) {
    pTria    pt;
    pPoint   ppt;
    int      it,k,ier,base,nm,ns,nnm,list[LMAX+2],ilist;
    char     i;

    if ( abs(info.imprim) > 5 || info.ddebug )
        fprintf(stdout,"  ** OPTIMIZING MESH\n");

    base = 1;
    for (k=1; k<=mesh->np; k++)  mesh->point[k].flag = base;

    it = nnm = 0;
    do {
        base++;
        nm = ns = 0;
        for (k=1; k<=mesh->nt; k++) {
            pt = &mesh->tria[k];
            if ( !MS_EOK(pt) || pt->ref < 0 )   continue;

            for (i=0; i<3; i++) {
                ppt = &mesh->point[pt->v[i]];

                if ( ppt->flag == base || MS_SIN(ppt->tag) || ppt->tag & MS_NOM )
                    continue;
                ier = 0;
                ilist = boulet(mesh,k,i,list);

                if ( MS_EDG(ppt->tag) ) {
                    ier = movridpt(mesh,met,list,ilist);
                    if ( ier )  ns++;
                }
                else
                    ier = movintpt(mesh,met,list,ilist);
                if ( ier ) {
                    nm++;
                    ppt->flag = base;
                }
            }
        }
        nnm += nm;
        if ( info.ddebug )  fprintf(stdout,"     %8d moved, %d geometry\n",nm,ns);
    }
    while ( ++it < maxit && nm > 0);

    if ( (abs(info.imprim) > 5 || info.ddebug) && nnm > 0 )
        fprintf(stdout,"     %8d vertices moved, %d iter.\n",nnm,it);

    return(nnm);
}

/* analyze triangles and split if needed */
static int anaelt(pMesh mesh,pSol met,char typchk) {
    pTria    pt;
    pPoint   ppt,p1,p2;
    Hash     hash;
    Bezier   pb;
    pGeom    go;
    double   s,o[3],no[3],to[3],dd,len;
    int      vx[3],i,j,ip,ip1,ip2,ier,k,ns,nc,nt;
    char     i1,i2;
    static double uv[3][2] = { {0.5,0.5}, {0.,0.5}, {0.5,0.} };

    hashNew(&hash,mesh->np);
    ns = 0;
    s  = 0.5;
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) || pt->ref < 0 )  continue;
        if ( MS_SIN(pt->tag[0]) || MS_SIN(pt->tag[1]) || MS_SIN(pt->tag[2]) )  continue;

        /* check element cut */
        pt->flag = 0;
        if ( typchk == 1 ) {
            if ( !chkedg(mesh,k) )  continue;
        }
        else if ( typchk == 2 ) {
            for (i=0; i<3; i++) {
                i1 = inxt[i];
                i2 = iprv[i];
                len = lenedg(mesh,met,pt->v[i1],pt->v[i2],0);
                if ( len > LLONG )  MS_SET(pt->flag,i);
            }
            if ( !pt->flag )  continue;
        }
        ns++;

        /* geometric support */
        ier = bezierCP(mesh,k,&pb);
        assert(ier);

        /* scan edges to split */
        for (i=0; i<3; i++) {
            if ( !MS_GET(pt->flag,i) )  continue;
            i1  = inxt[i];
            i2  = iprv[i];
            ip1 = pt->v[i1];
            ip2 = pt->v[i2];
            ip = hashGet(&hash,ip1,ip2);
            if ( !MS_EDG(pt->tag[i]) && ip > 0 )  continue;

            /* new point along edge */
            ier = bezierInt(&pb,uv[i],o,no,to);
            if ( !ip ) {
                ip = newPt(mesh,o,MS_EDG(pt->tag[i]) ? to : no);
                assert(ip);
                hashEdge(&hash,ip1,ip2,ip);
                p1  = &mesh->point[ip1];
                p2  = &mesh->point[ip2];
                ppt = &mesh->point[ip];

                if ( MS_EDG(pt->tag[i]) ) {
                    ++mesh->ng;
                    assert(mesh->ng < mesh->ngmax);
                    ppt->tag = pt->tag[i];
                    if ( p1->ref == pt->edg[i] || p2->ref == pt->edg[i] )
                        ppt->ref = pt->edg[i];
                    ppt->ig  = mesh->ng;
                    go = &mesh->geom[mesh->ng];
                    memcpy(go->n1,no,3*sizeof(double));

                    dd = go->n1[0]*ppt->n[0] + go->n1[1]*ppt->n[1] + go->n1[2]*ppt->n[2];
                    ppt->n[0] -= dd*go->n1[0];
                    ppt->n[1] -= dd*go->n1[1];
                    ppt->n[2] -= dd*go->n1[2];
                    dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
                    if ( dd > EPSD2 ) {
                        dd = 1.0 / sqrt(dd);
                        ppt->n[0] *= dd;
                        ppt->n[1] *= dd;
                        ppt->n[2] *= dd;
                    }
                }
                if ( met->m ) {
                    if ( typchk == 1 )
                        intmet33(mesh,met,ip1,ip2,ip,s);
                    else
                        intmet(mesh,met,k,i,ip,s);
                }
            }
            else if ( pt->tag[i] & MS_GEO ) {
                ppt = &mesh->point[ip];
                go  = &mesh->geom[ppt->ig];
                memcpy(go->n2,no,3*sizeof(double));

                /* a computation of the tangent with respect to these two normals is possible */
                ppt->n[0] = go->n1[1]*go->n2[2] - go->n1[2]*go->n2[1];
                ppt->n[1] = go->n1[2]*go->n2[0] - go->n1[0]*go->n2[2];
                ppt->n[2] = go->n1[0]*go->n2[1] - go->n1[1]*go->n2[0];
                dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
                if ( dd > EPSD2 ) {
                    dd = 1.0 / sqrt(dd);
                    ppt->n[0] *= dd;
                    ppt->n[1] *= dd;
                    ppt->n[2] *= dd;
                }
            }
        }
    }
    if ( !ns ) {
        free(hash.item);
        return(ns);
    }

    /* step 2. checking if split by adjacent */
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) || pt->ref < 0 )  continue;
        else if ( pt->flag == 7 )  continue;

        /* geometric support */
        ier = bezierCP(mesh,k,&pb);
        assert(ier);
        nc = 0;

        for (i=0; i<3; i++) {
            i1 = inxt[i];
            i2 = inxt[i1];
            if ( !MS_GET(pt->flag,i) && !MS_SIN(pt->tag[i]) ) {
                ip = hashGet(&hash,pt->v[i1],pt->v[i2]);
                if ( ip > 0 ) {
                    MS_SET(pt->flag,i);
                    nc++;
                    if ( pt->tag[i] & MS_GEO ) {
                        /* new point along edge */
                        ier = bezierInt(&pb,uv[i],o,no,to);
                        assert(ier);

                        ppt = &mesh->point[ip];
                        go  = &mesh->geom[ppt->ig];
                        memcpy(go->n2,no,3*sizeof(double));

                        /* a computation of the tangent with respect to these two normals is possible */
                        ppt->n[0] = go->n1[1]*go->n2[2] - go->n1[2]*go->n2[1];
                        ppt->n[1] = go->n1[2]*go->n2[0] - go->n1[0]*go->n2[2];
                        ppt->n[2] = go->n1[0]*go->n2[1] - go->n1[1]*go->n2[0];
                        dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
                        if ( dd > EPSD2 ) {
                            dd = 1.0 / sqrt(dd);
                            ppt->n[0] *= dd;
                            ppt->n[1] *= dd;
                            ppt->n[2] *= dd;
                        }
                    }
                }
            }
        }
        if ( nc > 0 )  ++ns;
    }
    if ( info.ddebug && ns ) {
        fprintf(stdout,"     %d analyzed  %d proposed\n",mesh->nt,ns);
        fflush(stdout);
    }

    /* step 3. splitting */
    ns = 0;
    nt = mesh->nt;
    for (k=1; k<=nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) || pt->ref < 0 )  continue;
        else if ( pt->flag == 0 )  continue;

        j  = -1;
        vx[0] = vx[1] = vx[2] = 0;
        for (i=0; i<3; i++) {
            i1 = inxt[i];
            i2 = inxt[i1];
            if ( MS_GET(pt->flag,i) ) {
                vx[i] = hashGet(&hash,pt->v[i1],pt->v[i2]);
                assert(vx[i]);
                j = i;
            }
        }
        if ( pt->flag == 1 || pt->flag == 2 || pt->flag == 4 ) {
            ier = split1(mesh,met,k,j,vx);
            assert(ier);
            ns++;
        }
        else if ( pt->flag == 7 ) {
            ier = split3(mesh,met,k,vx);
            assert(ier);
            ns++;
        }
        else {
            ier = split2(mesh,met,k,vx);
            assert(ier);
            ns++;
        }
    }
    if ( (info.ddebug || abs(info.imprim) > 5) && ns > 0 )
        fprintf(stdout,"     %7d splitted\n",ns);
    free(hash.item);

    return(ns);
}

/* check if splitting edge i of k is ok */
int chkspl(pMesh mesh,pSol met,int k,int i) {
    pTria    pt,pt1;
    pPoint   ppt;
    pGeom    go;
    Bezier   b;
    double   s,uv[2],o[3],no[3],to[3];
    int     *adja,jel,ip,ier;
    char     i1,i2,j,jj,j2;

    if ( mesh->ng > mesh->ngmax-2 )  return(0);
    pt = &mesh->tria[k];
    i1 = inxt[i];
    i2 = iprv[i];
    if ( MS_SIN(pt->tag[i1]) || MS_SIN(pt->tag[i2]) )  return(0);
    adja = &mesh->adja[3*(k-1)+1];
    jel  = adja[i] / 3;
    if ( jel ) {
        j   = adja[i] % 3;
        jj  = inxt[j];
        j2  = iprv[j];
        pt1 = &mesh->tria[jel];
        if ( MS_SIN(pt1->tag[jj]) || MS_SIN(pt1->tag[j2]) )  return(0);
    }

    ier = bezierCP(mesh,k,&b);
    assert(ier);

    /* create midedge point */
    uv[0] = 0.5;
    uv[1] = 0.5;
    if (i == 1)         uv[0] = 0.0;
    else if ( i == 2 )  uv[1] = 0.0;

    ier = bezierInt(&b,uv,o,no,to);
    assert(ier);
    ip = newPt(mesh,o,MS_EDG(pt->tag[i]) ? to : no);
    assert(ip);

    if ( MS_EDG(pt->tag[i]) ) {
        ++mesh->ng;
        ppt = &mesh->point[ip];
        ppt->ig  = mesh->ng;
        go = &mesh->geom[mesh->ng];
        memcpy(go->n1,no,3*sizeof(double));
    }
    s = 0.5;

    intmet(mesh,met,k,i,ip,s);

    return(ip);
}

/* attempt to collapse small edges */
static int colelt(pMesh mesh,pSol met,char typchk) {
    pTria    pt;
    pPoint   p1,p2;
    double   ll,ux,uy,uz;
    int      ier,list[LMAX+2],ilist,k,nc;
    char     i,i1,i2;

    nc = 0;
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) || pt->ref < 0 )   continue;

        /* check edge length */
        pt->flag = 0;
        ier = 0;
        for (i=0; i<3; i++) {
            if ( MS_SIN(pt->tag[i]) )  continue;

            i1 = inxt[i];
            i2 = iprv[i];
            p1 = &mesh->point[pt->v[i1]];
            p2 = &mesh->point[pt->v[i2]];
            if ( p1->tag & MS_NOM || p2->tag & MS_NOM )  continue;
            else if ( MS_SIN(p1->tag) )   continue;
            else if ( p1->tag > p2->tag || p1->tag > pt->tag[i] )  continue;

            /* check length */
            if ( typchk == 1 ) {
                ux = p2->c[0] - p1->c[0];
                uy = p2->c[1] - p1->c[1];
                uz = p2->c[2] - p1->c[2];
                ll = ux*ux + uy*uy + uz*uz;
                if ( ll > info.hmin*info.hmin )  continue;
            }
            else {
                ll = lenedg(mesh,met,pt->v[i1],pt->v[i2],0);
                if ( ll > LSHRT )  continue;
            }

            /* check if geometry preserved */
            ilist = chkcol(mesh,met,k,i,list,typchk);
            if ( ilist > 3 ) {
                nc += colver(mesh,list,ilist);
                break;
            }
            else if ( ilist == 3 ) {
                nc += colver3(mesh,list);
                break;
            }
            else if ( ilist == 2 ) {
                nc += colver2(mesh,list);
                break;
            }
        }
    }
    if ( nc > 0 && (abs(info.imprim) > 5 || info.ddebug) )
        fprintf(stdout,"     %8d vertices removed\n",nc);

    return(nc);
}

static int adpspl(pMesh mesh,pSol met) {
    pTria    pt;
    pPoint   p1,p2;
    double   len,lmax;
    int      ip,k,ns;
    char     i,i1,i2,imax;

    ns = 0;
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) || pt->ref < 0 )   continue;

        /* check edge length */
        pt->flag = 0;
        imax = -1;
        lmax = -1.0;
        for (i=0; i<3; i++) {
            i1  = inxt[i];
            i2  = iprv[i];
            len = lenedg(mesh,met,pt->v[i1],pt->v[i2],0);
            if ( len > lmax ) {
                lmax = len;
                imax = i;
            }
        }
        if ( lmax < LOPTL )  continue;
        else if ( MS_SIN(pt->tag[imax]) )  continue;

        /* check length */
        i1 = inxt[imax];
        i2 = iprv[imax];
        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt->v[i2]];
        if ( p1->tag & MS_NOM || p2->tag & MS_NOM )  continue;

        ip = chkspl(mesh,met,k,imax);
        if ( ip > 0 )
            ns += split1b(mesh,k,imax,ip);
    }
    return(ns);
}

/* analyze triangles and split or collapse to match gradation */
static int adpcol(pMesh mesh,pSol met) {
    pTria    pt;
    pPoint   p1,p2;
    double   len;
    int      k,list[LMAX+2],ilist,nc;
    char     i,i1,i2;

    nc = 0;
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) || pt->ref < 0 )   continue;

        /* check edge length */
        pt->flag = 0;
        for (i=0; i<3; i++) {
            if ( MS_SIN(pt->tag[i]) )  continue;

            /* check length */
            i1 = inxt[i];
            i2 = iprv[i];
            p1 = &mesh->point[pt->v[i1]];
            p2 = &mesh->point[pt->v[i2]];
            if ( p1->tag & MS_NOM || p2->tag & MS_NOM )  continue;

            len = lenedg(mesh,met,pt->v[i1],pt->v[i2],0);
            if ( len > LOPTS )  continue;

            p1 = &mesh->point[pt->v[i1]];
            p2 = &mesh->point[pt->v[i2]];
            if ( MS_SIN(p1->tag) )  continue;
            else if ( p1->tag > p2->tag || p1->tag > pt->tag[i] )  continue;

            /* check if geometry preserved */
            ilist = chkcol(mesh,met,k,i,list,2);
            if ( ilist > 3 ) {
                nc += colver(mesh,list,ilist);
                break;
            }
            else if ( ilist == 3 ) {
                nc += colver3(mesh,list);
                break;
            }
            else if ( ilist == 2 ) {
                nc += colver2(mesh,list);
                break;
            }
        }
    }
    return(nc);
}


/* analyze triangles and split or collapse to match gradation */
static int adptri(pMesh mesh,pSol met) {
    int        it,nnc,nns,nnf,nnm,maxit,nc,ns,nf,nm;

    /* iterative mesh modifications */
    it = nnc = nns = nnf = nnm = 0;
    maxit = 10;
    do {
        ns = adpspl(mesh,met);
        if ( ns < 0 ) {
            fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
            return(0);
        }

        nc = adpcol(mesh,met);
        if ( nc < 0 ) {
            fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
            return(0);
        }
        nf = nm = 0;

        nf = swpmsh(mesh,met,2);
        if ( nf < 0 ) {
            fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
            return(0);
        }

        nm = movtri(mesh,met,1);
        if ( nm < 0 ) {
            fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
            return(0);
        }

        nnc += nc;
        nns += ns;
        nnf += nf;
        nnm += nm;
        if ( (abs(info.imprim) > 4 || info.ddebug) && ns+nc+nf+nm > 0 )
            fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved\n",ns,nc,nf,nm);
        if ( ns < 10 && abs(nc-ns) < 3 )  break;
        else if ( it > 3 && abs(nc-ns) < 0.3 * MS_MAX(nc,ns) )  break;
    }
    while( ++it < maxit && nc+ns > 0 );

    nm = 0;
    nm = movtri(mesh,met,5);
    if ( nm < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh.\n");
        return(0);
    }
    nnm += nm;

    if ( abs(info.imprim) < 5 && (nnc > 0 || nns > 0) )
        fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved, %d iter. \n",nns,nnc,nnf,nnm,it);
    return(1);
}

/* analyze tetrahedra and split if needed */
static int anatri(pMesh mesh,pSol met,char typchk) {
    int     nc,ns,nf,nnc,nns,nnf,it,maxit;

    /* analyze tetras : initial splitting */
    nns = nnc = nnf = it = 0;
    maxit = 5;
    do {
        /* memory free */
        free(mesh->adja);
        mesh->adja = 0;

        /* analyze surface */
        ns = anaelt(mesh,met,typchk);
        if ( ns < 0 ) {
            fprintf(stdout,"  ## Unable to complete surface mesh. Exit program.\n");
            return(0);
        }

        if ( !hashTria(mesh) ) {
            fprintf(stdout,"  ## Hashing problem. Exit program.\n");
            return(0);
        }

        /* collapse short edges */
        nc = colelt(mesh,met,typchk);
        if ( nc < 0 ) {
            fprintf(stdout,"  ## Unable to collapse mesh. Exiting.\n");
            return(0);
        }

        /* attempt to swap */
        nf = swpmsh(mesh,met,typchk);
        if ( nf < 0 ) {
            fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
            return(0);
        }

        nnc += nc;
        nns += ns;
        nnf += nf;
        if ( (abs(info.imprim) > 4 || info.ddebug) && ns+nc > 0 )
            fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped\n",ns,nc,nf);
        if ( it > 3 && abs(nc-ns) < 0.1 * MS_MAX(nc,ns) )  break;
    }
    while ( ++it < maxit && ns+nc+nf > 0 );

    if ( (abs(info.imprim) < 5 || info.ddebug ) && nns+nnc > 0 )
        fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %d iter.\n",nns,nnc,nnf,it);

    return(1);
}

int mmgs1(pMesh mesh,pSol met) {
    if ( abs(info.imprim) > 4 )
        fprintf(stdout,"  ** MESH ANALYSIS\n");

    /*delref(mesh);
      setref(mesh,1566,6,1);
      setref(mesh,1016,6,1);
      setref(mesh,5120,4,1);
      return(1);*/

    /*--- stage 1: geometric mesh */
    if ( abs(info.imprim) > 4 || info.ddebug )
        fprintf(stdout,"  ** GEOMETRIC MESH\n");

    if ( !anatri(mesh,met,1) ) {
        fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
        return(0);
    }

    /*--- stage 2: computational mesh */
    if ( abs(info.imprim) > 4 || info.ddebug )
        fprintf(stdout,"  ** COMPUTATIONAL MESH\n");

    /* define metric map */
    if ( !defsiz(mesh,met) ) {
        fprintf(stdout,"  ## Metric undefined. Exit program.\n");
        return(0);
    }
    if ( info.hgrad > 0. && !gradsiz(mesh,met) ) {
        fprintf(stdout,"  ## Gradation problem. Exit program.\n");
        return(0);
    }

    if ( !anatri(mesh,met,2) ) {
        fprintf(stdout,"  ## Unable to proceed adaptation. Exit program.\n");
        return(0);
    }

    if ( !adptri(mesh,met) ) {
        fprintf(stdout,"  ## Unable to adapt. Exit program.\n");
        return(0);
    }

    return(1);
}

