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
 * \file mmgs/colver.c
 * \brief Functions for vertices collapsing.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

extern Info info;

/* check if geometry preserved by collapsing edge i */
int chkcol(pMesh mesh,pSol met,int k,char i,int *list,char typchk) {
    pTria     pt,pt0,pt1,pt2;
    pPoint    p1,p2;
    double    len,lon,ps,cosnold,cosnnew,kal,n0old[3],n1old[3],n00old[3];
    double    n0new[3],n1new[3],n00new[3];
    int      *adja,jel,kel,ip1,ip2,l,ll,ilist;
    char      i1,i2,j,jj,j2,lj,open,voy;

    pt0 = &mesh->tria[0];
    pt  = &mesh->tria[k];
    i1  = inxt[i];
    i2  = iprv[i];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];
    if ( typchk == 2 && met->m ) {
        lon = lenedg(mesh,met,ip1,ip2,0);
        lon = MS_MIN(lon,LSHRT);
        lon = MS_MAX(1.0/lon,LLONG);
    }

    /* collect all triangles around vertex i1 */
    ilist = boulechknm(mesh,k,i1,list);
    if ( ilist <= 0 )  return(0);

    /* check for open ball */
    adja = &mesh->adja[3*(k-1)+1];
    open = adja[i] == 0;

    if ( ilist > 3 ) {
        /* check references */
        if ( MS_EDG(pt->tag[i2]) ) {
            jel = list[1] / 3;
            pt1 = &mesh->tria[jel];
            if ( abs(pt->ref) != abs(pt1->ref) )  return(0);
        }

        /* analyze ball */
        for (l=1; l<ilist-1+open; l++) {
            jel = list[l] / 3;
            j   = list[l] % 3;
            jj  = inxt[j];
            j2  = iprv[j];
            pt1 = &mesh->tria[jel];

            /* check length */
            if ( typchk == 2 && met->m && !MS_EDG(mesh->point[ip2].tag) ) {
                ip1 = pt1->v[j2];
                len = lenedg(mesh,met,ip1,ip2,0);
                if ( len > lon )  return(0);
            }

            /* check normal flipping */
            if ( !nortri(mesh,pt1,n1old) )  return(0);
            memcpy(pt0,pt1,sizeof(Tria));
            pt0->v[j] = ip2;
            if ( !nortri(mesh,pt0,n1new) )  return(0);

            ps = n1new[0]*n1old[0] + n1new[1]*n1old[1]  + n1new[2]*n1old[2];
            if ( ps < 0.0 )  return(0);

            /* keep normals at 1st triangles */
            if ( l == 1 && !open ) {
                memcpy(n00old,n1old,3*sizeof(double));
                memcpy(n00new,n1new,3*sizeof(double));
            }

            /* check normals deviation */
            if ( !(pt1->tag[j2] & MS_GEO) ) {
                if ( l > 1 ) {
                    cosnold = n0old[0]*n1old[0] + n0old[1]*n1old[1] + n0old[2]*n1old[2];
                    cosnnew = n0new[0]*n1new[0] + n0new[1]*n1new[1] + n0new[2]*n1new[2];
                    if ( cosnold < ANGEDG ) {
                        if ( cosnnew < cosnold )  return(0);
                    }
                    else if ( cosnnew < ANGEDG )  return(0);
                }
                memcpy(n0old,n1old,3*sizeof(double));
                memcpy(n0new,n1new,3*sizeof(double));
            }

            /* check geometric support */
            if ( l == 1 ) {
                pt0->tag[j2] = MS_MAX(pt0->tag[j2],pt->tag[i1]);
            }
            else if ( l == ilist-3+open ) {
                ll = list[ilist-2+open] / 3;
                if ( ll > mesh->nt )  return(0);
                lj = list[ilist-2+open] % 3;
                pt0->tag[jj] = MS_MAX(pt0->tag[jj],mesh->tria[ll].tag[lj]);
            }
            if ( chkedg(mesh,0) )  return(0);

            /* check quality */
            if ( typchk == 2 && met->m )
                kal = ALPHAD*calelt(mesh,met,0);
            else
                kal = ALPHAD*calelt_iso(mesh,0,0);
            if ( kal < NULKAL )  return(0);
        }

        /* check angle between 1st and last triangles */
        if ( !open && !(pt->tag[i] & MS_GEO) ) {
            cosnold = n00old[0]*n1old[0] + n00old[1]*n1old[1] + n00old[2]*n1old[2];
            cosnnew = n00new[0]*n1new[0] + n00new[1]*n1new[1] + n00new[2]*n1new[2];
            if ( cosnold < ANGEDG ) {
                if ( cosnnew < cosnold )  return(0);
            }
            else if ( cosnnew < ANGEDG )  return(0);

            /* other checks for reference collapse */
            jel = list[ilist-1] / 3;
            j   = list[ilist-1] % 3;
            j   = iprv[j];
            pt  = &mesh->tria[jel];
            if ( MS_EDG(pt->tag[j]) ) {
                jel = list[ilist-2] / 3;
                pt1 = &mesh->tria[jel];
                if ( abs(pt->ref) != abs(pt1->ref) )  return(0);
            }
        }
    }

    /* specific test: no collapse if any interior edge is EDG */
    else if ( ilist == 3 ) {
        p1 = &mesh->point[pt->v[i1]];
        if ( MS_SIN(p1->tag) )  return(0);
        else if ( MS_EDG(pt->tag[i2]) && !MS_EDG(pt->tag[i]) )  return(0);
        else if ( !MS_EDG(pt->tag[i2]) && MS_EDG(pt->tag[i]) )  return(0);
        else if ( MS_EDG(pt->tag[i2]) && MS_EDG(pt->tag[i]) && MS_EDG(pt->tag[i1]) )  return(0);

        /* Check geometric approximation */
        jel = list[1] / 3;
        j   = list[1] % 3;
        j2  = iprv[j];
        pt0 = &mesh->tria[0];
        pt1 = &mesh->tria[jel];
        memcpy(pt0,pt,sizeof(Tria));
        pt0->v[i1] = pt1->v[j2];
        if ( chkedg(mesh,0) )  return(0);
    }

    /* for specific configurations along open ridge */
    else if ( ilist == 2 ) {
        if ( !open )  return(0);
        jel = list[1] / 3;
        j   = list[1] % 3;

        /* Topological test */
        adja = &mesh->adja[3*(jel-1)+1];
        kel = adja[j] / 3;
        voy = adja[j] % 3;
        pt2 = &mesh->tria[kel];

        if ( pt2->v[voy] == ip2) return(0);

        jj  = inxt[j];
        pt1 = &mesh->tria[jel];
        if ( abs(pt->ref) != abs(pt1->ref) )  return(0);
        else if ( !(pt1->tag[jj] & MS_GEO) )  return(0);

        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt1->v[jj]];
        if ( p2->tag > p1->tag || p2->ref != p1->ref )  return(0);

        /* Check geometric approximation */
        j2  = iprv[j];
        pt0 = &mesh->tria[0];
        memcpy(pt0,pt,sizeof(Tria));
        pt0->v[i1] = pt1->v[j2];
        if ( chkedg(mesh,0) )  return(0);
    }

    return(ilist);
}

/* collapse edge i of k, i1->i2 */
int colver(pMesh mesh,int *list,int ilist) {
    pTria    pt,pt1,pt2;
    int     *adja,k,iel,jel,kel,ip1,ip2;
    char     i,i1,i2,j,jj,open;

    iel = list[0] / 3;
    i1  = list[0] % 3;
    i   = iprv[i1];
    i2  = inxt[i1];
    pt  = &mesh->tria[iel];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];

    /* check for open ball */
    adja = &mesh->adja[3*(iel-1)+1];
    open = adja[i] == 0;

    /* update vertex ip1 -> ip2 */
    for (k=1; k<ilist-1+open; k++) {
        jel = list[k] / 3;
        jj  = list[k] % 3;
        pt1 = &mesh->tria[jel];
        pt1->v[jj] = ip2;
        pt1->base  = mesh->base;
    }

    /* update adjacent with 1st elt */
    jel = list[1] / 3;
    jj  = list[1] % 3;
    j   = iprv[jj];
    pt1 = &mesh->tria[jel];
    pt1->tag[j] = MS_MAX(pt->tag[i1],pt1->tag[j]);
    pt1->edg[j] = MS_MAX(pt->edg[i1],pt1->edg[j]);
    if ( adja[i1] ) {
        kel = adja[i1] / 3;
        k   = adja[i1] % 3;
        mesh->adja[3*(kel-1)+1+k] = 3*jel + j;
        mesh->adja[3*(jel-1)+1+j] = 3*kel + k;
        pt2 = &mesh->tria[kel];
        pt2->tag[k] = MS_MAX(pt1->tag[j],pt2->tag[k]);
        pt2->edg[k] = MS_MAX(pt1->edg[j],pt2->edg[k]);
    }
    else
        mesh->adja[3*(jel-1)+1+j] = 0;

    /* adjacent with last elt */
    if ( !open ) {
        iel = list[ilist-1] / 3;
        i1  = list[ilist-1] % 3;
        pt  = &mesh->tria[iel];

        jel = list[ilist-2] / 3;
        jj  = list[ilist-2] % 3;
        j   = inxt[jj];
        pt1 = &mesh->tria[jel];
        pt1->tag[j] = MS_MAX(pt->tag[i1],pt1->tag[j]);
        pt1->edg[j] = MS_MAX(pt->edg[i1],pt1->edg[j]);
        adja = &mesh->adja[3*(iel-1)+1];
        if ( adja[i1] ) {
            kel = adja[i1] / 3;
            k   = adja[i1] % 3;
            mesh->adja[3*(kel-1)+1+k] = 3*jel + j;
            mesh->adja[3*(jel-1)+1+j] = 3*kel + k;
            pt2 = &mesh->tria[kel];
            pt2->tag[k] = MS_MAX(pt1->tag[j],pt2->tag[k]);
            pt2->edg[k] = MS_MAX(pt1->edg[j],pt2->edg[k]);
        }
        else
            mesh->adja[3*(jel-1)+1+j] = 0;
    }

    delPt(mesh,ip1);
    delElt(mesh,list[0] / 3);
    if ( !open )  delElt(mesh,list[ilist-1] / 3);

    return(1);
}


/* collapse point list[0]/3 in case ilist = 3 : point is removed */
int colver3(pMesh mesh,int* list) {
    pTria   pt,pt1,pt2;
    int    *adja,iel,jel,kel,mel,ip;
    char    i,i1,i2,j,j2,k,m;

    /* update of new point for triangle list[0] */
    iel = list[0] / 3;
    i   = list[0] % 3;
    i1  = inxt[i];
    i2  = iprv[i];
    pt  = &mesh->tria[iel];
    ip  = pt->v[i];

    jel = list[1] / 3;
    j   = list[1] % 3;
    j2  = iprv[j];
    pt1 = &mesh->tria[jel];

    kel = list[2] / 3;
    k   = list[2] % 3;
    pt2 = &mesh->tria[kel];

    /* update info */
    pt->v[i]    = pt1->v[j2];
    pt->tag[i]  = MS_MAX(pt->tag[i],pt->tag[i1]);
    pt->edg[i]  = MS_MAX(pt->edg[i],pt->edg[i1]);
    pt->tag[i1] = pt1->tag[j];
    pt->edg[i1] = pt1->edg[j];
    pt->tag[i2] = pt2->tag[k];
    pt->edg[i2] = pt2->edg[k];
    pt->base    = mesh->base;

    /* update neighbours of new triangle */
    adja = &mesh->adja[3*(iel-1)+1];
    adja[i1] = mesh->adja[3*(jel-1)+1+j];
    adja[i2] = mesh->adja[3*(kel-1)+1+k];
    mel = adja[i] / 3;
    m   = adja[i] % 3;
    pt1 = &mesh->tria[mel];
    pt1->tag[m] = pt->tag[i];
    pt1->edg[m] = pt->edg[i];

    /* update of neighbours of old neighbours */
    adja = &mesh->adja[3*(jel-1)+1];
    mel   = adja[j] / 3;
    m     = adja[j] % 3;
    mesh->adja[3*(mel-1)+1+m] = 3*iel + i1;

    adja = &mesh->adja[3*(kel-1)+1];
    mel  = adja[k] / 3;
    m    = adja[k] % 3;
    mesh->adja[3*(mel-1)+1+m] = 3*iel + i2;

    /* remove vertex + elements */
    delPt(mesh,ip);
    delElt(mesh,jel);
    delElt(mesh,kel);

    return(1);
}


/* collapse point along open ridge */
int colver2(pMesh mesh,int* list) {
    pTria   pt,pt1;
    int    *adja,iel,jel,kel,ip;
    char    i,i1,i2,jj,j2,k;

    /* update of new point for triangle list[0] */
    iel = list[0] / 3;
    i1  = list[0] % 3;
    i   = iprv[i1];
    i2  = inxt[i1];
    pt  = &mesh->tria[iel];
    ip  = pt->v[i1];

    jel = list[1] / 3;
    j2  = list[1] % 3;
    jj  = iprv[j2];
    pt1 = &mesh->tria[jel];

    /* update info */
    pt->v[i1] = pt1->v[jj];
    pt->tag[i2] = pt1->tag[j2];
    pt->edg[i2] = pt1->edg[j2];
    pt->base  = mesh->base;

    /* update neighbours of new triangle */
    adja = &mesh->adja[3*(iel-1)+1];
    adja[i2] = mesh->adja[3*(jel-1)+1+j2];
    adja = &mesh->adja[3*(jel-1)+1];
    kel  = adja[j2] / 3;
    k    = adja[j2] % 3;
    mesh->adja[3*(kel-1)+1+k] = 3*iel + i2;

    /* remove vertex + element */
    delPt(mesh,ip);
    delElt(mesh,jel);

    return(1);
}

/* collapse edge i of k, i1->i2 */
int litcol(pMesh mesh,int k,char i,double kali) {
    pTria     pt,pt0,pt1;
    pPoint    p1,p2;
    double    kal,ps,cosnold,cosnnew,n0old[3],n0new[3],n1old[3],n1new[3],n00old[3],n00new[3];
    int      *adja,list[LMAX+2],jel,ip1,ip2,l,ilist;
    char      i1,i2,j,jj,j2,open;

    pt0 = &mesh->tria[0];
    pt  = &mesh->tria[k];
    i1  = inxt[i];
    i2  = iprv[i];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];

    /* collect all triangles around vertex i1 */
    ilist = boulet(mesh,k,i1,list);

    /* check for open ball */
    adja = &mesh->adja[3*(k-1)+1];
    open = adja[i] == 0;

    if ( ilist > 3 ) {
        /* check references */
        jel = list[1] / 3;
        pt1 = &mesh->tria[jel];
        if ( abs(pt->ref) != abs(pt1->ref) )  return(0);

        /* analyze ball */
        for (l=1; l<ilist-1+open; l++) {
            jel = list[l] / 3;
            j   = list[l] % 3;
            jj  = inxt[j];
            j2  = iprv[j];
            pt1 = &mesh->tria[jel];

            /* check normal flipping */
            if ( !nortri(mesh,pt1,n1old) )  return(0);
            memcpy(pt0,pt1,sizeof(Tria));
            pt0->v[j] = ip2;
            if ( !nortri(mesh,pt0,n1new) )  return(0);
            ps = n1new[0]*n1old[0] + n1new[1]*n1old[1]  + n1new[2]*n1old[2];
            if ( ps < 0.0 )  return(0);

            /* keep normals at 1st triangles */
            if ( l == 1 && !open ) {
                memcpy(n00old,n1old,3*sizeof(double));
                memcpy(n00new,n1new,3*sizeof(double));
            }

            /* check normals deviation */
            if ( !(pt1->tag[j2] & MS_GEO) ) {
                if ( l > 1 ) {
                    cosnold = n0old[0]*n1old[0] + n0old[1]*n1old[1] + n0old[2]*n1old[2];
                    cosnnew = n0new[0]*n1new[0] + n0new[1]*n1new[1] + n0new[2]*n1new[2];
                    if ( cosnold < ANGEDG ) {
                        if ( cosnnew < MS_MIN(0.0,cosnold) )  return(0);
                    }
                    else if ( cosnnew < ANGEDG )  return(0);
                }

                memcpy(n0old,n1old,3*sizeof(double));
                memcpy(n0new,n1new,3*sizeof(double));
            }
            /* check quality */
            kal = ALPHAD*calelt_iso(mesh,0,0);
            if ( kal < NULKAL )  return(0);
        }

        /* check angle between 1st and last triangles */
        if ( !open ) {
            cosnold = n00old[0]*n1old[0] + n00old[1]*n1old[1] + n00old[2]*n1old[2];
            cosnnew = n00new[0]*n1new[0] + n00new[1]*n1new[1] + n00new[2]*n1new[2];
            if ( cosnold < ANGEDG ) {
                if ( cosnnew < MS_MIN(0.0,cosnold) )  return(0);
            }
            else if ( cosnnew < ANGEDG )  return(0);

            /* other reference checks */
            jel = list[ilist-1] / 3;
            pt  = &mesh->tria[jel];
            jel = list[ilist-2] / 3;
            pt1 = &mesh->tria[jel];
            if ( abs(pt->ref) != abs(pt1->ref) )  return(0);
        }

        return(colver(mesh,list,ilist));
    }

    /* specific test: no collapse if any interior edge is EDG */
    else if ( ilist == 3 ) {
        p1 = &mesh->point[pt->v[i1]];
        if ( MS_SIN(p1->tag) )  return(0);
        else if (  MS_EDG(pt->tag[i2]) && !MS_EDG(pt->tag[i]) )  return(0);
        else if ( !MS_EDG(pt->tag[i2]) &&  MS_EDG(pt->tag[i]) )  return(0);
        else if (  MS_EDG(pt->tag[i2]) &&  MS_EDG(pt->tag[i]) && MS_EDG(pt->tag[i1]) )  return(0);

        return(colver3(mesh,list));
    }

    /* for specific configurations along open ridge */
    else if ( ilist == 2 ) {
        if ( !open )  return(0);
        jel = list[1] / 3;
        j   = list[1] % 3;
        jj  = inxt[j];
        pt1 = &mesh->tria[jel];
        if ( abs(pt->ref) != abs(pt1->ref) )  return(0);
        else if ( !(pt1->tag[jj] & MS_GEO) )  return(0);

        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt1->v[jj]];
        if ( p2->tag > p1->tag || p2->ref != p1->ref )  return(0);

        return(colver2(mesh,list));
    }

    return(0);
}



