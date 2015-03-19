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
 * \file mmg3d/bucket.c
 * \brief Functions for bucket computations in Delaunay mode.
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \remark Delaunay mode only (\a PATTERN flag set to \a OFF).
 * \todo doxygen documentation.
 */

#include "mmg3d.h"
#define PRECI 1
#define LFILT    0.2//0.7

/* create bucket structure and store initial vertices */
_MMG5_pBucket _MMG5_newBucket(MMG5_pMesh mesh,int nmax) {
    MMG5_pPoint        ppt;
    _MMG5_pBucket       bucket;
    double        dd;
    int           k,ic,ii,jj,kk;

    /* memory alloc */
    _MMG5_ADD_MEM(mesh,sizeof(_MMG5_Bucket),"bucket",return(NULL));
    _MMG5_SAFE_MALLOC(bucket,1,_MMG5_Bucket);
    bucket->size = nmax;
    _MMG5_ADD_MEM(mesh,(nmax*nmax*nmax+1)*sizeof(int),"bucket->head",return(NULL));
    _MMG5_SAFE_CALLOC(bucket->head,nmax*nmax*nmax+1,int);
    _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(int),"bucket->link",return(NULL));
    _MMG5_SAFE_CALLOC(bucket->link,mesh->npmax+1,int);

    /* insert vertices */
    dd = nmax / (double)PRECI;
    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) )  continue;
        if (ppt->tag & MG_BDY) continue;
        ii = MG_MAX(0,(int)(dd * ppt->c[0])-1);
        jj = MG_MAX(0,(int)(dd * ppt->c[1])-1);
        kk = MG_MAX(0,(int)(dd * ppt->c[2])-1);
        ic = (kk*nmax + jj)*nmax + ii;

        if ( !bucket->head[ic] )
            bucket->head[ic] = k;
        else {
            bucket->link[k]  = bucket->head[ic];
            bucket->head[ic] = k;
        }
    }

    return(bucket);
}


/* check and eventually insert vertex */
int _MMG5_buckin_ani(MMG5_pMesh mesh,MMG5_pSol sol,_MMG5_pBucket bucket,int ip) {
    MMG5_pPoint        ppt,pp1;
    double        dd,d2,det,ux,uy,uz,dmi,m1,m2,m3,dx,dy,dz;
    double        *ma,*mb;
    int           i,j,k,ii,jj,kk,ic,icc,siz,ip1;
    int           iadr,imin,imax,jmin,jmax,kmin,kmax;

    ppt = &mesh->point[ip];
    siz = bucket->size;
    dd  = siz / (double)PRECI;

    iadr = (ip-1)*sol->size + 1;
    ma   = &sol->m[iadr];
    dmi  = LFILT*LFILT;

    ii = MG_MAX(0,(int)(dd * ppt->c[0])-1);
    jj = MG_MAX(0,(int)(dd * ppt->c[1])-1);
    kk = MG_MAX(0,(int)(dd * ppt->c[2])-1);
    ic = (kk*siz + jj)*siz + ii;

    /* check current cell */
    if ( bucket->head[ic] ) {
        ip1 = bucket->head[ic];
        pp1 = &mesh->point[ip1];
        ux = pp1->c[0] - ppt->c[0];
        uy = pp1->c[1] - ppt->c[1];
        uz = pp1->c[2] - ppt->c[2];
        d2 =      ma[0]*ux*ux + ma[3]*uy*uy + ma[5]*uz*uz \
            + 2.0*(ma[1]*ux*uy + ma[2]*ux*uz + ma[4]*uy*uz);
        if ( d2 < dmi ) {
            iadr = (ip1-1)*sol->size + 1;
            mb = &sol->m[iadr];
            d2 =      mb[0]*ux*ux + mb[3]*uy*uy + mb[5]*uz*uz \
                + 2.0*(mb[1]*ux*uy + mb[2]*ux*uz + mb[4]*uy*uz);
            if ( d2 < dmi )  return(0);
        }

        while ( bucket->link[ip1] ) {
            ip1 = bucket->link[ip1];
            pp1 = &mesh->point[ip1];
            ux = pp1->c[0] - ppt->c[0];
            uy = pp1->c[1] - ppt->c[1];
            uz = pp1->c[2] - ppt->c[2];
            d2 =      ma[0]*ux*ux + ma[3]*uy*uy + ma[5]*uz*uz \
                + 2.0*(ma[1]*ux*uy + ma[2]*ux*uz + ma[4]*uy*uz);
            if ( d2 < dmi ) {
                iadr = (ip1-1)*sol->size + 1;
                mb = &sol->m[iadr];
                d2 =      mb[0]*ux*ux + mb[3]*uy*uy + mb[5]*uz*uz \
                    + 2.0*(mb[1]*ux*uy + mb[2]*ux*uz + mb[4]*uy*uz);
                if ( d2 < dmi )  return(0);
            }
        }
    }

    /* bounding box */
    det = ma[0] * (ma[3]*ma[5] - ma[4]*ma[4]) \
        - ma[1] * (ma[1]*ma[5] - ma[2]*ma[4]) \
        + ma[2] * (ma[1]*ma[4] - ma[3]*ma[2]);
    det = 1.0 / det;
    m1 = ma[3]*ma[5] - ma[4]*ma[4];
    m2 = ma[0]*ma[5] - ma[2]*ma[2];
    m3 = ma[0]*ma[3] - ma[1]*ma[1];
    if ( det < 0.0 || m1 < 0.0 )
        return(1);
    else {
        dx = LFILT * sqrt(m1 * det) ;
        dy = LFILT * sqrt(m2 * det) ;
        dz = LFILT * sqrt(m3 * det) ;
    }

    imin = (int)(dd * (ppt->c[0]-dx))-1;
    jmin = (int)(dd * (ppt->c[1]-dy))-1;
    kmin = (int)(dd * (ppt->c[2]-dz))-1;
    imax = (int)(dd * (ppt->c[0]+dx))-1;
    jmax = (int)(dd * (ppt->c[1]+dy))-1;
    kmax = (int)(dd * (ppt->c[2]+dz))-1;

    imin = MG_MAX(0,MG_MIN(imin,siz-1));
    imax = MG_MIN(siz-1,MG_MAX(0,imax));
    jmin = MG_MAX(0,MG_MIN(jmin,siz-1));
    jmax = MG_MIN(siz-1,MG_MAX(0,jmax));
    kmin = MG_MAX(0,MG_MIN(kmin,siz-1));
    kmax = MG_MIN(siz-1,MG_MAX(0,kmax));
    if ( imin == imax && jmin == jmax && kmin == kmax )  return(1);

    /* explore neighbours */
    for (k=kmin; k<=kmax; k++)
        for (j=jmin; j<=jmax; j++)
            for (i=imin; i<=imax; i++) {
                icc = (k*siz + j)*siz + i;
                ip1 = bucket->head[icc];
                if ( !ip1 )  continue;
                pp1 = &mesh->point[ip1];
                ux = pp1->c[0] - ppt->c[0];
                uy = pp1->c[1] - ppt->c[1];
                uz = pp1->c[2] - ppt->c[2];
                d2 =      ma[0]*ux*ux + ma[3]*uy*uy + ma[5]*uz*uz \
                    + 2.0*(ma[1]*ux*uy + ma[2]*ux*uz + ma[4]*uy*uz);
                if ( d2 < dmi ) {
                    iadr = (ip1-1)*sol->size + 1;
                    mb = &sol->m[iadr];
                    d2 =      mb[0]*ux*ux + mb[3]*uy*uy + mb[5]*uz*uz \
                        + 2.0*(mb[1]*ux*uy + mb[2]*ux*uz + mb[4]*uy*uz);
                    if ( d2 < dmi )  return(0);
                }

                while ( bucket->link[ip1] ) {
                    ip1 = bucket->link[ip1];
                    pp1 = &mesh->point[ip1];
                    ux = pp1->c[0] - ppt->c[0];
                    uy = pp1->c[1] - ppt->c[1];
                    uz = pp1->c[2] - ppt->c[2];
                    d2 =      ma[0]*ux*ux + ma[3]*uy*uy + ma[5]*uz*uz \
                        + 2.0*(ma[1]*ux*uy + ma[2]*ux*uz + ma[4]*uy*uz);
                    if ( d2 < dmi ) {
                        iadr = (ip1-1)*sol->size + 1;
                        mb = &sol->m[iadr];
                        d2 =      mb[0]*ux*ux + mb[3]*uy*uy + mb[5]*uz*uz \
                            + 2.0*(mb[1]*ux*uy + mb[2]*ux*uz + mb[4]*uy*uz);
                        if ( d2 < dmi )  return(0);
                    }
                }
            }

    return(1);
}


int _MMG5_buckin_iso(MMG5_pMesh mesh,MMG5_pSol sol,_MMG5_pBucket bucket,int ip) {
    MMG5_pPoint        ppt,pp1;
    double        dd,d2,ux,uy,uz,hpi,hp1,hp2;
    int           i,j,k,ii,jj,kk,ic,icc,siz,ip1;
    int           imin,imax,jmin,jmax,kmin,kmax;

    ppt = &mesh->point[ip];
    siz = bucket->size;
    dd  = siz / (double)PRECI;
    hpi = LFILT * sol->m[ip];
    hp1 = hpi*hpi;

    ii = MG_MAX(0,(int)(dd * ppt->c[0])-1);
    jj = MG_MAX(0,(int)(dd * ppt->c[1])-1);
    kk = MG_MAX(0,(int)(dd * ppt->c[2])-1);
    ic = (kk*siz + jj)*siz + ii;

    /* check current cell */
    if ( bucket->head[ic] ) {
        ip1 = bucket->head[ic];
        pp1 = &mesh->point[ip1];
        hp2 = LFILT * sol->m[ip1];
        ux = pp1->c[0] - ppt->c[0];
        uy = pp1->c[1] - ppt->c[1];
        uz = pp1->c[2] - ppt->c[2];
        d2 = ux*ux + uy*uy + uz*uz;
        if ( d2 < hp1 || d2 < hp2*hp2 )  {
            //printf("filtre current %d : %e %e %e %e\n",ip1,d2,hp1,d2,hp2*hp2);
            return(0);
        }

        while ( bucket->link[ip1] ) {
            ip1 = bucket->link[ip1];
            pp1 = &mesh->point[ip1];
            hp2 = LFILT * sol->m[ip1];
            ux = pp1->c[0] - ppt->c[0];
            uy = pp1->c[1] - ppt->c[1];
            uz = pp1->c[2] - ppt->c[2];
            d2 = ux*ux + uy*uy + uz*uz;
            if ( d2 < hp1 || d2 < hp2*hp2 )  {
                //printf("filtre link %d : %e %e %e %e\n",ip1,d2,hp1,d2,hp2*hp2);
                return(0);
            }
        }
    }

    /* explore neighbors */
    imin = (int)(dd * (ppt->c[0]-hpi))-1;
    jmin = (int)(dd * (ppt->c[1]-hpi))-1;
    kmin = (int)(dd * (ppt->c[2]-hpi))-1;
    imax = (int)(dd * (ppt->c[0]+hpi))-1;
    jmax = (int)(dd * (ppt->c[1]+hpi))-1;
    kmax = (int)(dd * (ppt->c[2]+hpi))-1;

    imin = MG_MAX(0,MG_MIN(imin,siz-1));
    imax = MG_MIN(siz-1,MG_MAX(0,imax));
    jmin = MG_MAX(0,MG_MIN(jmin,siz-1));
    jmax = MG_MIN(siz-1,MG_MAX(0,jmax));
    kmin = MG_MAX(0,MG_MIN(kmin,siz-1));
    kmax = MG_MIN(siz-1,MG_MAX(0,kmax));
    if ( imin == imax && jmin == jmax && kmin == kmax )  return(1);

    for (k=kmin; k<=kmax; k++)
        for (j=jmin; j<=jmax; j++)
            for (i=imin; i<=imax; i++) {
                icc = (k*siz + j)*siz + i;
                ip1 = bucket->head[icc];
                if ( !ip1 )  continue;
                pp1 = &mesh->point[ip1];
                hp2 = LFILT * sol->m[ip1];
                ux = pp1->c[0] - ppt->c[0];
                uy = pp1->c[1] - ppt->c[1];
                uz = pp1->c[2] - ppt->c[2];
                d2 = ux*ux + uy*uy + uz*uz;
                if ( d2 < hp1 || d2 < hp2*hp2 ) {
                    /*     printf("other cell %d %e < %e -- %e < %e \n",ip1,d2,MMG_length(mesh,sol,ip,ip1),d2,hp2*hp2);
                           printf("on filtre avec %d : %e %e %e\n",ip1,pp1->c[0],pp1->c[1],pp1->c[2]);
                    */ return(0);
                }

                while ( bucket->link[ip1] ) {
                    ip1 = bucket->link[ip1];
                    pp1 = &mesh->point[ip1];
                    hp2 = LFILT * sol->m[ip1];
                    ux = pp1->c[0] - ppt->c[0];
                    uy = pp1->c[1] - ppt->c[1];
                    uz = pp1->c[2] - ppt->c[2];
                    d2 = ux*ux + uy*uy + uz*uz;
                    if ( d2 < hp1 || d2 < hp2*hp2 )  {
                        //      printf("link cell %d %e < %e -- %e < %e \n",ip1,d2,hp1,d2,hp2*hp2);
                        return(0);
                    }
                }
            }

    return(1);
}


int _MMG5_addBucket(MMG5_pMesh mesh,_MMG5_pBucket bucket,int ip) {
    MMG5_pPoint        ppt;
    double        dd;
    int           ic,ii,jj,kk,siz;

    ppt = &mesh->point[ip];
    siz = bucket->size;
    dd  = siz / (double)PRECI;

    ii = MG_MAX(0,(int)(dd * ppt->c[0])-1);
    jj = MG_MAX(0,(int)(dd * ppt->c[1])-1);
    kk = MG_MAX(0,(int)(dd * ppt->c[2])-1);
    ic = (kk*siz + jj)*siz + ii;

    /* store new point */
    if ( !bucket->head[ic] ) {
        bucket->head[ic] = ip;
        bucket->link[ip] = 0;
    }
    else {
        //  assert(!bucket->link[ip]);
        bucket->link[ip] = bucket->head[ic];
        bucket->head[ic] = ip;
        assert(ip!=bucket->link[ip]);
    }

    return(1);
}


int _MMG5_delBucket(MMG5_pMesh mesh,_MMG5_pBucket bucket,int ip) {
    MMG5_pPoint        ppt;
    double        dd;
    int           ic,ii,jj,kk,siz,ip1;

    ppt = &mesh->point[ip];
    siz = bucket->size;
    dd  = siz / (double)PRECI;

    ii = MG_MAX(0,(int)(dd * ppt->c[0])-1);
    jj = MG_MAX(0,(int)(dd * ppt->c[1])-1);
    kk = MG_MAX(0,(int)(dd * ppt->c[2])-1);
    ic = (kk*siz + jj)*siz + ii;

    /* remove vertex from cell */
    if ( bucket->head[ic] ) {
        if ( bucket->head[ic] == ip ) {
            bucket->head[ic] = bucket->link[ip];
            bucket->link[ip] = 0;
        }
        else {
            ip1 = bucket->head[ic];
            while ( ip1 && bucket->link[ip1] != ip ) {
                ip1 = bucket->link[ip1];
            }
            if ( bucket->link[ip1] == ip ) {
                bucket->link[ip1] = bucket->link[ip];
                bucket->link[ip] = 0;
            } else {
                printf("point non trouve %d %c -- %d\n",ip,mesh->point[ip].tag,ic);
            }
        }
    }

    return(1);
}
