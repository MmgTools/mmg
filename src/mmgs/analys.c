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
 * \file mmgs/analys.c
 * \brief Mesh analysis.
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

/* topology: set adjacent, detect Moebius, flip faces, count connected comp. */
static int setadj(pMesh mesh){
    pTria   pt,pt1;
    int    *adja,*adjb,adji1,adji2,*pile,iad,ipil,ip1,ip2,gen;
    int     k,kk,iel,jel,nf,nr,nt,nre,ncc,ned,ref;
    char    i,ii,i1,i2,ii1,ii2,tag,voy;

    if ( abs(info.imprim) > 5  || info.ddebug )
        fprintf(stdout,"  ** SETTING TOPOLOGY\n");

    pile = (int*)malloc((mesh->nt+1)*sizeof(int));
    assert(pile);
    pile[1] = 1;
    ipil    = 1;
    nre = nr = nf = nt = ncc = ned = 0;

    while ( ipil > 0 ) {
        ncc++;

        do {
            k  = pile[ipil--];
            pt = &mesh->tria[k];
            pt->flag = 0;
            if ( !MS_EOK(pt) )  continue;

            pt->cc = ncc;
            adja = &mesh->adja[3*(k-1)+1];
            nt++;
            for (i=0; i<3; i++) {
                i1  = inxt[i];
                i2  = iprv[i];
                ip1 = pt->v[i1];
                ip2 = pt->v[i2];

                if ( MS_EDG(pt->tag[i]) ) {
                    mesh->point[ip1].tag |= pt->tag[i];
                    mesh->point[ip2].tag |= pt->tag[i];
                }

                /* open boundary */
                if ( !adja[i] ) {
                    pt->tag[i] |= MS_GEO;
                    mesh->point[ip1].tag |= MS_GEO;
                    mesh->point[ip2].tag |= MS_GEO;
                    nr++;
                    ned++;
                    continue;
                }
 
                kk = adja[i] / 3;
                ii = adja[i] % 3;
                if ( kk > k )  ned++;

                /* correct edge tag */
                pt1 = &mesh->tria[kk];
                if ( pt->tag[i] & MS_NOM && !(pt1->tag[ii] & MS_NOM) ) {
                    pt1->tag[ii] = pt->tag[i];
                    pt1->edg[ii] = pt->edg[i];
                    mesh->point[ip1].tag |= MS_NOM;
                    mesh->point[ip2].tag |= MS_NOM;
                }
                if ( pt1->cc > 0 )  continue;

                if ( abs(pt1->ref) != abs(pt->ref) ) {
                    pt->tag[i]   |= MS_REF;
                    pt1->tag[ii] |= MS_REF;
                    mesh->point[ip1].tag |= MS_REF;
                    mesh->point[ip2].tag |= MS_REF;
                    nre++;
                }

                /* store adjacent */ 
                if ( !pt1->flag ) {
                    pt1->flag    = 1;
                    pile[++ipil] = kk;
                }

                /* check orientation */
                ii1 = inxt[ii];
                ii2 = iprv[ii];
                if ( pt1->v[ii1] == ip1 ) {
                    /* Moebius strip */
                    if ( pt1->base < 0 ) {
                        pt1->ref      = -1;
                        pt->tag[i]   |= MS_REF;
                        pt1->tag[ii] |= MS_REF;
                        nre++;
                    }
                    /* flip orientation */
                    else if ( !(pt->tag[i] & MS_NOM) ) {
                        pt1->base   = -pt1->base;
                        pt1->v[ii1] = ip2;
                        pt1->v[ii2] = ip1;

                        /* update adj */
                        iad   = 3*(kk-1)+1;
                        adjb  = &mesh->adja[iad];
                        adji1 = mesh->adja[iad+ii1];
                        adji2 = mesh->adja[iad+ii2];
                        adjb[ii1] = adji2;
                        adjb[ii2] = adji1;

                        /* modif tag + ref */
                        tag = pt1->tag[ii1];
                        pt1->tag[ii1] = pt1->tag[ii2];
                        pt1->tag[ii2] = tag;
                        ref = pt1->edg[ii1];
                        pt1->edg[ii1] = pt1->edg[ii2];
                        pt1->edg[ii2] = ref;

                        /* modif voyeurs */
                        if ( adjb[ii1] ) {
                            iel = adjb[ii1] / 3;
                            voy = adjb[ii1] % 3;
                            mesh->adja[3*(iel-1)+1+voy] = 3*kk + ii1;
                        }
                        if ( adjb[ii2] ) {
                            iel = adjb[ii2] / 3;
                            voy = adjb[ii2] % 3;
                            mesh->adja[3*(iel-1)+1+voy] = 3*kk + ii2;
                        }
                        nf++;
                    }
                }
            }
        }
        while ( ipil > 0 );

        /* find next triangle */
        ipil = 0;
        for (kk=1; kk<=mesh->nt; kk++) {
            pt = &mesh->tria[kk];
            if ( pt->v[0] && (pt->cc == 0) ) {
                ipil = 1;
                pile[ipil] = kk;
                pt->flag   = 1;
                break;
            }
        }
    }

    /* bilan */
    nr = nre = 0;
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) )  continue;
        for (i=0; i<3; i++) {
            if ( !MS_EDG(pt->tag[i]) )  continue;

            adja = &mesh->adja[3*(k-1)+1];
            jel  = adja[i] / 3;
            if ( !jel || jel > k ) {
                if ( pt->tag[i] & MS_GEO )  nr++;
                if ( pt->tag[i] & MS_REF )  nre++;
            }
        }
    }

    info.ncc = ncc;
    if ( info.ddebug ) {
        fprintf(stdout,"  a- ridges: %d found.\n",nr);
        fprintf(stdout,"  a- connex: %d connected component(s)\n",ncc);
        fprintf(stdout,"  a- orient: %d flipped\n",nf);
    }
    else if ( abs(info.imprim) > 4 ) {
        gen = (2 - mesh->np + ned - nt) / 2;
        if ( !info.mani )  fprintf(stdout,"  ## [non-manifold model]\n");
        fprintf(stdout,"     Connected component: %d,  genus: %d,   reoriented: %d\n",ncc,gen,nf);
        fprintf(stdout,"     Edges: %d,  tagged: %d,  ridges: %d,  refs: %d\n",ned,nr+nre,nr,nre);
    }

    free(pile);
    return(1);
}

/* Detect non manifold points */
static void nmpoints(pMesh mesh) {
    pTria      pt;
    pPoint     p0;
    int        k,np,numt,iel,jel,nmp,*adja;
    char       i0,i1,i,jp;
  
    nmp = 0;
    /* Initialize point flags */
    for (k=1; k<=mesh->np; k++)
        mesh->point[k].s = 0;
  
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) )  continue; 

        for (i=0; i<3; i++) {
            np = pt->v[i];
            p0 = &mesh->point[np];
            if ( !p0->s )  p0->s = k;
        }  
    }
          
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) )  continue; 
    
        for (i=0; i<3; i++) {
            np = pt->v[i];
            p0 = &mesh->point[np];
            numt = p0->s;
            if ( k == numt )  continue;
            jel = k;
            jp  = i;
            
            /* Unfold ball of np, until numt is reached */
            do {
                iel = jel;
                i0  =  jp;
                i1  = inxt[i0];
                adja = &mesh->adja[3*(iel-1)+1];
                jel = adja[i1] / 3;
                jp  = adja[i1] % 3;
                jp  = inxt[jp];
            }
            while ( jel && (jel != numt) && (jel !=k) );

            /* Ball has been completely traveled without meeting numt */
            if ( jel == k ) {
                if ( !(p0->tag & MS_CRN) || !(p0->tag & MS_REQ) ) {
                    nmp++;
                    p0->tag |= MS_CRN + MS_REQ;
                }
                continue;
            }
            else if ( jel == numt ) 
                continue;
  
            jel = iel;
            jp = i0;
              
            /* At this point, jel =0, i.e. an open boundary has been hit : travel in the opposite sense */
            do {
                iel = jel;
                i0  =  jp;
                i1  = iprv[i0];
                adja = &mesh->adja[3*(iel-1)+1];
                jel = adja[i1] / 3;
                jp  = adja[i1] % 3;
                jp  = iprv[jp];
            }
            while ( jel && (jel != numt));  
          
            if ( jel != numt) {
                if ( !(p0->tag & MS_CRN) || !(p0->tag & MS_REQ) ) {
                    p0->tag |= MS_CRN + MS_REQ; 
                    nmp++;
                }
            }
        }
    }  

    /* reset point flags */
    for (k=1; k<=mesh->np; k++)
        mesh->point[k].s = 0;

    if ( nmp && abs(info.imprim) > 4 )
        fprintf(stdout,"  ## %d non manifold points detected\n",nmp);
}

/* improve badly shaped elts */
static int delbad(pMesh mesh) {
    pTria    pt;
    pPoint   p[3];
    double   s,kal,declic,ux,uy,uz,vx,vy,vz;
    int     *adja,k,iel,nd,ndd,it;
    char     i,ia,i1,i2,j,typ;

    it = ndd = 0;
    declic = BADKAL / ALPHAD;

    do {
        nd = 0;
        for (k=1; k<=mesh->nt; k++) {
            pt = &mesh->tria[k];
            if ( !MS_EOK(pt) )  continue;

            kal = calelt(mesh,NULL,k);
            if ( kal > declic )  continue;

            p[0] = &mesh->point[pt->v[0]];
            p[1] = &mesh->point[pt->v[1]];
            p[2] = &mesh->point[pt->v[2]];
            adja = &mesh->adja[3*(k-1)+1];
            typ  = typelt(p,&ia);
      
            /* needle */
            if ( typ == 1 ) {
                if ( litcol(mesh,k,ia,kal) ) {
                    nd++;
                    continue;
                }
            }
            /* obtuse */
            else if ( typ == 2 ) {
                /* delete boundary elt */
                if ( !adja[ia] ) {
                    /* update point coordinates on ridge */
                    i1 = inxt[ia];
                    i2 = iprv[ia];
                    p[0] = &mesh->point[pt->v[ia]];
                    p[1] = &mesh->point[pt->v[i1]];
                    p[2] = &mesh->point[pt->v[i2]];
                    ux = p[2]->c[0] - p[1]->c[0];
                    uy = p[2]->c[1] - p[1]->c[1];
                    uz = p[2]->c[2] - p[1]->c[2];
                    vx = p[0]->c[0] - p[1]->c[0];
                    vy = p[0]->c[1] - p[1]->c[1];
                    vz = p[0]->c[2] - p[1]->c[2];
                    s  = (ux*vx + uy*vy + uz*vz) / sqrt(ux*ux + uy*uy + uz*uz);
                    p[0]->c[0] = vx - s*ux;
                    p[0]->c[1] = vy - s*uy;
                    p[0]->c[2] = vz - s*uz;
          
                    delElt(mesh,k);
                    nd++;
                    continue;
                }
                if ( litswp(mesh,k,ia,kal) || litcol(mesh,k,ia,kal) ) {
                    nd++;
                    continue;
                }
            }

            /* brute force to improve */
            for (i=0; i<3; i++) {
                if ( litswp(mesh,k,i,kal) || litcol(mesh,k,i,kal) ) {
                    nd++;
                    break;
                }
                else if ( adja[i] ) {
                    iel = adja[i] / 3;
                    j   = adja[i] % 3;
                    if ( litcol(mesh,iel,j,kal) ) {
                        nd++;
                        break;
                    }
                }
            }
        }
        ndd += nd;
        if ( nd && (info.ddebug || info.imprim < 0) )  fprintf(stdout,"     %d improved\n",nd);
    }
    while ( nd > 0 && ++it < 5 );

    if ( abs(info.imprim) > 4 )
        fprintf(stdout,"     %d bad elements improved\n",ndd);

    return(1);
}


/* check for ridges: dihedral angle */
static int setdhd(pMesh mesh) {
    pTria    pt,pt1;
    double   n1[3],n2[3],dhd;
    int     *adja,k,kk,nr;
    char     i,ii,i1,i2;

    nr = 0;
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) )  continue;

        nortri(mesh,pt,n1);
        adja = &mesh->adja[3*(k-1)+1];
        for (i=0; i<3; i++) {
            if ( pt->tag[i] & MS_GEO )  continue;
            kk = adja[i] / 3;
            ii = adja[i] % 3;

            /* check angle w. neighbor */
            if ( k < kk ) {
                pt1 = &mesh->tria[kk];
                nortri(mesh,pt1,n2);
                dhd = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
                if ( dhd <= info.dhd ) {
                    pt->tag[i]   |= MS_GEO;
                    pt1->tag[ii] |= MS_GEO;
                    i1 = inxt[i];
                    i2 = inxt[i1];
                    mesh->point[pt->v[i1]].tag |= MS_GEO;
                    mesh->point[pt->v[i2]].tag |= MS_GEO;
                    nr++;
                }
            }
        }
    }

    if ( abs(info.imprim) > 4 && nr > 0 )
        fprintf(stdout,"     %d ridges updated\n",nr);

    return(1);
}

/* check for singularities (corners) */
int singul(pMesh mesh) {
    pTria     pt;
    pPoint    ppt,p1,p2;
    double    ux,uy,uz,vx,vy,vz,dd;
    int       list[LMAX+2],k,nc,ng,nr,ns,nre;
    char      i;

    nre = nc = 0;
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) )  continue;

        for (i=0; i<3; i++) {
            ppt = &mesh->point[pt->v[i]];
            ppt->s++;
            if ( !MS_VOK(ppt) || MS_SIN(ppt->tag) )  continue;
            else if ( MS_EDG(ppt->tag) ) {
                ns = bouler(mesh,k,i,list,&ng,&nr);
                if ( !ns ) continue;
        
                if ( (ng + nr) > 2 ) {
                    ppt->tag |= MS_CRN + MS_REQ; 
                    nre++;   
                    nc++;
                } 
                else if ( ng == 1 && nr == 1 ) {
                    ppt->tag |= MS_REQ; 
                    nre++;
                } 
                else if ( ng == 1 && !nr ) {
                    ppt->tag |= MS_CRN + MS_REQ; 
                    nre++;
                    nc++;
                }
                else if ( nr == 1 && !ng ) {
                    ppt->tag |= MS_CRN + MS_REQ; 
                    nre++;
                    nc++;
                }
     
                /* check ridge angle */
                else {
                    p1 = &mesh->point[list[1]];
                    p2 = &mesh->point[list[2]];
                    ux = p1->c[0] - ppt->c[0];
                    uy = p1->c[1] - ppt->c[1];
                    uz = p1->c[2] - ppt->c[2];
                    vx = p2->c[0] - ppt->c[0];
                    vy = p2->c[1] - ppt->c[1];
                    vz = p2->c[2] - ppt->c[2];
                    dd = (ux*ux + uy*uy + uz*uz) * (vx*vx + vy*vy + vz*vz);
                    if ( fabs(dd) > EPSD ) {
                        dd = (ux*vx + uy*vy + uz*vz) / sqrt(dd);
                        if ( dd > -info.dhd ) {
                            ppt->tag |= MS_CRN;
                            nc++;
                        }
                    }
                }
            }
        }
    }

    /* check for handle */
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) )  continue;
        for (i=0; i<3; i++) {
            ppt = &mesh->point[pt->v[i]];
            if ( !ppt->s )  continue;
            nr = boulet(mesh,k,i,list);
            if ( nr != ppt->s ) {
                ppt->tag |= MS_CRN + MS_REQ;
                ppt->s = 0;
                nc++;
            }
        }
    }

    if ( abs(info.imprim) > 4 && nre > 0 )
        fprintf(stdout,"     %d corners, %d singular points detected\n",nc,nre);
    return(1);
}


/* compute normals at C1 vertices, for C0: tangents */
static int norver(pMesh mesh) {
    pTria     pt;
    pPoint    ppt;
    pGeom     go;
    double    n[3],dd;
    int      *adja,k,kk,ier,ng,nn,nt,nf;
    char      i,ii,i1,i2;

    if ( abs(info.imprim) > 4 || info.ddebug )
        fprintf(stdout,"  ** DEFINING GEOMETRY\n"); 

    /* 1. process C1 vertices, normals */
    nn = ng = nt = nf = 0;
    ++mesh->base;
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) )  continue;

        for (i=0; i<3; i++) {
            ppt = &mesh->point[pt->v[i]];
            if ( MS_SIN(ppt->tag) )  continue;
            else if ( MS_EDG(ppt->tag) ) {
                ng++;
                continue;
            }
            else if ( ppt->flag == mesh->base )  continue;
            else if ( mesh->nc1 )  continue;

            ier = boulen(mesh,k,i,ppt->n);
            if ( ier ) {
                ppt->flag = mesh->base;
                nn++;
            }
            else
                nf++;
        }
    }

    /* memory to store normals on both sides of ridges */
    if(!ng) {
        mesh->ngmax = NGMAX;
        mesh->geom  = (pGeom)calloc(mesh->ngmax+1,sizeof(Geom));
        assert(mesh->geom);

    } else if ( ng ) {
        mesh->ngmax = MS_MAX(1.5*ng,NGMAX);
        mesh->geom  = (pGeom)calloc(mesh->ngmax+1,sizeof(Geom));
        assert(mesh->geom);

        /* 2. process C0 vertices on curves, tangents */
        for (k=1; k<=mesh->nt; k++) {
            pt = &mesh->tria[k];
            if ( !MS_EOK(pt) )  continue;
      
            adja = &mesh->adja[3*(k-1)+1];
            for (i=0; i<3; i++) {
                i1  = inxt[i];
                i2  = inxt[i1];
                ppt = &mesh->point[pt->v[i]];

                if ( ppt->tag & MS_CRN || ppt->flag == mesh->base )  continue;
                else if ( !MS_EDG(pt->tag[i1]) )  continue;

                ier = boulen(mesh,k,i,n);
                if ( !ier )  continue;

                ++mesh->ng;
                assert(mesh->ng < mesh->ngmax);
                ppt->ig = mesh->ng;
                go = &mesh->geom[mesh->ng];
                memcpy(go->n1,n,3*sizeof(double));         

                /* compute n2 along ridge */
                if ( pt->tag[i1] & MS_GEO ) {
                    if ( adja[i1] ) {
                        kk  = adja[i1] / 3;
                        ii  = adja[i1] % 3;
                        ii  = inxt[ii];

                        ier = boulen(mesh,kk,ii,n);
                        if ( !ier )  continue;
                        memcpy(go->n2,n,3*sizeof(double));

                        /* compute tangent as intersection of n1 + n2 */
                        ppt->n[0] = go->n1[1]*go->n2[2] - go->n1[2]*go->n2[1];
                        ppt->n[1] = go->n1[2]*go->n2[0] - go->n1[0]*go->n2[2];
                        ppt->n[2] = go->n1[0]*go->n2[1] - go->n1[1]*go->n2[0];
                        ppt->flag = mesh->base;
                        dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
                        if ( dd > EPSD2 ) {
                            dd = 1.0 / sqrt(dd);
                            ppt->n[0] *= dd;
                            ppt->n[1] *= dd;
                            ppt->n[2] *= dd;
                        }
                        ++nt;
                        continue;
                    }
                }

                /* compute tgte */
                ier = boulec(mesh,k,i,ppt->n);
                if ( !ier )  continue;
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
                    ppt->flag = mesh->base;
                    ++nt;
                }
            }
        }
    }

    if ( abs(info.imprim) > 4 && nn+nt > 0 )
        fprintf(stdout,"     %d normals,  %d tangents updated  (%d failed)\n",nn,nt,nf);

    return(1);
}

/* regularization procedure for derivatives, dual Laplacian */
static int regnor(pMesh mesh) {
    pTria    pt;
    pPoint   ppt,p0;
    double  *tabl,n[3],lm1,lm2,dd,nx,ny,nz,res0,res;
    int      i,k,iad,it,nn,nit,iel,ilist,list[LMAX];

    /* assign seed to vertex */
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MS_EOK(pt) )  continue;
        for (i=0; i<3; i++) {
            ppt = &mesh->point[pt->v[i]];
            if ( !ppt->s )  ppt->s = k;
        }
    }

    /* allocate memory for normals */
    tabl = (double*)calloc(3*mesh->np+1,sizeof(double));
    assert(tabl);

    it   = 0;
    nit  = 2;
    res0 = 0.0;
    lm1  = 0.4;
    lm2  = 0.399;
    while ( it++ < nit ) {
        /* step 1: laplacian */
        for (k=1; k<=mesh->np; k++) {
            ppt = &mesh->point[k];
            if ( !MS_VOK(ppt) || ppt->tag > MS_REF )  continue;

            iel = ppt->s;
            assert(iel);
            pt = &mesh->tria[iel];
            i  = 0;
            if ( pt->v[1] == k )  i = 1;
            else if ( pt->v[2] == k ) i = 2;

            ilist = boulep(mesh,iel,i,list);  

            /* average normal */
            nx = ny = nz = 0.0;
            for (i=1; i<=ilist; i++) {
                p0  = &mesh->point[list[i]];
                if ( p0->tag > MS_REF )  continue;
                nx += p0->n[0];
                ny += p0->n[1];
                nz += p0->n[2];
            }
            dd  = nx*nx + ny*ny + nz*nz;
            if ( dd > EPSD2 ) {
                dd = 1.0 / sqrt(dd);
                nx *= dd;
                ny *= dd;
                nz *= dd;
            }

            /* Laplacian */
            iad = 3*(k-1)+1;
            tabl[iad+0] = ppt->n[0] + lm1 * (nx - ppt->n[0]);
            tabl[iad+1] = ppt->n[1] + lm1 * (ny - ppt->n[1]);
            tabl[iad+2] = ppt->n[2] + lm1 * (nz - ppt->n[2]);
        }

        /* step 2: anti-laplacian */
        res = 0;
        nn  = 0;
        for (k=1; k<=mesh->np; k++) {
            ppt = &mesh->point[k];
            if ( !MS_VOK(ppt) || ppt->tag > MS_REF )  continue;

            iel = ppt->s;
            assert(iel);
            pt = &mesh->tria[iel];
            i = 0;
            if ( pt->v[1] == k )  i = 1;
            else if ( pt->v[2] == k ) i = 2;

            ilist = boulep(mesh,iel,i,list);  

            /* average normal */
            nx = ny = nz = 0.0;
            for (i=1; i<=ilist; i++) {
                iad = 3*(list[i]-1) + 1;
                nx += tabl[iad+0];
                ny += tabl[iad+1];
                nz += tabl[iad+2];
            }
            dd  = nx*nx + ny*ny + nz*nz;
            if ( dd > EPSD2 ) {
                dd = 1.0 / sqrt(dd);
                nx *= dd;
                ny *= dd;
                nz *= dd;
            }

            /* antiLaplacian */
            iad = 3*(k-1)+1;
            n[0] = tabl[iad+0] - lm2 * (nx - tabl[iad+0]);
            n[1] = tabl[iad+1] - lm2 * (ny - tabl[iad+1]);
            n[2] = tabl[iad+2] - lm2 * (nz - tabl[iad+2]);
            nn++;
            res += (ppt->n[0]-n[0])*(ppt->n[0]*n[0]) + (ppt->n[1]-n[1])*(ppt->n[1]*n[1]) + (ppt->n[2]-n[2])*(ppt->n[2]*n[2]); 
        }

        if ( it == 1 )  res0 = res;
        if ( res0 > EPSD )  res  = res / res0;
        if ( info.imprim < 0 || info.ddebug ) {
            fprintf(stdout,"     iter %5d  res %.3E\r",it,res); 
            fflush(stdout);
        }
        if ( it > 1 && res < EPS )  break;
    }
    if ( info.imprim < 0 || info.ddebug )  fprintf(stdout,"\n");

    if ( abs(info.imprim) > 4 )
        fprintf(stdout,"     %d normals regularized: %.3e\n",nn,res);

    free(tabl);
    return(1);
}


/* preprocessing stage: mesh analysis */
int analys(pMesh mesh) {

    /* create adjacency */
    if ( !hashTria(mesh) ) {
        fprintf(stdout,"  ## Hashing problem. Exit program.\n");
        return(0);
    }

    /* delete badly shaped elts */
    /*if ( info.badkal && !delbad(mesh) ) {
      fprintf(stdout,"  ## Geometry trouble. Exit program.\n");
      return(0);
      }*/

    /* identify connexity */
    if ( !setadj(mesh) ) {
        fprintf(stdout,"  ## Topology problem. Exit program.\n");
        return(0);
    }

    /* check for nomanifold point */
    nmpoints(mesh);

    /* check for ridges */
    if ( info.dhd > ANGLIM && !setdhd(mesh) ) {
        fprintf(stdout,"  ## Geometry problem. Exit program.\n");
        return(0);
    }

    /* identify singularities */
    if ( !singul(mesh) ) {
        fprintf(stdout,"  ## Singularity problem. Exit program.\n");
        return(0);
    }

    /* define normals */
    if ( !mesh->ng ) {
        if ( !norver(mesh) ) {
            fprintf(stdout,"  ## Normal problem. Exit program.\n");
            return(0);
        }
        /* regularize normals */
        if ( info.nreg && !regnor(mesh) ) {
            fprintf(stdout,"  ## Normal regularization problem. Exit program.\n");
            return(0);
        }
    }

    return(1);
}

