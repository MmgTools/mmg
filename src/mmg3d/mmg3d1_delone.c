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
 * \file mmg3d/mmg3d1_delone.c
 * \brief Perform volume and surface mesh adaptation in delaunay mode.
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 *
 * Perform volume and surface mesh adaptation in delaunay mode (\a
 * PATTERN preprocessor flag set to OFF).
 *
 * \todo Clean the boucle for (code copy...)
 */
#include "mmg3d.h"

char  ddb;

#define _MMG5_LOPTL_MMG5_DEL     1.41
#define _MMG5_LOPTS_MMG5_DEL     0.6

int MMG_npuiss,MMG_nvol,MMG_npres,MMG_npd;


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param bucket pointer toward the bucket structure.
 * \param ne number of elements.
 * \param ifilt pointer to store the number of vertices filtered by the bucket.
 * \param ns pointer to store the number of vertices insertions.
 * \param nc pointer to store the number of collapse.
 * \param warn pointer to store a flag that warn the user in case of
 * reallocation difficulty.
 * \param it iteration index.
 * \return -1 if fail and we don't save the mesh, 0 if fail but we try to save
 * the mesh, 1 otherwise.
 *
 * \ref adpsplcol loop: split edges longer than \ref _MMG5_LOPTL_MMG5_DEL and
 * collapse edges shorter than \ref _MMG5_LOPTS_MMG5_DEL.
 *
 */
static inline int
_MMG5_boucle_for(MMG5_pMesh mesh, MMG5_pSol met,_MMG5_pBucket bucket,int ne,
                 int* ifilt,int* ns,int* nc,int* warn,int it) {
    MMG5_pTetra     pt;
    MMG5_pxTetra    pxt;
    MMG5_Tria       ptt;
    MMG5_pPoint     p0,p1,ppt;
    MMG5_pxPoint    pxp;
    double     dd,len,lmax,o[3],to[3],ro[3],no1[3],no2[3],v[3];
    int        k,ip,ip1,ip2,list[_MMG5_LMAX+2],ilist,ref;
    char       imax,tag,j,i,i1,i2,ifa0,ifa1;
    int        lon,ret,ier;
    double     lmin;
    int        imin,iq;
    int        ii;
    double     lmaxtet,lmintet;
    int        imaxtet,imintet;

    for (k=1; k<=ne; k++) {
        pt = &mesh->tetra[k];
        if ( !MG_EOK(pt)  || (pt->tag & MG_REQ) )   continue;

        pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

        /* 1) find longest and shortest edge  and try to manage it*/
        imax = -1; lmax = 0.0;
        imin = -1; lmin = DBL_MAX;
        for (ii=0; ii<6; ii++) {
            if ( pt->xt && (pxt->tag[ii] & MG_REQ) )  continue;
            ip1  = _MMG5_iare[ii][0];
            ip2  = _MMG5_iare[ii][1];
            len = _MMG5_lenedg(mesh,met,pt->v[ip1],pt->v[ip2]);
            if ( len > lmax ) {
                lmax = len;
                imax = ii;
            }
            if ( len < lmin ) {
                lmin = len;
                imin = ii;
            }
        }
        if ( imax==-1 )
            fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are boundary and required\n",
                    __FILE__,__LINE__,k);
        if ( imin==-1 )
            fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are boundary and required\n",
                    __FILE__,__LINE__,k);

        if ( lmax >= _MMG5_LOPTL_MMG5_DEL )  {
            /* proceed edges according to lengths */
            ifa0 = _MMG5_ifar[imax][0];
            ifa1 = _MMG5_ifar[imax][1];
            i  = (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
            j  = _MMG5_iarfinv[i][imax];
            i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
            i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
            ip1 = pt->v[i1];
            ip2 = pt->v[i2];
            p0  = &mesh->point[ip1];
            p1  = &mesh->point[ip2];

            /* Case of a boundary face */
            if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
                if ( !(MG_GET(pxt->ori,i)) ) continue;
                ref = pxt->edg[_MMG5_iarf[i][j]];
                tag = pxt->tag[_MMG5_iarf[i][j]];
                if ( tag & MG_REQ )  continue;
                tag |= MG_BDY;
                ilist = _MMG5_coquil(mesh,k,imax,list);
                if ( !ilist )  continue;
                else if ( ilist<0 ) return(-1);
                if ( tag & MG_NOM ){
                    if( !_MMG5_BezierNom(mesh,ip1,ip2,0.5,o,no1,to) )
                        continue;
                    else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
                        _MMG5_tet2tri(mesh,k,i,&ptt);
                        _MMG5_nortri(mesh,&ptt,no1);
                        if ( !MG_GET(pxt->ori,i) ) {
                            no1[0] *= -1.0;
                            no1[1] *= -1.0;
                            no1[2] *= -1.0;
                        }
                    }
                }
                else if ( tag & MG_GEO ) {
                    if ( !_MMG5_BezierRidge(mesh,ip1,ip2,0.5,o,no1,no2,to) )
                        continue;
                    if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
                        _MMG5_tet2tri(mesh,k,i,&ptt);
                        _MMG5_nortri(mesh,&ptt,no1);
                        no2[0] = to[1]*no1[2] - to[2]*no1[1];
                        no2[1] = to[2]*no1[0] - to[0]*no1[2];
                        no2[2] = to[0]*no1[1] - to[1]*no1[0];
                        dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
                        if ( dd > _MMG5_EPSD2 ) {
                            dd = 1.0 / sqrt(dd);
                            no2[0] *= dd;
                            no2[1] *= dd;
                            no2[2] *= dd;
                        }
                    }
                }
                else if ( tag & MG_REF ) {
                    if ( !_MMG5_BezierRef(mesh,ip1,ip2,0.5,o,no1,to) )
                        goto collapse;
                }
                else {
                    if ( !_MMG5_norface(mesh,k,i,v) )  goto collapse;
                    if ( !_MMG5_BezierReg(mesh,ip1,ip2,0.5,v,o,no1) ) goto collapse;

                }
                ier = _MMG5_simbulgept(mesh,list,ilist,o);
                if ( !ier ) {
                    ier = _MMG5_dichoto1b(mesh,list,ilist,o,ro);
                    memcpy(o,ro,3*sizeof(double));
                }
                ip = _MMG5_newPt(mesh,o,tag);

                if ( !ip ) {
                    /* reallocation of point table */
                    if ( bucket ) {
                        _MMG5_POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
                                                 *warn=1;
                                                 goto collapse,
                                                 o,tag);
                    }
                    else {
                        printf("ERROR: function not available in delaunay mode. Exiting\n");
                        exit(EXIT_FAILURE);
                    }
                }

                if ( met->m )
                    met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

                ier = _MMG5_split1b(mesh,met,list,ilist,ip,1);

                /* if we realloc memory in _MMG5_split1b pt and pxt pointers are not valid */
                pt = &mesh->tetra[k];
                pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

                if ( ier < 0 ) {
                    fprintf(stdout,"  ## Error: unable to split.\n");
                    _MMG5_delPt(mesh,ip);
                    return(-1);
                }
                else if ( !ier ) {
                    _MMG5_delPt(mesh,ip);
                    goto collapse;
                } else {
                    (*ns)++;

                    ppt = &mesh->point[ip];
                    if ( MG_EDG(tag) || (tag & MG_NOM) )
                        ppt->ref = ref;
                    else
                        ppt->ref = pxt->ref[i];
                    ppt->tag = tag;
                    if ( met->m )
                        met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

                    pxp = &mesh->xpoint[ppt->xp];
                    if ( tag & MG_NOM ){
                        memcpy(pxp->n1,no1,3*sizeof(double));
                        memcpy(pxp->t,to,3*sizeof(double));
                    }
                    else if ( tag & MG_GEO ) {
                        memcpy(pxp->n1,no1,3*sizeof(double));
                        memcpy(pxp->n2,no2,3*sizeof(double));
                        memcpy(pxp->t,to,3*sizeof(double));
                    }
                    else if ( tag & MG_REF ) {
                        memcpy(pxp->n1,no1,3*sizeof(double));
                        memcpy(pxp->t,to,3*sizeof(double));
                    }
                    else
                        memcpy(pxp->n1,no1,3*sizeof(double));
                }
                continue;
            }
            else if(pt->xt){
                if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) {
                    continue;
                }
                ilist = _MMG5_coquil(mesh,k,imax,list);
                if ( !ilist )    continue;
                else if ( ilist<0 ) return(-1);
                o[0] = 0.5*(p0->c[0] + p1->c[0]);
                o[1] = 0.5*(p0->c[1] + p1->c[1]);
                o[2] = 0.5*(p0->c[2] + p1->c[2]);
                ip = _MMG5_newPt(mesh,o,MG_NOTAG);

                if ( !ip )  {
                    /* reallocation of point table */
                    if ( bucket ) {
                        _MMG5_POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
                                                 *warn=1;
                                                 goto collapse,
                                                 o,MG_NOTAG);
                    }
                    else {
                        printf("ERROR: function not available in delaunay mode. Exiting\n");
                        exit(EXIT_FAILURE);
                    }

                }

                if ( met->m )
                    met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

                ier = _MMG5_split1b(mesh,met,list,ilist,ip,1);
                if ( ier < 0 ) {
                    fprintf(stdout,"  ## Error: unable to split.\n");
                    _MMG5_delPt(mesh,ip);
                    return(-1);
                }
                else if ( !ier ) {
                    _MMG5_delPt(mesh,ip);
                    goto collapse;
                }
                else {
                    ppt = &mesh->point[ip];
                    met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
                    _MMG5_addBucket(mesh,bucket,ip);
                    (*ns)++;
                    continue;
                }

                /* Case of an internal face */
            } else {
                ilist = _MMG5_coquil(mesh,k,imax,list);
                if ( !ilist )    continue;
                else if ( ilist<0 ) return(-1);
                else if(ilist%2) goto collapse; //bdry edge
                o[0] = 0.5*(p0->c[0] + p1->c[0]);
                o[1] = 0.5*(p0->c[1] + p1->c[1]);
                o[2] = 0.5*(p0->c[2] + p1->c[2]);
                ip = _MMG5_newPt(mesh,o,MG_NOTAG);

                if ( !ip )  {
                    /* reallocation of point table */
                    if ( bucket ) {
                        _MMG5_POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
                                                 *warn=1;
                                                 goto collapse,
                                                 o,MG_NOTAG);
                    }
                    else {
                        printf("ERROR: function not available in delaunay mode. Exiting\n");
                        exit(EXIT_FAILURE);
                    }
                }

                if ( met->m )
                    met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

                /* Delaunay */
                if ( !_MMG5_buckin_iso(mesh,met,bucket,ip) ) {
                    _MMG5_delPt(mesh,ip);
                    (*ifilt)++;
                    goto collapse;
                } else {
                    lon = _MMG5_cavity(mesh,met,k,ip,list,ilist/2);
                    if ( lon < 1 ) {
                        MMG_npd++;
                        _MMG5_delPt(mesh,ip);
                        goto collapse;
                    } else {
                        ret = _MMG5_delone(mesh,met,ip,list,lon);
                        if ( ret > 0 ) {
                            ppt = &mesh->point[ip];
                            met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);

                            _MMG5_addBucket(mesh,bucket,ip);
                            (*ns)++;
                            continue;
                        }
                        else if ( ret == 0 ) {
                            MMG_npd++;
                            _MMG5_delPt(mesh,ip);
                            goto collapse;//continue;
                        }
                        else { /*allocation problem ==> saveMesh*/
                            _MMG5_delPt(mesh,ip);
                            return(0);
                        }
                    }
                }
            }
        }
    collapse:
        if(lmin <= _MMG5_LOPTS_MMG5_DEL) {
            ifa0 = _MMG5_ifar[imin][0];
            ifa1 = _MMG5_ifar[imin][1];
            i  =  (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
            j  = _MMG5_iarfinv[i][imin];
            i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
            i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
            ip = pt->v[i1];
            iq = pt->v[i2];
            p0 = &mesh->point[ip];
            p1 = &mesh->point[iq];

            if ( (p0->tag > p1->tag) || (p0->tag & MG_REQ) )  continue;

            /* Case of a boundary face */
            ilist = 0;
            if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
                tag = pxt->tag[_MMG5_iarf[i][j]];
                if ( tag & MG_REQ )  continue;
                tag |= MG_BDY;
                if ( p0->tag > tag )   continue;
                if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) continue;
                ilist = _MMG5_chkcol_bdy(mesh,k,i,j,list);
                if ( ilist > 0 ) {
                    ier = _MMG5_colver(mesh,list,ilist,i2);

                    if ( ier < 0 ) return(-1);
                    else if(ier) {
                        _MMG5_delPt(mesh,ier);
                        (*nc)++;
                        continue;
                    }
                }
                else if (ilist < 0 )  return(-1);
            }
            /* Case of an internal face */
            else {
                if ( p0->tag & MG_BDY )  continue;
                ilist = _MMG5_chkcol_int(mesh,met,k,i,j,list,2);
                if ( ilist > 0 ) {
                    ier = _MMG5_colver(mesh,list,ilist,i2);
                    if ( ilist < 0 ) continue;
                    if ( ier < 0 ) return(-1);
                    else if(ier) {
                        _MMG5_delBucket(mesh,bucket,ier);
                        _MMG5_delPt(mesh,ier);
                        (*nc)++;
                        continue;
                    }
                }
                else if (ilist < 0 )  return(-1);
            }
        } //end if lmin < _MMG5_LOPTS_MMG5_DEL

        /*2) longest and shortest edges are stucked => try another edges*/
        imaxtet = imax;
        imintet = imin;
        lmaxtet = lmax;
        lmintet = lmin;
        for (ii=0; ii<6; ii++) {
            if ( pt->xt && (pxt->tag[ii] & MG_REQ) )  continue;
            if ( (ii==imintet) && (lmintet < _MMG5_LOPTS_MMG5_DEL)) continue;
            if ( (ii==imaxtet) && (lmaxtet > _MMG5_LOPTL_MMG5_DEL) ) continue;

            ip1  = _MMG5_iare[ii][0];
            ip2  = _MMG5_iare[ii][1];
            len = _MMG5_lenedg(mesh,met,pt->v[ip1],pt->v[ip2]);

            imax = ii;
            lmax = len;
            imin = ii;
            lmin = len;
            if ( lmax >= _MMG5_LOPTL_MMG5_DEL )  {
                /* proceed edges according to lengths */
                ifa0 = _MMG5_ifar[imax][0];
                ifa1 = _MMG5_ifar[imax][1];
                i  = (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
                j  = _MMG5_iarfinv[i][imax];
                i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
                i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
                ip1 = pt->v[i1];
                ip2 = pt->v[i2];
                p0  = &mesh->point[ip1];
                p1  = &mesh->point[ip2];

                /* Case of a boundary face */
                if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
                    if ( !(MG_GET(pxt->ori,i)) ) continue;
                    ref = pxt->edg[_MMG5_iarf[i][j]];
                    tag = pxt->tag[_MMG5_iarf[i][j]];
                    if ( tag & MG_REQ )  continue;
                    tag |= MG_BDY;
                    ilist = _MMG5_coquil(mesh,k,imax,list);
                    if ( !ilist )  continue;
                    else if ( ilist<0 ) return(-1);
                    if ( tag & MG_NOM ){
                        if( !_MMG5_BezierNom(mesh,ip1,ip2,0.5,o,no1,to) )
                            continue;
                        else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
                            _MMG5_tet2tri(mesh,k,i,&ptt);
                            _MMG5_nortri(mesh,&ptt,no1);
                            if ( !MG_GET(pxt->ori,i) ) {
                                no1[0] *= -1.0;
                                no1[1] *= -1.0;
                                no1[2] *= -1.0;
                            }
                        }
                    }
                    else if ( tag & MG_GEO ) {
                        if ( !_MMG5_BezierRidge(mesh,ip1,ip2,0.5,o,no1,no2,to) )
                            continue;
                        if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
                            _MMG5_tet2tri(mesh,k,i,&ptt);
                            _MMG5_nortri(mesh,&ptt,no1);
                            no2[0] = to[1]*no1[2] - to[2]*no1[1];
                            no2[1] = to[2]*no1[0] - to[0]*no1[2];
                            no2[2] = to[0]*no1[1] - to[1]*no1[0];
                            dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
                            if ( dd > _MMG5_EPSD2 ) {
                                dd = 1.0 / sqrt(dd);
                                no2[0] *= dd;
                                no2[1] *= dd;
                                no2[2] *= dd;
                            }
                        }
                    }
                    else if ( tag & MG_REF ) {
                        if ( !_MMG5_BezierRef(mesh,ip1,ip2,0.5,o,no1,to) )
                            goto collapse2;
                    }
                    else {
                        if ( !_MMG5_norface(mesh,k,i,v) )  goto collapse2;
                        if ( !_MMG5_BezierReg(mesh,ip1,ip2,0.5,v,o,no1) ) goto collapse2;

                    }
                    ier = _MMG5_simbulgept(mesh,list,ilist,o);
                    if ( !ier ) {
                        ier = _MMG5_dichoto1b(mesh,list,ilist,o,ro);
                        memcpy(o,ro,3*sizeof(double));
                    }
                    ip = _MMG5_newPt(mesh,o,tag);

                    if ( !ip ){
                        /* reallocation of point table */
                        if ( bucket ) {
                            _MMG5_POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
                                                     *warn=1;
                                                     goto collapse2//break
                                                     ,o,tag);
                        }
                        else {
                            printf("ERROR: function not available in delaunay mode. Exiting\n");
                            exit(EXIT_FAILURE);
                        }
                    }

                    if ( met->m )
                        met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

                    ier = _MMG5_split1b(mesh,met,list,ilist,ip,1);
                    /* if we realloc memory in _MMG5_split1b pt and pxt pointers are not valid */
                    pt = &mesh->tetra[k];
                    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

                    if ( ier < 0 ) {
                        fprintf(stdout,"  ## Error: unable to split.\n");
                        return(-1);
                    }
                    else if ( !ier ) {
                        _MMG5_delPt(mesh,ip);
                        goto collapse2;//continue;
                    } else {
                        (*ns)++;
                        //_MMG5_addBucket(mesh,bucket,ip);

                        ppt = &mesh->point[ip];
                        if ( MG_EDG(tag) || (tag & MG_NOM) )
                            ppt->ref = ref;
                        else
                            ppt->ref = pxt->ref[i];
                        ppt->tag = tag;
                        if ( met->m )
                            met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

                        pxp = &mesh->xpoint[ppt->xp];
                        if ( tag & MG_NOM ){
                            memcpy(pxp->n1,no1,3*sizeof(double));
                            memcpy(pxp->t,to,3*sizeof(double));
                        }
                        else if ( tag & MG_GEO ) {
                            memcpy(pxp->n1,no1,3*sizeof(double));
                            memcpy(pxp->n2,no2,3*sizeof(double));
                            memcpy(pxp->t,to,3*sizeof(double));
                        }
                        else if ( tag & MG_REF ) {
                            memcpy(pxp->n1,no1,3*sizeof(double));
                            memcpy(pxp->t,to,3*sizeof(double));
                        }
                        else
                            memcpy(pxp->n1,no1,3*sizeof(double));
                    }
                    break;//imax continue;
                }
                else if(pt->xt){
                    if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) {
                        continue;
                    }
                    ilist = _MMG5_coquil(mesh,k,imax,list);
                    if ( !ilist )    continue;
                    else if ( ilist<0 ) return(-1);
                    o[0] = 0.5*(p0->c[0] + p1->c[0]);
                    o[1] = 0.5*(p0->c[1] + p1->c[1]);
                    o[2] = 0.5*(p0->c[2] + p1->c[2]);
                    ip = _MMG5_newPt(mesh,o,MG_NOTAG);

                    if ( !ip )  {
                        /* reallocation of point table */
                        if ( bucket ) {
                            _MMG5_POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
                                                     *warn=1;
                                                     goto collapse2
                                                     ,o,MG_NOTAG);
                        }
                        else {
                            printf("ERROR: function not available in delaunay mode. Exiting\n");
                            exit(EXIT_FAILURE);
                        }
                    }

                    if ( met->m )
                        met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

                    ier = _MMG5_split1b(mesh,met,list,ilist,ip,1);
                    if ( ier < 0 ) {
                        fprintf(stdout,"  ## Error: unable to split.\n");
                        _MMG5_delPt(mesh,ip);
                        return(-1);
                    }
                    else if ( !ier ) {
                        _MMG5_delPt(mesh,ip);
                        goto collapse2;
                    }
                    else {
                        ppt = &mesh->point[ip];
                        met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
                        _MMG5_addBucket(mesh,bucket,ip);
                        (*ns)++;
                        break;//imax continue;
                    }
                    printf("on doit pas passer la\n");
                    /* Case of an internal face */
                } else {
                    ilist = _MMG5_coquil(mesh,k,imax,list);
                    if ( !ilist )    continue;
                    else if ( ilist<0 ) return(-1);
                    else if(ilist%2) goto collapse2; //bdry edge
                    o[0] = 0.5*(p0->c[0] + p1->c[0]);
                    o[1] = 0.5*(p0->c[1] + p1->c[1]);
                    o[2] = 0.5*(p0->c[2] + p1->c[2]);
                    ip = _MMG5_newPt(mesh,o,MG_NOTAG);

                    if ( !ip )  {
                        /* reallocation of point table */
                        if ( bucket ) {
                            _MMG5_POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
                                                     *warn=1;
                                                     goto collapse2,
                                                     o,MG_NOTAG);
                        } else {
                            printf("ERROR: function not available in delaunay mode. Exiting\n");
                            exit(EXIT_FAILURE);
                        }
                    }

                    if ( met->m )
                        met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

                    if ( /*lmax>4 &&*/ /*it &&*/  !_MMG5_buckin_iso(mesh,met,bucket,ip) ) {
                        _MMG5_delPt(mesh,ip);
                        (*ifilt)++;
                        goto collapse2;
                    } else {
                        lon = _MMG5_cavity(mesh,met,k,ip,list,ilist/2);
                        if ( lon < 1 ) {
                            MMG_npd++;
                            _MMG5_delPt(mesh,ip);
                            goto collapse2;
                        } else {
                            ret = _MMG5_delone(mesh,met,ip,list,lon);
                            if ( ret > 0 ) {
                                ppt = &mesh->point[ip];
                                met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);

                                _MMG5_addBucket(mesh,bucket,ip);
                                (*ns)++;
                                break;//imax continue;
                            }
                            else if ( ret == 0 ) {
                                MMG_npd++;
                                _MMG5_delPt(mesh,ip);
                                goto collapse2;//continue;
                            }
                            else { /*allocation problem ==> savemesh*/
                                _MMG5_delPt(mesh,ip);
                                return(0);
                            }
                        }
                    }
                }
            }
        collapse2:
            if(lmin > _MMG5_LOPTS_MMG5_DEL) continue;
            ifa0 = _MMG5_ifar[imin][0];
            ifa1 = _MMG5_ifar[imin][1];
            i  =  (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
            j  = _MMG5_iarfinv[i][imin];
            i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
            i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
            ip = pt->v[i1];
            iq = pt->v[i2];
            p0 = &mesh->point[ip];
            p1 = &mesh->point[iq];

            if ( (p0->tag > p1->tag) || (p0->tag & MG_REQ) )  continue;

            /* Case of a boundary face */
            ilist = 0;
            if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
                tag = pxt->tag[_MMG5_iarf[i][j]];
                if ( tag & MG_REQ )  continue;
                tag |= MG_BDY;
                if ( p0->tag > tag )   continue;
                if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) continue;
                ilist = _MMG5_chkcol_bdy(mesh,k,i,j,list);
                if ( ilist > 0 ) {
                    ier = _MMG5_colver(mesh,list,ilist,i2);
                    if ( ier < 0 ) return(-1);
                    else if(ier) {
                        _MMG5_delPt(mesh,ier);
                        (*nc)++;
                        break;
                    }
                }
                else if (ilist < 0 )  return(-1);
            }
            /* Case of an internal face */
            else {
                if ( p0->tag & MG_BDY )  continue;
                ilist = _MMG5_chkcol_int(mesh,met,k,i,j,list,2);
                if ( ilist > 0 ) {
                    ier = _MMG5_colver(mesh,list,ilist,i2);
                    if ( ilist < 0 ) continue;
                    if ( ier < 0 ) return(-1);
                    else if(ier) {
                        _MMG5_delBucket(mesh,bucket,ier);
                        _MMG5_delPt(mesh,ier);
                        (*nc)++;
                        break;
                    }
                }
                else if (ilist < 0 )  return(-1);
            }
        }//end for ii
    }

    return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param bucket pointer toward the bucket structure.
 * \return -1 if fail and we dont try to end the remesh process,
 * 0 if fail but we try to end the remesh process and 1 if success.
 *
 * Split edges longer than \ref _MMG5_LOPTL_MMG5_DEL and collapse edges shorter
 * than \ref _MMG5_LOPTS_MMG5_DEL.
 *
 */
static int
_MMG5_adpsplcol(MMG5_pMesh mesh,MMG5_pSol met,_MMG5_pBucket bucket, int* warn) {
    int        ifilt,ne,ier;
    int        ns,nc,it,nnc,nns,nnf,nnm,maxit,nf,nm;
    double     maxgap;

    /* Iterative mesh modifications */
    it = nnc = nns = nnf = nnm = 0;
    maxit = 10;
    mesh->gap = maxgap = 0.5;
    MMG_npuiss = MMG_nvol = MMG_npres = MMG_npd = 0;
    do {
        if ( !mesh->info.noinsert ) {
            *warn=0;
            ns = nc = 0;
            nf = nm = 0;
            ifilt = 0;
            ne = mesh->ne;
            ier = _MMG5_boucle_for(mesh,met,bucket,ne,&ifilt,&ns,&nc,warn,it);
            if(ier<0) exit(EXIT_FAILURE);
            else if(!ier) return(-1);
        } /* End conditional loop on mesh->info.noinsert */
        else  ns = nc = ifilt = 0;

        if ( !mesh->info.noswap ) {
            nf = _MMG5_swpmsh(mesh,met,bucket);
            if ( nf < 0 ) {
                fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
                return(0);
            }
            nnf += nf;
            if(it==2 || it==6/*&& it==1 || it==3 || it==5 || it > 8*/) {
                nf += _MMG5_swptet(mesh,met,1.053,bucket);
            } else {
                nf += 0;
            }
            if ( nf < 0 ) {
                fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
                return(0);
            }
        }
        else  nf = 0;
        nnf+=nf;

        if ( !mesh->info.nomove ) {
            nm = _MMG5_movtet(mesh,met,-1);
            if ( nm < 0 ) {
                fprintf(stdout,"  ## Unable to improve mesh.\n");
                return(0);
            }
        }
        else  nm = 0;
        nnm += nm;
        nnc += nc;
        nns += ns;

        /* decrease size of gap for reallocation */

        if ( mesh->gap > maxgap/(double)maxit )
            mesh->gap -= maxgap/(double)maxit;
        else
            mesh->gap -= mesh->gap/(double)maxit;


        if ( ((abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && ns+nc > 0) )
            fprintf(stdout,"     %8d filtered %8d splitted, %8d collapsed, %8d swapped, %8d moved\n",ifilt,ns,nc,nf,nm);
        if ( ns < 10 && abs(nc-ns) < 3 )  break;
        else if ( it > 3 && abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;

    }
    while( ++it < maxit && nc+ns > 0 );

    return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param bucket pointer toward the bucket structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Mesh optimization using egde swapping and point relocation.
 *
 */
static int
_MMG5_optet(MMG5_pMesh mesh, MMG5_pSol met,_MMG5_pBucket bucket) {
    int it,nnm,nnf,maxit,nm,nf;
    double declic;

    /* shape optim */
    it = nnm = nnf = 0;
    maxit = 4;
    declic = 1.053;
    do {
        /* badly shaped process */
        if ( !mesh->info.noswap ) {
            nf = _MMG5_swpmsh(mesh,met,bucket);
            if ( nf < 0 ) {
                fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
                return(0);
            }
            nnf += nf;

            nf = _MMG5_swptet(mesh,met,declic,bucket);
            if ( nf < 0 ) {
                fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
                return(0);
            }
        }
        else  nf = 0;

        if ( !mesh->info.nomove ) {
            nm = _MMG5_movtet(mesh,met,0);
            if ( nm < 0 ) {
                fprintf(stdout,"  ## Unable to improve mesh.\n");
                return(0);
            }
        }
        else  nm = 0;
        nnm += nm;

        /* ier = _MMG5_badelt(mesh,met);
           if ( ier < 0 ) {
           fprintf(stdout,"  ## Unable to remove bad elements.\n");
           return(0);
           }*/

        if ( (abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && nf+nm > 0 ){
            fprintf(stdout,"                                            ");
            fprintf(stdout,"%8d swapped, %8d moved\n",nf,nm);
        }
    }
    while( ++it < maxit && nm+nf > 0 );

    if ( !mesh->info.nomove ) {
        nm = _MMG5_movtet(mesh,met,3);
        if ( nm < 0 ) {
            fprintf(stdout,"  ## Unable to improve mesh.\n");
            return(0);
        }
    }
    else  nm = 0;
    nnm += nm;
    if ( (abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && nm > 0 )
        fprintf(stdout,"                                            ");
    fprintf(stdout,"                  %8d moved\n",nm);

    return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param bucket pointer toward the bucket structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Analyze tetrahedra and split long / collapse short, according to
 * prescribed metric.
 *
 */
static int
_MMG5_adptet_delone(MMG5_pMesh mesh,MMG5_pSol met,_MMG5_pBucket bucket) {
    int      nnf,ns,nf;
    int      warn;

    /*initial swap*/
    if ( !mesh->info.noswap ) {
        nf = _MMG5_swpmsh(mesh,met,bucket);
        if ( nf < 0 ) {
            fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
            return(0);
        }
        nnf = nf;
        nf = _MMG5_swptet(mesh,met,1.053,bucket);
        if ( nf < 0 ) {
            fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
            return(0);
        }
        nnf+=nf;
    } else  nnf = nf = 0;

#ifdef DEBUG
    fprintf(stdout,"$$$$$$$$$$$$$$$$$$ INITIAL SWAP %7d\n",nnf);
    _MMG5_outqua(mesh,met);
#endif

    /* Iterative mesh modifications */
    warn = 0;

    ns = _MMG5_adpsplcol(mesh,met,bucket,&warn);

    if ( ns < 0 ) {
        fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
        return(0);
    }

    if ( warn ) {
        fprintf(stdout,"  ## Error:");
        fprintf(stdout," unable to allocate a new point in last call of adpspl.\n");
        fprintf(stdout,"  ## Check the mesh size or ");
        fprintf(stdout,"increase the maximal authorized memory with the -m option.\n");
        fprintf(stdout,"  ## Uncomplete mesh. Exiting\n" );
        return(0);
    }

    /* renumerotation if available */
    if ( !_MMG5_scotchCall(mesh,met) )
        return(0);

    if(!_MMG5_optet(mesh,met,bucket)) return(0);

    return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if failed, 1 if success.
 *
 * Main adaptation routine.
 *
 */
int _MMG5_mmg3d1_delone(MMG5_pMesh mesh,MMG5_pSol met) {
    _MMG5_pBucket bucket;

    if ( abs(mesh->info.imprim) > 3 )
        fprintf(stdout,"  ** MESH ANALYSIS\n");

    if ( mesh->info.iso && !_MMG5_chkmani(mesh) ) {
        fprintf(stdout,"  ## Non orientable implicit surface. Exit program.\n");
        return(0);
    }

    /**--- stage 1: geometric mesh */
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
        fprintf(stdout,"  ** GEOMETRIC MESH\n");

    if ( !_MMG5_anatet(mesh,met,1,0) ) {
        fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
        return(0);
    }

#ifdef DEBUG
    _MMG5_outqua(mesh,met);
#endif

    /**--- stage 2: computational mesh */
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
        fprintf(stdout,"  ** COMPUTATIONAL MESH\n");

    /* define metric map */
    if ( !_MMG5_defsiz(mesh,met) ) {
        fprintf(stdout,"  ## Metric undefined. Exit program.\n");
        return(0);
    }

    if ( mesh->info.hgrad > 0. && !_MMG5_gradsiz(mesh,met) ) {
        fprintf(stdout,"  ## Gradation problem. Exit program.\n");
        return(0);
    }
    if ( !_MMG5_anatet(mesh,met,2,0) ) {
        fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
        return(0);
    }

#ifdef DEBUG
    puts("---------------------------Fin anatet---------------------");
    _MMG5_outqua(mesh,met);
#endif

    /* renumerotation if available */
    if ( !_MMG5_scotchCall(mesh,met) )
        return(0);

    /* CEC : create filter */
    bucket = _MMG5_newBucket(mesh,mesh->info.bucket); //M_MAX(mesh->mesh->info.bucksiz,BUCKSIZ));
    if ( !bucket )  return(0);

    if ( !_MMG5_adptet_delone(mesh,met,bucket) ) {
        fprintf(stdout,"  ## Unable to adapt. Exit program.\n");
        return(0);
    }

#ifdef DEBUG
    puts("---------------------Fin adptet-----------------");
    _MMG5_outqua(mesh,met);
#endif
    /* in test phase: check if no element with 2 bdry faces */
    if ( !_MMG5_chkfemtopo(mesh) ) {
        fprintf(stdout,"  ## Topology of mesh unsuited for fem computations. Exit program.\n");
        return(0);
    }

    if ( mesh->info.iso && !_MMG5_chkmani(mesh) ) {
        fprintf(stdout,"  ## Non orientable implicit surface. Exit program.\n");
        return(0);
    }

    /*free bucket*/
    _MMG5_DEL_MEM(mesh,bucket->head,(bucket->size*bucket->size*bucket->size+1)*sizeof(int));
    _MMG5_DEL_MEM(mesh,bucket->link,(mesh->npmax+1)*sizeof(int));
    _MMG5_DEL_MEM(mesh,bucket,sizeof(_MMG5_Bucket));

    return(1);
}
