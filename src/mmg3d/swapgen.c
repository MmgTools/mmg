/* =============================================================================
**  This file is part of the Mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
**
**  Mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with Mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the Mmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file mmg3d/swapgen.c
 * \brief Functions for swapping process inside the mesh.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg3d.h"

/**
 * \param mesh pointer toward the mesh structure
 * \param start tetrahedra in which the swap should be performed
 * \param ia edge that we want to swap
 * \param ilist pointer to store the size of the shell of the edge
 * \param list pointer to store the shell of the edge
 * \param crit improvment coefficient
 * \return 0 if fail, the index of point corresponding to the swapped
 * configuration otherwise (\f$4*k+i\f$).
 *
 * Check whether swap of edge \a ia in \a start should be performed, and
 * return \f$4*k+i\f$ the index of point corresponding to the swapped
 * configuration. The shell of edge is built during the process.
 *
 */
int _MMG5_chkswpgen(MMG5_pMesh mesh,int start,int ia,int *ilist,int *list,double crit) {
    MMG5_pTetra    pt,pt0;
    MMG5_pPoint    p0;
    double    calold,calnew,caltmp;
    int       na,nb,np,adj,piv,npol,refdom,k,l,iel;
    int       *adja,pol[_MMG5_LMAX+2];
    char      i,ipa,ipb,ip,ier;

    pt  = &mesh->tetra[start];
    refdom = pt->ref;

    pt0 = &mesh->tetra[0];
    na  = pt->v[_MMG5_iare[ia][0]];
    nb  = pt->v[_MMG5_iare[ia][1]];
    calold = pt->qual;


    /* Store shell of ia in list, and associated pseudo polygon in pol */
    (*ilist) = 0;
    npol = 0;
    list[(*ilist)] = 6*start+ia;
    (*ilist)++;
    adja = &mesh->adja[4*(start-1)+1];
    adj  = adja[_MMG5_ifar[ia][0]] / 4;      // start travelling by face (ia,0)
    piv  = pt->v[_MMG5_ifar[ia][1]];
    pol[npol] = 4*start + _MMG5_ifar[ia][1];
    npol++;

    while ( adj && adj != start ) {
        pt = &mesh->tetra[adj];
        if ( pt->tag & MG_REQ ) return(0);

        /* Edge is on a boundary between two different domains */
        if ( pt->ref != refdom )  return(0);
        calold = MG_MIN(calold, pt->qual);
        /* identification of edge number in tetra adj */
        for (i=0; i<6; i++) {
            ipa = _MMG5_iare[i][0];
            ipb = _MMG5_iare[i][1];
            if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                 (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
        }
        assert(i<6);
        list[(*ilist)] = 6*adj +i;
        (*ilist)++;
        /* overflow */
        if ( (*ilist) > _MMG5_LMAX-3 )  return(0);

        /* set new triangle for travel */
        adja = &mesh->adja[4*(adj-1)+1];
        if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
            pol[npol] = 4*adj + _MMG5_ifar[i][1];
            npol++;
            adj = adja[ _MMG5_ifar[i][0] ] / 4;
            piv = pt->v[ _MMG5_ifar[i][1] ];
        }
        else {
            assert(pt->v[ _MMG5_ifar[i][1] ] == piv);
            pol[npol] = 4*adj + _MMG5_ifar[i][0];
            npol++;
            adj = adja[ _MMG5_ifar[i][1] ] /4;
            piv = pt->v[ _MMG5_ifar[i][0] ];
        }
    }
    //CECILE : je vois pas pourquoi ca ameliore de faire ce test
    //plus rapide mais du coup on elimine des swap...
    //4/01/14 commentaire
    //if ( calold*_MMG5_ALPHAD > 0.5 )  return(0);

    /* Prevent swap of an external boundary edge */
    if ( !adj )  return(0);
#ifdef SINGUL
    /* Prevent swap of an internal ridge or required edge */
    if ( pt->xt && ( (mesh->xtetra[pt->xt].tag[i] & MG_GEO) ||
                     (mesh->xtetra[pt->xt].tag[i] & MG_REQ) ) ) return(0);
#endif

    assert(npol == (*ilist)); // du coup, apres on pourra virer npol

    /* Find a configuration that enhances the worst quality within the shell */
    for (k=0; k<npol; k++) {
        iel = pol[k] / 4;
        ip  = pol[k] % 4;
        np  = mesh->tetra[iel].v[ip];
        calnew = 1.0;
        ier = 1;

        if ( mesh->info.fem ) {
            p0 = &mesh->point[np];
            if ( p0->tag & MG_BDY ) {
                for (l=0; l<npol;l++) {
                    if ( k < npol-1 ) {
                        if ( l == k || l == k+1 )  continue;
                    }
                    else {
                        if ( l == npol-1 || l == 0 )  continue;
                    }
                    iel = pol[l] / 4;
                    ip  = pol[l] % 4;
                    pt = &mesh->tetra[iel];
                    p0 = &mesh->point[pt->v[ip]];
                    if ( p0->tag & MG_BDY ) {
                        ier = 0;
                        break;
                    }
                }
            }
            if ( !ier )  continue;
            ier = 1;
        }

        for (l=0; l<(*ilist); l++) {
            /* Do not consider tets of the shell of collapsed edge */
            if ( k < npol-1 ) {
                if ( l == k || l == k+1 )  continue;
            }
            else {
                if ( l == npol-1 || l == 0 )  continue;
            }
            iel = list[l] / 6;
            i   = list[l] % 6;
            pt  = &mesh->tetra[iel];

            /* First tetra obtained from iel */
            memcpy(pt0,pt,sizeof(MMG5_Tetra));
            pt0->v[_MMG5_iare[i][0]] = np;
            caltmp = _MMG5_orcal(mesh,0);
            calnew = MG_MIN(calnew,caltmp);
            /* Second tetra obtained from iel */
            memcpy(pt0,pt,sizeof(MMG5_Tetra));
            pt0->v[_MMG5_iare[i][1]] = np;
            caltmp = _MMG5_orcal(mesh,0);
            calnew = MG_MIN(calnew,caltmp);
            ier = (calnew > crit*calold);
            if ( !ier )  break;
        }
        if ( ier )  return(pol[k]);
    }
    return(0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param nconf configuration.
 * \param ilist number of tetrahedra in the shell of the edge that we want
 *  to swap.
 * \param list pointer toward the shell of the edge that we want to swap.
 * \param bucket pointer toward the bucket structure in Delaunay mode,
 * NULL pointer in pattern mode.
 * \return -1 if lack of memory, 0 if fail to swap, 1 otherwise.
 *
 * Perform swap of edge whose shell is passed according to configuration nconf.
 *
 */
int _MMG5_swpgen(MMG5_pMesh mesh,MMG5_pSol met,int nconf,int ilist,int *list,_MMG5_pBucket bucket) {
    MMG5_pTetra    pt;
    MMG5_pPoint    p0,p1;
    int       iel,na,nb,np,nball,ret,start;
    double    m[3];
    char      ia,ip,iq;
    int       ier;

    iel = list[0] / 6;
    ia  = list[0] % 6;

    pt = &mesh->tetra[iel];
    na = pt->v[_MMG5_iare[ia][0]];
    nb = pt->v[_MMG5_iare[ia][1]];
    p0 = &mesh->point[na];
    p1 = &mesh->point[nb];

    /* Temporarily create midpoint at swapped edge */
    m[0] = 0.5*(p0->c[0] + p1->c[0]);
    m[1] = 0.5*(p0->c[1] + p1->c[1]);
    m[2] = 0.5*(p0->c[2] + p1->c[2]);

    np  = _MMG5_newPt(mesh,m,0);
    if(!np){
        if ( bucket ) {
            _MMG5_POINT_AND_BUCKET_REALLOC(mesh,met,np,mesh->gap,
                                     printf("  ## Error: unable to allocate a new point\n");
                                     _MMG5_INCREASE_MEM_MESSAGE();
                                     return(-1)
                                     ,m,0);
        }
        else {
            _MMG5_POINT_REALLOC(mesh,met,np,mesh->gap,
                          printf("  ## Error: unable to allocate a new point\n");
                          _MMG5_INCREASE_MEM_MESSAGE();
                          return(-1)
                          ,m,0);
        }
    }
    if ( met->m )  met->m[np] = 0.5*(met->m[na]+met->m[nb]);

    /** First step : split of edge (na,nb) */
    ret = 2*ilist + 0;
    ier = _MMG5_split1b(mesh,met,list,ret,np,0);
    /* pointer adress may change if we need to realloc memory during split */
    pt = &mesh->tetra[iel];

    if ( ier < 0 ) {
        fprintf(stdout,"  ## Warning: unable to swap internal edge.\n");
        return(-1);
    }
    else if ( !ier )  return(0);

    /** Second step : collapse of np towards enhancing configuration */
    start = nconf / 4;
    iq = nconf % 4;

    pt = &mesh->tetra[start];
    for (ip=0; ip<4; ip++) {
        if ( pt->v[ip] == np )  break;
    }
    assert(ip<4);

    memset(list,0,(_MMG5_LMAX+2)*sizeof(int));
    nball = _MMG5_boulevolp(mesh,start,ip,list);

    ier = _MMG5_colver(mesh,list,nball,iq);
    if ( ier < 0 ) {
        fprintf(stdout,"  ## Warning: unable to swap internal edge.\n");
        return(-1);
    }
    else if ( ier ) _MMG5_delPt(mesh,ier);

    return(1);
}

#ifdef SINGUL
/**
 * \param mesh pointer toward the mesh structure
 * \param k index of the tetrahedron in which perform the swap
 * \param ip index of point to perform the swap
 *
 * Perform swap 2->3. Return 0 if the swap lead to an invalid configuration.
 *
 */
int _MMG5_swap23(MMG5_pMesh mesh,int k,int ip) {
    MMG5_pTetra        pt[3],pt0;
    MMG5_xTetra        xt[3],zerotet;
    MMG5_pxTetra       pxt0;
    double        qual0,qual1,qual2;
    int           newtet[3],*adja,iq,np,nq,i;
    unsigned char tau1[4],tau2[4],*taued1,*taued2,isxt[3];

    /** Table that associates to each (even) permutation of the 4 vertices of a   *
     *  tetrahedron the corresponding permutation of its edges                    *
     *  (warning: this table is not exactly the same than the table in split.c ). *
     *  Labels :
     *  0  : [0,1,2,3]
     *  1  : [0,2,3,1]
     *  2  : [0,3,1,2]
     *  3  : [1,0,3,2]
     *  4  : [1,3,2,0]
     *  5  : [1,2,0,3]
     *  6  : [2,0,1,3]
     *  7  : [2,1,3,0]
     *  8  : [2,3,0,1]
     *  9  : [3,0,2,1]
     *  10 : [3,2,1,0]
     *  11 : [3,1,0,2]  */
    unsigned char permedg[12][6] = {
        {0,1,2,3,4,5}, {1,2,0,5,3,4}, {2,0,1,4,5,3}, {0,4,3,2,1,5},
        {4,3,0,5,2,1}, {3,0,4,1,5,2}, {1,3,5,0,2,4}, {3,5,1,4,0,2},
        {5,1,3,2,4,0}, {2,5,4,1,0,3}, {5,4,2,3,1,0}, {4,2,5,0,3,1} };


    adja      = &mesh->adja[4*(k-1)+1];
    newtet[0] = k;
    newtet[1] = adja[ip]/4;
    iq        = adja[ip]%4;

    pt0       = &mesh->tetra[0];
    pt[0]     = &mesh->tetra[k];
    pt[1]     = &mesh->tetra[newtet[1]];
    np        = pt[0]->v[ip];
    nq        = pt[1]->v[iq];

    memset(&zerotet,0,sizeof(MMG5_xTetra));
    zerotet.ori = 15;

    if ( pt[0]->xt ) {
        pxt0 = &mesh->xtetra[pt[0]->xt];
        memcpy(&xt[0],pxt0,sizeof(MMG5_xTetra));
        memcpy(&xt[1],pxt0,sizeof(MMG5_xTetra));
        memcpy(&xt[2],pxt0,sizeof(MMG5_xTetra));
    }
    else {
        memset(&xt[0],0,sizeof(MMG5_xTetra));
        memset(&xt[1],0,sizeof(MMG5_xTetra));
        memset(&xt[2],0,sizeof(MMG5_xTetra));
    }

    if ( pt[1]->xt )
        pxt0 = &mesh->xtetra[pt[1]->xt];
    else
        pxt0 = &zerotet;

    if ( (pt[0]->tag & MG_REQ) || (pt[1]->tag & MG_REQ) ||
         (xt[0].ftag[ip] & MG_REQ) || (pxt0->ftag[iq] & MG_REQ) )
        return(0);

    /* For now, we don't deal with 2 domains in this swap */
    if ( pt[1]->ref != pt[0]->ref ) {
        fprintf(stdout,"  ## Warning: ");
        fprintf(stdout,"unable to perform swap 2->3 with 2 domains.\n");
        return(0);
    }
    if ( pt[1]->tag != pt[0]->tag ) {
        fprintf(stdout,"  ## Warning: ");
        fprintf(stdout,"unable to perform swap 2->3 with 2 tetras");
        fprintf(stdout," with different tags.\n");
        return(0);
    }


    tau1[0] = ip;
    tau1[1] = _MMG5_idir[ip][0];
    tau1[2] = _MMG5_idir[ip][1];
    tau1[3] = _MMG5_idir[ip][2];
    taued1  = &permedg[3*ip][0];

    /* Check if the first new tetra is valid */
    pt0->v[tau1[0]] = np;
    pt0->v[tau1[1]] = nq;
    pt0->v[tau1[2]] = pt[0]->v[tau1[2]];
    pt0->v[tau1[3]] = pt[0]->v[tau1[3]];
    qual0 = _MMG5_orcal(mesh,0);
    if ( qual0 < _MMG5_NULKAL )  return(0);

    /* Check if the second new tetra is valid */
    pt0->v[tau1[1]] = pt[0]->v[tau1[1]];
    pt0->v[tau1[2]] = nq;
    qual1 = _MMG5_orcal(mesh,0);
    if ( qual1 < _MMG5_NULKAL )  return(0);

    /* Check if the third new tetra is valid */
    pt0->v[tau1[2]] = pt[0]->v[tau1[2]];
    pt0->v[tau1[3]] = nq;
    qual2 = _MMG5_orcal(mesh,0);
    if ( qual2 < _MMG5_NULKAL )  return(0);

    /* All the new tetras are valid, we can swap */
    /* Update vertices and xt fields */
    memcpy(pt0,pt[0],sizeof(MMG5_Tetra));
    if ( pt[1]->v[_MMG5_idir[iq][0]] ==  pt0->v[tau1[1]] ) {
        tau2[0] = iq;
        tau2[1] = _MMG5_idir[iq][0];
        tau2[2] = _MMG5_idir[iq][1];
        tau2[3] = _MMG5_idir[iq][2];
        taued2  = &permedg[3*iq][0];
    }
    else if ( pt[1]->v[_MMG5_idir[iq][0]] == pt0->v[tau1[2]] ) {
        tau2[0] = iq;
        tau2[1] = _MMG5_idir[iq][1];
        tau2[2] = _MMG5_idir[iq][2];
        tau2[3] = _MMG5_idir[iq][0];
        taued2  = &permedg[3*iq+1][0];
    }
    else if ( pt[1]->v[_MMG5_idir[iq][0]] == pt0->v[tau1[3]] ) {
        tau2[0] = iq;
        tau2[1] = _MMG5_idir[iq][2];
        tau2[2] = _MMG5_idir[iq][0];
        tau2[3] = _MMG5_idir[iq][1];
        taued2  = &permedg[3*iq+2][0];
    }
    else {
        fprintf(stdout,"%s:%d: ",__FILE__,__LINE__);
        fprintf(stdout,"Error: we should not passed here. Exit.\n");
        exit(0);
    }

    /* First tet*/
    pt[0]->v[tau1[1]] = nq;
    pt[0]->base = pt[0]->mark = pt[0]->flag = 0;

    xt[0].tag[taued1[0]] = 0;  xt[0].edg[taued1[0]] = 0;
    xt[0].tag[taued1[3]] = pxt0->tag[taued2[1]];
    xt[0].edg[taued1[3]] = pxt0->edg[taued2[1]];
    xt[0].tag[taued1[4]] = pxt0->tag[taued2[2]];
    xt[0].edg[taued1[4]] = pxt0->edg[taued2[2]];
    xt[0].ref  [tau1[2]] = 0;  xt[0].ref  [tau1[3]] = 0;
    xt[0].ftag [tau1[2]] = 0;  xt[0].ftag [tau1[3]] = 0;
    xt[0].ref  [tau1[0]] = pxt0->ref [tau2[1]];
    xt[0].ftag [tau1[0]] = pxt0->ftag[tau2[1]];

    MG_SET(xt[0].ori, tau1[2]);  MG_SET(xt[0].ori, tau1[3]);
    if ( MG_GET(pxt0->ori, tau2[1]) )  MG_SET(xt[0].ori, tau1[0]);
    else  MG_CLR(xt[0].ori, tau1[0]);

    /* Second tet*/
    pt[1]->v[tau1[0]] = pt0->v[tau1[0]]; pt[1]->v[tau1[1]] = pt0->v[tau1[1]];
    pt[1]->v[tau1[2]] = nq             ; pt[1]->v[tau1[3]] = pt0->v[tau1[3]];
    pt[1]->base = pt[1]->mark = pt[1]->flag = 0;

    xt[1].tag[taued1[1]] = 0;  xt[1].edg[taued1[1]] = 0;
    xt[1].tag[taued1[3]] = pxt0->tag[taued2[0]];
    xt[1].edg[taued1[3]] = pxt0->edg[taued2[0]];
    xt[1].tag[taued1[5]] = pxt0->tag[taued2[1]];
    xt[1].edg[taued1[5]] = pxt0->edg[taued2[1]];
    xt[1].ref  [tau1[1]] = 0;  xt[1].ref [tau1[3]] = 0;
    xt[1].ftag [tau1[1]] = 0;  xt[1].ftag[tau1[3]] = 0;

    MG_SET(xt[1].ori, tau1[1]);  MG_SET(xt[1].ori, tau1[3]);
    if ( MG_GET(pxt0->ori, tau2[3]) )  MG_SET(xt[1].ori, tau1[0]);
    else  MG_CLR(xt[1].ori, tau1[0]);


    /* Third tet */
    newtet[2] =  _MMG5_newElt(mesh);
    if ( !newtet[2] ) {
        fprintf(stdout,"%s:%d: Error: unable to allocate a new element\n"
                ,__FILE__,__LINE__);
        return(0);
    }
    pt[2] = &mesh->tetra[newtet[2]];
    pt[2]->v[tau1[0]] = pt0->v[tau1[0]] ; pt[2]->v[tau1[1]] = pt0->v[tau1[1]];
    pt[2]->v[tau1[2]] = pt0->v[tau1[2]] ; pt[2]->v[tau1[3]] = nq;
    pt[2]->ref  = pt0->ref;
    pt[2]->tag  = pt0->tag;
    pt[2]->base = pt[2]->mark = pt[2]->flag = 0;

    xt[2].tag[taued1[2]] = 0;  xt[2].edg[taued1[2]] = 0;
    xt[2].tag[taued1[4]] = pxt0->tag[taued2[0]];
    xt[2].edg[taued1[4]] = pxt0->edg[taued2[0]];
    xt[2].tag[taued1[5]] = pxt0->tag[taued2[2]];
    xt[2].edg[taued1[5]] = pxt0->edg[taued2[2]];
    xt[2].ref  [tau1[1]] = 0;  xt[2].ref [tau1[2]] = 0;
    xt[2].ftag [tau1[1]] = 0;  xt[2].ftag[tau1[2]] = 0;

    MG_SET(xt[2].ori, tau1[1]);  MG_SET(xt[2].ori, tau1[2]);
    if ( MG_GET(pxt0->ori, tau2[2]) )  MG_SET(xt[2].ori, tau1[0]);
    else  MG_CLR(xt[2].ori, tau1[0]);

    /* Assignation of the xt newtet[2] to the appropriate tets */
    memset(isxt,0,3*sizeof(char));
    for (i=0; i<4; i++) {
        if ( xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
        if ( xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
        if ( xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
    }
    if ( mesh->info.sing ) {
        for (i=0; i<6; i++) {
            if ( xt[0].edg[i] || xt[0].tag[i] ) isxt[0] = 1;
            if ( xt[1].edg[i] || xt[1].tag[i] ) isxt[1] = 1;
            if ( xt[2].edg[i] || xt[2].tag[i] ) isxt[2] = 1;
        }
    }

    if ( pt[0]->xt ) {
        if ( isxt[0] ) {
            memcpy(&mesh->xtetra[pt[0]->xt],&xt[0],sizeof(MMG5_xTetra));
            if ( pt[1]->xt ) {
                if ( isxt[1] ) {
                    memcpy(pxt0,&xt[1],sizeof(MMG5_xTetra));
                    if ( isxt[2] ) {
                        mesh->xt++;
                        if ( mesh->xt > mesh->xtmax ) {
                            /* realloc of xtetras table */
                            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                         "larger xtetra table",
                                         mesh->xt--;
                                         return(0));
                        }
                        pt[2]->xt = mesh->xt;
                        memcpy(&mesh->xtetra[pt[2]->xt],&xt[2],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[2] */
                        pt[2]->xt = 0;
                }
                else /* !isxt[1] */ {
                    if ( isxt[2] ) {
                        pt[2]->xt = pt[1]->xt;
                        memcpy(&mesh->xtetra[pt[2]->xt],&xt[2],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[2] */
                        pt[2]->xt = 0;
                    pt[1]->xt = 0;
                }
            }
            else /* !pt[1]->xt */ {
                for (i=1; i<3; i++) {
                    if ( isxt[i] ) {
                        mesh->xt++;
                        if ( mesh->xt > mesh->xtmax ) {
                            /* realloc of xtetras table */
                            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                         "larger xtetra table",
                                         mesh->xt--;
                                         return(0));
                        }
                        pt[i]->xt = mesh->xt;
                        memcpy(&mesh->xtetra[pt[i]->xt],&xt[i],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[i] */
                        pt[i]->xt = 0;
                }
            }
        }
        else /* !isxt[0] */ {
            if ( pt[1]->xt ) {
                if ( isxt[1] ) {
                    memcpy(&mesh->xtetra[pt[1]->xt],&xt[1],sizeof(MMG5_xTetra));
                    if ( isxt[2] ) {
                        pt[2]->xt = pt[0]->xt;
                        memcpy(&mesh->xtetra[pt[2]->xt],&xt[2],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[2] */
                        pt[2]->xt = 0;
                }
                else /* !isxt[1] */{
                    if ( isxt[2] ) {
                        pt[2]->xt = pt[0]->xt;
                        memcpy(&mesh->xtetra[pt[2]->xt],&xt[2],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[2] */
                        pt[2]->xt = 0;
                    pt[1]->xt = 0;
                }
            }
            else /* !pt[1]->xt */ {
                if ( isxt[1] ) {
                    pt[1]->xt = pt[0]->xt;
                    memcpy(&mesh->xtetra[pt[1]->xt],&xt[1],sizeof(MMG5_xTetra));
                    if ( isxt[2] ) {
                        mesh->xt++;
                        if ( mesh->xt > mesh->xtmax ) {
                            /* realloc of xtetras table */
                            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                         "larger xtetra table",
                                         mesh->xt--;
                                         return(0));
                        }
                        pt[2]->xt = mesh->xt;
                        memcpy(&mesh->xtetra[pt[2]->xt],&xt[2],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[2] */
                        pt[2]->xt = 0;
                }
                else /* !isxt[1] */ {
                    if ( isxt[2] ) {
                        pt[2]->xt = pt[0]->xt;
                        memcpy(&mesh->xtetra[pt[2]->xt],&xt[2],sizeof(MMG5_xTetra));
                    }
                    else
                        pt[2]->xt = 0;
                }
            }
            pt[0]->xt = 0;
        }
    }
    else /* !pt[0]->xt */ {
        if ( isxt[0] ) {
            if ( pt[1]->xt ) {
                if ( isxt[1] ) {
                    memcpy(&mesh->xtetra[pt[1]->xt],&xt[1],sizeof(MMG5_xTetra));
                    mesh->xt++;
                    if ( mesh->xt > mesh->xtmax ) {
                        /* realloc of xtetras table */
                        _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                     "larger xtetra table",
                                     mesh->xt--;
                                     return(0));
                    }
                    pt[0]->xt = mesh->xt;
                    memcpy(&mesh->xtetra[pt[0]->xt],&xt[0],sizeof(MMG5_xTetra));
                    if ( isxt[2] ) {
                        mesh->xt++;
                        if ( mesh->xt > mesh->xtmax ) {
                            /* realloc of xtetras table */
                            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                         "larger xtetra table",
                                         mesh->xt--;
                                         return(0));
                        }
                        pt[2]->xt = mesh->xt;
                        memcpy(&mesh->xtetra[pt[2]->xt],&xt[2],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[2] */
                        pt[2]->xt = 0;
                }
                else /* !isxt[1] */ {
                    pt[0]->xt = pt[1]->xt;
                    pt[1]->xt = 0;
                    memcpy(&mesh->xtetra[pt[0]->xt],&xt[0],sizeof(MMG5_xTetra));
                    if ( isxt[2] ) {
                        mesh->xt++;
                        if ( mesh->xt > mesh->xtmax ) {
                            /* realloc of xtetras table */
                            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                         "larger xtetra table",
                                         mesh->xt--;
                                         return(0));
                        }
                        pt[2]->xt = mesh->xt;
                        memcpy(&mesh->xtetra[pt[2]->xt],&xt[2],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[2] */
                        pt[2]->xt = 0;
                }
            }
            else /* !pt[1]->xt */ {
                mesh->xt++;
                if ( mesh->xt > mesh->xtmax ) {
                    /* realloc of xtetras table */
                    _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 return(0));
                }
                pt[0]->xt = mesh->xt;
                memcpy(&mesh->xtetra[pt[0]->xt],&xt[0],sizeof(MMG5_xTetra));
                for (i=1; i<3; i++) {
                    if ( isxt[i] ) {
                        mesh->xt++;
                        if ( mesh->xt > mesh->xtmax ) {
                            /* realloc of xtetras table */
                            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                         "larger xtetra table",
                                         mesh->xt--;
                                         return(0));
                        }
                        pt[i]->xt = mesh->xt;
                        memcpy(&mesh->xtetra[pt[i]->xt],&xt[i],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[i] */
                        pt[i]->xt = 0;
                }
            }
        }
        else /* !isxt[0] */ {
            if ( pt[1]->xt ) {
                if ( isxt[1] ) {
                    memcpy(&mesh->xtetra[pt[1]->xt],&xt[1],sizeof(MMG5_xTetra));
                    if ( isxt[2] ) {
                        mesh->xt++;
                        if ( mesh->xt > mesh->xtmax ) {
                            /* realloc of xtetras table */
                            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                         "larger xtetra table",
                                         mesh->xt--;
                                         return(0));
                        }
                        pt[2]->xt = mesh->xt;
                        memcpy(&mesh->xtetra[pt[2]->xt],&xt[2],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[2] */
                        pt[2]->xt = 0;
                }
                else /* !isxt[1] */ {
                    if ( isxt[2] ) {
                        pt[2]->xt = pt[1]->xt;
                        memcpy(&mesh->xtetra[pt[2]->xt],&xt[2],sizeof(MMG5_xTetra));
                    }
                    else /* !isxt[2] */
                        pt[1]->xt = 0;
                }
            }
        }
    }

    /* Quality update */
    pt[0]->qual=qual0;
    pt[1]->qual=qual1;
    pt[2]->qual=qual2;

    /* Update of adjacency relations */
    /* Third tet */
    adja = &mesh->adja[4*(newtet[2]-1)+1];
    adja[tau1[0]] = mesh->adja[4*(newtet[1]-1)+tau2[2]+1];
    adja[tau1[1]] = 4*newtet[0] + tau1[3];
    adja[tau1[2]] = 4*newtet[1] + tau1[3];
    adja[tau1[3]] = mesh->adja[4*(newtet[0]    -1)+tau1[3]+1];
    if ( adja[tau1[0]] )
        mesh->adja[4*(adja[tau1[0]]/4-1)+adja[tau1[0]]%4+1] = 4*newtet[2] + tau1[0];
    if ( adja[tau1[3]] )
        mesh->adja[4*(adja[tau1[3]]/4-1)+adja[tau1[3]]%4+1] = 4*newtet[2] + tau1[3];

    /* First faces of first and second tet */
    adja = &mesh->adja[4*(newtet[0]-1)+1];
    adja[tau1[0]] = mesh->adja[4*(newtet[1]-1)+tau2[1]+1];
    if ( adja[tau1[0]] )
        mesh->adja[4*(adja[tau1[0]]/4-1)+adja[tau1[0]]%4+1] = 4*newtet[0] + tau1[0];

    adja = &mesh->adja[4*(newtet[1]-1)+1];
    adja[tau1[0]] = mesh->adja[4*(newtet[1]-1)+tau2[3]+1];
    if ( adja[tau1[0]] )
        mesh->adja[4*(adja[tau1[0]]/4-1)+adja[tau1[0]]%4+1] = 4*newtet[1] + tau1[0];

    /* Remaining faces of second tet */
    adja[tau1[1]] = 4*newtet[0] + tau1[2];
    adja[tau1[2]] = mesh->adja[4*(newtet[0]-1)+tau1[2]+1];
    adja[tau1[3]] = 4*newtet[2] + tau1[2];
    if ( adja[tau1[2]] )
        mesh->adja[4*(adja[tau1[2]]/4-1)+adja[tau1[2]]%4+1] = 4*newtet[1] + tau1[2];

    /* Remaining faces of first tet */
    adja = &mesh->adja[4*(newtet[0]-1)+1];
    adja[tau1[2]] = 4*newtet[1] + tau1[1];
    adja[tau1[3]] = 4*newtet[2] + tau1[1];
    if ( adja[tau1[1]] )
        mesh->adja[4*(adja[tau1[1]]/4-1)+adja[tau1[1]]%4+1] = 4*newtet[0] + tau1[1];

    return(1);
}
#endif
