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
 * \file mmgs/split.c
 * \brief Functions to create new points.
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


/* split element k along edge i */
int split1(pMesh mesh,pSol met,int k,int i,int *vx) {
    pTria      pt,pt1;
    pPoint     ppt;
    int        iel;
    char       i1,i2;

    iel = newElt(mesh);
    assert(iel);

    pt  = &mesh->tria[k];
    pt->flag = 0;
    pt1 = &mesh->tria[iel];
    pt1 = memcpy(pt1,pt,sizeof(Tria));

    i1 = inxt[i];
    i2 = inxt[i1];


    if ( pt->edg[i] > 0 ) {
        ppt = &mesh->point[vx[i]];
        ppt->ref = pt->edg[i];
    }

    pt->v[i2]   = pt1->v[i1]   = vx[i];
    pt->tag[i1] = pt1->tag[i2] = MS_NOTAG;
    pt->edg[i1] = pt1->edg[i2] = 0;

    return(1);
}


/* Split tria k, along edge i, inserting point ip, updating adjacency relations */
int split1b(pMesh mesh,int k,char i,int ip) {
    pTria     pt,pt1;
    pPoint    ppt;
    Bezier    b;
    pGeom     go;
    double    uv[2],o[3],no[3],to[3];
    int      *adja,iel,jel,kel,mel,ier;
    char      i1,i2,j,j1,j2,m;

    pt = &mesh->tria[k];
    pt->flag = 0;
    pt->base = mesh->base;

    iel = newElt(mesh);
    if ( !iel )  return(0);

    pt1 = &mesh->tria[iel];
    memcpy(pt1,pt,sizeof(Tria));
    memcpy(&mesh->adja[3*(iel-1)+1],&mesh->adja[3*(k-1)+1],3*sizeof(int));

    ppt = &mesh->point[ip];
    if ( pt->edg[i] )  ppt->ref = pt->edg[i];
    if ( pt->tag[i] )  ppt->tag = pt->tag[i];

    adja = &mesh->adja[3*(k-1)+1];
    jel = adja[i] / 3;
    j   = adja[i] % 3;

    /* update normal n2 if need be */
    if ( jel && pt->tag[i] & MS_GEO ) {
        ier = bezierCP(mesh,jel,&b);
        assert(ier);
        uv[0] = 0.5;
        uv[1] = 0.5;
        if ( j == 1 )       uv[0] = 0.0;
        else if ( j == 2 )  uv[1] = 0.0;

        ier = bezierInt(&b,uv,o,no,to);
        assert(ier);
        go = &mesh->geom[ppt->ig];
        memcpy(go->n2,no,3*sizeof(double));
    }

    /* update two triangles */
    i1  = inxt[i];
    i2  = iprv[i];
    pt->v[i2]   = ip;
    pt->tag[i1] = MS_NOTAG;
    pt->edg[i1] = 0;
    pt1->v[i1]   = ip;
    pt1->tag[i2] = MS_NOTAG;
    pt1->edg[i2] = 0;

    /* update adjacency relations */
    mel = adja[i1] / 3;
    m   = adja[i1] % 3;
    mesh->adja[3*(k-1)+1+i1]   = 3*iel+i2;
    mesh->adja[3*(iel-1)+1+i2] = 3*k+i1;
    mesh->adja[3*(iel-1)+1+i1] = 3*mel+m;
    if(mel)
        mesh->adja[3*(mel-1)+1+m]  = 3*iel+i1;

    if ( jel ) {
        kel = newElt(mesh);
        assert(kel);
        pt  = &mesh->tria[jel];
        pt1 = &mesh->tria[kel];
        pt->flag = 0;
        pt->base = mesh->base;
        memcpy(pt1,pt,sizeof(Tria));
        memcpy(&mesh->adja[3*(kel-1)+1],&mesh->adja[3*(jel-1)+1],3*sizeof(int));

        j1 = inxt[j];
        j2 = iprv[j];
        pt->v[j1]    = ip;
        pt->tag[j2]  = MS_NOTAG;
        pt->edg[j2]  = 0;
        pt1->v[j2]   = ip;
        pt1->tag[j1] = MS_NOTAG;
        pt1->edg[j1] = 0;

        /* update adjacency */
        adja = &mesh->adja[3*(jel-1)+1];
        mel  = adja[j2] / 3;
        m    = adja[j2] % 3;
        mesh->adja[3*(jel-1)+1+j2] = 3*kel+j1;
        mesh->adja[3*(kel-1)+1+j1] = 3*jel+j2;
        mesh->adja[3*(kel-1)+1+j2] = 3*mel+m;
        if(mel)
            mesh->adja[3*(mel-1)+1+m]  = 3*kel+j2;

        mesh->adja[3*(iel-1)+1+i]  = 3*kel+j;
        mesh->adja[3*(kel-1)+1+j]  = 3*iel+i;
    }

    return(1);
}


/* split element k along 2 edges i1 and i2 */
int split2(pMesh mesh,pSol met,int k,int *vx) {
    pTria    pt,pt1,pt2;
    pPoint   p0,p1,p2,p3,p4;
    int      iel,jel;
    char     i,i1,i2;

    /* create 2 elements */
    iel = newElt(mesh);
    assert(iel);
    jel = newElt(mesh);
    assert(jel);

    pt  = &mesh->tria[k];
    pt->flag = 0;
    pt1 = &mesh->tria[iel];
    pt2 = &mesh->tria[jel];
    pt1 = memcpy(pt1,pt,sizeof(Tria));
    pt2 = memcpy(pt2,pt,sizeof(Tria));

    i = 0;
    if ( !vx[0] )  i = 1;
    else if ( !vx[1] )  i = 2;
    i1 = inxt[i];
    i2 = inxt[i1];

    p0 = &mesh->point[pt->v[i]];
    p1 = &mesh->point[pt->v[i1]];
    p2 = &mesh->point[pt->v[i2]];
    p3 = &mesh->point[vx[i]];
    p4 = &mesh->point[vx[i1]];

    /* update refs */
    if ( pt->edg[i] > 0 )   p3->ref = pt->edg[i];
    if ( pt->edg[i1] > 0 )  p4->ref = pt->edg[i1];

    /* check alternate configs */

    if ( 1 ) {
        pt->v[i1] = pt1->v[i2] = pt2->v[i1] = vx[i];
        pt->v[i2] = pt2->v[i]  = vx[i1];

        pt->tag[i] = pt->tag[i2] = pt1->tag[i1] = pt2->tag[i2] = MS_NOTAG;
        pt->edg[i] = pt->edg[i2] = pt1->edg[i1] = pt2->edg[i2] = 0;
    }
    else {
        pt->v[i2]  = pt1->v[i]  = pt2->v[i] = vx[i1];
        pt1->v[i2] = pt2->v[i1] = vx[i];

        pt->tag[i] = pt1->tag[i1] = pt1->tag[i2] = pt2->tag[i2] = MS_NOTAG;
        pt->edg[i] = pt1->edg[i1] = pt1->edg[i2] = pt2->edg[i2] = 0;
    }

    return(1);
}


/* split all 3 edges of element k */
int split3(pMesh mesh,pSol met,int k,int *vx) {
    pTria    pt,pt1,pt2,pt3;
    pPoint   p0,p1,p2,p3,p4,p5;
    int      iel,jel,kel;

    /* create 3 elements */
    iel = newElt(mesh);
    assert(iel);
    jel = newElt(mesh);
    assert(jel);
    kel = newElt(mesh);
    assert(jel);

    pt  = &mesh->tria[k];
    pt->flag = 0;
    pt1 = &mesh->tria[iel];
    pt2 = &mesh->tria[jel];
    pt3 = &mesh->tria[kel];
    pt1 = memcpy(pt1,pt,sizeof(Tria));
    pt2 = memcpy(pt2,pt,sizeof(Tria));
    pt3 = memcpy(pt3,pt,sizeof(Tria));

    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    p3 = &mesh->point[vx[0]];
    p4 = &mesh->point[vx[1]];
    p5 = &mesh->point[vx[2]];

    /* update refs */
    if ( pt->edg[0] > 0 )  p3->ref = pt->edg[0];
    if ( pt->edg[1] > 0 )  p4->ref = pt->edg[1];
    if ( pt->edg[2] > 0 )  p5->ref = pt->edg[2];

    /* update topo */
    pt->v[1]  = pt1->v[0] = pt3->v[0] = vx[2];
    pt->v[2]  = pt2->v[0] = pt3->v[2] = vx[1];
    pt1->v[2] = pt2->v[1] = pt3->v[1] = vx[0];

    pt->tag[0]  = pt1->tag[1] = pt2->tag[2] = MS_NOTAG;
    pt->edg[0]  = pt1->edg[1] = pt2->edg[2] = 0;

    pt3->tag[0] = pt3->tag[1] = pt3->tag[2] = MS_NOTAG;
    pt3->edg[0] = pt3->edg[1] = pt3->edg[2] = 0;

    return(1);
}

