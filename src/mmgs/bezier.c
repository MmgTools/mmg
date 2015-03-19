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
 * \file mmgs/bezier.c
 * \brief Functions for Bezier surface computation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

extern char ddb;

/* return Bezier control points on triangle iel (cf. Vlachos) */
int bezierCP(pMesh mesh,int iel,pBezier pb) {
    pTria     pt;
    pPoint    p[3];
    double   *n1,*n2,nt[3],ps,ps2,dd,ux,uy,uz;
    char      i,i1,i2;

    pt   = &mesh->tria[iel];
    p[0] = &mesh->point[pt->v[0]];
    p[1] = &mesh->point[pt->v[1]];
    p[2] = &mesh->point[pt->v[2]];

    memset(pb,0,sizeof(Bezier));

    /* first 3 CP = vertices, normals */
    for (i=0; i<3; i++) {
        memcpy(&pb->b[i],p[i]->c,3*sizeof(double));
        pb->p[i] = p[i];

        if ( MS_SIN(p[i]->tag) ) {
            nortri(mesh,pt,pb->n[i]);
        }
        else if ( MS_EDG(p[i]->tag) ) {
            nortri(mesh,pt,nt);
            n1 = &mesh->geom[p[i]->ig].n1[0];
            n2 = &mesh->geom[p[i]->ig].n2[0];

            ps  = n1[0]*nt[0] + n1[1]*nt[1] + n1[2]*nt[2];
            ps2 = n2[0]*nt[0] + n2[1]*nt[1] + n2[2]*nt[2];
            if ( fabs(ps) > fabs(ps2) )
                memcpy(&pb->n[i],n1,3*sizeof(double));
            else
                memcpy(&pb->n[i],n2,3*sizeof(double));
            memcpy(&pb->t[i],p[i]->n,3*sizeof(double));
        }
        else
            memcpy(&pb->n[i],p[i]->n,3*sizeof(double));
    }

    /* compute control points along edges */
    for (i=0; i<3; i++) {
        i1 = inxt[i];
        i2 = inxt[i1];

        ux = p[i2]->c[0] - p[i1]->c[0];
        uy = p[i2]->c[1] - p[i1]->c[1];
        uz = p[i2]->c[2] - p[i1]->c[2];

        /* choose normals */
        n1 = pb->n[i1];
        n2 = pb->n[i2];

        /* check for boundary curve */
        if ( MS_EDG(pt->tag[i]) ) {
            if ( MS_SIN(p[i1]->tag) ) {
                dd = 1.0 / 3.0;
                pb->b[2*i+3][0] = p[i1]->c[0] + dd * ux;
                pb->b[2*i+3][1] = p[i1]->c[1] + dd * uy;
                pb->b[2*i+3][2] = p[i1]->c[2] + dd * uz;
            }
            else {
                dd = (ux*pb->t[i1][0] + uy*pb->t[i1][1] + uz*pb->t[i1][2]) / 3.0;
                pb->b[2*i+3][0] = p[i1]->c[0] + dd * pb->t[i1][0];
                pb->b[2*i+3][1] = p[i1]->c[1] + dd * pb->t[i1][1];
                pb->b[2*i+3][2] = p[i1]->c[2] + dd * pb->t[i1][2];
            }
            if ( MS_SIN(p[i2]->tag) ) {
                dd = 1.0 / 3.0;
                pb->b[2*i+4][0] = p[i2]->c[0] - dd * ux;
                pb->b[2*i+4][1] = p[i2]->c[1] - dd * uy;
                pb->b[2*i+4][2] = p[i2]->c[2] - dd * uz;
            }
            else {
                dd = -(ux*pb->t[i2][0] + uy*pb->t[i2][1] + uz*pb->t[i2][2]) / 3.0;
                pb->b[2*i+4][0] = p[i2]->c[0] + dd * pb->t[i2][0];
                pb->b[2*i+4][1] = p[i2]->c[1] + dd * pb->t[i2][1];
                pb->b[2*i+4][2] = p[i2]->c[2] + dd * pb->t[i2][2];
            }

            /* tangent evaluation */
            ps = ux*(pb->t[i1][0]+pb->t[i2][0]) + uy*(pb->t[i1][1]+pb->t[i2][1]) + uz*(pb->t[i1][2]+pb->t[i2][2]);
            ps = 2.0*ps / (ux*ux + uy*uy + uz*uz);
            pb->t[i+3][0] = pb->t[i1][0] + pb->t[i2][0] - ps*ux;
            pb->t[i+3][1] = pb->t[i1][1] + pb->t[i2][1] - ps*uy;
            pb->t[i+3][2] = pb->t[i1][2] + pb->t[i2][2] - ps*uz;
            dd = pb->t[i+3][0]*pb->t[i+3][0] + pb->t[i+3][1]*pb->t[i+3][1] + pb->t[i+3][2]*pb->t[i+3][2];
            if ( dd > EPSD2 ) {
                dd = 1.0 / sqrt(dd);
                pb->t[i+3][0] *= dd;
                pb->t[i+3][1] *= dd;
                pb->t[i+3][2] *= dd;
            }
        }
        else { /* internal edge */
            ps = ux*n1[0] + uy*n1[1] + uz*n1[2];
            pb->b[2*i+3][0] = (2.0*p[i1]->c[0] + p[i2]->c[0] - ps*n1[0]) / 3.0;
            pb->b[2*i+3][1] = (2.0*p[i1]->c[1] + p[i2]->c[1] - ps*n1[1]) / 3.0;
            pb->b[2*i+3][2] = (2.0*p[i1]->c[2] + p[i2]->c[2] - ps*n1[2]) / 3.0;

            ps = -(ux*n2[0] + uy*n2[1] + uz*n2[2]);
            pb->b[2*i+4][0] = (2.0*p[i2]->c[0] + p[i1]->c[0] - ps*n2[0]) / 3.0;
            pb->b[2*i+4][1] = (2.0*p[i2]->c[1] + p[i1]->c[1] - ps*n2[1]) / 3.0;
            pb->b[2*i+4][2] = (2.0*p[i2]->c[2] + p[i1]->c[2] - ps*n2[2]) / 3.0;
        }

        /* normal evaluation */
        ps = ux*(n1[0]+n2[0]) + uy*(n1[1]+n2[1]) + uz*(n1[2]+n2[2]);
        ps = 2.0*ps / (ux*ux + uy*uy + uz*uz);
        pb->n[i+3][0] = n1[0] + n2[0] - ps*ux;
        pb->n[i+3][1] = n1[1] + n2[1] - ps*uy;
        pb->n[i+3][2] = n1[2] + n2[2] - ps*uz;
        dd = pb->n[i+3][0]*pb->n[i+3][0] + pb->n[i+3][1]*pb->n[i+3][1] + pb->n[i+3][2]*pb->n[i+3][2];
        if ( dd > EPSD2 ) {
            dd = 1.0 / sqrt(dd);
            pb->n[i+3][0] *= dd;
            pb->n[i+3][1] *= dd;
            pb->n[i+3][2] *= dd;
        }
    }

    /* Central Bezier coefficient */
    for (i=0; i<3; i++) {
        pb->b[9][0] -= (0.5*ATHIRD*pb->b[i][0]);
        pb->b[9][1] -= (0.5*ATHIRD*pb->b[i][1]);
        pb->b[9][2] -= (0.5*ATHIRD*pb->b[i][2]);
    }

    for (i=0; i<3; i++) {
        pb->b[9][0] += 0.25 * (pb->b[2*i+3][0] + pb->b[2*i+4][0]);
        pb->b[9][1] += 0.25 * (pb->b[2*i+3][1] + pb->b[2*i+4][1]);
        pb->b[9][2] += 0.25 * (pb->b[2*i+3][2] + pb->b[2*i+4][2]);
    }

    return(1);
}


/* return point o at (u,v) in Bezier patch and normal */
int bezierInt(pBezier pb,double uv[2],double o[3],double no[3],double to[3]) {
    double    dd,u,v,w,ps,ux,uy,uz;
    char      i;

    memset(to,0,3*sizeof(double));
    u = uv[0];
    v = uv[1];
    w = 1 - u - v;

    /* coordinates + normals */
    for (i=0; i<3; i++) {
        o[i]  = pb->b[0][i]*w*w*w + pb->b[1][i]*u*u*u + pb->b[2][i]*v*v*v \
            + 3.0 * (pb->b[3][i]*u*u*v + pb->b[4][i]*u*v*v + pb->b[5][i]*w*v*v \
                     + pb->b[6][i]*w*w*v + pb->b[7][i]*w*w*u + pb->b[8][i]*w*u*u)\
            + 6.0 * pb->b[9][i]*u*v*w;

        /* quadratic interpolation of normals */
        no[i] = pb->n[0][i]*w*w + pb->n[1][i]*u*u + pb->n[2][i]*v*v \
            + 2.0*(pb->n[3][i]*u*v + pb->n[4][i]*v*w + pb->n[5][i]*u*w);

        /* linear interpolation, not used here
           no[i] = pb->n[0][i]*w + pb->n[1][i]*u + pb->n[2][i]*v; */
    }

    /* tangent */
    if ( w < EPSD2 ) {
        ux = pb->b[2][0] - pb->b[1][0];
        uy = pb->b[2][1] - pb->b[1][1];
        uz = pb->b[2][2] - pb->b[1][2];
        dd = ux*ux + uy*uy + uz*uz;
        if ( dd > EPSD2 ) {
            dd = 1.0 / sqrt(dd);
            ux *= dd;
            uy *= dd;
            uz *= dd;
        }

        if ( MS_SIN(pb->p[1]->tag) ) {
            pb->t[1][0] = ux;
            pb->t[1][1] = uy;
            pb->t[1][2] = uz;
        }
        if ( MS_SIN(pb->p[2]->tag) ) {
            pb->t[2][0] = ux;
            pb->t[2][1] = uy;
            pb->t[2][2] = uz;
        }

        ps = pb->t[1][0]* pb->t[2][0] + pb->t[1][1]* pb->t[2][1] + pb->t[1][2]* pb->t[2][2];
        if ( ps > 0.0 ) {
            to[0] = pb->t[1][0]*u + pb->t[2][0]*v;
            to[1] = pb->t[1][1]*u + pb->t[2][1]*v;
            to[2] = pb->t[1][2]*u + pb->t[2][2]*v;
        }
        else {
            to[0] = -pb->t[1][0]*u + pb->t[2][0]*v;
            to[1] = -pb->t[1][1]*u + pb->t[2][1]*v;
            to[2] = -pb->t[1][2]*u + pb->t[2][2]*v;
        }
    }

    if ( u < EPSD2 ) {
        ux = pb->b[2][0] - pb->b[0][0];
        uy = pb->b[2][1] - pb->b[0][1];
        uz = pb->b[2][2] - pb->b[0][2];
        dd = ux*ux + uy*uy + uz*uz;
        if ( dd > EPSD2 ) {
            dd = 1.0 / sqrt(dd);
            ux *= dd;
            uy *= dd;
            uz *= dd;
        }

        if ( MS_SIN(pb->p[0]->tag) ) {
            pb->t[0][0] = ux;
            pb->t[0][1] = uy;
            pb->t[0][2] = uz;
        }
        if ( MS_SIN(pb->p[2]->tag) ) {
            pb->t[2][0] = ux;
            pb->t[2][1] = uy;
            pb->t[2][2] = uz;
        }

        ps = pb->t[0][0]* pb->t[2][0] + pb->t[0][1]* pb->t[2][1] + pb->t[0][2]* pb->t[2][2];
        if ( ps > 0.0 ) {
            to[0] = pb->t[0][0]*w + pb->t[2][0]*v;
            to[1] = pb->t[0][1]*w + pb->t[2][1]*v;
            to[2] = pb->t[0][2]*w + pb->t[2][2]*v;
        }
        else {
            to[0] = -pb->t[0][0]*w + pb->t[2][0]*v;
            to[1] = -pb->t[0][1]*w + pb->t[2][1]*v;
            to[2] = -pb->t[0][2]*w + pb->t[2][2]*v;
        }
    }

    if ( v < EPSD2 ) {
        ux = pb->b[1][0] - pb->b[0][0];
        uy = pb->b[1][1] - pb->b[0][1];
        uz = pb->b[1][2] - pb->b[0][2];
        dd = ux*ux + uy*uy + uz*uz;
        if ( dd > EPSD2 ) {
            dd = 1.0 / sqrt(dd);
            ux *= dd;
            uy *= dd;
            uz *= dd;
        }

        if ( MS_SIN(pb->p[0]->tag) ) {
            pb->t[0][0] = ux;
            pb->t[0][1] = uy;
            pb->t[0][2] = uz;
        }
        if ( MS_SIN(pb->p[1]->tag) ) {
            pb->t[1][0] = ux;
            pb->t[1][1] = uy;
            pb->t[1][2] = uz;
        }

        ps = pb->t[0][0]* pb->t[1][0] + pb->t[0][1]* pb->t[1][1] + pb->t[0][2]* pb->t[1][2];
        if ( ps > 0.0 ) {
            to[0] = pb->t[0][0]*w + pb->t[1][0]*u;
            to[1] = pb->t[0][1]*w + pb->t[1][1]*u;
            to[2] = pb->t[0][2]*w + pb->t[1][2]*u;
        }
        else {
            to[0] = -pb->t[0][0]*w + pb->t[1][0]*u;
            to[1] = -pb->t[0][1]*w + pb->t[1][1]*u;
            to[2] = -pb->t[0][2]*w + pb->t[1][2]*u;
        }
    }

    dd = no[0]*no[0] + no[1]*no[1] + no[2]*no[2];
    if ( dd > EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        no[0] *= dd;
        no[1] *= dd;
        no[2] *= dd;
    }

    dd = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];
    if ( dd > EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        to[0] *= dd;
        to[1] *= dd;
        to[2] *= dd;
    }

    return(1);
}

/* Computes the Bezier coefficients associated to the underlying curve to [p0p1]
   isrid = 0 if p0p1 is not a special edge, 1 otherwise */
inline void bezierEdge(pMesh mesh,int i0,int i1,double b0[3],double b1[3],char isrid,double v[3]) {
    pPoint    p0,p1;
    double    ux,uy,uz,*n1,*n2,*t,ps1,ps2;

    p0 = &mesh->point[i0];
    p1 = &mesh->point[i1];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    if ( isrid ) {
        if ( MS_SIN(p0->tag) ) {
            b0[0] = p0->c[0] + ATHIRD*ux;
            b0[1] = p0->c[1] + ATHIRD*uy;
            b0[2] = p0->c[2] + ATHIRD*uz;
        }
        else {
            t = &p0->n[0];
            ps1 = t[0]*ux + t[1]*uy + t[2]*uz;
            b0[0] = p0->c[0] + ATHIRD*ps1*t[0];
            b0[1] = p0->c[1] + ATHIRD*ps1*t[1];
            b0[2] = p0->c[2] + ATHIRD*ps1*t[2];
        }

        if (MS_SIN(p1->tag) ) {
            b1[0] = p1->c[0] - ATHIRD*ux;
            b1[1] = p1->c[1] - ATHIRD*uy;
            b1[2] = p1->c[2] - ATHIRD*uz;
        }
        else {
            t = &p1->n[0];
            ps1 = -(t[0]*ux + t[1]*uy + t[2]*uz);
            b1[0] = p1->c[0] + ATHIRD*ps1*t[0];
            b1[1] = p1->c[1] + ATHIRD*ps1*t[1];
            b1[2] = p1->c[2] + ATHIRD*ps1*t[2];
        }
    }

    /* regular edge */
    else {
        if ( MS_SIN(p0->tag) ) {
            b0[0] = p0->c[0] + ATHIRD*ux;
            b0[1] = p0->c[1] + ATHIRD*uy;
            b0[2] = p0->c[2] + ATHIRD*uz;
        }
        else {
            if ( MS_GEO & p0->tag ) {
                n1 = &mesh->geom[p0->ig].n1[0];
                n2 = &mesh->geom[p0->ig].n2[0];
                ps1 = v[0]*n1[0] + v[1]*n1[1] + v[2]*n1[2];
                ps2 = v[0]*n2[0] + v[1]*n2[1] + v[2]*n2[2];
                if ( ps1 < ps2 ) {
                    n1 = &mesh->geom[p0->ig].n2[0];
                    ps1 = ps2;
                }
            }
            else if ( MS_REF & p0->tag ) {
                n1 = &mesh->geom[p0->ig].n1[0];
                ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
            }
            else {
                n1 = &p0->n[0];
                ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
            }
            b0[0] = ATHIRD*(2.0*p0->c[0] + p1->c[0] - ps1*n1[0]);
            b0[1] = ATHIRD*(2.0*p0->c[1] + p1->c[1] - ps1*n1[1]);
            b0[2] = ATHIRD*(2.0*p0->c[2] + p1->c[2] - ps1*n1[2]);
        }

        if ( MS_SIN(p1->tag) ) {
            b1[0] = p1->c[0] - ATHIRD*ux;
            b1[1] = p1->c[1] - ATHIRD*uy;
            b1[2] = p1->c[2] - ATHIRD*uz;
        }
        else {
            if ( MS_GEO & p1->tag ) {
                n1 = &mesh->geom[p1->ig].n1[0];
                n2 = &mesh->geom[p1->ig].n2[0];
                ps1 = -(v[0]*n1[0] + v[1]*n1[1] + v[2]*n1[2]);
                ps2 = -(v[0]*n2[0] + v[1]*n2[1] + v[2]*n2[2]);
                if ( fabs(ps2) < fabs(ps1) ) {
                    n1 = &mesh->geom[p1->ig].n2[0];
                    ps1 = ps2;
                }
            }
            else if ( MS_REF & p1->tag ) {
                n1 = &mesh->geom[p1->ig].n1[0];
                ps1 = -(ux*n1[0] + uy*n1[1] + uz*n1[2]);
            }
            else {
                n1 = &p1->n[0];
                ps1 = -(ux*n1[0] + uy*n1[1] + uz*n1[2]);
            }
            b1[0] = ATHIRD*(2.0*p1->c[0] + p0->c[0] - ps1*n1[0]);
            b1[1] = ATHIRD*(2.0*p1->c[1] + p0->c[1] - ps1*n1[1]);
            b1[2] = ATHIRD*(2.0*p1->c[2] + p0->c[2] - ps1*n1[2]);
        }
    }
}
