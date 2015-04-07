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
 * \file mmgs/scalem.c
 * \brief Scale and unscale mesh and solution.
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


int scaleMesh(pMesh mesh,pSol met) {
    pPoint    ppt;
    pPar      par;
    double    dd,d1;
    int       i,k;

    /* compute bounding box */
    for (i=0; i<3; i++) {
        info.min[i] =  DBL_MAX;
        info.max[i] = -DBL_MAX;
    }
    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MS_VOK(ppt) )  continue;
        for (i=0; i<3; i++) {
            if ( ppt->c[i] > info.max[i] )  info.max[i] = ppt->c[i];
            if ( ppt->c[i] < info.min[i] )  info.min[i] = ppt->c[i];
        }
    }
    info.delta = 0.0;
    for (i=0; i<3; i++) {
        dd = info.max[i] - info.min[i];
        if ( dd > info.delta )  info.delta = dd;
    }
    if ( info.delta < EPSD ) {
        fprintf(stdout,"  ## Unable to scale mesh.\n");
        return(0);
    }

    /* normalize coordinates */
    dd = 1.0 / info.delta;
    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MS_VOK(ppt) )  continue;
        ppt->c[0] = dd * (ppt->c[0] - info.min[0]);
        ppt->c[1] = dd * (ppt->c[1] - info.min[1]);
        ppt->c[2] = dd * (ppt->c[2] - info.min[2]);
    }

    /* normalize values */
    info.hmin  *= dd;
    info.hmax  *= dd;
    info.hausd *= dd;

    /* normalize sizes */
    if ( met->m ) {
        if ( met->size == 1 ) {
            for (k=1; k<=mesh->np; k++)    met->m[k] *= dd;
        }
        else {
            d1 = 1.0 / (dd*dd);
            for (k=1; k<=6*mesh->np; k++)  met->m[k] *= d1;
        }
    }

    /* normalize local parameters */
    for (k=0; k<info.npar; k++) {
        par = &info.par[k];
        par->hmin  *= dd;
        par->hmax  *= dd;
        par->hausd *= dd;
    }
    return(1);
}

int unscaleMesh(pMesh mesh,pSol met) {
    pPoint     ppt;
    double     dd;
    int        k,i;

    /* de-normalize coordinates */
    dd = info.delta;
    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MS_VOK(ppt) )  continue;
        ppt->c[0] = ppt->c[0] * dd + info.min[0];
        ppt->c[1] = ppt->c[1] * dd + info.min[1];
        ppt->c[2] = ppt->c[2] * dd + info.min[2];
    }

    /* unscale sizes */
    if ( met->m ) {
        if ( met->size == 6 ) {
            dd = 1.0 / (dd*dd);
            for (k=1; k<=mesh->np; k++) {
                ppt = &mesh->point[k];
                if ( !MS_VOK(ppt) )  continue;
                for (i=0; i<6; i++)  met->m[6*(k)+1+i] *= dd;
            }
        }
        else {
            dd = 1.0 / dd;
            for (k=1; k<=mesh->np ; k++) {
                ppt = &mesh->point[k];
                if ( MS_VOK(ppt) )  met->m[k] *= dd;
            }
        }
    }
    return(1);
}
