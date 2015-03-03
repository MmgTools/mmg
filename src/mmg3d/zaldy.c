/* =============================================================================
**  This file is part of the MMG3D 5 software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
**
**  MMG3D 5 is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  MMG3D 5 is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with MMG3D 5 (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the MMG3D 5 distribution only if you accept them.
** =============================================================================
*/

/**
 * \file mmg3d/zaldy.c
 * \brief Memory management.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmg3d.h"

/** get new point address */
int newPt(pMesh mesh,double c[3],char tag) {
    pPoint  ppt;
    int     curpt;

    if ( !mesh->npnil )  return(0);
    curpt = mesh->npnil;
    if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
    ppt   = &mesh->point[curpt];
    memcpy(ppt->c,c,3*sizeof(double));
    mesh->npnil = ppt->tmp;
    ppt->tmp    = 0;

    ppt->ref = 0;
    ppt->xp = 0;
    ppt->flag = 0;
    /* point on geometry */
    if ( tag & MG_BDY ) {
        mesh->xp++;
        if(mesh->xp > mesh->xpmax){
            /* reallocation of xpoint table */
            TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,xPoint,
                         "larger xpoint table",
                         return(0));
        }
        ppt->xp  = mesh->xp;
    }
    assert(tag < 127);
    assert(tag >= 0);
    ppt->tag = tag;
    ppt->tagdel = 0;
    return(curpt);
}

void delPt(pMesh mesh,int ip) {
    pPoint   ppt;
    xPoint  *pxp;

    ppt = &mesh->point[ip];
    if ( ppt->xp ) {
        pxp = &mesh->xpoint[ppt->xp];
        memset(pxp,0,sizeof(xPoint));
    }
    memset(ppt,0,sizeof(Point));
    ppt->tag    = MG_NUL;
    ppt->tmp    = mesh->npnil;
    mesh->npnil = ip;
    if ( ip == mesh->np ) {
        while ( !MG_VOK((&mesh->point[mesh->np])) )  mesh->np--;
    }
}

/** get new elt address */
int newElt(pMesh mesh) {
    int     curiel;

    if ( !mesh->nenil )  return(0);
    curiel = mesh->nenil;

    if ( mesh->nenil > mesh->ne )  mesh->ne = mesh->nenil;
    mesh->nenil = mesh->tetra[curiel].v[3];
    mesh->tetra[curiel].v[3] = 0;
    mesh->tetra[curiel].mark=0;

    return(curiel);
}


void delElt(pMesh mesh,int iel) {
    pTetra   pt;
    int      iadr;

    pt = &mesh->tetra[iel];
    if ( !MG_EOK(pt) ) {
        fprintf(stdout,"  ## INVALID ELEMENT %d.\n",iel);
        exit(EXIT_FAILURE);
    }
    memset(pt,0,sizeof(Tetra));
    pt->v[3] = mesh->nenil;
    iadr = 4*(iel-1) + 1;
    if ( mesh->adja )
        memset(&mesh->adja[iadr],0,4*sizeof(int));
    mesh->nenil = iel;
    if ( iel == mesh->ne ) {
        while ( !MG_EOK((&mesh->tetra[mesh->ne])) )  mesh->ne--;
    }
}

long long memSize (void) {
    long long mem;

#if (defined(__APPLE__) && defined(__MACH__))
    size_t size;

    size = sizeof(mem);
    if ( sysctlbyname("hw.memsize",&mem,&size,NULL,0) == -1)
        return(0);

#elif defined(__unix__) || defined(__unix) || defined(unix)
    mem = ((long long)sysconf(_SC_PHYS_PAGES))*
        ((long long)sysconf(_SC_PAGESIZE));
#else
    printf("  ## WARNING: UNKNOWN SYSTEM, RECOVER OF MAXIMAL MEMORY NOT AVAILABLE.\n");
    return(0);
#endif

    return(mem);
}
/** memory repartition for the -m option */
void memOption(pMesh mesh) {
    long long  million = 1048576L;
    int        ctri,npask,bytes,memtmp;

    mesh->memMax = memSize();

    mesh->npmax = MG_MAX(1.5*mesh->np,NPMAX);
    mesh->nemax = MG_MAX(1.5*mesh->ne,NEMAX);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,NTMAX);

    if ( mesh->info.mem <= 0 ) {
        if ( mesh->memMax )
            /* maximal memory = 50% of total physical memory */
            mesh->memMax = mesh->memMax*50/100;
        else {
            /* default value = 800 Mo */
            printf("  Maximum memory set to default value: %d Mo.\n",MEMMAX);
            mesh->memMax = MEMMAX*million;
        }
    }
    else {
        /* memory asked by user if possible, otherwise total physical memory */
        if ( (long long)mesh->info.mem*million > mesh->memMax && mesh->memMax ) {
            fprintf(stdout,"  ## Warning: asking for %d Mo of memory ",mesh->info.mem);
            fprintf(stdout,"when only %lld available.\n",mesh->memMax/million);
        }
        else {
            mesh->memMax= (long long)(mesh->info.mem)*million;
        }

        /* if asked memory is lower than default NPMAX/NEMAX/NTMAX we take lower values */
#ifdef SINGUL
        /* Remarks:
         * 1-- in insertion part, we have memory allocated to store *
         * edges and singular points (in Singul) but we don't need to take this *
         * into account because xpoints and xtetra are free and need more memory *
         * 2-- we need more xtetra so we increase the memory to save for triangles */
        if ( mesh->info.sing )
            ctri = 4;
        else
            ctri = 2;
#else
        ctri = 2;
#endif

        /* Euler-poincare: ne = 6*np; nt = 2*np; na = np/5 *
         * point+tria+tets+adja+adjt+sol+item *
         * warning: we exceed memory in saveMesh when we call hNew */
        bytes = sizeof(Point) + sizeof(xPoint) +
            6*sizeof(Tetra) + ctri*sizeof(xTetra) +
            4*6*sizeof(int) + ctri*3*sizeof(int) +
            sizeof(Sol)+4*sizeof(hedge);
#ifdef USE_SCOTCH
        /* bytes = bytes + vertTab + edgeTab + PermVrtTab *
         * + vertOldTab + sortPartTab - adja */
        bytes = bytes + 3*6*sizeof(int);
#endif

        /*init allocation need 38Mo*/
        npask = (double)(mesh->info.mem-38) / bytes * (int)million;
        mesh->npmax = MG_MIN(npask,mesh->npmax);
        mesh->ntmax = MG_MIN(ctri*npask,mesh->ntmax);
        mesh->nemax = MG_MIN(6*npask,mesh->nemax);
        /*check if the memory asked is enough to load the mesh*/
        if(mesh->np &&
           (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne)) {
            memtmp = mesh->np * bytes /(int)million + 38;
            memtmp = MG_MAX(memtmp, mesh->nt * bytes /(ctri* (int)million) + 38);
            memtmp = MG_MAX(memtmp, mesh->ne * bytes /(6*(int)million) + 38);
            mesh->memMax = (long long) memtmp+1;
            fprintf(stdout,"  ## ERROR: asking for %d Mo of memory ",mesh->info.mem);
            fprintf(stdout,"is not enough to load mesh. You need to ask %d Mo minimum\n",
                    memtmp+1);
        }
        if(mesh->info.mem < 39) {
            mesh->memMax = (long long) 39;
            fprintf(stdout,"  ## ERROR: asking for %d Mo of memory ",mesh->info.mem);
            fprintf(stdout,"is not enough to load mesh. You need to ask %d Mo minimum\n",
                    39);
        }
    }

    if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
        fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED (Mo)    %lld\n",
                mesh->memMax/million);

    if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
        fprintf(stdout,"  NPMAX    %d\n",mesh->npmax);
        fprintf(stdout,"  NTMAX    %d\n",mesh->ntmax);
        fprintf(stdout,"  NEMAX    %d\n",mesh->nemax);
    }

    return;
}

/** allocate main structure */
int zaldy(pMesh mesh) {
    int     k;

    memOption(mesh);

    ADD_MEM(mesh,(mesh->npmax+1)*sizeof(Point),"initial vertices",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(mesh->point,mesh->npmax+1,Point);

    ADD_MEM(mesh,(mesh->nemax+1)*sizeof(Tetra),"initial tetrahedra",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(mesh->tetra,mesh->nemax+1,Tetra);

    if ( mesh->nt ) {
        ADD_MEM(mesh,(mesh->nt+1)*sizeof(Tria),"initial triangles",return(0));
        SAFE_CALLOC(mesh->tria,mesh->nt+1,Tria);
        memset(&mesh->tria[0],0,sizeof(Tria));
    }
    if ( mesh->na ) {
        ADD_MEM(mesh,(mesh->na+1)*sizeof(Edge),"initial edges",return(0));
        SAFE_CALLOC(mesh->edge,(mesh->na+1),Edge);
    }

    /* keep track of empty links */
    mesh->npnil = mesh->np + 1;
    mesh->nenil = mesh->ne + 1;

    for (k=mesh->npnil; k<mesh->npmax-1; k++)
        mesh->point[k].tmp  = k+1;

    for (k=mesh->nenil; k<mesh->nemax-1; k++)
        mesh->tetra[k].v[3] = k+1;

    return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Free xtetra structure.
 *
 */
void freeXTets(pMesh mesh) {
    pTetra pt;
    int    k;

    for (k=1; k<=mesh->ne; k++) {
        pt     = &mesh->tetra[k];
        pt->xt = 0;
    }
    if ( mesh->xtetra )
        DEL_MEM(mesh,mesh->xtetra,(mesh->xtmax+1)*sizeof(xTetra));
    mesh->xt = 0;
}
