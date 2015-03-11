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
 * \file mmg3d/libmmg3d.c
 * \brief Private API functions for MMG3D library.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \todo documentation doxygen
 *
 * Private API functions for MMG3D library: incompatible functions
 * with the main binary.
 *
 */

#include "mmg3d.h"
#include "shared_func.h"

/**
 * Pack the mesh \a mesh and its associated metric \a met and return \a val.
 */
#define _MMG5_RETURN_AND_PACK(mesh,met,val)do         \
    {                                           \
        MMG5_packMesh(mesh,met);                     \
        return(val);                            \
    }while(0)

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param sing pointer toward the sing structure (only for insertion of
 * singularities mode).
 *
 * Deallocations before return.
 *
 */
void MMG5_Free_all(MMG5_pMesh mesh,MMG5_pSol met
#ifdef SINGUL
                   ,MMG5_pSingul singul
#endif
    ){

#ifdef SINGUL
    MMG5_Free_structures(mesh,met,singul);
#else
    MMG5_Free_structures(mesh,met);
#endif

#ifdef SINGUL
    _MMG5_SAFE_FREE(singul);
#endif
    _MMG5_SAFE_FREE(met);
    _MMG5_SAFE_FREE(mesh);
}

/**
 * \param mesh pointer toward the mesh structure (unused).
 *
 * Set pointer for MMG5_saveMesh function.
 *
 */
void MMG5_Set_saveFunc(MMG5_pMesh mesh) {
    MMG5_saveMesh = _MMG5_saveLibraryMesh;
}

/** Free adja, xtetra and xpoint tables */
static inline
void MMG5_Free_topoTables(MMG5_pMesh mesh) {
    int k;

    mesh->xp = 0;
    if ( mesh->adja )
        _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

    _MMG5_freeXTets(mesh);

    if ( mesh->xpoint )
        _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
    for(k=1; k <=mesh->np; k++) {
        mesh->point[k].xp = 0;
    }

    return;
}

/**
 * \param mesh pointer toward the mesh structure (unused).
 * \param met pointer toward the solution (metric) structure.
 *
 * Pack the sparse mesh and create triangles and edges before getting
 * out of library
 *
 */
static inline
int MMG5_packMesh(MMG5_pMesh mesh,MMG5_pSol met) {
    MMG5_pTetra   pt,ptnew;
    MMG5_pPoint   ppt,pptnew;
    MMG5_hgeom   *ph;
    int     np,nc,nr, k,ne,nbl,imet,imetnew,i;
    int     iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;

    /* compact vertices */
    np = nc = nr = 0;
    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) )  continue;
        ppt->tmp = ++np;
        if ( ppt->tag & MG_CRN )  nc++;
        ppt->ref = abs(ppt->ref);
    }

    /* compact tetrahedra */
    ne  = 0;
    nbl = 1;
    for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( !MG_EOK(pt) )  continue;

        pt->v[0] = mesh->point[pt->v[0]].tmp;
        pt->v[1] = mesh->point[pt->v[1]].tmp;
        pt->v[2] = mesh->point[pt->v[2]].tmp;
        pt->v[3] = mesh->point[pt->v[3]].tmp;
        ne++;
        if ( k!=nbl ) {
            ptnew = &mesh->tetra[nbl];
            memcpy(ptnew,pt,sizeof(MMG5_Tetra));

            iadr = 4*(k-1) + 1;
            adja = &mesh->adja[iadr];
            iadrnew = 4*(nbl-1) + 1;
            adjanew = &mesh->adja[iadrnew];
            for(i=0 ; i<4 ; i++) {
                adjanew[i] = adja[i];
                if(!adja[i]) continue;
                iadrv = 4*(adja[i]/4-1) +1;
                adjav = &mesh->adja[iadrv];
                voy = i;
                adjav[adja[i]%4] = 4*nbl + voy;
            }
        }
        nbl++;
    }
    mesh->ne = ne;

    /* compact metric */
    nbl = 1;
    if ( met->m ) {
        for (k=1; k<=mesh->np; k++) {
            ppt = &mesh->point[k];
            if ( !MG_VOK(ppt) )  continue;
            imet    = (k-1) * met->size + 1;
            imetnew = (nbl-1) * met->size + 1;

            for (i=0; i<met->size; i++)
                met->m[imetnew + i] = met->m[imet + i];
            ++nbl;
        }
    }

    /*compact vertices*/
    np  = 0;
    nbl = 1;
    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) )  continue;
        np++;
        if ( k!=nbl ) {
            pptnew = &mesh->point[nbl];
            memmove(pptnew,ppt,sizeof(MMG5_Point));
            memset(ppt,0,sizeof(MMG5_Point));
            ppt->tag    = MG_NUL;
        }
        nbl++;
    }
    mesh->np = np;
    met->np  = np;

    /* rebuild triangles*/
    mesh->nt = 0;
    _MMG5_chkNumberOfTri(mesh);
    if ( !_MMG5_bdryTria(mesh) ) {
        fprintf(stdout," ## Error: unable to rebuild triangles\n");
        return(0);
    }

    /* build hash table for edges */
    if ( mesh->htab.geom )
        _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));

    mesh->na = 0;
    /* in the wost case (all edges are marked), we will have around 1 edge per *
     * triangle (we count edges only one time) */
    mesh->memCur += (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));
    if ( (mesh->memCur) > (mesh->memMax) ) {
        mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));
        fprintf(stdout,"  ## Warning:");
        fprintf(stdout," unable to allocate htab.\n");
    }
    else if ( _MMG5_hNew(&mesh->htab,mesh->nt,3*(mesh->nt),0) ) {
        for (k=1; k<=mesh->ne; k++) {
            pt   = &mesh->tetra[k];
            if ( MG_EOK(pt) &&  pt->xt ) {
                for (i=0; i<6; i++) {
                    if ( mesh->xtetra[pt->xt].edg[i] ||
                         ( MG_EDG(mesh->xtetra[pt->xt].tag[i] ) ||
                           (mesh->xtetra[pt->xt].tag[i] & MG_REQ) ) )
                        _MMG5_hEdge(mesh,pt->v[_MMG5_iare[i][0]],pt->v[_MMG5_iare[i][1]],
                              mesh->xtetra[pt->xt].edg[i],mesh->xtetra[pt->xt].tag[i]);
                }
            }
        }

        /* edges + ridges + required edges */
        for (k=0; k<=mesh->htab.max; k++) {
            ph = &mesh->htab.geom[k];
            if ( !(ph->a) )  continue;
            mesh->na++;
        }
        if ( mesh->na ) {
            _MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"edges",);
            _MMG5_SAFE_CALLOC(mesh->edge,mesh->na+1,MMG5_Edge);

            mesh->na = 0;
            for (k=0; k<=mesh->htab.max; k++) {
                ph = &mesh->htab.geom[k];
                if ( !ph->a )  continue;
                mesh->na++;
                mesh->edge[mesh->na ].a  = mesh->point[ph->a].tmp;
                mesh->edge[mesh->na ].b  = mesh->point[ph->b].tmp;
                mesh->edge[mesh->na].tag = ( ph->tag | MG_REF ) ;
                mesh->edge[mesh->na].ref = ph->ref;
                if ( MG_GEO & ph->tag ) nr++;
            }
        }
        _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));
    }
    else
        mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));

    for(k=1 ; k<=mesh->np ; k++)
        mesh->point[k].tmp = 0;

    mesh->npnil = mesh->np + 1;
    for(k=mesh->npnil; k<mesh->npmax-1; k++)
        mesh->point[k].tmp  = k+1;

    mesh->nenil = mesh->ne + 1;
    for(k=mesh->nenil; k<mesh->nemax-1; k++)
        mesh->tetra[k].v[3] = k+1;

    /* to could save the mesh, the adjacency have to be correct */
    if ( mesh->info.ddebug && (!_MMG5_chkmsh(mesh,1,1) ) ) {
        fprintf(stdout,"  ##  Problem. Invalid mesh.\n");
        return(0);
    }

    MMG5_Free_topoTables(mesh);

    if ( mesh->info.imprim ) {
        fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",mesh->np,nc);
        if ( mesh->na )
            fprintf(stdout,"     NUMBER OF EDGES      %8d   RIDGES  %8d\n",mesh->na,nr);
        if ( mesh->nt )
            fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",mesh->nt);
        fprintf(stdout,"     NUMBER OF ELEMENTS   %8d\n",mesh->ne);
    }
    return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param sing pointer toward the sing structure (only for insertion of
 * singularities mode).
 * \return Return \ref MMG5_SUCCESS if success.
 * \return Return \ref MMG5_LOWFAILURE if failed but a conform mesh is saved.
 * \return Return \ref MMG5_STRONGFAILURE if failed and we can't save the mesh.
 *
 * Main program for the library .
 *
 */
int MMG5_mmg3dlib(MMG5_pMesh mesh,MMG5_pSol met
#ifdef SINGUL
                  ,MMG5_pSingul sing
#endif
    ) {
    mytime    ctim[TIMEMAX];
    char      stim[32];
#ifdef SINGUL
    int       ier;
#else
    /* sing is not used but must be declared */
    MMG5_pSingul   sing;
    MMG5_Singul    singul;
    sing = &singul;
    memset(sing,0,sizeof(MMG5_Singul));
#endif

    fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",MG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

    signal(SIGABRT,_MMG5_excfun);
    signal(SIGFPE,_MMG5_excfun);
    signal(SIGILL,_MMG5_excfun);
    signal(SIGSEGV,_MMG5_excfun);
    signal(SIGTERM,_MMG5_excfun);
    signal(SIGINT,_MMG5_excfun);

    tminit(ctim,TIMEMAX);
    chrono(ON,&(ctim[0]));

#ifdef USE_SCOTCH
    _MMG5_warnScotch(mesh);
#endif

    fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");
    /* load data */
    chrono(ON,&(ctim[1]));
    _MMG5_warnOrientation(mesh);

    if ( met->np && (met->np != mesh->np) ) {
        fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
        _MMG5_DEL_MEM(mesh,met->m,(met->size*met->npmax+1)*sizeof(double));
        met->np = 0;
    }
    else if ( met->size!=1 ) {
        fprintf(stdout,"  ## ERROR: ANISOTROPIC METRIC NOT IMPLEMENTED.\n");
        return(MMG5_STRONGFAILURE);
    }
#ifdef SINGUL
    if ( mesh->info.sing ) {
        if ( !mesh->info.iso ) {
            if ( !sing->namein )
                fprintf(stdout,"  ## WARNING: NO SINGULARITIES PROVIDED.\n");
        }
        else if ( sing->namein ) {
            fprintf(stdout,"  ## WARNING: SINGULARITIES MUST BE INSERTED IN");
            fprintf(stdout," A PRE-REMESHING PROCESS.\n");
            fprintf(stdout,"              FILE %s IGNORED\n",sing->namein);
        }
    }
#endif

    chrono(OFF,&(ctim[1]));
    printim(ctim[1].gdif,stim);
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

    /* analysis */
    chrono(ON,&(ctim[2]));
    _MMG5_setfunc(mesh,met);
    MMG5_Set_saveFunc(mesh);
    if ( abs(mesh->info.imprim) > 0 )  _MMG5_outqua(mesh,met);
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    if ( mesh->info.imprim )  fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");

    if ( !_MMG5_scaleMesh(mesh,met,sing) ) return(MMG5_STRONGFAILURE);
    if ( mesh->info.iso ) {
        if ( !met->np ) {
            fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
            return(MMG5_STRONGFAILURE);
        }
        if ( !_MMG5_mmg3d2(mesh,met) ) return(MMG5_STRONGFAILURE);
    }

#ifdef SINGUL
    if ( mesh->info.sing ) {
        if ( !mesh->info.iso ) {
            if ( !met->np && !_MMG5_DoSol(mesh,met) )
                _MMG5_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
            if ( !( ier=_MMG5_inserSingul(mesh,met,sing) ) )
                return(MMG5_STRONGFAILURE);
            else if (ier > 0 ) {
                chrono(OFF,&ctim[2]);
                printim(ctim[2].gdif,stim);
                fprintf(stdout,"  -- INSERTION OF SINGULARITIES COMPLETED.     %s\n\n",stim);
                chrono(ON,&ctim[2]);
            }
        }
    }
#endif

#ifdef DEBUG
    if ( !met->np && !_MMG5_DoSol(mesh,met,&mesh->info) ) {
        if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
        _MMG5_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
    }
#endif
    if ( !_MMG5_analys(mesh) ) {
        if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
        _MMG5_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
    }

    if ( mesh->info.imprim > 4 && !mesh->info.iso && met->m ) _MMG5_prilen(mesh,met);

    chrono(OFF,&(ctim[2]));
    printim(ctim[2].gdif,stim);
    if ( mesh->info.imprim )
        fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

    /* mesh adaptation */
    chrono(ON,&(ctim[3]));
    if ( mesh->info.imprim )
        fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");

    /* renumerotation if available */
    if ( !_MMG5_scotchCall(mesh,met) )
    {
        if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
        _MMG5_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
    }

#ifdef SINGUL
    if ( mesh->info.sing && (!mesh->info.iso) ) {
        if ( _MMG5_colSing(mesh,met)<0 ) {
            fprintf(stdout,"  ## Collapse of singularities problem.\n");
            // return(MMG5_STRONGFAILURE);
        }
    }
#endif


#ifdef PATTERN
    if ( !_MMG5_mmg3d1_pattern(mesh,met) ) {
        if ( !(mesh->adja) && !_MMG5_hashTetra(mesh,1) ) {
            fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
            return(MMG5_STRONGFAILURE);
        }
        if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
        _MMG5_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
    }
#else
    /** Patterns in iso mode, delauney otherwise */
    if ( !mesh->info.iso ) {
        if ( !_MMG5_mmg3d1_delone(mesh,met) ) {
            if ( !(mesh->adja) && !_MMG5_hashTetra(mesh,1) ) {
                fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
                return(MMG5_STRONGFAILURE);
            }
            if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
            _MMG5_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
        }
    }
    else {
        if ( !_MMG5_mmg3d1_pattern(mesh,met) ) {
            if ( !(mesh->adja) && !_MMG5_hashTetra(mesh,1) ) {
                fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
                return(MMG5_STRONGFAILURE);
            }
            if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
            _MMG5_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
        }
    }

#endif

#ifdef SINGUL
    if ( mesh->info.sing && (!mesh->info.iso) ) {
        if ( !_MMG5_solveUnsignedTet(mesh,met) ) {
            fprintf(stdout,"  ## Solve of undetermined tetrahedra problem.\n");
            if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
            _MMG5_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
        }
    }
#endif

    chrono(OFF,&(ctim[3]));
    printim(ctim[3].gdif,stim);
    if ( mesh->info.imprim )
        fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG3d: IMB-LJLL \n  %s\n",MG_STR,MG_STR);

    /* save file */
    _MMG5_outqua(mesh,met);
    if ( mesh->info.imprim > 4 && !mesh->info.iso )
        _MMG5_prilen(mesh,met);

    chrono(ON,&(ctim[1]));
    if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");
    if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
    if ( !MMG5_packMesh(mesh,met) )     return(MMG5_STRONGFAILURE);
    chrono(OFF,&(ctim[1]));

    chrono(OFF,&ctim[0]);
    printim(ctim[0].gdif,stim);
    fprintf(stdout,"\n   MMG3DLIB: ELAPSED TIME  %s\n",stim);
    return(MMG5_SUCCESS);
}
