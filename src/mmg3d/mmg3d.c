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
 * \file mmg3d/mmg3d.c
 * \brief Main file for MMG3D executable: perform 3d mesh adaptation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg3d.h"
#include "shared_func.h"

mytime         MMG5_ctim[TIMEMAX];

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

}

/**
 * \param mesh pointer toward the mesh structure (unused).
 *
 * Set pointer for MMG5_saveMesh function.
 *
 */
void MMG5_Set_saveFunc(MMG5_pMesh mesh) {
    MMG5_saveMesh = _MMG5_saveAllMesh;
}

/**
 *
 * Messages at end of the code.
 *
 */
static void _MMG5_endcod() {
    char    stim[32];

    chrono(OFF,&MMG5_ctim[0]);
    printim(MMG5_ctim[0].gdif,stim);
    fprintf(stdout,"\n   ELAPSED TIME  %s\n",stim);
}

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \return \ref MMG5_SUCCESS if success.
 * \return \ref MMG5_LOWFAILURE if failed but a conform mesh is saved.
 * \return \ref MMG5_STRONGFAILURE if failed and we can't save the mesh.
 *
 * Main program for MMG3D executable: perform mesh adaptation.
 *
 */
int main(int argc,char *argv[]) {
    MMG5_Mesh      mesh;
    MMG5_Sol       met;
    MMG5_Singul    sing;
    int       ier;
    char      stim[32];

    fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",MG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

    signal(SIGABRT,_MMG5_excfun);
    signal(SIGFPE,_MMG5_excfun);
    signal(SIGILL,_MMG5_excfun);
    signal(SIGSEGV,_MMG5_excfun);
    signal(SIGTERM,_MMG5_excfun);
    signal(SIGINT,_MMG5_excfun);
    atexit(_MMG5_endcod);

    tminit(MMG5_ctim,TIMEMAX);
    chrono(ON,&MMG5_ctim[0]);

    /* assign default values */
    memset(&mesh,0,sizeof(MMG5_Mesh));
    memset(&met,0,sizeof(MMG5_Sol));
#ifdef SINGUL
    memset(&sing,0,sizeof(MMG5_Singul));
#endif

    MMG5_Init_parameters(&mesh);

    met.size      = 1;

    /* command line */
#ifdef SINGUL
    if ( !MMG5_parsar(argc,argv,&mesh,&met,&sing) )  return(MMG5_STRONGFAILURE);
#else
    if ( !MMG5_parsar(argc,argv,&mesh,&met) )  return(MMG5_STRONGFAILURE);
#endif
#ifdef USE_SCOTCH
    _MMG5_warnScotch(&mesh);
#endif

    /* load data */
    fprintf(stdout,"\n  -- INPUT DATA\n");
    chrono(ON,&MMG5_ctim[1]);
    _MMG5_warnOrientation(&mesh);
    /* read mesh file */
    if ( !MMG5_loadMesh(&mesh) ) _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);

    /* read metric if any */
    ier = MMG5_loadMet(&mesh,&met);
    if ( !ier )
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
    else if ( ier > 0 && met.np != mesh.np ) {
        fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
        _MMG5_DEL_MEM(&mesh,met.m,(met.size*met.npmax+1)*sizeof(double));
        met.np = 0;
    } else if ( met.size!=1 ) {
        fprintf(stdout,"  ## ERROR: ANISOTROPIC METRIC NOT IMPLEMENTED.\n");
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
    }
#ifdef SINGUL
    if ( mesh.info.sing ) {
        if ( !mesh.info.iso ) {
            if ( !sing.namein )
                fprintf(stdout,"  ## WARNING: NO SINGULARITIES PROVIDED.\n");
            else
                if ( !MMG5_loadSingul(&mesh,&sing) )
                    _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
        }
        else if ( sing.namein ) {
            fprintf(stdout,"  ## WARNING: SINGULARITIES MUST BE INSERTED IN");
            fprintf(stdout," A PRE-REMESHING PROCESS.\n");
            fprintf(stdout,"              FILE %s IGNORED\n",sing.namein);
        }
    }
#endif
    if ( !MMG5_parsop(&mesh,&met) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_LOWFAILURE);

    chrono(OFF,&MMG5_ctim[1]);
    printim(MMG5_ctim[1].gdif,stim);
    fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

    /* analysis */
    chrono(ON,&MMG5_ctim[2]);
    _MMG5_setfunc(&mesh,&met);
    MMG5_Set_saveFunc(&mesh);

    if ( abs(mesh.info.imprim) > 0 )  _MMG5_outqua(&mesh,&met);
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",
            MG_STR,MG_VER,MG_REL,MG_STR);
    if ( mesh.info.imprim )  fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");

    if ( !_MMG5_scaleMesh(&mesh,&met,&sing) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
    if ( mesh.info.iso ) {
        if ( !met.np ) {
            fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
            _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
        }
        if ( !_MMG5_mmg3d2(&mesh,&met) )
            _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
    }

#ifdef SINGUL
    if ( mesh.info.sing ) {
        if ( !mesh.info.iso ) {
            if ( !met.np && !_MMG5_DoSol(&mesh,&met) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_LOWFAILURE);
            if ( !( ier=_MMG5_inserSingul(&mesh,&met,&sing) ) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            else if (ier > 0 ) {
                chrono(OFF,&MMG5_ctim[2]);
                printim(MMG5_ctim[2].gdif,stim);
                fprintf(stdout,"  -- INSERTION OF SINGULARITIES COMPLETED.     %s\n\n",stim);
                chrono(ON,&MMG5_ctim[2]);
            }
        }
    }
#endif

    if ( !mesh.info.iso && !met.np && !_MMG5_DoSol(&mesh,&met) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_LOWFAILURE);

    if ( !_MMG5_analys(&mesh) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_LOWFAILURE);

    if ( mesh.info.imprim > 3 && !mesh.info.iso && met.m ) _MMG5_prilen(&mesh,&met);

    chrono(OFF,&MMG5_ctim[2]);
    printim(MMG5_ctim[2].gdif,stim);
    if ( mesh.info.imprim )
        fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

    /* mesh adaptation */
    chrono(ON,&MMG5_ctim[3]);
    if ( mesh.info.imprim )
        fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met.size < 6 ? "ISOTROPIC" : "ANISOTROPIC");

    /* renumerotation if available */
    if ( !_MMG5_scotchCall(&mesh,&met) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);

#ifdef SINGUL
    if ( mesh.info.sing && (!mesh.info.iso) ) {
        if ( _MMG5_colSing(&mesh,&met)<0 ) {
            fprintf(stdout,"  ## Collapse of singularities problem.\n");
            // _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
        }
    }
#endif


#ifdef PATTERN
    if ( !_MMG5_mmg3d1_pattern(&mesh,&met) ) {
        if ( !(mesh.adja) && !_MMG5_hashTetra(&mesh,1) ) {
            fprintf(stdout,"  ## Hashing problem. Unable to save mesh.\n");
            _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
        }
        if ( !_MMG5_unscaleMesh(&mesh,&met) )
            _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
        if ( !MMG5_saveMesh(&mesh) )
            _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
        if ( met.m && !MMG5_saveMet(&mesh,&met) )
            _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_LOWFAILURE);
    }
#else
    /* Pattern in iso mode, delauney otherwise */
    if ( !mesh.info.iso ) {
        if( !_MMG5_mmg3d1_delone(&mesh,&met) ) {
            if ( !(mesh.adja) && !_MMG5_hashTetra(&mesh,1) ) {
                fprintf(stdout,"  ## Hashing problem. Unable to save mesh.\n");
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            }
            if ( !_MMG5_unscaleMesh(&mesh,&met) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            if ( !MMG5_saveMesh(&mesh) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            if ( met.m && !MMG5_saveMet(&mesh,&met) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_LOWFAILURE);
        }
    }
    else {
        if( !_MMG5_mmg3d1_pattern(&mesh,&met) ) {
            if ( !(mesh.adja) && !_MMG5_hashTetra(&mesh,1) ) {
                fprintf(stdout,"  ## Hashing problem. Unable to save mesh.\n");
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            }
            if ( !_MMG5_unscaleMesh(&mesh,&met) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            if ( !MMG5_saveMesh(&mesh) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            if ( met.m && !MMG5_saveMet(&mesh,&met) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_LOWFAILURE);
        }
    }
#endif


#ifdef SINGUL
    if ( mesh.info.sing && (!mesh.info.iso) ) {
        if ( !_MMG5_solveUnsignedTet(&mesh,&met) ) {
            fprintf(stdout,"  ## Solve of undetermined tetrahedra problem.\n");
            if ( !_MMG5_unscaleMesh(&mesh,&met) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            if ( !MMG5_saveMesh(&mesh) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            if ( met.m && !MMG5_saveMet(&mesh,&met) )
                _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
            _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_LOWFAILURE);
        }
    }
#endif

    chrono(OFF,&MMG5_ctim[3]);
    printim(MMG5_ctim[3].gdif,stim);
    if ( mesh.info.imprim )
        fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG3d: IMB-LJLL \n  %s\n",MG_STR,MG_STR);

    /* save file */
    _MMG5_outqua(&mesh,&met);

    if ( mesh.info.imprim > 3 && !mesh.info.iso )
        _MMG5_prilen(&mesh,&met);

    chrono(ON,&MMG5_ctim[1]);
    if ( mesh.info.imprim )  fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh.nameout);
    if ( !_MMG5_unscaleMesh(&mesh,&met) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);

    if ( !MMG5_saveMesh(&mesh) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);

    if ( !MMG5_saveMet(&mesh,&met) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_STRONGFAILURE);
    chrono(OFF,&MMG5_ctim[1]);
    if ( mesh.info.imprim )  fprintf(stdout,"  -- WRITING COMPLETED\n");

    /* free mem */
    _MMG5_RETURN_AND_FREE(&mesh,&met,&sing,MMG5_SUCCESS);
}
