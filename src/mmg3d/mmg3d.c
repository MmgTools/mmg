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
 *
 * Deallocations before return.
 *
 */
void MMG5_Free_all(MMG5_pMesh mesh,MMG5_pSol met
  ){

  MMG5_Free_structures(mesh,met);
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
  MMG5_Sol       met,disp;
  int       ier;
  char      stim[32];

  fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  _MMG5_Set_APIFunc();

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
  memset(&disp,0,sizeof(MMG5_Sol));

  MMG5_Init_parameters(&mesh);

  met.size      = 1;
  disp.size     = 2;

  /* command line */
  if ( !MMG5_parsar(argc,argv,&mesh,&met) )  return(MMG5_STRONGFAILURE);

#ifdef USE_SCOTCH
  _MMG5_warnScotch(&mesh);
#endif

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&MMG5_ctim[1]);
  _MMG5_warnOrientation(&mesh);
  /* read mesh file */
  if ( MMG5_loadMesh(&mesh) < 1 )  _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);

  /* read displacement if any */
  if ( mesh.info.lag > -1 ) {
    if ( !MMG5_Set_inputSolName(&mesh,&disp,met.namein) )
      exit(EXIT_FAILURE);
    ier = MMG5_loadMet(&mesh,&disp);
    if ( ier == 0 ) {
      fprintf(stdout,"  ## ERROR: NO DISPLACEMENT FOUND.\n");
      _MMG5_RETURN_AND_FREE(&mesh,&disp,MMG5_STRONGFAILURE);
    }
    else if ( ier == -1 ) {
      fprintf(stdout,"  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
      _MMG5_RETURN_AND_FREE(&mesh,&disp,MMG5_STRONGFAILURE);
    }
  }
  /* read metric if any */
  else {
    ier = MMG5_loadMet(&mesh,&met);
    if ( ier == -1 ) {
      fprintf(stdout,"  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
    }
    if ( mesh.info.iso && !ier ) {
      fprintf(stdout,"  ## ERROR: NO ISOVALUE DATA.\n");
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
    }
    if ( !MMG5_parsop(&mesh,&met) )
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_LOWFAILURE);
  }
  chrono(OFF,&MMG5_ctim[1]);
  printim(MMG5_ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&MMG5_ctim[2]);
  _MMG5_setfunc(&mesh,&met);
  MMG5_Set_saveFunc(&mesh);

  _MMG5_outqua(&mesh,&met);
  fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",
          MG_STR,MG_VER,MG_REL,MG_STR);
  if ( mesh.info.imprim )  fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");

  /* scaling mesh */
  if ( mesh.info.lag == -1 ) {
    if ( !_MMG5_scaleMesh(&mesh,&met) )
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
  }
  else {
    if ( !_MMG5_scaleMesh(&mesh,&disp) )
      _MMG5_RETURN_AND_FREE(&mesh,&disp,MMG5_STRONGFAILURE);
  }

  /* specific meshing */
  if ( mesh.info.iso ) {
    if ( !_MMG5_mmg3d2(&mesh,&met) )
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
  }
  else if ( mesh.info.lag >= 0 ) {
    if ( !_MMG5_mmg3d3(&mesh,&disp) )
      _MMG5_RETURN_AND_FREE(&mesh,&disp,MMG5_STRONGFAILURE);
  }
  else {
    if ( !met.np && !_MMG5_DoSol(&mesh,&met) )
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_LOWFAILURE);
  }

  /* mesh analysis */
  if ( !_MMG5_analys(&mesh) )
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_LOWFAILURE);

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
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);


#ifdef PATTERN
  if ( !_MMG5_mmg3d1_pattern(&mesh,&met) ) {
    if ( !(mesh.adja) && !_MMG5_hashTetra(&mesh,1) ) {
      fprintf(stdout,"  ## Hashing problem. Unable to save mesh.\n");
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(&mesh,&met) )
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
    if ( !MMG5_saveMesh(&mesh) )
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
    if ( met.m && !MMG5_saveMet(&mesh,&met) )
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_LOWFAILURE);
  }
#else
  /* Pattern in iso mode, delauney otherwise */
  if ( !mesh.info.iso ) {
    if( !_MMG5_mmg3d1_delone(&mesh,&met) ) {
      if ( !(mesh.adja) && !_MMG5_hashTetra(&mesh,1) ) {
        fprintf(stdout,"  ## Hashing problem. Unable to save mesh.\n");
        _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
      }
      if ( !_MMG5_unscaleMesh(&mesh,&met) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
      if ( !MMG5_saveMesh(&mesh) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
      if ( met.m && !MMG5_saveMet(&mesh,&met) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_LOWFAILURE);
    }
  }
  else {
    if( !_MMG5_mmg3d1_pattern(&mesh,&met) ) {
      if ( !(mesh.adja) && !_MMG5_hashTetra(&mesh,1) ) {
        fprintf(stdout,"  ## Hashing problem. Unable to save mesh.\n");
        _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
      }
      if ( !_MMG5_unscaleMesh(&mesh,&met) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
      if ( !MMG5_saveMesh(&mesh) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
      if ( met.m && !MMG5_saveMet(&mesh,&met) )
        _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_LOWFAILURE);
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
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);

  if ( !MMG5_saveMesh(&mesh) )
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);

  if ( !MMG5_saveMet(&mesh,&met) )
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
  chrono(OFF,&MMG5_ctim[1]);
  if ( mesh.info.imprim )  fprintf(stdout,"  -- WRITING COMPLETED\n");

  /* free mem */
  _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_SUCCESS);
}
