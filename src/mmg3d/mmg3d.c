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

mytime         MMG5_ctim[TIMEMAX];


/**
 * Print elapsed time at end of process.
 */
static void _MMG5_endcod() {
  char   stim[32];

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

  MMG5_pMesh      mesh;
  MMG5_pSol       met,disp;
  int             ier;
  char            stim[32];

  fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  _MMG3D_Set_commonFunc();

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
  mesh = NULL;
  met  = NULL;
  disp = NULL;

  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                  MMG5_ARG_ppDisp,&disp,
                  MMG5_ARG_end);
  /* reset default values for file names */
  MMG3D_Free_names(MMG5_ARG_start,
                   MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                   MMG5_ARG_ppDisp,&disp,
                   MMG5_ARG_end);

  /* command line */
  if ( !MMG3D_parsar(argc,argv,mesh,met) )  return(MMG5_STRONGFAILURE);

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&MMG5_ctim[1]);
  /* read mesh file */
  if ( MMG3D_loadMesh(mesh,mesh->namein)<1 )
    _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_STRONGFAILURE);

  /* Set default metric size */
  if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Scalar) )
    _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_STRONGFAILURE);

  /* read displacement if any */
  if ( mesh->info.lag > -1 ) {
    if ( !MMG3D_Set_solSize(mesh,disp,MMG5_Vertex,0,MMG5_Vector) )
      _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_STRONGFAILURE);

    if ( !MMG3D_Set_inputSolName(mesh,disp,met->namein) )
      _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_STRONGFAILURE);

    ier = MMG3D_loadSol(mesh,disp,disp->namein);
    if ( ier == 0 ) {
      fprintf(stdout,"  ## ERROR: NO DISPLACEMENT FOUND.\n");
      _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_STRONGFAILURE);
    }
    else if ( ier == -1 ) {
      fprintf(stdout,"  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
      _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_STRONGFAILURE);
    }
  }
  /* read metric if any */
  else {
    ier = MMG3D_loadSol(mesh,met,met->namein);
    if ( ier == -1 ) {
      fprintf(stdout,"  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
      _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_STRONGFAILURE);
    }
    if ( mesh->info.iso && !ier ) {
      fprintf(stdout,"  ## ERROR: NO ISOVALUE DATA.\n");
      _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_STRONGFAILURE);
    }
    if ( !MMG3D_parsop(mesh,met) )
      _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_LOWFAILURE);
  }

  chrono(OFF,&MMG5_ctim[1]);
  printim(MMG5_ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  if ( mesh->info.lag > -1 ) {
    ier = MMG3D_mmg3dmov(mesh,met,disp);
  }
  else if ( mesh->info.iso ) {
    ier = MMG3D_mmg3dls(mesh,met);
  }
  else {
    ier = MMG3D_mmg3dlib(mesh,met);
  }

  if ( ier != MMG5_STRONGFAILURE ) {
    chrono(ON,&MMG5_ctim[1]);
    if ( mesh->info.imprim )
      fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh->nameout);
    if ( !MMG3D_saveMesh(mesh,mesh->nameout) )
      _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_STRONGFAILURE);

    if ( !MMG3D_saveSol(mesh,met,met->nameout) )
      _MMG5_RETURN_AND_FREE(mesh,met,disp,MMG5_STRONGFAILURE);

    chrono(OFF,&MMG5_ctim[1]);
    if ( mesh->info.imprim )
      fprintf(stdout,"  -- WRITING COMPLETED\n");
  }

  /* free mem */
  _MMG5_RETURN_AND_FREE(mesh,met,disp,ier);
}
