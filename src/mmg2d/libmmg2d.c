/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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
#include "mmg2d.h"

/**
 * Pack the mesh \a mesh and its associated metric \a met and return \a val.
 */
#define _MMG2D_RETURN_AND_PACK(mesh,met,val)do                          \
  {                                                                     \
  if ( !MMG2_pack(mesh,met) ) {                                         \
    mesh->npi = mesh->np;                                               \
    mesh->nti = mesh->nt;                                               \
    mesh->nai = mesh->na;                                               \
    mesh->nei = mesh->ne;                                               \
    met->npi  = met->np;                                                \
    return(MMG5_LOWFAILURE);                                            \
  }                                                                     \
    _LIBMMG5_RETURN(mesh,met,val);                                      \
    }while(0)


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the solution structure.
 *
 * Truncate a scalar metric by hmax and hmin values.
 *
 */
static inline
void _MMG2D_scalarSolTruncature(MMG5_pMesh mesh, MMG5_pSol met) {
  int         k;

  /* vertex size */
  for (k=1; k<=mesh->np; k++) {
    met->m[k] = MG_MIN(mesh->info.hmax,MG_MAX(mesh->info.hmin,met->m[k]));
  }
  return;
}

int MMG2D_mmg2dlib(MMG5_pMesh mesh,MMG5_pSol sol)
{
  mytime    ctim[TIMEMAX];
  char      stim[32];

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  /*uncomment to callback*/
  //MMG2D_callbackinsert = titi;
  /* interrupts */
  signal(SIGABRT,_MMG2_excfun);
  signal(SIGFPE,_MMG2_excfun);
  signal(SIGILL,_MMG2_excfun);
  signal(SIGSEGV,_MMG2_excfun);
  signal(SIGTERM,_MMG2_excfun);
  signal(SIGINT,_MMG2_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));
  
  /* Check options */
  if ( !mesh->nt ) {
    fprintf(stdout,"\n  ## ERROR: NO TRIANGLES IN THE MESH. \n");
    fprintf(stdout,"          To generate a mesh from boundaries call the"
            " MMG2D_mmg2dmesh function\n.");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.iso ) {
    fprintf(stdout,"  ## Error: level-set discretisation unavailable"
            " (MMG2D_IPARAM_iso):\n"
            "          You must call the MMG2D_mmg2dls function to use this option.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.lag >= 0 ) {
    fprintf(stdout,"  ## Error: lagrangian mode unavailable (MMG2D_IPARAM_lag):\n"
            "            You must call the MMG2D_mmg2dmov function to move a rigidbody.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG2DLIB: INPUT DATA\n");

  /* Load data */
  chrono(ON,&(ctim[1]));

  sol->np   = mesh->np;
  sol->ver  = mesh->ver;

  /* Allocate memory if no metric is supplied */
  if ( !sol->m ) {
    _MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax+1))*sizeof(double),
                  "initial solution",return(MMG5_STRONGFAILURE));
    _MMG5_SAFE_CALLOC(sol->m,sol->size*(mesh->npmax+1),double,MMG5_STRONGFAILURE);
    sol->np = 0;
  }
  else if ( sol->np && ( sol->np != mesh->np ) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER : %d != %d\n",sol->np,mesh->np);
    //exit(EXIT_FAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);

  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* /\* default values *\/ */
  /* mesh->info.imprim = opt[5]; */
  /* mesh->info.mem    = 0; */
  /* mesh->info.ddebug = opt[1]; */
  /* mesh->info.iso = 0; */
  /* mesh->info.lag = -1; */
  /* mesh->info.hmin = -1; */
  /* mesh->info.hmax = -1; */
  /* mesh->info.hausd = 0.01; */
  /* switch(opt[0]) { */
  /* case 0: */
  /* case 1: */
  /* case 2: */
  /*   break; */
  /* case 6: */
  /*   mesh->info.iso = 1; */
  /*   break; */
  /* case 9: */
  /* case 99: */
  /*   mesh->info.lag = 0; */
  /*   break; */
  /* default: */
  /*   fprintf(stdout,"option not recognized %d\n",opt[0]); */
  /*   exit(EXIT_FAILURE); */
  /* } */
  /* mesh->info.noswap = opt[2]; */
  /* mesh->info.nomove = opt[4]; */
  /* mesh->info.noinsert = opt[3]; */
  /* mesh->info.hgrad  = optdbl[0]; */
  /* if(opt[6]) */
  /*   mesh->info.dhd  = -1; */
  /* else */
  /*   mesh->info.dhd  = 180.-optdbl[1]; */
  /* /\*this options are not used inside library version*\/ */
  /* //qdegrad[0] = 10./_MMG2D_ALPHA; */
  /* //qdegrad[1] = 1.3; */
  /* mesh->info.renum = 0; */

  /* sol->type = 1; */

  /* Set function pointers */
  MMG2D_setfunc(mesh,sol);
  _MMG2D_Set_commonFunc();

  fprintf(stdout,"\n  %s\n   MODULE MMG2D-IMB/LJLL : %s (%s) %s\n  %s\n",
          MG_STR,MG_VER,MG_REL,sol->size == 1 ? "ISO" : "ANISO",MG_STR);
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM NUMBER OF POINTS    (NPMAX) : %8d\n",mesh->npmax);
    fprintf(stdout,"  MAXIMUM NUMBER OF TRIANGLES (NTMAX) : %8d\n",mesh->ntmax);
  }

  /* Data analysis */
  chrono(ON,&ctim[2]);
  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : DATA ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");

  /* Scale input mesh */
  if ( !MMG2_scaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  if ( !sol->np && mesh->info.optim ) {
    if ( !MMG2D_doSol(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
    _MMG2D_scalarSolTruncature(mesh,sol);
  }

  /* Create adjacency relations in the mesh */
  if ( !MMG2_hashTria(mesh) )
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  if ( mesh->info.ddebug && !_MMG5_chkmsh(mesh,1,0) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /* Print initial quality history */
  MMG2_outqua(mesh,sol);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* remeshing */
  chrono(ON,&ctim[3]);

  if ( mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 2 : MESH ADAPTATION\n");

  /* Mesh analysis */
  if (! _MMG2_analys(mesh) )
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);


  /* Mesh improvement */
  if ( !MMG2_mmg2d1n(mesh,sol) ) {
    if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 3 COMPLETED.     %s\n",stim);
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG2D: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  /* Unscale mesh */
  if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /* Print output quality history */
  MMG2_outqua(mesh,sol);

  /* Print edge length histories */
  if ( abs(mesh->info.imprim) > 1 )  {
    MMG2_prilen(mesh,sol);
  }

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  if (!MMG2_pack(mesh,sol) ) _LIBMMG5_RETURN(mesh,sol,MMG5_LOWFAILURE);

  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMG2DLIB: ELAPSED TIME  %s\n",stim);
  _LIBMMG5_RETURN(mesh,sol,MMG5_SUCCESS);

}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail (lack of memory), 1 otherwise.
 *
 * Clean the mesh structure when we just call the MMG2D_Free_Triangles and
 * MMG2D_Free_Edges functions between 2 call of the MMG2D_mmg2dmesh function:
 *   - Allocate the tria and edge structures if needed;
 *   - Reset the tags at vertices.
 *
 */
static inline
int _MMG2D_restart(MMG5_pMesh mesh){
  int k;

 /** If needed, reallocate the missing structures */
  if ( !mesh->tria ) {
    /* If we call the library more than one time and if we free the triangles
     * using the MMG2D_Free_triangles function we need to reallocate it */
    _MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),
                  "initial triangles",return(0));
    _MMG5_SAFE_CALLOC(mesh->tria,mesh->ntmax+1,MMG5_Tria,0);
    mesh->nenil = mesh->nt + 1;
    for ( k=mesh->nenil; k<mesh->ntmax-1; k++) {
      mesh->tria[k].v[2] = k+1;
    }
  }
  if ( mesh->na && !mesh->edge ) {
    /* If we call the library more than one time and if we free the triangles
     * using the MMG2D_Free_triangles function we need to reallocate it */
    _MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),
                  "initial edges",return(0));
    _MMG5_SAFE_CALLOC(mesh->edge,mesh->namax+1,MMG5_Edge,0);
    mesh->nanil = mesh->na + 1;
  }

  for ( k=1; k<=mesh->np;  ++k ) {
    mesh->point[k].tag = 0;
  }
  return 1;
}


int MMG2D_mmg2dmesh(MMG5_pMesh mesh,MMG5_pSol sol) {
  mytime    ctim[TIMEMAX];
  char      stim[32];

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  /*uncomment for callback*/
  //MMG2D_callbackinsert = titi;

  /* interrupts */
  signal(SIGABRT,_MMG2_excfun);
  signal(SIGFPE,_MMG2_excfun);
  signal(SIGILL,_MMG2_excfun);
  signal(SIGSEGV,_MMG2_excfun);
  signal(SIGTERM,_MMG2_excfun);
  signal(SIGINT,_MMG2_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->nt ) {
    fprintf(stdout,"  ## Error: your mesh contains already triangles.\n"
            " The mesh generation option is unavailable.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }

  else if ( mesh->info.iso ) {
    fprintf(stdout,"  ## Error: level-set discretisation unavailable"
            " (MMG2D_IPARAM_iso):\n"
            "          You must call the MMG2D_mmg2dls function to use this option.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }

  else if ( mesh->info.lag >= 0 ) {
    fprintf(stdout,"  ## Error: lagrangian mode unavailable (MMG2D_IPARAM_lag):\n"
            "            You must call the MMG2D_mmg2dmov function to move a rigidbody.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG2DMESH: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  sol->ver  = mesh->ver;

  if ( !sol->m ) {
    /* mem alloc */
    _MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax+1))*sizeof(double),
                  "initial solution",return(MMG5_STRONGFAILURE));
    _MMG5_SAFE_CALLOC(sol->m,sol->size*(mesh->npmax+1),double,MMG5_STRONGFAILURE);
    sol->np = 0;
  }

  else   if ( sol->np && (sol->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER : %d != %d\n",sol->np,mesh->np);
    //exit(EXIT_FAILURE);
  }
  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* Create function pointers */
  MMG2D_setfunc(mesh,sol);
  _MMG2D_Set_commonFunc();

  fprintf(stdout,"\n  %s\n   MODULE MMG2D-IMB/LJLL : %s (%s) %s\n  %s\n",
          MG_STR,MG_VER,MG_REL,sol->size == 1 ? "ISO" : "ANISO",MG_STR);
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM NUMBER OF POINTS    (NPMAX) : %8d\n",mesh->npmax);
    fprintf(stdout,"  MAXIMUM NUMBER OF TRIANGLES (NTMAX) : %8d\n",mesh->ntmax);
  }

  /* analysis */
  chrono(ON,&ctim[2]);

  if ( !_MMG2D_restart(mesh) ) {
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }

  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : DATA ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");
  if ( !MMG2_scaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  if ( mesh->info.ddebug && !_MMG5_chkmsh(mesh,1,0) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* Specific meshing */
  chrono(ON,&ctim[3]);

  /* Memory alloc */
  _MMG5_ADD_MEM(mesh,(3*mesh->ntmax+5)*sizeof(int),"adjacency table",
                printf("  Exit program.\n");
                return(MMG5_STRONGFAILURE));
  _MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,int,MMG5_STRONGFAILURE);

  /* Delaunay triangulation of the set of points contained in the mesh, enforcing the edges of the mesh */
  if ( mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 2 : MESH GENERATION\n");

  if ( !MMG2_mmg2d2(mesh,sol) )  {
    if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);

  /* remeshing */
  chrono(ON,&ctim[4]);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 3 : MESH IMPROVEMENT\n");
  }

  /* Mesh analysis */
  if (! _MMG2_analys(mesh) )
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /* Mesh improvement - call new version of mmg2d1 */
  if ( !MMG2_mmg2d1n(mesh,sol) ) {
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[4]));
  printim(ctim[5].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 3 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG2D: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  /* Unscale mesh */
  if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /* Print quality histories */
  MMG2_outqua(mesh,sol);

  /* Print edge length histories */
  if ( abs(mesh->info.imprim) > 1 )  {
    MMG2_prilen(mesh,sol);
  }

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  if (!MMG2_pack(mesh,sol) ) _LIBMMG5_RETURN(mesh,sol,MMG5_LOWFAILURE);

  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMG2DMESH: ELAPSED TIME  %s\n",stim);
  _LIBMMG5_RETURN(mesh,sol,MMG5_SUCCESS);

}

int MMG2D_mmg2dls(MMG5_pMesh mesh,MMG5_pSol sol)
{
  mytime    ctim[TIMEMAX];
  char      stim[32];

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  /*uncomment for callback*/
  //MMG2D_callbackinsert = titi;

  /* interrupts */
  signal(SIGABRT,_MMG2_excfun);
  signal(SIGFPE,_MMG2_excfun);
  signal(SIGILL,_MMG2_excfun);
  signal(SIGSEGV,_MMG2_excfun);
  signal(SIGTERM,_MMG2_excfun);
  signal(SIGINT,_MMG2_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->info.lag >= 0 ) {
    fprintf(stdout,"  ## Error: lagrangian mode unavailable (MMG2D_IPARAM_lag):\n"
            "            You must call the MMG2D_mmg2dmov function to move a rigidbody.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }

  /* Load data */
  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG2DLS: INPUT DATA\n");
  chrono(ON,&(ctim[1]));

  sol->ver  = mesh->ver;

  if ( !mesh->nt ) {
    fprintf(stdout,"\n  ## ERROR: NO TRIANGLES IN THE MESH \n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }
  else if ( !sol->m ) {
    fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  } else   if ( sol->size != 1 ) {
    fprintf(stdout,"  ## ERROR: WRONG DATA TYPE.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  } else if ( sol->np && (sol->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,sol->m,(sol->size*(sol->npmax+1))*sizeof(double));
    sol->np = 0;
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);


  /* Set pointers */
  MMG2D_setfunc(mesh,sol);
  _MMG2D_Set_commonFunc();

  fprintf(stdout,"\n  %s\n   MODULE MMG2D-IMB/LJLL : %s (%s) %s\n  %s\n",
          MG_STR,MG_VER,MG_REL,sol->size == 1 ? "ISO" : "ANISO",MG_STR);
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM NUMBER OF POINTS    (NPMAX) : %8d\n",mesh->npmax);
    fprintf(stdout,"  MAXIMUM NUMBER OF TRIANGLES (NTMAX) : %8d\n",mesh->ntmax);
  }

  /* analysis */
  chrono(ON,&ctim[2]);
  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : DATA ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");
  if ( !MMG2_scaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  if ( mesh->nt && !MMG2_hashTria(mesh) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  if ( mesh->info.ddebug && !_MMG5_chkmsh(mesh,1,0) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /* Print initial quality */
  MMG2_outqua(mesh,sol);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  chrono(ON,&ctim[3]);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 2 : LEVEL-SET DISCRETIZATION\n");
  }

  /* Discretization of the mesh->info.ls isovalue of sol in the mesh */
  if (! MMG2_mmg2d6(mesh,sol) )
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /* Mesh analysis */
  if (! _MMG2_analys(mesh) )
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);

  /* Mesh improvement - call new version of mmg2d1 */
  chrono(ON,&ctim[4]);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 3 : MESH IMPROVEMENT\n");
  }

  if ( (!mesh->info.noinsert) && !MMG2_mmg2d1n(mesh,sol) ) {
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  /* End of mmg2dls */
  chrono(OFF,&(ctim[4]));
  printim(ctim[4].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 3 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG2D: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  /* Unscale mesh */
  if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /* Print quality histories */
  MMG2_outqua(mesh,sol);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  /* Pack mesh */
  if (!MMG2_pack(mesh,sol) ) _LIBMMG5_RETURN(mesh,sol,MMG5_LOWFAILURE);

  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMG2DLIB: ELAPSED TIME  %s\n",stim);
  _LIBMMG5_RETURN(mesh,sol,MMG5_SUCCESS);

}

int MMG2D_mmg2dmov(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol disp) {
  mytime    ctim[TIMEMAX];
  char      stim[32];

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  /*uncomment for callback*/
  //MMG2D_callbackinsert = titi;

  /* interrupts */
  signal(SIGABRT,_MMG2_excfun);
  signal(SIGFPE,_MMG2_excfun);
  signal(SIGILL,_MMG2_excfun);
  signal(SIGSEGV,_MMG2_excfun);
  signal(SIGTERM,_MMG2_excfun);
  signal(SIGINT,_MMG2_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check data compatibility */
  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG2DMOV: INPUT DATA\n");
  chrono(ON,&(ctim[1]));

  disp->ver  = mesh->ver;

  if ( !mesh->nt ) {
    fprintf(stdout,"\n  ## ERROR: NO TRIANGLES IN THE MESH \n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  else if ( !disp->m ) {
    fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  else if ( disp->size != 2 ) {
    fprintf(stdout,"  ## ERROR: LAGRANGIAN MOTION OPTION NEED A VECTOR DISPLACEMENT.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  else if ( disp->np && (disp->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,disp->m,(disp->size*(disp->npmax+1))*sizeof(double));
    disp->np = 0;
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* Set pointers */
  MMG2D_setfunc(mesh,met);
  _MMG2D_Set_commonFunc();

  fprintf(stdout,"\n  %s\n   MODULE MMG2D-IMB/LJLL : %s (%s) %s\n  %s\n",
          MG_STR,MG_VER,MG_REL,met->size == 1 ? "ISO" : "ANISO",MG_STR);
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM NUMBER OF POINTS    (NPMAX) : %8d\n",mesh->npmax);
    fprintf(stdout,"  MAXIMUM NUMBER OF TRIANGLES (NTMAX) : %8d\n",mesh->ntmax);
  }

  /* Analysis */
  chrono(ON,&ctim[2]);
  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : DATA ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");
  if ( !MMG2_scaleMesh(mesh,disp) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  if ( mesh->nt && !MMG2_hashTria(mesh) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  if ( mesh->info.ddebug && !_MMG5_chkmsh(mesh,1,0) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  /* Print initial quality */
  MMG2_outqua(mesh,met);

  /* Mesh analysis */
  if (! _MMG2_analys(mesh) )
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* Lagrangian motion */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 2 : LAGRANGIAN MOTION\n");
  }

#ifdef USE_ELAS

  /* Lagrangian mode */
  if ( !MMG2_mmg2d9(mesh,disp,met) ) {
    disp->npi = disp->np;
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
#endif

  /* End with a classical remeshing stage, provided mesh->info.lag > 1 */
  if ( (mesh->info.lag >= 1) && !MMG2_mmg2d1n(mesh,met) ) {
    _MMG2D_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG2D: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  /* Unscale mesh */
  if ( !MMG2_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  /* Print quality histories */
  MMG2_outqua(mesh,met);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  /* Pack mesh */
  if (!MMG2_pack(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);

  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMG2DMOV: ELAPSED TIME  %s\n",stim);

  _LIBMMG5_RETURN(mesh,met,MMG5_SUCCESS);
}
