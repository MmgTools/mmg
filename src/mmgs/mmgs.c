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
 * \file mmgs/mmgs.c
 * \brief Main file for MMGS executable: perform surface mesh adaptation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgs.h"
#include "shared_func.h"

#include <math.h>

mytime         MMG5_ctim[TIMEMAX];

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (metric or solution).
 * \param ... optional arguments: not used for now. To end by the NULL value.
 *
 * Deallocations before return.
 *
 * \Remark To call with NULL as last argument.
 *
 */
void MMGS_Free_all(MMG5_pMesh mesh,MMG5_pSol met,... ){

  MMGS_Free_structures(mesh,met);
}

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
 * \param mesh pointer toward the mesh structure (unused).
 *
 * Set pointer for MMGS_saveMesh function.
 *
 */
void MMGS_Set_saveFunc(MMG5_pMesh mesh) {
  _MMGS_saveMeshinternal = _MMGS_saveAllMesh;
}

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 1.
 *
 * Store command line arguments.
 *
 */
static
int _MMG5_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met) {
  int    i;
  char   namein[128];

  /* First step: search if user want to see the default parameters values. */
  for ( i=1; i< argc; ++i ) {
    if ( !strcmp(argv[i],"-val") ) {
      _MMG5_defaultValues(mesh);
    }
  }

  /* Second step: read all other arguments. */
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        _MMG5_usage(argv[0]);
        break;
      case 'a': /* ridge angle */
        if ( !strcmp(argv[i],"-ar") && ++i < argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_angleDetection,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        break;
      case 'A': /* anisotropy */
        if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor) )
          exit(EXIT_FAILURE);
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hmin,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hmax,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hausd,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hgrad,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else
          _MMG5_usage(argv[0]);
        break;
      case 'd':
        if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_debug,1) )
          exit(EXIT_FAILURE);
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMGS_Set_inputMeshName(mesh, argv[i]) )
              exit(EXIT_FAILURE);

            if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,5) )
              exit(EXIT_FAILURE);
          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            _MMG5_usage(argv[0]);
          }
        }
        break;
      case 'm':
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_mem,atoi(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          _MMG5_usage(argv[0]);
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_angle,0) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_noswap,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_noinsert,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_nomove,1) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-nreg") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_nreg,1) )
            exit(EXIT_FAILURE);
        }
        break;
      case 'o':
        if ( !strcmp(argv[i],"-out") ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            if ( !MMGS_Set_outputMeshName(mesh,argv[i]) )
              exit(EXIT_FAILURE);
          }else{
            fprintf(stderr,"Missing filname for %c%c%c\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            _MMG5_usage(argv[0]);
          }
        }
        break;
#ifdef USE_SCOTCH
      case 'r':
        if ( !strcmp(argv[i],"-rn") ) {
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {
              if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_renum,atoi(argv[i])) )
                exit(EXIT_FAILURE);
            }
            else {
              fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
              _MMG5_usage(argv[0]);
            }
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            _MMG5_usage(argv[0]);
          }
        }
        break;
#endif
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMGS_Set_inputSolName(mesh,met,argv[i]) )
              exit(EXIT_FAILURE);
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            _MMG5_usage(argv[0]);
          }
        }
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) ) {
            if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          _MMG5_usage(argv[0]);
        }
        break;
      default:
        fprintf(stderr,"Unrecognized option %s\n",argv[i]);
        _MMG5_usage(argv[0]);
      }
    }
    else {
      if ( mesh->namein == NULL ) {
        if ( !MMGS_Set_inputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
        if ( mesh->info.imprim == -99 ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,5) )
            exit(EXIT_FAILURE);
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMGS_Set_outputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
      }
      else {
        fprintf(stdout,"Argument %s ignored\n",argv[i]);
        _MMG5_usage(argv[0]);
      }
    }
    i++;
  }

  /* check file names */
  if ( mesh->info.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    fscanf(stdin,"%d",&i);
    if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,i) )
      exit(EXIT_FAILURE);
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",namein);
    if ( !MMGS_Set_inputMeshName(mesh,namein) )
      exit(EXIT_FAILURE);
  }

  if ( mesh->nameout == NULL ) {
    if ( !MMGS_Set_outputMeshName(mesh,"") )
      exit(EXIT_FAILURE);
  }

  if ( met->namein == NULL ) {
    if ( !MMGS_Set_inputSolName(mesh,met,"") )
      exit(EXIT_FAILURE);
  }
  if ( met->nameout == NULL ) {
    if ( !MMGS_Set_outputSolName(mesh,met,"") )
      exit(EXIT_FAILURE);
  }
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 1.
 *
 * Read local parameters file. This file must have the same name as
 * the mesh with the \a .mmg3d5 extension or must be named \a
 * DEFAULT.mmg3d5.
 *
 */
static int _MMG5_parsop(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_Par   *par;
  float      fp1,fp2;
  int        ref,i,j,ret,npar;
  char       *ptr,buf[256],data[256];
  FILE       *in;

  /* check for parameter file */
  strcpy(data,mesh->namein);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".mmgs5");
  in = fopen(data,"r");
  if ( !in ) {
    sprintf(data,"%s","DEFAULT.mmgs5");
    in = fopen(data,"r");
    if ( !in )  return(1);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* read parameters */
  mesh->info.npar = 0;
  while ( !feof(in) ) {
    /* scan line */
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* check for condition type */
    if ( !strcmp(data,"parameters") ) {
      fscanf(in,"%d",&npar);
      if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_numberOfLocalParam,npar) )
        exit(EXIT_FAILURE);

      for (i=0; i<mesh->info.npar; i++) {
        par = &mesh->info.par[i];
        fscanf(in,"%d %s ",&ref,buf);
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          par->elt = MS_Ver;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  par->elt = MS_Tri;
        else {
          fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
          continue;
        }
        ret = fscanf(in,"%f %f",&fp1,&fp2);
        if ( !MMGS_Set_localParameter(mesh,met,MMG5_Triangle,ref,
                                      fp1,fp2,mesh->info.hausd) )
          exit(EXIT_FAILURE);
      }
    }
  }
  fclose(in);
  return(1);
}

int main(int argc,char *argv[]) {
  MMG5_Mesh mesh;
  MMG5_Sol  met;
  int       ier;
  char      stim[32];

  fprintf(stdout,"  -- MMGS, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  _MMG5_Set_commonFunc();

  /* trap exceptions */
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

  MMGS_Init_parameters(&mesh);

  met.size    = 1;

  /* command line */
  if ( !_MMG5_parsar(argc,argv,&mesh,&met) )  return(MMG5_STRONGFAILURE);

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&MMG5_ctim[1]);

  if ( !MMGS_loadMesh(&mesh) )
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);

  ier = MMGS_loadMet(&mesh,&met);
  if ( ier==-1 ) {
      fprintf(stdout,"  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
  }

  if ( !_MMG5_parsop(&mesh,&met) )
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_LOWFAILURE);

  if ( !_MMG5_scaleMesh(&mesh,&met) )
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);

  chrono(OFF,&MMG5_ctim[1]);
  printim(MMG5_ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&MMG5_ctim[2]);
  _MMG5_setfunc(&mesh,&met);
  MMGS_Set_saveFunc(&mesh);

  fprintf(stdout,"\n  %s\n   MODULE MMGS-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
  if ( mesh.info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  if ( !_MMG5_analys(&mesh) )
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_LOWFAILURE);

  if ( abs(mesh.info.imprim) > 0 ) _MMG5_inqua(&mesh,&met);

  if ( mesh.info.imprim > 1 && met.m ) _MMG5_prilen(&mesh,&met,0);

  chrono(OFF,&MMG5_ctim[2]);
  if ( mesh.info.imprim ) {
    printim(MMG5_ctim[2].gdif,stim);
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  /* solve */
  chrono(ON,&MMG5_ctim[3]);
  if ( mesh.info.imprim )
    fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met.size < 6 ? "ISOTROPIC" : "ANISOTROPIC");
  if ( !_MMG5_mmgs1(&mesh,&met) )  {
    if ( (!mesh.adja) && !_MMG5_hashTria(&mesh) ) {
      fprintf(stdout,"  ## Hashing problem. Unable to save mesh.\n");
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(&mesh,&met) )
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
    if ( !MMGS_saveMesh(&mesh) )
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
    if ( met.m && !MMGS_saveMet(&mesh,&met) )
      _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_LOWFAILURE);
  }
  chrono(OFF,&MMG5_ctim[3]);
  if ( mesh.info.imprim ) {
    printim(MMG5_ctim[3].gdif,stim);
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }
  fprintf(stdout,"\n  %s\n   END OF MODULE MMGS-LJLL \n  %s\n",MG_STR,MG_STR);

  /* save file */
  _MMG5_outqua(&mesh,&met);
  if ( mesh.info.imprim > 1 )  _MMG5_prilen(&mesh,&met,1);

  chrono(ON,&MMG5_ctim[1]);
  if ( mesh.info.imprim )  fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh.nameout);
  if ( !_MMG5_unscaleMesh(&mesh,&met) )
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
  if ( !MMGS_saveMesh(&mesh) )
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
  if ( !MMGS_saveMet(&mesh,&met) )
    _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_STRONGFAILURE);
  chrono(OFF,&MMG5_ctim[1]);
  if ( mesh.info.imprim )  fprintf(stdout,"  -- WRITING COMPLETED\n");

  /* release memory */
  /* free mem */
  _MMG5_RETURN_AND_FREE(&mesh,&met,MMG5_SUCCESS);

  return(0);
}
