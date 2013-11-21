/** \include Example for using of mmg3dlib
 * mmg3d: 3d mesh adaptation
 *
 * Written by Cecile Dobrzynski (IMB), Charles Dapogny and Pascal Frey (LJLL)
 * Copyright (c) 2004- IMB/LJLL.
 * All rights reserved.
 */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "libmmg3d5.h"
#include "mmg3d.h"

int       opt_i[10];
double    opt_d[6];
mytime    ctim[TIMEMAX];

static void usage(char *prog) {
  fprintf(stdout,"\nUsage: %s [-v [n]] [opts..] filein [fileout]\n",prog);

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-h      Print this message\n");
  fprintf(stdout,"-v [n]  Tune level of verbosity, [-10..10]\n");
  fprintf(stdout,"-m [n]  Set memory size to n Mbytes\n");
  fprintf(stdout,"-d      Turn on debug mode\n");

  fprintf(stdout,"\n**  File specifications\n");
  fprintf(stdout,"-in  file  input triangulation\n");
  fprintf(stdout,"-out file  output triangulation\n");
  fprintf(stdout,"-sol file  load solution file\n");
#ifdef SINGUL
  fprintf(stdout,"-sing file load file containing singularities\n");
#endif

  fprintf(stdout,"\n**  Parameters\n");
  fprintf(stdout,"-ar val    angle detection\n");
  fprintf(stdout,"-nr        no angle detection\n");
  fprintf(stdout,"-hmin val  minimal mesh size\n");
  fprintf(stdout,"-hmax val  maximal mesh size\n");
  fprintf(stdout,"-hausd val control Hausdorff distance\n");
  fprintf(stdout,"-hgrad val control gradation\n");
  fprintf(stdout,"-ls        levelset meshing \n");
  fprintf(stdout,"-noswap    no edge or face flipping\n");
  fprintf(stdout,"-nomove    no point relocation\n");
  fprintf(stdout,"-noinsert  no point insertion/deletion \n");

#ifdef USE_SCOTCH
  fprintf(stdout,"-rn [n]    Turn on or off the renumbering using SCOTCH (0/1) \n");
#endif

  exit(MMG5_LOWFAILURE);
}


static int parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,pSingul sing) {
  int     i;
  char   *ptr;

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        usage(argv[0]);
        break;

      case 'a':
        if ( !strcmp(argv[i],"-ar") && ++i < argc )
          opt_d[MMG5_dhd] = atof(argv[i]);
        break;
      case 'd':  /* debug */
        opt_i[MMG5_debug] = 1;
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc )
          opt_d[MMG5_hmin] = atof(argv[i]);
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc )
          opt_d[MMG5_hmax] = atof(argv[i]);
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc )
          opt_d[MMG5_hausd] = atof(argv[i]);
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc )
          opt_d[MMG5_hgrad] = atof(argv[i]);
        else
          usage(argv[0]);
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            mesh->namein = (char*) calloc(strlen(argv[i])+1,sizeof(char));
            strcpy(mesh->namein,argv[i]);
            opt_d[MMG5_imprim] = 5;
          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            usage(argv[0]);
          }
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-ls") ) {
          opt_i[MMG5_iso] = 1;
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            opt_d[MMG5_ls] = atof(argv[i]);
          }
          else if ( i == argc ) {
            fprintf(stderr,"Missing argument option %c%c\n",argv[i-1][1],argv[i-1][2]);
            usage(argv[0]);
          }
          else i--;
        }
        break;
      case 'm':  /* memory */
        if ( ++i < argc && isdigit(argv[i][0]) )
          opt_i[MMG5_mem] = atoi(argv[i]);
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          usage(argv[0]);
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") ) {
          opt_i[MMG5_angle] =  0;
          opt_d[  MMG5_dhd] = -1.0;
        }
        else if ( !strcmp(argv[i],"-noswap") )
          opt_i[MMG5_noswap] = 1;
        else if( !strcmp(argv[i],"-noinsert") )
          opt_i[MMG5_noinsert] = 1;
        else if( !strcmp(argv[i],"-nomove") )
          opt_i[MMG5_nomove] = 1;
        break;
      case 'o':
        if ( !strcmp(argv[i],"-out") ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            mesh->nameout = (char*) calloc(strlen(argv[i])+1,sizeof(char));
            strcpy(mesh->nameout,argv[i]);
          }else{
            fprintf(stderr,"Missing filname for %c%c%c\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            usage(argv[0]);
          }
        }
        break;
#ifdef USE_SCOTCH
      case 'r':
        if ( !strcmp(argv[i],"-rn") ) {
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) )
              opt_i[MMG5_renum] = atoi(argv[i]);
            else {
              fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
              usage(argv[0]);
            }
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            usage(argv[0]);
          }
        }
        break;
#endif
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            met->namein = (char*) calloc(strlen(argv[i])+1,sizeof(char));
            strcpy(met->namein,argv[i]);
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            usage(argv[0]);
          }
        }
#ifdef SINGUL
        else if ( !strcmp(argv[i],"-sf") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            sing->namein = (char*) calloc(strlen(argv[i])+1,sizeof(char));
            strcpy(sing->namein,argv[i]);
            opt_i[MMG5_sing] = 1;
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            usage(argv[0]);
          }
        }
        else if ( !strcmp(argv[i],"-sing") )
          opt_i[MMG5_sing] = 1;
#endif
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) )
            opt_i[MMG5_imprim] = atoi(argv[i]);
          else
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          usage(argv[0]);
        }
        break;
      default:
        fprintf(stderr,"Unrecognized option %s\n",argv[i]);
        usage(argv[0]);
      }
    }
    else {
      if ( mesh->namein == NULL ) {
        mesh->namein = (char*) calloc(strlen(argv[i])+1,sizeof(char));
        strcpy(mesh->namein,argv[i]);
        if ( opt_i[MMG5_imprim] == -99 )  opt_i[MMG5_imprim] = 5;
      }
      else if ( mesh->nameout == NULL ){
        mesh->nameout = (char*) calloc(strlen(argv[i])+1,sizeof(char));
        strcpy(mesh->nameout,argv[i]);
      }
      else {
        fprintf(stdout,"Argument %s ignored\n",argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* check file names */
  if ( opt_i[MMG5_imprim] == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    fscanf(stdin,"%d",&i);
    opt_i[MMG5_imprim] = i;
  }

  if ( mesh->namein == NULL ) {
    mesh->namein = (char *)calloc(128,sizeof(char));
    assert(mesh->namein);
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",mesh->namein);
  }
  if ( mesh->nameout == NULL ) {
    mesh->nameout = (char *)calloc(128,sizeof(char));
    assert(mesh->nameout);
    strcpy(mesh->nameout,mesh->namein);
    ptr = strstr(mesh->nameout,".mesh");
    if ( ptr ) *ptr = '\0';
    strcat(mesh->nameout,".o.mesh");
    ptr = strstr(mesh->nameout,".meshb");
    if ( ptr )  strcat(mesh->nameout,"b");
  }

  if ( met->namein == NULL ) {
    met->namein = (char *)calloc(128,sizeof(char));
    assert(met->namein);
    strcpy(met->namein,mesh->namein);
    ptr = strstr(met->namein,".mesh");
    if ( ptr ) *ptr = '\0';
    strcat(met->namein,".sol");
  }
  if ( met->nameout == NULL ) {
    met->nameout = (char *)calloc(128,sizeof(char));
    assert(met->nameout);
    strcpy(met->nameout,mesh->nameout);
    ptr = strstr(met->nameout,".mesh");
    if ( ptr ) *ptr = '\0';
    strcat(met->nameout,".sol");
  }
  return(1);
}

/** Deallocations of names */
static inline
void freeName(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSingul singul) {

  free(mesh->nameout);
  mesh->nameout = NULL;
  free(mesh->namein);
  mesh->namein = NULL;
  if ( met->namein ) {
    free(met->namein);
    met->namein = NULL;
  }
  if ( met->nameout ) {
    free(met->nameout);
    met->nameout = NULL;
  }
#ifdef SINGUL
  if ( info.sing && singul->namein ) {
    free(singul->namein);
    singul->namein = NULL;
  }
#endif
}

static inline void endcod() {
  char    stim[32];

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3DLIB: ELAPSED TIME  %s\n",stim);
}

int main(int argc,char *argv[]) {
  MMG5_Mesh      mesh;
  MMG5_Sol       met;
  MMG5_Singul    sing;
  int            ier;
  char           stim[32];

  atexit(endcod);

  tminit(ctim,TIMEMAX);
  chrono(ON,&ctim[0]);

  /* assign default values */
  memset(&mesh,0,sizeof(MMG5_Mesh));
  memset(&met,0,sizeof(MMG5_Sol));
#ifdef SINGUL
  memset(&sing,0,sizeof(MMG5_Singul));
#endif

  MMG5_mmg3dinit(opt_i,opt_d);

  /* command line */
  if ( !parsar(argc,argv,&mesh,&met,&sing) )  return(1);

  /* infos needed to call MMG5_loadMesh: info.iso, info.imprim, info.mem */
  info.iso    = opt_i[MMG5_iso];
  info.imprim = opt_i[MMG5_imprim];
  info.mem    = opt_i[MMG5_mem];

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&ctim[1]);
  /* read mesh file */
  if ( !MMG5_loadMesh(&mesh) ) {
    MMG5_freeAll(&mesh,&met
#ifdef SINGUL
                 ,&sing
#endif
                 );
    freeName(&mesh,&met,&sing);
    return(MMG5_STRONGFAILURE);
  }

  /* read metric if any */
  /* infos needed to call MMG5_loadMet: met.size, met.npmax */
  met.size = 1;
  met.npmax = mesh.npmax;
  ier = MMG5_loadMet(&met);
  if ( !ier ) {
    MMG5_freeAll(&mesh,&met
#ifdef SINGUL
                 ,&sing
#endif
                 );
    freeName(&mesh,&met,&sing);
    return(MMG5_STRONGFAILURE);
  }
#ifdef SINGUL
  if ( opt_i[MMG5_sing] && sing.namein ) {
    ier = MMG5_loadSingul(&sing);
    if ( !ier ) {
      MMG5_freeAll(&mesh,&met,&sing);
      freeName(&mesh,&met,&sing);
      return(MMG5_STRONGFAILURE);
    }
  }
#endif
  chrono(OFF,&ctim[1]);
  printim(ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  ier = MMG5_mmg3dlib(opt_i,opt_d,&mesh,&met
#ifdef SINGUL
                      ,&sing
#endif
                      );

  if ( ier != MMG5_STRONGFAILURE ) {
    chrono(ON,&ctim[1]);
    if ( opt_i[MMG5_imprim] )  fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh.nameout);
    if ( !MMG5_saveMesh(&mesh) )         {
      MMG5_freeAll(&mesh,&met
#ifdef SINGUL
                   ,&sing
#endif
                   );
      freeName(&mesh,&met,&sing);
      return(EXIT_FAILURE);
    }
    if ( !MMG5_saveMet(&mesh,&met) )     {
      MMG5_freeAll(&mesh,&met
#ifdef SINGUL
                   ,&sing
#endif
                   );
      freeName(&mesh,&met,&sing);
      return(EXIT_FAILURE);
    }
    chrono(OFF,&ctim[1]);
    if ( opt_i[MMG5_imprim] )  fprintf(stdout,"  -- WRITING COMPLETED\n");
  }

  /* free mem */
  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3D: ELAPSED TIME  %s\n",stim);
  MMG5_freeAll(&mesh,&met
#ifdef SINGUL
               ,&sing
#endif
               );
  freeName(&mesh,&met,&sing);
  return(ier);
}
