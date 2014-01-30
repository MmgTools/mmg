/** \include Example for using of mmg3dlib
 * mmg3d: 3d mesh adaptation
 *
 * Written by Cecile Dobrzynski (IMB), Charles Dapogny and Pascal Frey (LJLL)
 * Copyright (c) 2004- IMB/LJLL.
 * All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <ctype.h>

#include "libmmg3d5.h"

mytime    ctim[TIMEMAX];

#ifdef SINGUL
#define RETURN_AND_FREE(mesh,met,sing,val)do    \
    {                                                 \
      MMG5_Free_all(mesh,met,sing);                   \
      return(val);                                    \
    }while(0)
#else
#define RETURN_AND_FREE(mesh,met,sing,val)do    \
    {                                                 \
      MMG5_Free_all(mesh,met);                        \
      return(val);                                    \
    }while(0)
#endif

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
  fprintf(stdout,"-sf  file load file containing singularities\n");
#endif
#ifdef USE_SCOTCH
  fprintf(stdout,"-rn [n]    Turn on or off the renumbering using SCOTCH [1/0] \n");
#endif
#ifdef SINGUL
  fprintf(stdout,"-sing      Preserve internal singularities\n");
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

  exit(MMG5_LOWFAILURE);
}


static int parsar(int argc,char *argv[],MMG5_pMesh mesh,
                  MMG5_pSol met,MMG5_pSingul sing) {
  int     i;

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        usage(argv[0]);
        break;

      case 'a':
        if ( !strcmp(argv[i],"-ar") && ++i < argc )
          if ( !MMG5_Set_dparameters(mesh,met,MMG5_DPARAM_angleDetection,
                                     atof(argv[i])) )
            exit(EXIT_FAILURE);
        break;
      case 'd':  /* debug */
        if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_debug,1) )
          exit(EXIT_FAILURE);
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !MMG5_Set_dparameters(mesh,met,MMG5_DPARAM_hmin,
                                     atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !MMG5_Set_dparameters(mesh,met,MMG5_DPARAM_hmax,
                                     atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc ) {
          if ( !MMG5_Set_dparameters(mesh,met,MMG5_DPARAM_hausd,
                                     atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc ) {
          if ( !MMG5_Set_dparameters(mesh,met,MMG5_DPARAM_hgrad,
                                     atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else
          usage(argv[0]);
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMG5_Set_inputMeshName(mesh, argv[i]) )
              exit(EXIT_FAILURE);

            if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_verbose,5) )
              exit(EXIT_FAILURE);
          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            usage(argv[0]);
          }
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-ls") ) {
          if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_iso,1) )
            exit(EXIT_FAILURE);
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG5_Set_dparameters(mesh,met,MMG5_DPARAM_ls,atof(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else if ( i == argc ) {
            fprintf(stderr,"Missing argument option %c%c\n",argv[i-1][1],argv[i-1][2]);
            usage(argv[0]);
          }
          else i--;
        }
        break;
      case 'm':  /* memory */
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_mem,atoi(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          usage(argv[0]);
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_angle,0) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_noswap,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_noinsert,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_nomove,1) )
            exit(EXIT_FAILURE);
        }
        break;
      case 'o':
        if ( !strcmp(argv[i],"-out") ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            if ( !MMG5_Set_outputMeshName(mesh,argv[i]) )
              exit(EXIT_FAILURE);
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
            if ( isdigit(argv[i][0]) ) {
              if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_renum,1) )
                exit(EXIT_FAILURE);
            }
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
            if ( !MMG5_Set_inputSolName(mesh,met,argv[i]) )
              exit(EXIT_FAILURE);
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            usage(argv[0]);
          }
        }
#ifdef SINGUL
        else if ( !strcmp(argv[i],"-sf") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG5_Set_inputSingulName(mesh,sing,argv[i]) )
              exit(EXIT_FAILURE);
            if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_sing,1) )
              exit(EXIT_FAILURE);
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            usage(argv[0]);
          }
        }
        else if ( !strcmp(argv[i],"-sing") )
          if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_sing,1) )
            exit(EXIT_FAILURE);
#endif
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) ) {
            if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_verbose,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
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
        if ( !MMG5_Set_inputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
        if ( mesh->info.imprim == -99 ) {
          if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_verbose,5) )
              exit(EXIT_FAILURE);
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMG5_Set_outputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
      }
      else {
        printf("ah?? %s %s\n",mesh->namein,mesh->nameout );
        fprintf(stdout,"Argument %s ignored\n",argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* check file names */
  if ( mesh->info.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    fscanf(stdin,"%d",&i);
    if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_verbose,i) )
      exit(EXIT_FAILURE);
  }

  if ( mesh->namein == NULL ) {
    mesh->namein = (char *)calloc(128,sizeof(char));
    if ( !mesh->namein ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }

    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",mesh->namein);
  }
  if ( mesh->nameout == NULL ) {
    if ( !MMG5_Set_outputMeshName(mesh,"") )
      exit(EXIT_FAILURE);
  }

  if ( met->namein == NULL ) {
    if ( !MMG5_Set_inputSolName(mesh,met,"") )
      exit(EXIT_FAILURE);
  }
  if ( met->nameout == NULL ) {
    if ( !MMG5_Set_outputSolName(mesh,met,"") )
      exit(EXIT_FAILURE);
  }
  return(1);
}

/** Read parammeter file */
static inline
int parsop(MMG5_pMesh mesh,MMG5_pSol met) {
  float       fp1;
  int         ref,i,j,ret,npar;
  char       *ptr,buf[256],data[256];
  FILE       *in;

  /* check for parameter file */
  strcpy(data,mesh->namein);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".mmg3d5");
  in = fopen(data,"r");
  if ( !in ) {
    sprintf(data,"%s","DEFAULT.mmg3d5");
    in = fopen(data,"r");
    if ( !in )  return(1);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* read parameters */
  while ( !feof(in) ) {
    /* scan line */
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* check for condition type */
    if ( !strcmp(data,"parameters") ) {
      fscanf(in,"%d",&npar);
      if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_numberOfLocalParam,npar) )
        exit(EXIT_FAILURE);

      for (i=0; i<mesh->info.npar; i++) {
        ret = fscanf(in,"%d %s ",&ref,buf);
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        if ( strcmp(buf,"triangles") && strcmp(buf,"triangle") ) {
          fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
          continue;
        }
        ret = fscanf(in,"%f",&fp1);
        if ( !MMG5_Set_localParameters(mesh,met,MMG5_Triangle,ref,fp1) )
          exit(EXIT_FAILURE);
      }
    }
  }
  fclose(in);
  return(1);
}

static inline void endcod() {
  char    stim[32];

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3DLIB: ELAPSED TIME  %s\n",stim);
}

int main(int argc,char *argv[]) {
  MMG5_pMesh      mesh;
  MMG5_pSol       met;
  MMG5_pSingul    sing;
  int             ier;
  char            stim[32];

  atexit(endcod);

  tminit(ctim,TIMEMAX);
  chrono(ON,&ctim[0]);

  /* assign default values */
  mesh = NULL;
  met  = NULL;
  sing = NULL;
#ifndef SINGUL
#ifndef SINGUL
  MMG5_Init_mesh(&mesh,&met);
  //  memset(&sing,0,sizeof(MMG5_Singul));
#else
  MMG5_Init_mesh(&mesh,&met,&sing);
#endif
  /* reset default values for file names */
#ifndef SINGUL
  MMG5_Free_names(mesh,met);
#else
  MMG5_Free_names(mesh,met,sing);
#endif

  /* command line */
  if ( !parsar(argc,argv,mesh,met,sing) )  return(1);

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&ctim[1]);
  /* read mesh file */
  if ( !MMG5_loadMesh(mesh) ) {
    MMG5_Free_all(mesh,met
#ifdef SINGUL
                 ,sing
#endif
                 );
    return(MMG5_STRONGFAILURE);
  }
  if ( !MMG5_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Scalar) ) {
    MMG5_Free_all(mesh,met
#ifdef SINGUL
                  ,sing
#endif
                  );
    return(MMG5_STRONGFAILURE);
  }

  /* read metric if any */
  ier = MMG5_loadMet(mesh,met);
  if ( !ier ) {
    MMG5_Free_all(mesh,met
#ifdef SINGUL
                 ,sing
#endif
                 );
    return(MMG5_STRONGFAILURE);
  }
#ifdef SINGUL
  if ( mesh->info.sing && sing->namein ) {
    ier = MMG5_loadSingul(sing);
    if ( !ier ) {
      MMG5_Free_all(mesh,met,sing);
      return(MMG5_STRONGFAILURE);
    }
  }
#endif
  if ( !parsop(mesh,met) )
    RETURN_AND_FREE(mesh,met,sing,MMG5_LOWFAILURE);

  chrono(OFF,&ctim[1]);
  printim(ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  ier = MMG5_mmg3dlib(mesh,met
#ifdef SINGUL
                      ,sing
#endif
                      );

  if ( ier != MMG5_STRONGFAILURE ) {
    chrono(ON,&ctim[1]);
    if ( mesh->info.imprim )
      fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh->nameout);
    if ( !MMG5_saveMesh(mesh) )         {
      MMG5_Free_all(mesh,met
#ifdef SINGUL
                   ,sing
#endif
                   );
      return(EXIT_FAILURE);
    }
    if ( !MMG5_saveMet(mesh,met) )     {
      MMG5_Free_all(mesh,met
#ifdef SINGUL
                   ,sing
#endif
                   );
      return(EXIT_FAILURE);
    }
    chrono(OFF,&ctim[1]);
    if ( mesh->info.imprim )
      fprintf(stdout,"  -- WRITING COMPLETED\n");
  }

  /* free mem */
  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3D: ELAPSED TIME  %s\n",stim);
  MMG5_Free_all(mesh,met
#ifdef SINGUL
               ,sing
#endif
               );
  return(ier);
}
