/** mmg3d: 3d mesh adaptation
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

int       opt_i[9];
double    opt_d[5];
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

  exit(MG_LOWFAILURE);
}


static int parsar(int argc,char *argv[],pMesh mesh,pSol met) {
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
          opt_d[dhd] = atof(argv[i]);
        break;
      case 'd':  /* debug */
        opt_i[debug] = 1;
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc )
          opt_d[hmin] = atof(argv[i]);
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc )
          opt_d[hmax] = atof(argv[i]);
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc )
          opt_d[hausd] = atof(argv[i]);
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc )
          opt_d[hgrad] = atof(argv[i]);
        else
          usage(argv[0]);
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            mesh->namein = (char*) calloc(strlen(argv[i])+1,sizeof(char));
            strcpy(mesh->namein,argv[i]);
            opt_d[imprim] = 5;
          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            usage(argv[0]);
          }
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-ls") ) {
          opt_i[iso] = 1;
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            opt_d[ls] = atof(argv[i]);
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
          opt_i[mem] = atoi(argv[i]);
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          usage(argv[0]);
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") ) {
          opt_i[angle] =  0;
          opt_d[  dhd] = -1.0;
        }
        else if ( !strcmp(argv[i],"-noswap") )
          opt_i[noswap] = 1;
        else if( !strcmp(argv[i],"-noinsert") )
          opt_i[noinsert] = 1;
        else if( !strcmp(argv[i],"-nomove") )
          opt_i[nomove] = 1;
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
              opt_i[renum] = atoi(argv[i]);
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
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) )
            opt_i[imprim] = atoi(argv[i]);
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
        if ( opt_i[imprim] == -99 )  opt_i[imprim] = 5;
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
  if ( opt_i[imprim] == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    fscanf(stdin,"%d",&i);
    opt_i[imprim] = i;
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

static inline void endcod() {
  char    stim[32];

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3DLIB: ELAPSED TIME  %s\n",stim);
}

/** Deallocations before return */
static void freeAll(pMesh mesh,pSol met){
  free(mesh->point);
  mesh->point=NULL;
  free(mesh->tetra);
  mesh->tetra=NULL;
  free(mesh->adja);
  mesh->adja=NULL;
  free(mesh->nameout);
  mesh->nameout=NULL;
  free(mesh->namein);
  mesh->namein=NULL;
  if ( mesh->xpoint ) {
    free(mesh->xpoint);
    mesh->xpoint=NULL;
  }
  if ( mesh->htab.geom ) {
    free(mesh->htab.geom);
    mesh->htab.geom=NULL;
  }
  if ( mesh->tria ) {
    free(mesh->tria);
    mesh->tria=NULL;
  }
  if ( mesh->xtetra ) {
    free(mesh->xtetra);
    mesh->xtetra=NULL;
  }

  /* met */
  if ( met->namein ) {
    free(met->namein);
    met->namein=NULL;
  }
  if ( met->nameout ) {
    free(met->nameout);
    met->nameout=NULL;
  }
  if ( !opt_i[iso] && met->m ) {
    free(met->m);
    met->m = NULL;
  }
  //AJETER?? memset(&met,0,sizeof(Sol));
}


int main(int argc,char *argv[]) {
  Mesh      mesh;
  Sol       met;
  int       ier;
  char      stim[32];

  atexit(endcod);

  tminit(ctim,TIMEMAX);
  chrono(ON,&ctim[0]);

  /* assign default values */
  memset(&mesh,0,sizeof(Mesh));
  memset(&met,0,sizeof(Sol));

  mmg3dinit(opt_i,opt_d);

  /* command line */
  if ( !parsar(argc,argv,&mesh,&met) )  return(1);

  info.iso    = opt_i[iso];
  info.imprim = opt_i[imprim];
  info.mem    = opt_i[mem];

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&ctim[1]);
  /* read mesh file */
  if ( !loadMesh(&mesh) ) RETURN_AND_FREE(&mesh,&met,MG_STRONGFAILURE);

  /* read metric if any */
  met.size = 1;
  met.npmax = mesh.npmax;
  ier = loadMet(&met);
  if ( !ier )
    RETURN_AND_FREE(&mesh,&met,MG_STRONGFAILURE);

  chrono(OFF,&ctim[1]);
  printim(ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  ier = mmg3dlib(opt_i,opt_d,&mesh,&met);

  if ( ier != MG_STRONGFAILURE ) {
    chrono(ON,&ctim[1]);
    if ( opt_i[imprim] )  fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh.nameout);
    if ( !unscaleMesh(&mesh,&met) ) { RETURN_AND_FREE(&mesh,&met,EXIT_FAILURE);}
    if ( !saveMesh(&mesh) )         { RETURN_AND_FREE(&mesh,&met,EXIT_FAILURE);}
    if ( !saveMet(&mesh,&met) )     { RETURN_AND_FREE(&mesh,&met,EXIT_FAILURE);}
    chrono(OFF,&ctim[1]);
    if ( opt_i[imprim] )  fprintf(stdout,"  -- WRITING COMPLETED\n");
  }

  /* free mem */
  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3D: ELAPSED TIME  %s\n",stim);
  RETURN_AND_FREE(&mesh,&met,ier);
}
