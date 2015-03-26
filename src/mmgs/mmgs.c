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
#include <math.h>

/* globals */
mytime         MMG5_ctim[TIMEMAX];

unsigned char inxt[3] = {1,2,0};
unsigned char iprv[3] = {2,0,1};


static void excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
  case SIGABRT:
    fprintf(stdout,"  Abnormal stop\n");  break;
  case SIGFPE:
    fprintf(stdout,"  Floating-point exception\n"); break;
  case SIGILL:
    fprintf(stdout,"  Illegal instruction\n"); break;
  case SIGSEGV:
    fprintf(stdout,"  Segmentation fault\n");  break;
  case SIGTERM:
  case SIGINT:
    fprintf(stdout,"  Program killed\n");  break;
  }
  exit(1);
}

static void usage(char *name) {
  fprintf(stdout,"\nUsage: %s [-v [n]] [-h] filein [fileout]\n",name);

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-h      Print this message\n");
  fprintf(stdout,"-v [n]  Turn on numerical information, [-10..10]\n");
  fprintf(stdout,"-d      Debbuging mesh->info.\n");

  fprintf(stdout,"\n**  File specifications\n");
  fprintf(stdout,"-in  file  input triangulation\n");
  fprintf(stdout,"-out file  output triangulation\n");
  fprintf(stdout,"-sol file  load solution or metric file\n");

  fprintf(stdout,"\n**  Parameters\n");
  fprintf(stdout,"-ar val    angle detection\n");
  fprintf(stdout,"-no        no mesh optimisation\n");
  fprintf(stdout,"-nr        no angle detection\n");
  fprintf(stdout,"-nreg      nornal regul.\n");
  fprintf(stdout,"-hmin val  minimal mesh size\n");
  fprintf(stdout,"-hmax val  maximal mesh size\n");
  fprintf(stdout,"-hausd val control Hausdorff distance\n");
  fprintf(stdout,"-hgrad val control gradation\n");
  fprintf(stdout,"-A         enable anisotropy\n");

  exit(1);
}


static int parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met) {
  int    i;
  char   namein[128];

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        usage(argv[0]);
        break;
      case 'a': /* ridge angle */
        if ( !strcmp(argv[i],"-ar") && ++i < argc ) {
          mesh->info.dhd = atof(argv[i]);
          mesh->info.dhd = MG_MAX(0.0, MG_MIN(180.0,mesh->info.dhd));
          mesh->info.dhd = cos(mesh->info.dhd*M_PI/180.0);
        }
        break;
      case 'A': /* anisotropy */
        met->size = 6;
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc )
          mesh->info.hmin = atof(argv[i]);
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc )
          mesh->info.hmax = atof(argv[i]);
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc ) {
          mesh->info.hausd = atof(argv[i]);
        }
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc ) {
          mesh->info.hgrad = atof(argv[i]);
          if ( mesh->info.hgrad < 0.0 )
            mesh->info.hgrad = -1.0;
          else
            mesh->info.hgrad = log(mesh->info.hgrad);
        }
        else
          usage(argv[0]);
        break;
      case 'd':
        mesh->info.ddebug = 1;
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMG5_Set_inputMeshName(mesh, argv[i]) )
              exit(EXIT_FAILURE);
            mesh->info.imprim = 5;
          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            usage(argv[0]);
          }
        }
        break;
      case 'm':
        if ( !strcmp(argv[i],"-m") ) {
          ++i;
          mesh->info.mem = atoi(argv[i]);
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") )
          mesh->info.dhd = -1.0;
        else if ( !strcmp(argv[i],"-nreg") )
          mesh->info.nreg = 1;
        else if ( !strcmp(argv[i],"-no") )
          mesh->info.opt = 0;
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
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG5_Set_inputSolName(mesh,met,argv[i]) )
              exit(EXIT_FAILURE);
            mesh->info.imprim = 5;
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
            mesh->info.imprim = atoi(argv[i]);
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
        if ( mesh->info.imprim == -99 )  mesh->info.imprim = 5;
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMG5_Set_outputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
      }
      else if ( met->namein == NULL ) {
        if ( !MMG5_Set_inputSolName(mesh,met,argv[i]) )
          exit(EXIT_FAILURE);
      }
      else if ( met->nameout == NULL ) {
        if ( !MMG5_Set_outputSolName(mesh,met,argv[i]) )
          exit(EXIT_FAILURE);
      }
      else {
        fprintf(stderr,"Argument %s ignored\n",argv[i]);
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
    mesh->info.imprim = i;
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",namein);
    if ( !MMG5_Set_inputMeshName(mesh,namein) )
      exit(EXIT_FAILURE);
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

static int parsop(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_Par   *par;
  float       fp1,fp2;
  int         i,j,ret;
  char       *ptr,buf[256],data[256];
  FILE       *in;

  /* check for parameter file */
  strcpy(data,mesh->namein);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".mmgs");
  in = fopen(data,"r");
  if ( !in ) {
    sprintf(data,"%s","DEFAULT.mmgs");
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
      fscanf(in,"%d",&mesh->info.npar);
      _MMG5_ADD_MEM(mesh,mesh->info.npar*sizeof(MMG5_Par),
                    "parameters",
                    printf("  Exit program.\n");
                    exit(EXIT_FAILURE));
      _MMG5_SAFE_CALLOC(mesh->info.par,mesh->info.npar,MMG5_Par);

      for (i=0; i<mesh->info.npar; i++) {
        par = &mesh->info.par[i];
        fscanf(in,"%d %s ",&par->ref,buf);
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          par->elt = MS_Ver;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  par->elt = MS_Tri;
        else {
          fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
          continue;
        }
        ret = fscanf(in,"%f %f",&fp1,&fp2);
        par->hmin  = fp1;
        par->hmax  = fp2;
        par->hausd = mesh->info.hausd;
      }
    }
  }
  fclose(in);
  return(1);
}

static void endcod() {
  char   stim[32];

  chrono(OFF,&MMG5_ctim[0]);
  printim(MMG5_ctim[0].gdif,stim);
  fprintf(stdout,"\n   ELAPSED TIME  %s\n",stim);
}

/* set function pointers w/r iso/aniso */
static void setfunc(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( met->size < 6 ) {
    calelt  = calelt_iso;
    defsiz  = defsiz_iso;
    gradsiz = gradsiz_iso;
    _MMG5_lenedg  = _MMG5_lenedg_iso;
    intmet  = intmet_iso;
    movintpt= movintpt_iso;
    movridpt= movridpt_iso;
  }
  else {
    calelt  = calelt_ani;
    defsiz  = defsiz_ani;
    gradsiz = gradsiz_ani;
    _MMG5_lenedg  = _MMG5_lenedg_ani;
    intmet  = intmet_ani;
    movintpt= movintpt_ani;
    movridpt= movridpt_ani;
  }
  MMG5_Set_outputMeshName = MMG5_MMGS_Set_outputMeshName;
  MMG5_Init_parameters = MMG5_MMGS_Init_parameters;
}

int main(int argc,char *argv[]) {
  MMG5_Mesh mesh;
  MMG5_Sol  met;
  int       ier;
  char      stim[32];

  fprintf(stdout,"  -- MMGS, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  setfunc(&mesh,&met);
  
  /* trap exceptions */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  atexit(endcod);

  tminit(MMG5_ctim,TIMEMAX);
  chrono(ON,&MMG5_ctim[0]);

  /* assign default values */
  memset(&mesh,0,sizeof(MMG5_Mesh));
  memset(&met,0,sizeof(MMG5_Sol));

  memset(&mesh,0,sizeof(MMG5_Mesh));
  memset(&met,0,sizeof(MMG5_Sol));
  mesh.info.imprim = -99;
  mesh.info.ddebug = 0;
  mesh.info.mem    = -1;
  mesh.info.dhd    = _MMG5_ANGEDG;
  mesh.info.hmin   = 0.0;
  mesh.info.hmax   = FLT_MAX;
  mesh.info.hausd  = 0.01;
  mesh.info.hgrad  = 0.1;
  mesh.info.badkal = 0;
  mesh.info.nreg   = 0;
  mesh.info.opt    = 1;
  mesh.info.mani   = 1;
  met.size    = 1;

  /* initial value for memMax and gap */
  mesh.gap = 0.2;
  mesh.memMax = _MMG5_memSize();
  if ( mesh.memMax )
    /* maximal memory = 50% of total physical memory */
    mesh.memMax = mesh.memMax*50/100;
  else {
    /* default value = 800 Mo */
    printf("  Maximum memory set to default value: %d Mo.\n",_MMG5_MEMMAX);
    mesh.memMax = _MMG5_MEMMAX << 20;
  }

  /* command line */
  if ( !parsar(argc,argv,&mesh,&met) )  return(1);

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&MMG5_ctim[1]);
  if ( !loadMesh(&mesh) )  return(1);
  met.npmax = mesh.npmax;
  met.dim   = 3;
  ier = MMG5_loadMet(&mesh,&met);
  if ( !ier )
    return(1);
  else if ( ier > 0 && met.np != mesh.np ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(&mesh,met.m,(met.size*met.npmax+1)*sizeof(double));
    memset(&met,0,sizeof(MMG5_Sol));
  }
  if ( !parsop(&mesh,&met) )     return(1);
  if ( !_MMG5_scaleMesh(&mesh,&met) )  return(1);
  chrono(OFF,&MMG5_ctim[1]);
  printim(MMG5_ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&MMG5_ctim[2]);

  inqua(&mesh,&met);
  fprintf(stdout,"\n  %s\n   MODULE MMGS-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
  if ( mesh.info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  if ( !analys(&mesh) )  return(1);
  chrono(OFF,&MMG5_ctim[2]);
  if ( mesh.info.imprim ) {
    printim(MMG5_ctim[2].gdif,stim);
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  }
  /* solve */
  chrono(ON,&MMG5_ctim[3]);
  if ( mesh.info.imprim )
    fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met.size < 6 ? "ISOTROPIC" : "ANISOTROPIC");
  if ( !mmgs1(&mesh,&met) )  return(1);
  chrono(OFF,&MMG5_ctim[3]);
  if ( mesh.info.imprim ) {
    printim(MMG5_ctim[3].gdif,stim);
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }
  fprintf(stdout,"\n  %s\n   END OF MODULE MMGS-LJLL \n  %s\n",MG_STR,MG_STR);

  /* save file */
  outqua(&mesh,&met);
  chrono(ON,&MMG5_ctim[1]);
  if ( mesh.info.imprim )  fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh.nameout);
  if ( !_MMG5_unscaleMesh(&mesh,&met) )  return(1);
  if ( !saveMesh(&mesh) )      return(1);
  if ( !saveMet(&mesh,&met) )  return(1);
  chrono(OFF,&MMG5_ctim[1]);
  if ( mesh.info.imprim )  fprintf(stdout,"  -- WRITING COMPLETED\n");

  /* release memory */
  if ( mesh.point )
    _MMG5_DEL_MEM(&mesh,mesh.point,(mesh.npmax+1)*sizeof(MMG5_Point));
  if ( mesh.adja )
    _MMG5_DEL_MEM(&mesh,mesh.adja,(3*mesh.ntmax+5)*sizeof(int));
  if ( mesh.tria )
    _MMG5_DEL_MEM(&mesh,mesh.tria,(mesh.ntmax+1)*sizeof(MMG5_Tria));
  if ( mesh.edge )
    _MMG5_DEL_MEM(&mesh,mesh.tria,(mesh.na+1)*sizeof(MMG5_Tria));
  if ( met.m )
    _MMG5_DEL_MEM(&mesh,met.m,(met.size*met.npmax+1)*sizeof(double));
  if ( mesh.info.par )
    _MMG5_DEL_MEM(&mesh,mesh.info.par,mesh.info.npar*sizeof(MMG5_Par));
  if ( mesh.xpoint )
    _MMG5_DEL_MEM(&mesh,mesh.xpoint,(mesh.xpmax+1)*sizeof(MMG5_xPoint));

  MMG5_Free_names(&mesh,&met);

  return(0);
}
