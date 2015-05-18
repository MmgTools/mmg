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

#include "mmg2d.h"

mytime   ctim[TIMEMAX];

int MMG2_iare[3][2] = {{1,2},{2,0},{0,1}};
int MMG2_iopp[3][2] = {{1,2},{0,2},{0,1}};
unsigned int MMG2_idir[5] = {0,1,2,0,1};

unsigned char _MMG5_iprv2[3] = {2,0,1};
unsigned char _MMG5_inxt2[3] = {1,2,0};
/**
 * Set common pointer functions between mmgs and mmg3d to the matching mmg3d
 * functions.
 */
void _MMG5_Set_commonFunc() {
  fprintf(stdout,"UNUSED FUNCTION IN MMG2D\n");
  return;
}


static void excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
  case SIGABRT:
    fprintf(stdout,"  Abnormal stop\n");  exit(1);
  case SIGFPE:
    fprintf(stdout,"  Floating-point exception\n"); exit(1);
  case SIGILL:
    fprintf(stdout,"  Illegal instruction\n"); exit(1);
  case SIGSEGV:
    fprintf(stdout,"  Segmentation fault\n");  exit(1);
  case SIGTERM:
  case SIGINT:
    fprintf(stdout,"  Program killed\n");  exit(1);
  }
  exit(1);
}


static void usage(char *name) {
  fprintf(stdout,"\nUsage: %s [-v [n]] [-h] [-m [n]] [opts..] -in filein [-out fileout]\n",name);

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-h      Print this message\n");
  fprintf(stdout,"-v [n]  Turn on numerical information, [-10..10]\n");
  fprintf(stdout,"-m [n]  Set memory size to n Mbytes\n");

  fprintf(stdout,"\n");


  fprintf(stdout,"-ar    val    angle detection (default 45.)\n");
  fprintf(stdout,"-nr           no ridge/corners detection \n");
  fprintf(stdout,"-hmin  val    minimal mesh size\n");
  fprintf(stdout,"-hmax  val    maximal mesh size\n");
  fprintf(stdout,"-hausd val    angle used for boundary collapses (default 135 degrees)\n");
  fprintf(stdout,"-hgrad val    mesh gradation (-1 = no gradation)\n");

  fprintf(stdout,"-noinsert     no insertion/suppression point\n");
  fprintf(stdout,"-noswap       no edge flipping\n");
  fprintf(stdout,"-nomove       no point relocation\n");

  fprintf(stdout,"\n");
  fprintf(stdout,"-ls             create mesh of isovalue 0\n");
  fprintf(stdout,"-lag [0/1/2]    Lagrangian mesh displacement according to mode 0/1/2\n");
  fprintf(stdout,"-mov filedep    (with -lag option)\n");
  fprintf(stdout,"-nsd val        only if no given triangle, save the subdomain nb (0==all subdomain)\n");
  fprintf(stdout,"-msh val        read and write to gmsh visu if val = 1 (out) if val=2 (in and out)\n");
  fprintf(stdout,"-degrad Qw Qdeg (with -lag option) : threshold for optimisation (default 10. 1.3)\n");
 
  /* fprintf(stdout,"-per          obsolete : to deal with periodic mesh on a square\n");*/

  fprintf(stdout,"\n\n");
  /* 
  
     fprintf(stdout,"-optim       mesh optimization\n");
     fprintf(stdout,"-nsurf       no surfacic modifications\n");
  */
  exit(1);
}


int parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol sol,double *qdegrad) {
  int     i;
  char   *ptr;

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
      case 'h':
        if ( !strcmp(argv[i],"-hgrad") ) {
          ++i;
          mesh->info.hgrad   = atof(argv[i]);
        }
        else if ( !strcmp(argv[i],"-hausd") ) {
          ++i;
          mesh->info.hausd   = atof(argv[i]);
        } else if ( !strcmp(argv[i],"-hmin") ) {
          ++i;
          mesh->info.hmin   = atof(argv[i]);
        } else if ( !strcmp(argv[i],"-hmax") ) {
          ++i;
          mesh->info.hmax   = atof(argv[i]);
        } else {
          fprintf(stderr,"  Unrecognized option %s\n",argv[i]);
          usage(argv[0]);
        }
      case 'a':
        if ( !strcmp(argv[i],"-ar") ) { 
          ++i;
          mesh->info.dhd = 180. - atof(argv[i]);
        }
        break;
      case 'd':
        if ( !strcmp(argv[i],"-degrad") ) {
          ++i;
          qdegrad[0] = atof(argv[i++])/ALPHA;
          qdegrad[1] = atof(argv[i]);
        } else {
          mesh->info.ddebug = 1;
          break;
        }
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          ++i;
          mesh->namein = argv[i];
          mesh->info.imprim   = 3;
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-lag") ) {
          ++i;
          mesh->info.lag  = atoi(argv[i]);
        } else if ( !strcmp(argv[i],"-ls") ) {
          ++i;
          mesh->info.iso  = atoi(argv[i]);
        }
        break;
      case 'm':
        if ( !strcmp(argv[i],"-mov") ) {
          ++i;
#warning todo potion 9
          // mesh->namedep = argv[i];
        } else if ( !strcmp(argv[i],"-msh") ) {
          ++i;
          mesh->info.nreg = atoi(argv[i]);
        }else if ( !strcmp(argv[i],"-m") ) {
          ++i;
          mesh->info.mem = atoi(argv[i]);
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-noswap") ) {
          mesh->info.noswap  = 1;
        } else if ( !strcmp(argv[i],"-nomove") ) {
          mesh->info.nomove  = 1;
        } else if ( !strcmp(argv[i],"-nr") ) {
          mesh->info.dhd  = -1;
        } else if ( !strcmp(argv[i],"-noinsert") ) {
          mesh->info.noinsert  = 1;
        } else if ( !strcmp(argv[i],"-nsubdomain") ) {
          ++i;
          mesh->info.renum = atoi(argv[i]);		    
        } else {
          printf("wrong option %s \n",argv[i]);
          usage(argv[0]);
        }
        break;
      case 'o':
        if ( !strcmp(argv[i],"-out") ) {
          ++i;
          mesh->nameout = argv[i];
          mesh->info.imprim = 5;
        }
        break;
      case 'p':
        if ( !strcmp(argv[i],"-per") ) {
          fprintf(stdout,"WARNING OBSOLETE OPTION\n");
          mesh->info.renum = 1;
        }
        break;
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          ++i;
          sol->namein = argv[i];
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
        fprintf(stderr,"  Unrecognized option %s\n",argv[i]);
        usage(argv[0]);
      }
    }

    else {
      if ( mesh->namein == NULL ) {
        mesh->namein = argv[i];
        if ( mesh->info.imprim == -99 )  mesh->info.imprim = 3;
      }
      else if ( mesh->nameout == NULL ) {
        mesh->nameout = argv[i];
      }
      else {
        fprintf(stdout,"  Argument %s ignored\n",argv[i]);
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
    mesh->namein = (char *)calloc(128,sizeof(char));
    assert(mesh->namein);
    fprintf(stdout,"  -- FILE BASENAME ?\n");
    fflush(stdin); 
    fscanf(stdin,"%s",mesh->namein);
  }
  if ( mesh->nameout == NULL ) {
    mesh->nameout = (char *)calloc(128,sizeof(char));
    assert(mesh->nameout);
    strcpy(mesh->nameout,mesh->namein);
    ptr = strstr(mesh->nameout,".mesh");
    if ( ptr ) {
      *ptr = '\0';  
      strcat(mesh->nameout,".o.mesh");
    } else {
      ptr = strstr(mesh->nameout,".mesh");
      //ptr = strstr(mesh->nameout,".meshb");
      if ( ptr ) {
        *ptr = '\0';  
        //strcat(mesh->nameout,".o.meshb");        
        strcat(mesh->nameout,".o.mesh");        
      } else { 
        strcat(mesh->nameout,".o.mesh");     
      }
    }
  }
  if ( sol->namein == NULL ) {
    sol->namein = (char *)calloc(128,sizeof(char));
    assert(sol->namein);
    strcpy(sol->namein,mesh->namein);
    ptr = strstr(sol->namein,".mesh");
    if ( ptr ) *ptr = '\0';
  }
  /* if ( mesh->namedep == NULL ) { */
  /*   mesh->namedep = (char *)calloc(128,sizeof(char)); */
  /*   assert(mesh->namedep); */
  /*   strcpy(mesh->namedep,mesh->namein); */
  /*   ptr = strstr(mesh->namedep,".mesh"); */
  /*   if ( ptr ) *ptr = '\0'; */
  /* } */

  return(1);
}


static void endcod() {
  double   ttot,ttim[TIMEMAX];
  int      k,call[TIMEMAX];

  //chrono(OFF,&ctim[0]);
#warning endcod : comment for merge
  /* for (k=0; k<TIMEMAX; k++) { */
  /*   call[k] = ctim[k].call; */
  /*   ttim[k] = ctim[k].call ? gttime(ctim[k]) : 0.0; */
  /* } */
  /* ttot    = ttim[1]+ttim[2]+ttim[3]+ttim[4]; */
  /* ttim[0] = M_MAX(ttim[0],ttot); */

  /* fprintf(stdout,"\n  -- CPU REQUIREMENTS\n"); */
  /* fprintf(stdout,"  in/out	%8.2f %%    %3d. calls,   %7.2f sec/call\n", */
  /*         100.*ttim[1]/ttim[0],call[1],ttim[1]/(float)call[1]); */
  /* fprintf(stdout,"  analysis	%8.2f %%    %3d. calls,   %7.2f sec/call\n", */
  /*         100.*ttim[2]/ttim[0],call[2],ttim[2]/(float)call[2]); */
  /* fprintf(stdout,"  remeshing	%8.2f %%    %3d. calls,   %7.2f sec/call\n", */
  /*         100.*ttim[3]/ttim[0],call[3],ttim[3]/(float)call[3]); */
  /* fprintf(stdout,"  total     %8.2f %%    %3d. calls,	%7.2f sec/call\n", */
  /*         100.*ttot/ttim[0],call[0],ttot/(float)call[0]); */

  /* fprintf(stdout,"\n   ELAPSED TIME  %.2f SEC.  (%.2f)\n",ttim[0],ttot); */
}


/* set function pointers */
int setfunc(int type) {
  if ( type == 3 ) {
    MMG2_length = long_ani;
    MMG2_caltri = caltri_ani;
    MMG2_buckin = buckin_ani;
    MMG2_lissmet = lissmet_ani;
    MMG2_optlen     = optlen_ani;
/*    interp     = interp_ani;
 */
  }
  else {
    MMG2_length  = long_iso;
    MMG2_caltri  = caltri_iso;
    MMG2_buckin  = buckin_iso;
    MMG2_lissmet = lissmet_iso;

    MMG2_optlen     = optlen_iso;
/*    interp     = interp_iso;
 */
  }

  return(1);
}


int main(int argc,char *argv[]) {
  MMG5_Mesh      mesh;
  MMG5_Sol	    sol;
  double declic,qdegrad[2];  
  int    nsiter,k;

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",M_VER,M_REL);
  fprintf(stdout,"     %s\n",M_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  /* interrupts */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  atexit(endcod);

  //tminit(ctim,TIMEMAX);
  //chrono(ON,&ctim[0]);

  /* default values */
  memset(&mesh,0,sizeof(MMG5_Mesh));
  memset(&sol,0,sizeof(MMG5_Sol));
  mesh.info.imprim = -99;
  mesh.info.mem    = 500;
  mesh.info.ddebug = 0;
  mesh.info.lag    = -1;
  mesh.info.iso    = 0;
  mesh.info.noswap = 0;
  mesh.info.nomove = 0;
  mesh.info.noinsert = 0;
  mesh.info.hgrad  = 1.3;
  mesh.info.renum  = 0;
  mesh.info.hausd  = 135.;
  mesh.info.dhd  = 45.;
  mesh.info.hmin     = -1.;    
  mesh.info.hmax     = -1.;    
  
  qdegrad[0] = 10./ALPHA;
  qdegrad[1] = 1.3;   
  mesh.info.renum = 0;
  mesh.info.nreg = 0;
  sol.type = 1;

  if ( !parsar(argc,argv,&mesh,&sol,qdegrad) )  return(1);

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  //chrono(ON,&ctim[1]);
  if ( !MMG2_loadMesh(&mesh,mesh.namein) )  return(1); 
  if ( !MMG2_loadSol(&sol,sol.namein,mesh.npmax) )  {    
    sol.np = mesh.np;
    sol.size = 1;
    sol.ver  = mesh.ver;
    /* mem alloc */
    sol.m = (double*)M_calloc((sol.size*mesh.npmax)+1,sizeof(double),"mmg");
    assert(sol.m);  
    sol.np = 0;
  } else   if ( sol.np && (sol.np != mesh.np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER : %d != %d\n",sol.np,mesh.np);
    //exit(1);
  }
  setfunc(sol.size);//setfunc(sol.type);
  
 
  if(mesh.info.lag>=0 ) {
#warning option 9
    /* /\*alloc Disp*\/ */
    /* mesh.disp.mv = (double*)M_calloc(2*(mesh.npmax + 1),sizeof(double),"displ.mv"); */
    /* assert(mesh.disp.mv); */
    /* mesh.disp.alpha = (short*)M_calloc(mesh.npmax+1,sizeof(short),"displ.alpha"); */
    /* assert(mesh.disp.alpha); */
    /* if ( !MMG2_loadVect(&mesh,mesh.namedep) )  return(1);   */    
  } 
  //chrono(OFF,&ctim[1]);
  // fprintf(stdout,"  -- DATA READING COMPLETED.     %.2f sec.\n",
  //      gttime(ctim[1]));


  fprintf(stdout,"\n  %s\n   MODULE MMG2D-IMB/LJLL : %s (%s) %s\n  %s\n",
          M_STR,M_VER,M_REL,sol.size == 1 ? "ISO" : "ANISO",M_STR);
  fprintf(stdout,"  MAXIMUM NUMBER OF POINTS    (NPMAX) : %8d\n",mesh.npmax);
  fprintf(stdout,"  MAXIMUM NUMBER OF TRIANGLES (NTMAX) : %8d\n",mesh.ntmax);

  /* analysis */
  //chrono(ON,&ctim[2]);
  if ( mesh.info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : DATA ANALYSIS\n");
  if ( abs(mesh.info.imprim) > 4 )
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");
  if ( !MMG2_scaleMesh(&mesh,&sol) )  return(1);
  if ( !sol.np && !MMG2_doSol(&mesh,&sol) )  return(1);

  if ( mesh.nt && !MMG2_hashel(&mesh) )  return(1);
  if ( !mesh.info.renum && !MMG2_chkmsh(&mesh,1) )        return(1);
  /*geom : corner detection*/    
  if ( mesh.info.dhd>0. )
    if( !MMG2_evalgeom(&mesh) ) return(1);

  /*mesh gradation*/
  if( mesh.nt && mesh.info.hgrad > 0 ) {
    if ( mesh.info.imprim )   fprintf(stdout,"\n  -- GRADATION : %8f\n",mesh.info.hgrad);
    MMG2_lissmet(&mesh,&sol); 
  }
  if ( mesh.nt && abs(mesh.info.imprim) > 1 )  MMG2_outqua(&mesh,&sol);
  
  if ( mesh.nt && abs(mesh.info.imprim) > 1 )  {
    MMG2_prilen(&mesh,&sol);
  }                       
  
  //chrono(OFF,&ctim[2]);
  // fprintf(stdout,"  -- PHASE 1 COMPLETED.     %.2f sec.\n",gttime(ctim[2]));

  /* remeshing */
  //chrono(ON,&ctim[3]);
  /* specific meshing */
  if ( mesh.info.iso ) {
    fprintf(stdout,"Fit an embedded mesh\n");
    MMG2_mmg2d6(&mesh,&sol);
    MMG2_saveMesh(&mesh,mesh.nameout);
    return(0);
  } else if ( mesh.info.lag >= 0 ) {
    /* if ( !MMG2_mmg2d9(&mesh,&sol,&qdegrad) ) { */
    /*   M_free(mesh.disp.mv); */
    /*   M_free(mesh.disp.alpha); */
        
    /*   if ( abs(mesh.info.imprim) > 3 )  { */
    /*       MMG2_outqua(&mesh,&sol); */
    /*       MMG2_prilen(&mesh,&sol);  */
    /*   } */
    /*   chrono(ON,&ctim[1]); */
    /*   fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh.nameout); */
    /*   MMG2_saveMesh(&mesh,mesh.nameout); */
    /*   MMG2_saveSol(&mesh,&sol,mesh.nameout); */
    /*   fprintf(stdout,"  -- WRITING COMPLETED\n"); */
    /*   chrono(OFF,&ctim[1]); */
        
    /*   /\* free memory*\/ */
    /*   M_free(mesh.point); */
    /*   M_free(mesh.tria); */
    /*   M_free(mesh.edge); */
    /*   //M_free(mesh.adja); */
    /*   M_free(sol.met);     */
    /*   return(1);    */
    /* }  */
    /* M_free(mesh.disp.mv); */
    /* M_free(mesh.disp.alpha); */
#warning option 9
    printf("exit option 9 not implemented\n");
    exit(1);
  } else {
 
    if(!mesh.nt) {
      fprintf(stdout,"\n  -- PHASE 2 : MESH GENERATION\n");
      if ( !MMG2_mmg2d2(&mesh,&sol) )  return(1);
    } else {
      fprintf(stdout,"\n  -- PHASE 2 : MESH ADAPTATION\n");
      if ( (!mesh.info.noinsert) && !MMG2_mmg2d1(&mesh,&sol) )  return(1);
    }
    
    /* optimisation */
    //chrono(ON,&ctim[4]);
    fprintf(stdout,"\n  -- PHASE 3 : MESH OPTIMISATION\n");
    //if ( !optlap(&mesh,&sol) ) return(1);
    if ( !MMG2_mmg2d0(&mesh,&sol) )  return(1);  
    // chrono(OFF,&ctim[4]);
    //fprintf(stdout,"  -- PHASE 3 COMPLETED.     %.2f sec.\n",gttime(ctim[4]));
      
  }
  // chrono(OFF,&ctim[3]);
  //fprintf(stdout,"  -- PHASE 2 COMPLETED.     %.2f sec.\n",gttime(ctim[3]));


  
  if ( !MMG2_unscaleMesh(&mesh,&sol) )  return(1);
 
  if ( abs(mesh.info.imprim) > 1 )  {
    MMG2_outqua(&mesh,&sol);
    MMG2_prilen(&mesh,&sol); 
  }
  //chrono(ON,&ctim[1]);
  fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh.nameout);
  MMG2_saveMesh(&mesh,mesh.nameout);
  if( sol.np )
    MMG2_saveSol(&mesh,&sol,mesh.nameout);
  fprintf(stdout,"  -- WRITING COMPLETED\n");
  //chrono(OFF,&ctim[1]);

  /* free memory*/
  M_free(mesh.point);
  M_free(mesh.tria);
  M_free(mesh.edge);
  free(mesh.adja);
  M_free(sol.m);    
  
  if ( mesh.info.imprim < -4 )  M_memDump();   
  return(0);
}
  

