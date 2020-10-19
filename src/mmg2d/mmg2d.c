/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
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

mytime   MMG5_ctim[TIMEMAX];

/**
 * Print elapsed time at end of process.
 */
static void MMG5_endcod(void) {
  char   stim[32];

  chrono(OFF,&MMG5_ctim[0]);
  printim(MMG5_ctim[0].gdif,stim);
  fprintf(stdout,"\n   ELAPSED TIME  %s\n",stim);
}

static int MMG2D_usage(char *name) {

  /* Common generic options, file options and mode options */
  MMG5_mmgUsage(name);

  /* Lagrangian option (only for mmg2d/3d) */
  MMG5_lagUsage();

  /* Common parameters (first section) */
  MMG5_paramUsage1();

  /* Parameters shared by mmg2d and 3d only*/
  MMG5_2d3dUsage();

  /* Specific parameters */
  fprintf(stdout,"-3dMedit val read and write to gmsh visu if val = 1 (out) if val=2 (in and out)\n");

  fprintf(stdout,"\n");

  fprintf(stdout,"-nosurf      no surface modifications\n");

  /* Common parameters (second section) */
  MMG5_paramUsage2();

  /* Common options for advanced users */
  MMG5_advancedUsage();

  fprintf(stdout,"\n\n");

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
static inline int MMG5_defaultValues(MMG5_pMesh mesh) {

  MMG5_mmgDefaultValues(mesh);

  fprintf(stdout,"\n\n");

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param bdyRefs pointer toward the list of the boundary references.
 * \return npar, the number of local parameters at edges if success,
 * 0 otherwise.
 *
 * Count the local default values at edges and fill the list of the boundary
 * references.
 *
 */
static inline
int MMG2D_countLocalParamAtEdg( MMG5_pMesh mesh,MMG5_iNode **bdyRefs) {
  int         npar,k,ier;

  /** Count the number of different boundary references and list it */
  (*bdyRefs) = NULL;

  k = mesh->na? mesh->edge[1].ref : 0;

  /* Try to alloc the first node */
  ier = MMG5_Add_inode( mesh, bdyRefs, k );
  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate the first boundary"
           " reference node.\n",__func__);
    return 0;
  }
  else {
    assert(ier);
    npar = 1;
  }

  for ( k=1; k<=mesh->na; ++k ) {
    ier = MMG5_Add_inode( mesh, bdyRefs, mesh->edge[k].ref );

    if ( ier < 0 ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to list the edge references.\n"
              "              Uncomplete parameters file.\n",__func__);
      break;
    }
    else if ( ier ) ++npar;
  }

  return npar;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param bdryRefs pointer toward the list of the boundary references.
 * \param out pointer toward the file in which to write.
 * \return 1 if success, 0 otherwise.
 *
 * Write the local default values at edges in the parameter file.
 *
 */
static inline
int MMG2D_writeLocalParamAtEdg( MMG5_pMesh mesh, MMG5_iNode *bdryRefs,
                                 FILE *out ) {
  MMG5_iNode *cur;

  cur = bdryRefs;
  while( cur ) {
    fprintf(out,"%d Edge %e %e %e \n",cur->val,
            mesh->info.hmin, mesh->info.hmax,mesh->info.hausd);
    cur = cur->nxt;
  }

  MMG5_Free_ilinkedList(mesh,bdryRefs);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 otherwise.
 *
 * Write a DEFAULT.mmg2d file containing the default values of parameters that
 * can be locally defined.
 *
 */
static inline
int MMG2D_writeLocalParam( MMG5_pMesh mesh ) {
  MMG5_iNode  *edgRefs,*triRefs;
  int          nparEdg,nparTri;
  char         *ptr,data[MMG5_FILESTR_LGTH];
  FILE         *out;

  /** Save the local parameters file */
  strcpy(data,mesh->namein);

  ptr = MMG5_Get_filenameExt(data);

  if ( ptr ) *ptr = '\0';
  strcat(data,".mmg2d");

  if ( !(out = fopen(data,"wb")) ) {
    fprintf(stderr,"\n  ** UNABLE TO OPEN %s.\n",data);
    return 0;
  }

  fprintf(stdout,"\n  %%%% %s OPENED\n",data);

  nparEdg = MMG2D_countLocalParamAtEdg( mesh, &edgRefs);
  if ( !nparEdg ) {
    fclose(out);
    return 0;
  }

  nparTri = MMG5_countLocalParamAtTri( mesh, &triRefs);
  if ( !nparTri ) {
    fclose(out);
    return 0;
  }

  fprintf(out,"parameters\n %d\n",nparTri+nparEdg);

  /** Write local param at triangles */
  if (! MMG2D_writeLocalParamAtEdg(mesh,edgRefs,out) ) {
    fclose(out);
    return 0;
  }

  /** Write local param at tetra */
  if (! MMG5_writeLocalParamAtTri(mesh,triRefs,out) ) {
    fclose(out);
    return 0;
  }

  fclose(out);
  fprintf(stdout,"  -- WRITING COMPLETED\n");

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (metric).
 * \param sol pointer toward a sol structure (metric).
 *
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 *
 * Program to save the local default parameter file: read the mesh and metric
 * (needed to compite the hmax/hmin parameters), scale the mesh and compute the
 * hmax/hmin param, unscale the mesh and write the default parameter file.
 *
 */
static inline
int MMG2D_defaultOption(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  mytime    ctim[TIMEMAX];
  double    hsiz;
  char      stim[32];

  signal(SIGABRT,MMG2D_excfun);
  signal(SIGFPE,MMG2D_excfun);
  signal(SIGILL,MMG2D_excfun);
  signal(SIGSEGV,MMG2D_excfun);
  signal(SIGTERM,MMG2D_excfun);
  signal(SIGINT,MMG2D_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  if ( mesh->info.npar ) {
    fprintf(stderr,"\n  ## Error: %s: "
            "unable to save of a local parameter file with"
            " the default parameters values because local parameters"
            " are provided.\n",__func__);
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
  }


  if ( mesh->info.imprim > 0 ) fprintf(stdout,"\n  -- INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));

  if ( met && met->np && (met->np != mesh->np) ) {
    fprintf(stderr,"\n  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    MMG5_DEL_MEM(mesh,met->m);
    met->np = 0;
  }

  if ( sol && sol->np && (sol->np != mesh->np) ) {
    fprintf(stderr,"\n  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    MMG5_DEL_MEM(mesh,sol->m);
    sol->np = 0;
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMG2D_setfunc(mesh,met);
  MMG2D_Set_commonFunc();

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  %s\n   MODULE MMG2D: IMB-LJLL : %s (%s)\n  %s\n",
            MG_STR,MMG_VERSION_RELEASE,MMG_RELEASE_DATE,MG_STR);
    fprintf(stdout,"\n  -- DEFAULT PARAMETERS COMPUTATION\n");
  }

  /* scaling mesh and hmin/hmax computation*/
  if ( !MMG2D_scaleMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);

  /* specific meshing + update hmin/hmax */
  if ( mesh->info.optim ) {
    if ( !MMG2D_doSol(mesh,met) ) {
      if ( !MMG5_unscaleMesh(mesh,met,sol) )
        _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
    }
    MMG2D_solTruncatureForOptim(mesh,met);
  }
  if ( mesh->info.hsiz > 0. ) {
    if ( !MMG5_Compute_constantSize(mesh,met,&hsiz) ) {
     if ( !MMG5_unscaleMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
     _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }

  /* unscaling mesh */
  if ( !MMG5_unscaleMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);

  /* Save the local parameters file */
  mesh->mark = 0;
  if ( !MMG2D_writeLocalParam(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: Unable to save the local parameters file.\n"
            "            Exit program.\n",__func__);
     _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
  }

  _LIBMMG5_RETURN(mesh,met,sol,MMG5_SUCCESS);
}

// In ls mode : metric must be provided using -met option (-sol or default is the ls).
// In adp mode : -sol or -met or default allow to store the metric.
int parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  MMG5_pSol tmp = NULL;
  int     i;
  char    namein[MMG5_FILESTR_LGTH];

  /* First step: search if user want to see the default parameters values. */
  for ( i=1; i< argc; ++i ) {
    if ( !strcmp(argv[i],"-val") ) {
      MMG5_defaultValues(mesh);
      return 0;
    }
  }

  /* Second step: read all other arguments. */
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        MMG2D_usage(argv[0]);
        return 0;
      case 'a':
        if ( !strcmp(argv[i],"-ar") && ++i < argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_angleDetection,
                                     atof(argv[i])) )
            return 0;
        }
        break;
      case 'A': /* anisotropy */
        if ( !MMG2D_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor) )
          return 0;
        break;
      case 'd':
        if ( !strcmp(argv[i],"-default") ) {
          mesh->mark=1;
        } else {  /* debug */
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_debug,1) )
            return 0;
        }
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_hmin,
                                     atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_hmax,
                                     atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hsiz") && ++i < argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_hsiz,
                                     atof(argv[i])) )
            return 0;

        }
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_hausd,
                                     atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hgradreq") && ++i <= argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_hgradreq,
                                     atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_hgrad,
                                     atof(argv[i])) )
            return 0;
        }
        else {
          MMG2D_usage(argv[0]);
          return 0;
        }
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMG2D_Set_inputMeshName(mesh, argv[i]) )
              return 0;

            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,5) )
              return 0;
          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-lag") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_lag,atoi(argv[i])) )
              return 0;
          }
          else if ( i == argc ) {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-ls") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_iso,1) )
            return 0;
          if ( ++i < argc && (isdigit(argv[i][0]) ||
                              (argv[i][0]=='-' && isdigit(argv[i][1])) ) ) {
            if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_ls,atof(argv[i])) )
              return 0;
          }
          else i--;
        }
        break;
      case 'm':  /* memory */
        if ( !strcmp(argv[i],"-met") ) {
          if ( !met ) {
            fprintf(stderr,"No metric structure allocated for %c%c%c option\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            return 0;
          }
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG2D_Set_inputSolName(mesh,met,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        else if (!strcmp(argv[i],"-m") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_mem,atoi(argv[i])) )
              return 0;
          }
          else {
            fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nreg") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nreg,1) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_angle,0) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nsd") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_numsubdomain,atoi(argv[i])) )
              return 0;
          }
          else {
            fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        } else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_noswap,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_noinsert,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nomove,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nosurf") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nosurf,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nosizreq") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nosizreq,1) ) {
            return 0;
          }
        }
        break;
      case 'o':
        if ( (!strcmp(argv[i],"-out")) || (!strcmp(argv[i],"-o")) ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            if ( !MMG2D_Set_outputMeshName(mesh,argv[i]) )
              return 0;
          }else{
            fprintf(stderr,"Missing filname for %c%c%c\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-opnbdy") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_opnbdy,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-optim") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_optim,1) )
            return 0;
        }
        break;
      case 'r':
        if ( !strcmp(argv[i],"-rmc") ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_rmc,0) )
            return 0;
          if ( ++i < argc && (isdigit(argv[i][0]) ) ) {
            if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_rmc,atof(argv[i])) )
              return 0;
          }
          else i--;
        }
        break;
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          /* For retrocompatibility, store the metric if no sol structure available */
          tmp = sol ? sol : met;

          assert(tmp);
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG2D_Set_inputSolName(mesh,tmp,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,atoi(argv[i])) )
              return 0;
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          MMG2D_usage(argv[0]);
          return 0;
        }
        break;
      case '3':
        if(!strcmp(argv[i],"-3dMedit") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_3dMedit,atoi(argv[i])) )
              return 0;
          }
          else {
            fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        break;
      default:
        fprintf(stderr,"Unrecognized option %s\n",argv[i]);
        MMG2D_usage(argv[0]);
        return 0;
      }

    }

    else {
      if ( mesh->namein == NULL ) {
        if ( !MMG2D_Set_inputMeshName(mesh,argv[i]) )
          return 0;
        if ( mesh->info.imprim == -99 )  {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,5) )
            return 0;
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMG2D_Set_outputMeshName(mesh,argv[i]) )
          return 0;
      }
      else {
        fprintf(stdout,"  Argument %s ignored\n",argv[i]);
        MMG2D_usage(argv[0]);
        return 0;
      }
    }
    i++;
  }

  /** check file names */
  if ( mesh->info.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    MMG_FSCANF(stdin,"%d",&i);
    if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,i) )
      return 0;
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    MMG_FSCANF(stdin,"%127s",namein);
    if ( !MMG2D_Set_inputMeshName(mesh,namein) )
      return 0;
  }
  if ( mesh->nameout == NULL ) {
    if ( !MMG2D_Set_outputMeshName(mesh,"") )
      return 0;
  }

  /* adp mode: if the metric name has been stored in sol, move it in met */
  if ( met->namein==NULL && sol && sol->namein && !(mesh->info.iso || mesh->info.lag>=0) ) {
    if ( !MMG2D_Set_inputSolName(mesh,met,sol->namein) )
      return 0;
    MMG5_DEL_MEM(mesh,sol->namein);
  }

  /* default : store solution (resp. displacement) name in iso
   * (resp. lagrangian) mode, metric name otherwise */
  tmp = ( mesh->info.iso || mesh->info.lag >=0 ) ? sol : met;
  assert ( tmp );
  if ( tmp->namein == NULL ) {
    if ( !MMG2D_Set_inputSolName(mesh,tmp,"") )
      return 0;
  }
  if ( met->nameout == NULL ) {
    if ( !MMG2D_Set_outputSolName(mesh,met,"") )
      return 0;
  }

  return 1;
}

int main(int argc,char *argv[]) {
  MMG5_pMesh    mesh;
  MMG5_pSol     sol,met,disp,ls;
  int           ier,ierSave,fmtin,fmtout;
  char          stim[32],*ptr;

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",MMG_VERSION_RELEASE,MMG_RELEASE_DATE);
  fprintf(stdout,"     %s\n",MMG_COPYRIGHT);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  /* Print timer at exit */
  atexit(MMG5_endcod);

  MMG2D_Set_commonFunc();
  tminit(MMG5_ctim,TIMEMAX);
  chrono(ON,&MMG5_ctim[0]);

  /* assign default values */
  mesh = NULL;
  met  = NULL;
  ls   = NULL;
  disp = NULL;

  if ( !MMG2D_Init_mesh(MMG5_ARG_start,
                        MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                        MMG5_ARG_ppLs,&ls,
                        MMG5_ARG_ppDisp,&disp,
                        MMG5_ARG_end) )
    return MMG5_STRONGFAILURE;

  /* reset default values for file names */
  if ( !MMG2D_Free_names(MMG5_ARG_start,
                         MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                         MMG5_ARG_ppLs,&ls,
                         MMG5_ARG_ppDisp,&disp,
                         MMG5_ARG_end) )
    return MMG5_STRONGFAILURE;

  /* Set default metric size */
  if ( !MMG2D_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Scalar) )
    MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);

  /* Read command line */
  if ( !parsar(argc,argv,mesh,met,ls) )  return MMG5_STRONGFAILURE;

  /* load data */
  if ( mesh->info.imprim >= 0 )
    fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&MMG5_ctim[1]);

  /* For each mode: pointer over the solution structure to load */
  if ( mesh->info.lag >= 0 ) {
    sol = disp;
  }
  else if ( mesh->info.iso ) {
    sol = ls;
  }
  else {
    sol = met;
  }

  /* read mesh/sol files */
  ptr   = MMG5_Get_filenameExt(mesh->namein);
  fmtin = MMG5_Get_format(ptr,MMG5_FMT_MeditASCII);

  switch ( fmtin ) {
  case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
    ier = MMG2D_loadMshMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_VtkVtp ):
    ier = MMG2D_loadVtpMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_VtkVtu ):
    ier = MMG2D_loadVtuMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_VtkVtk ):
    ier = MMG2D_loadVtkMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_MeditASCII ): case ( MMG5_FMT_MeditBinary ):
    ier = MMG2D_loadMesh(mesh,mesh->namein);
    if ( ier <  1 ) { break; }

    /* Read displacement in lag mode */
    if ( mesh->info.lag >= 0 ) {
      /* In Lagrangian mode, the name of the displacement file has been parsed in ls */
      if ( !MMG2D_Set_inputSolName(mesh,disp,ls->namein) ) {
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
      MMG5_DEL_MEM(mesh,ls->namein);
    }

    if ( mesh->info.lag >= 0 || mesh->info.iso ) {
      /* displacement or isovalue are mandatory */
      if (  MMG2D_loadSol(mesh,sol,sol->namein) < 1 ) {
        /* displacement or isovalue are mandatory */
        fprintf(stdout,"  ## ERROR: UNABLE TO LOAD SOLUTION.\n");
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
    }
    else {
      /* Facultative metric */
      if ( MMG2D_loadSol(mesh,met,met->namein) == -1 ) {
        fprintf(stdout,"\n  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
    }
    /* In iso mode: read metric if any */
    if ( mesh->info.iso && met->namein ) {
      if (  MMG2D_loadSol(mesh,met,met->namein) < 1 ) {
        fprintf(stdout,"  ## ERROR: UNABLE TO LOAD METRIC.\n");
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
    }
    break;
  default:
    fprintf(stderr,"  ** I/O AT FORMAT %s NOT IMPLEMENTED.\n",MMG5_Get_formatName(fmtin) );
    MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
  }

  if ( ier < 1) {
    if ( ier==0 ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",mesh->namein);
      fprintf(stderr,"  ** UNABLE TO OPEN INPUT FILE.\n");
    }
    MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
  }

  /* Check input data */
  if ( mesh->info.lag >= 0 ) {
    if ( met->namein ) {
      fprintf(stdout,"  ## WARNING: MESH ADAPTATION UNAVAILABLE IN"
              " LAGRANGIAN MODE. METRIC IGNORED.\n");
      MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
    }
  }
  else if ( mesh->info.iso ) {
    if ( ls->m == NULL ) {
      fprintf(stderr,"\n  ## ERROR: NO ISOVALUE DATA.\n");
      MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
    }
  }

  /* Read parameter file */
  if ( !MMG2D_parsop(mesh,met) )
    MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);

  chrono(OFF,&MMG5_ctim[1]);
  if ( mesh->info.imprim >= 0 ) {
    printim(MMG5_ctim[1].gdif,stim);
    fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);
  }

  if ( mesh->mark ) {
    /* Save a local parameters file containing the default parameters */
    ier = MMG2D_defaultOption(mesh,met,disp);
    MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,ier);
  }
  /* Lagrangian mode */
  else if ( mesh->info.lag > -1 ) {
    ier = MMG2D_mmg2dmov(mesh,met,disp);
  }
  /* Level Set mode */
  else if ( mesh->info.iso ) {
    ier = MMG2D_mmg2dls(mesh,ls,met);
  }
  /* Mesh generation mode */
  else if ( !mesh->nt ) {
    ier = MMG2D_mmg2dmesh(mesh,met);
  }
  /* Remeshing mode */
  else {
    if ( met && ls && met->namein && ls->namein ) {
      fprintf(stdout,"\n  ## ERROR: IMPOSSIBLE TO PROVIDE BOTH A METRIC"
              " AND A SOLUTION IN ADAPTATION MODE.\n");
      MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
    }

    ier = MMG2D_mmg2dlib(mesh,met);
  }

  if ( ier != MMG5_STRONGFAILURE ) {
    chrono(ON,&MMG5_ctim[1]);
    if ( mesh->info.imprim > 0 )
      fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh->nameout);

    ptr    = MMG5_Get_filenameExt(mesh->nameout);
    fmtout = MMG5_Get_format(ptr,fmtin);

    switch ( fmtout ) {
    case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
      ierSave = MMG2D_saveMshMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_VtkVtp ):
      ierSave = MMG2D_saveVtpMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_VtkVtu ):
      ierSave = MMG2D_saveVtuMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_VtkVtk ):
      ierSave = MMG2D_saveVtkMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_Tetgen ):
      ierSave = MMG2D_saveTetgenMesh(mesh,mesh->nameout);
      /* This format dont allow to save a solution: use a .sol file */
      if ( !ierSave ) {
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
      if ( met && met->np ) {
        ierSave = MMG2D_saveSol(mesh,met,mesh->nameout);
      }
      break;
    default:
      ierSave = MMG2D_saveMesh(mesh,mesh->nameout);
      if ( !ierSave ) {
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
      if ( met && met->np ) {
        ierSave = MMG2D_saveSol(mesh,met,mesh->nameout);
      }
      break;
    }

    if ( !ierSave )
      MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);

    chrono(OFF,&MMG5_ctim[1]);
    if ( mesh->info.imprim > 0 ) fprintf(stdout,"  -- WRITING COMPLETED\n");
  }

  /* free mem */
  MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,ier);
}
