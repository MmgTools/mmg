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
 * \file mmg3d/pampautils.c
 * \brief API functions only usefull for the PaMPA library.
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning Some functions are copying from \ref mmg3d/shared_func.c.
 **/

#include "mmg3d.h"

/* COPY OF PART OF shared_func.c */

/**
 * \param mesh pointer toward the mesh structure.
 * \warning Copy of the \ref warnOrientation function of the
 * \ref mmg3d/shared_func.h file.
 *
 * Warn user that some tetrahedra of the mesh have been reoriented.
 *
 */
static inline
void _MMG5_pampa_warnOrientation(MMG5_pMesh mesh) {
  if ( mesh->xt ) {
    if ( mesh->xt != mesh->ne ) {
      fprintf(stdout,"  ## Warning: %d tetra on %d reoriented.\n",
              mesh->xt,mesh->ne);
      fprintf(stdout,"  Your mesh may be non-conform.\n");
    }
    else {
      fprintf(stdout,"  ## Warning: all tetra reoriented.\n");
    }
  }
  mesh->xt = 0;
}

/**
 * \param sigid signal number.
 * \warning Copy of the \a excfun function of the \ref mmg3d/shared_func.h
 * file.
 *
 * Signal handling: specify error messages depending from catched signal.
 *
 */
static inline
void _MMG5_pampa_excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
  case SIGABRT:
    fprintf(stdout,"  *** potential lack of memory.\n");  break;
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
  exit(EXIT_FAILURE);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \warning Copy of the \a setfunc function of the \ref mmg3d/shared_func.h
 * file.
 *
 * Set function pointers for lenedgeCoor, _MMG5_hashTetra and saveMesh.
 *
 */
void MMG5_pampa_setfunc(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( met->size < 6 )
    MMG5_lenedgCoor = _MMG5_lenedgCoor_iso;
  else
    MMG5_lenedgCoor = _MMG5_lenedgCoor_ani;
  MMG5_hashTetra = _MMG5_hashTetra;
  MMG5_saveMesh = _MMG5_saveLibraryMesh;
}
/* END COPY */

/**
 * \brief Return adjacent elements of a tetrahedron.
 * \param mesh pointer toward the mesh structure.
 * \param kel tetrahedron index.
 * \param *v0 pointer toward the index of the adjacent element of \a kel through
 * its face number 0.
 * \param *v1 pointer toward the index of the adjacent element of \a kel through
 * its face number 1.
 * \param *v2 pointer toward the index of the adjacent element of \a kel through
 * its face number 2.
 * \param *v3 pointer toward the index of the adjacent element of \a kel through
 * its face number 3.
 * \return 1.
 *
 * Find the indices of the 4 adjacent elements of tetrahedron \a
 * kel. \f$v_i = 0\f$ if the \f$i^{th}\f$ face has no adjacent element
 * (so we are on a boundary face).
 *
 */
int MMG5_Get_adjaTet(MMG5_pMesh mesh, int kel, int *v0, int *v1, int *v2, int *v3) {

  if ( ! mesh->adja ) {
    if (! _MMG5_hashTetra(mesh, 0))
      return(0);
  }

  (*v0) = mesh->adja[4*(kel-1)+1]/4;
  (*v1) = mesh->adja[4*(kel-1)+2]/4;
  (*v2) = mesh->adja[4*(kel-1)+3]/4;
  (*v3) = mesh->adja[4*(kel-1)+4]/4;

  return(1);
}

/**
 * \param *prog pointer toward the program name.
 *
 * Print help for mmg3d5 options.
 *
 */
void _MMG5_usage(char *prog) {

  _MMG5_mmgUsage(prog);

  fprintf(stdout,"-lag [0/1/2] Lagrangian mesh displacement according to mode 0/1/2\n");
  fprintf(stdout,"-ls     val  create mesh of isovalue val\n");
  fprintf(stdout,"-optim       mesh optimization\n");
  fprintf(stdout,"-noinsert    no point insertion/deletion \n");
  fprintf(stdout,"-noswap      no edge or face flipping\n");
  fprintf(stdout,"-nomove      no point relocation\n");
  fprintf(stdout,"-nsurf       no surfacic modifications\n");
#ifndef PATTERN
  fprintf(stdout,"-bucket val  Specify the size of bucket per dimension \n");
#endif
#ifdef USE_SCOTCH
  fprintf(stdout,"-rn [n]      Turn on or off the renumbering using SCOTCH [1/0] \n");
#endif
  exit(EXIT_FAILURE);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
void _MMG5_defaultValues(MMG5_pMesh mesh) {

  _MMG5_mmgDefaultValues(mesh);

#ifndef PATTERN
  fprintf(stdout,"Bucket size per dimension (-bucket) : %d\n",
          mesh->info.bucket);
#endif

  exit(EXIT_FAILURE);
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
int MMG5_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met) {
  int     i;
  char    namein[128];

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

      case 'a':
        if ( !strcmp(argv[i],"-ar") && ++i < argc )
          if ( !MMG5_Set_dparameter(mesh,met,MMG5_DPARAM_angleDetection,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        break;
      case 'A': /* anisotropy */
        if ( !MMG5_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor) )
          exit(EXIT_FAILURE);
        break;
#ifndef PATTERN
      case 'b':
        if ( !strcmp(argv[i],"-bucket") && ++i < argc )
          if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_bucket,
                                    atoi(argv[i])) )
            exit(EXIT_FAILURE);
        break;
#endif
      case 'd':  /* debug */
        if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_debug,1) )
          exit(EXIT_FAILURE);
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !MMG5_Set_dparameter(mesh,met,MMG5_DPARAM_hmin,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !MMG5_Set_dparameter(mesh,met,MMG5_DPARAM_hmax,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc ) {
          if ( !MMG5_Set_dparameter(mesh,met,MMG5_DPARAM_hausd,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc ) {
          if ( !MMG5_Set_dparameter(mesh,met,MMG5_DPARAM_hgrad,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else
          _MMG5_usage(argv[0]);
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMG5_Set_inputMeshName(mesh, argv[i]) )
              exit(EXIT_FAILURE);

            if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_verbose,5) )
              exit(EXIT_FAILURE);
          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            _MMG5_usage(argv[0]);
          }
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-lag") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_lag,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else if ( i == argc ) {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            _MMG5_usage(argv[0]);
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            _MMG5_usage(argv[0]);
            i--;
          }
        }
        else if ( !strcmp(argv[i],"-ls") ) {
          if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_iso,1) )
            exit(EXIT_FAILURE);
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG5_Set_dparameter(mesh,met,MMG5_DPARAM_ls,atof(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else if ( i == argc ) {
            fprintf(stderr,"Missing argument option %c%c\n",argv[i-1][1],argv[i-1][2]);
            _MMG5_usage(argv[0]);
          }
          else i--;
        }
        break;
      case 'm':  /* memory */
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_mem,atoi(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          _MMG5_usage(argv[0]);
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_angle,0) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_noswap,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_noinsert,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_nomove,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-nosurf") ) {
          if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_nosurf,1) )
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
            _MMG5_usage(argv[0]);
          }
        }
        else if( !strcmp(argv[i],"-optim") ) {
          if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_optim,1) )
            exit(EXIT_FAILURE);
        }
        break;
#ifdef USE_SCOTCH
      case 'r':
        if ( !strcmp(argv[i],"-rn") ) {
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {
              if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_renum,atoi(argv[i])) )
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
            if ( !MMG5_Set_inputSolName(mesh,met,argv[i]) )
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
            if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_verbose,atoi(argv[i])) )
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
        if ( !MMG5_Set_inputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
        if ( mesh->info.imprim == -99 ) {
          if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_verbose,5) )
            exit(EXIT_FAILURE);
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMG5_Set_outputMeshName(mesh,argv[i]) )
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
    if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_verbose,i) )
      exit(EXIT_FAILURE);
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
int MMG5_parsop(MMG5_pMesh mesh,MMG5_pSol met) {
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
      if ( !MMG5_Set_iparameter(mesh,met,MMG5_IPARAM_numberOfLocalParam,npar) )
        exit(EXIT_FAILURE);

      for (i=0; i<mesh->info.npar; i++) {
        ret = fscanf(in,"%d %s ",&ref,buf);
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        if ( strcmp(buf,"triangles") && strcmp(buf,"triangle") ) {
          fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
          continue;
        }
        ret = fscanf(in,"%f",&fp1);
        if ( !MMG5_Set_localParameter(mesh,met,MMG5_Triangle,ref,fp1) )
          exit(EXIT_FAILURE);
      }
    }
  }
  fclose(in);
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param *info pointer toward the info structure.
 * \return 1.
 *
 * Store the info structure in the mesh structure.
 *
 */
int _MMG5_stockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(&mesh->info,info,sizeof(MMG5_Info));
  _MMG5_memOption(mesh);
  if( mesh->info.mem > 0) {
    if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne)) {
      return(0);
    } else if(mesh->info.mem < 39)
      return(0);
  }
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param *info pointer toward the info structure.
 *
 * Recover the info structure stored in the mesh structure.
 *
 */
void _MMG5_destockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(info,&mesh->info,sizeof(MMG5_Info));
  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param critmin minimum quality for elements.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param *eltab pointer toward the table of invalid elements.
 *
 * Search invalid elements (in term of quality or edge length).
 *
 */
int MMG5_mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol met,
                    double critmin, double lmin, double lmax, int *eltab) {

  mytime    ctim[TIMEMAX];
  char      stim[32];
  int       ier;

  fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"    %s %s\n",__DATE__,__TIME__);

  signal(SIGABRT,_MMG5_pampa_excfun);
  signal(SIGFPE,_MMG5_pampa_excfun);
  signal(SIGILL,_MMG5_pampa_excfun);
  signal(SIGSEGV,_MMG5_pampa_excfun);
  signal(SIGTERM,_MMG5_pampa_excfun);
  signal(SIGINT,_MMG5_pampa_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  _MMG5_pampa_warnOrientation(mesh);

  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));
    met->np = 0;
  }
  else if ( met->size!=1 ) {
    fprintf(stdout,"  ## ERROR: ANISOTROPIC METRIC NOT IMPLEMENTED.\n");
    return(MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMG5_pampa_setfunc(mesh,met);
  fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);

  if ( !_MMG5_scaleMesh(mesh,met) ) return(MMG5_STRONGFAILURE);
  if ( mesh->info.iso ) {
    if ( !met->np ) {
      fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
      return(MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_mmg3d2(mesh,met) ) return(MMG5_STRONGFAILURE);
  }

  MMG5_searchqua(mesh,met,critmin,eltab);
  ier = MMG5_searchlen(mesh,met,lmin,lmax,eltab);
  if ( !ier )
    return(MMG5_LOWFAILURE);

  return(MMG5_SUCCESS);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param critmin minimum quality for elements.
 * \param *eltab pointer toward the table of invalid elements.
 *
 * Store elements which have worse quality than \a critmin in \a eltab,
 * \a eltab is allocated and could contain \a mesh->ne elements.
 *
 */
void MMG5_searchqua(MMG5_pMesh mesh,MMG5_pSol met,double critmin, int *eltab) {
  MMG5_pTetra   pt;
  double   rap;
  int      k;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    if( !MG_EOK(pt) )
      continue;

    rap = _MMG5_ALPHAD *
      _MMG5_caltet(mesh,met,pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
    if ( rap == 0.0 || rap < critmin ) {
      eltab[k] = 1;
    }
  }
  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param *eltab table of invalid elements.
 *
 * Store in \a eltab elements which have edge lengths shorter than \a lmin
 * or longer than \a lmax, \a eltab is allocated and could contain \a mesh->ne
 * elements.
 *
 */
int MMG5_searchlen(MMG5_pMesh mesh, MMG5_pSol met, double lmin,
                   double lmax, int *eltab) {
  MMG5_pTetra          pt;
  MMG5_pxTetra    pxt;
  _MMG5_Hash           hash;
  double          len;
  int             k,np,nq;
  char            ia,i0,i1,ier;

  /* Hash all edges in the mesh */
  if ( !_MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return(0);

  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = _MMG5_iare[ia][0];
      i1 = _MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      if(!_MMG5_hashEdge(mesh,&hash,np,nq,0)){
        fprintf(stdout,"%s:%d: Error: function _MMG5_hashEdge return 0\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
      }
    }
  }

  /* Pop edges from hash table, and analyze their length */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    for(ia=0; ia<6; ia++) {
      i0 = _MMG5_iare[ia][0];
      i1 = _MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      /* Remove edge from hash ; ier = 1 if edge has been found */
      ier = _MMG5_hashPop(&hash,np,nq);
      if( ier ) {
        if ( pt->xt )
          len = _MMG5_lenedg(mesh,met,np,nq,(pxt->tag[ia] & MG_GEO));
        else
          len = _MMG5_lenedg(mesh,met,np,nq,0);

        if( (len < lmin) || (len > lmax) ) {
          eltab[k] = 1;
          break;
        }
      }
    }
  }
  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  return(1);
}

/**
 * \brief Compute edge length from edge's coordinates.
 * \param *ca pointer toward the coordinates of the first edge's extremity.
 * \param *cb pointer toward the coordinates of the second edge's extremity.
 * \param *ma pointer toward the metric associated to the first edge's extremity.
 * \param *mb pointer toward the metric associated to the second edge's extremity.
 * \return edge length.
 *
 * Compute length of edge \f$[ca,cb]\f$ (with \a ca and \a cb
 * coordinates of edge extremities) according to the isotropic size
 * prescription.
 *
 */
inline double _MMG5_lenedgCoor_iso(double *ca,double *cb,double *ma,double *mb) {
  double   h1,h2,l,r,len;

  h1 = *ma;
  h2 = *mb;
  l = (cb[0]-ca[0])*(cb[0]-ca[0]) + (cb[1]-ca[1])*(cb[1]-ca[1]) \
    + (cb[2]-ca[2])*(cb[2]-ca[2]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;
  len = fabs(r) < _MMG5_EPS ? l / h1 : l / (h2-h1) * log(r+1.0);

  return(len);
}

/**
 * \brief Compute edge length from edge's coordinates.
 * \param *ca pointer toward the coordinates of the first edge's extremity.
 * \param *cb pointer toward the coordinates of the second edge's extremity.
 * \param *ma pointer toward the metric associated to the first edge's extremity.
 * \param *mb pointer toward the metric associated to the second edge's extremity.
 * \return edge length.
 *
 * Compute length of edge \f$[ca,cb]\f$ (with \a ca and \a cb
 * coordinates of edge extremities) according to the anisotropic size
 * prescription.
 *
 */
inline double _MMG5_lenedgCoor_ani(double *ca,double *cb,double *sa,double *sb) {
  fprintf(stdout,"under develop : first thing to do\n");
  return(0.0);
}
