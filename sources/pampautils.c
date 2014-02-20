/** functions for the pampa library */
#include "mmg3d.h"

/* COPY OF PART OF shared_func.c */

/** Warn user that some tetra of the mesh have been reoriented */
static inline
void pampa_warnOrientation(MMG5_pMesh mesh) {
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

static inline
void pampa_excfun(int sigid) {
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

/** set function pointers */
static inline
void pampa_setfunc(pMesh mesh,pSol met) {
  if ( met->size < 6 )
    MMG5_lenedgCoor = lenedgCoor_iso;
  else
    MMG5_lenedgCoor = lenedgCoor_ani;
}
/* END COPY */

/** Return the index of the 4 adjacent elements of tetra kel.
    vi = 0 if the i th face has no adjacent (boundary face)  */
int Get_adjaTet(pMesh mesh, int kel, int *v0, int *v1, int *v2, int *v3) {

  if ( ! mesh->adja ) {
    return(0);
  }

  (*v0) = mesh->adja[4*(kel-1)+1]/4;
  (*v1) = mesh->adja[4*(kel-1)+2]/4;
  (*v2) = mesh->adja[4*(kel-1)+3]/4;
  (*v3) = mesh->adja[4*(kel-1)+4]/4;

  return(1);
}

/** help for mmg3d5 options */
void usage(char *prog) {
  fprintf(stdout,"\nUsage: %s [-v [n]] [opts..] filein [fileout]\n",prog);

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-h      Print this message\n");
  fprintf(stdout,"-v [n]  Tune level of verbosity, [-10..10]\n");
  fprintf(stdout,"-m [n]  Set maximal memory size to n Mbytes\n");
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
  fprintf(stdout,"-ar     val  angle detection\n");
  fprintf(stdout,"-nr          no angle detection\n");
  fprintf(stdout,"-hmin   val  minimal mesh size\n");
  fprintf(stdout,"-hmax   val  maximal mesh size\n");
  fprintf(stdout,"-hausd  val  control Hausdorff distance\n");
  fprintf(stdout,"-hgrad  val  control gradation\n");
  fprintf(stdout,"-ls          levelset meshing \n");
  fprintf(stdout,"-noswap      no edge or face flipping\n");
  fprintf(stdout,"-nomove      no point relocation\n");
  fprintf(stdout,"-noinsert    no point insertion/deletion \n");
#ifndef PATTERN
  fprintf(stdout,"-bucket val Specify the size of bucket per dimension \n");
#endif
  exit(EXIT_FAILURE);
}

/** Store command line arguments */
int parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met
#ifdef SINGUL
           ,MMG5_pSingul sing
#endif
           ) {
  int     i;
  char    namein[128];

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
#ifndef PATTERN
      case 'b':
        if ( !strcmp(argv[i],"-bucket") && ++i < argc )
          if ( !MMG5_Set_iparameters(mesh,met,MMG5_IPARAM_bucket,
                                     atoi(argv[i])) )
            exit(EXIT_FAILURE);
        break;
#endif
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

/** Read local parameters file */
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

/** Store the info structure in mesh */
void stockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(&mesh->info,info,sizeof(MMG5_Info));
  return;
}

/** Recover the info structure stored in mesh */
void destockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(info,&mesh->info,sizeof(MMG5_Info));
  return;
}

/** Search elements which are not correct (quality or edge length) */
int mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol met,
#ifdef SINGUL
               MMG5_pSingul sing,
#endif
               double critmin, double lmin, double lmax, int *eltab) {

  mytime    ctim[TIMEMAX];
  char      stim[32];
  int       ier;
#ifndef SINGUL
  /* sing is not used but must be declared */
  pSingul   sing;
  Singul    singul;
  sing = &singul;
  memset(sing,0,sizeof(Singul));
#endif

  fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"    %s %s\n",__DATE__,__TIME__);

  signal(SIGABRT,pampa_excfun);
  signal(SIGFPE,pampa_excfun);
  signal(SIGILL,pampa_excfun);
  signal(SIGSEGV,pampa_excfun);
  signal(SIGTERM,pampa_excfun);
  signal(SIGINT,pampa_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  pampa_warnOrientation(mesh);

  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    DEL_MEM(mesh,met->m,(met->size*met->npmax+1)*sizeof(double));
    met->np = 0;
  }
  else if ( met->size!=1 ) {
    fprintf(stdout,"  ## ERROR: ANISOTROPIC METRIC NOT IMPLEMENTED.\n");
    return(MMG5_STRONGFAILURE);
  }
#ifdef SINGUL
  if ( mesh->info.sing ) {
    if ( !mesh->info.iso ) {
      if ( !sing->namein )
        fprintf(stdout,"  ## WARNING: NO SINGULARITIES PROVIDED.\n");
    }
    else if ( sing->namein ) {
      fprintf(stdout,"  ## WARNING: SINGULARITIES MUST BE INSERTED IN");
      fprintf(stdout," A PRE-REMESHING PROCESS.\n");
      fprintf(stdout,"              FILE %s IGNORED\n",sing->namein);
    }
  }
#endif

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  pampa_setfunc(mesh,met);
  fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);

  if ( !scaleMesh(mesh,met,sing) ) return(MMG5_STRONGFAILURE);
  if ( mesh->info.iso ) {
    if ( !met->np ) {
      fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
      return(MMG5_STRONGFAILURE);
    }
    if ( !mmg3d2(mesh,met) ) return(MMG5_STRONGFAILURE);
  }

  MMG5_searchqua(mesh,met,critmin,eltab);
  ier = MMG5_searchlen(mesh,met,lmin,lmax,eltab);
  if ( !ier )
    return(MMG5_LOWFAILURE);

  return(MMG5_SUCCESS);
}

/** Store in eltab elements which have worse quality than critmin,
 * eltab is allocated and could contain mesh->ne elements */
void searchqua(MMG5_pMesh mesh,MMG5_pSol met,double critmin, int *eltab) {
  pTetra   pt;
  double   rap;
  int      k;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    if( !MG_EOK(pt) )
      continue;

    rap = ALPHAD * caltet(mesh,met,pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
    if ( rap == 0.0 || rap < critmin ) {
      eltab[k] = 1;
    }
  }
  return;
}

/** Store in eltab elements which have edge lengths shorter than lmin or longer than
 * lmax, eltab is allocated and could contain mesh->ne elements */
int searchlen(MMG5_pMesh mesh, MMG5_pSol met, double lmin, double lmax, int *eltab) {
  pTetra          pt;
  Hash            hash;
  double          len;
  int             k,np,nq;
  char            ia,i0,i1,ier;

  /* Hash all edges in the mesh */
  if ( !hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return(0);

  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = iare[ia][0];
      i1 = iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      if(!hashEdge(mesh,&hash,np,nq,0)){
        fprintf(stdout,"%s:%d: Error: function hashEdge return 0\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
      }
    }
  }

  /* Pop edges from hash table, and analyze their length */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = iare[ia][0];
      i1 = iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      /* Remove edge from hash ; ier = 1 if edge has been found */
      ier = hashPop(&hash,np,nq);
      if( ier ) {
        len = lenedg(mesh,met,np,nq);

        if( (len < lmin) || (len > lmax) ) {
          eltab[k] = 1;
          break;
        }
      }
    }
  }
  return(1);
}

/** Compute length of edge [ca,cb] (with ca and cb coordinates of edge extremities)
    according to the size prescription */
inline double lenedgCoor_iso(double *ca,double *cb,double *ma,double *mb) {
  double   h1,h2,l,r,len;

  h1 = *ma;
  h2 = *mb;
  l = (cb[0]-ca[0])*(cb[0]-ca[0]) + (cb[1]-ca[1])*(cb[1]-ca[1]) \
    + (cb[2]-ca[2])*(cb[2]-ca[2]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;
  len = fabs(r) < EPS ? l / h1 : l / (h2-h1) * log(r+1.0);

  return(len);
}

inline double lenedgCoor_ani(double *ca,double *cb,double *sa,double *sb) {
  return(0.0);
}
