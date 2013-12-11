/* global variables */
Info   info;
unsigned char inxt2[3] = {1,2,0};
unsigned char iprv2[3] = {2,0,1};
unsigned char idir[4][3] = { {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1} };
char idirinv[4][4] = {{-1,0,1,2},{0,-1,2,1},{0,1,-1,2},{0,2,1,-1}};
unsigned char iarf[4][3] = { {5,4,3}, {5,1,2}, {4,2,0}, {3,0,1} };
unsigned char iarfinv[4][6] = { {-1,-1,-1,2,1,0}, {-1,1,2,-1,-1,0},{2,-1,1,-1,0,-1},{1,2,-1,0,-1,-1}};
unsigned char inxt3[7] = { 1,2,3,0,1,2,3 };
unsigned char iprv3[7] = { 3,0,1,2,3,0,1 };
unsigned char iare[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
unsigned char ifar[6][2] = { {2,3}, {1,3}, {1,2}, {0,3}, {0,2}, {0,1} };
unsigned char isar[6][2] = { {2,3}, {3,1}, {1,2}, {0,3}, {2,0}, {0,1} };
unsigned char arpt[4][3] = { {0,1,2}, {0,4,3}, {1,3,5}, {2,5,4} };

/* shared functions */

#ifdef USE_SCOTCH
/** Warn user that we overflow asked memory during scotch call */
static inline
void warnScotch(int mem) {
  if ( info.imprim > 4 || info.ddebug ) {
    if ( info.mem >= 0 ) {
      fprintf(stdout,"  ## Warning: we will overflow the memory asked with \"-m\"");
      fprintf(stdout," option during Scotch call.\n" );
    }
  }
}
#endif

static inline
void excfun(int sigid) {
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
void setfunc(pMesh mesh,pSol met) {
  if ( met->size < 6 ) {
    caltet = caltet_iso;
    lenedg = lenedg_iso;
    defsiz = defsiz_iso;
    gradsiz = gradsiz_iso;
  }
  else {
    caltet = caltet_ani;
    lenedg = lenedg_ani;
    /*defsiz = defsiz_ani;
      gradsiz = gradsiz_ani;*/
  }
}

/** Read parammeter file */
static inline
int parsop(pMesh mesh,pSol met) {
  Par        *par;
  float       fp1;
  int         i,j,ret;
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
      fscanf(in,"%d",&info.npar);
      info.par = (Par*)calloc(info.npar,sizeof(Par));
      if ( !info.par ) {
        perror("  ## Memory problem: calloc");
        exit(EXIT_FAILURE);
      }

      for (i=0; i<info.npar; i++) {
        par = &info.par[i];
        fscanf(in,"%d %s ",&par->ref,buf);
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )
          par->elt = MMG5_Triangle;
        else {
          fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
          continue;
        }
        ret = fscanf(in,"%f",&fp1);
        par->hausd = fp1;
      }
    }
  }
  fclose(in);
  return(1);
}
