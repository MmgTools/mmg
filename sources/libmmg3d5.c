/**
 *
 * Written by Cecile Dobrzynski (IMB), Charles Dapogny and Pascal Frey (LJLL)
 * Copyright (c) 2004- IMB/LJLL.
 * All rights reserved.
 *
 * mmg3dlib(int options_i[10],double options_d[5]  ): to use mmg3d via a library
 *
 * option_i:
 *    option_i[   MMG5_imprim] = [-10..10] , Tune level of verbosity;
 *    option_i[      MMG5_mem] = [n/-1]    , Set memory size to n Mbytes/keep the default value;
 *    option_i[    MMG5_debug] = [1/0]     , Turn on/off debug mode;
 *    option_i[    MMG5_angle] = [1/0]     , Turn on/off angle detection;
 *    option_i[      MMG5_iso] = [1/0]     , Turn on/off levelset meshing;
 *    option_i[ MMG5_noinsert] = [1/0]     , avoid/allow point insertion/deletion;
 *    option_i[   MMG5_noswap] = [1/0]     , avoid/allow edge or face flipping;
 *    option_i[   MMG5_nomove] = [1/0]     , avoid/allow point relocation;
 *    option_i[    MMG5_renum] = [1/0]     , Turn on/off the renumbering using SCOTCH;
 *    option_i[    MMG5_sing ] = [1/0]     , Turn on/off the insertion of singularities
 *                                           (need to compile with -DSINGUL flag);
 *
 *    option_d[  MMG5_dhd] = [val]     , angle detection;
 *    option_d[ MMG5_hmin] = [val]     , minimal mesh size;
 *    option_d[ MMG5_hmax] = [val]     , maximal mesh size;
 *    option_d[MMG5_hausd] = [val]     , control Hausdorff distance;
 *    option_d[MMG5_hgrad] = [val]     , control gradation;
 *    option_d[   MMG5_ls] = [val]     , level set value;
 **/

#include "compil.date"
#include "mmg3d.h"

/* global */
/** Warning: if the library is run on multithread, global variables may be overwritten */
Info info;

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


static inline void excfun(int sigid) {
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

/** Deallocations before return */
void freeAll(pMesh mesh,pSol met
#ifdef SINGUL
             ,pSingul singul
#endif
             ){

  /* mesh */
  free(mesh->point);
  mesh->point = NULL;
  free(mesh->tetra);
  mesh->tetra = NULL;
  free(mesh->adja);
  mesh->adja = NULL;
  if ( mesh->xpoint ) {
    free(mesh->xpoint);
    mesh->xpoint = NULL;
  }
  if ( mesh->htab.geom ) {
    free(mesh->htab.geom);
    mesh->htab.geom = NULL;
  }
  if ( mesh->edge ) {
    free(mesh->edge);
    mesh->edge = NULL;
  }
  if ( mesh->tria ) {
    free(mesh->tria);
    mesh->tria = NULL;
  }
  if ( mesh->xtetra ) {
    free(mesh->xtetra);
    mesh->xtetra = NULL;
  }

  /* met */
  if ( !info.iso && met->m ) {
    free(met->m);
    met->m = NULL;
  }
#ifdef SINGUL
  /* singul */
  if ( info.sing ) {
    if ( singul->point ) {
      free(singul->point);
      singul->point=NULL;
    }
    if ( singul->edge ) {
      free(singul->edge);
      singul->edge=NULL;
    }
  }
#endif
}

/** Recover mesh data */
static inline
int inputdata(pMesh mesh,pSol met) {
  pPoint  ppt;
  int	  	k,i;

  /* Fill dimension and version data if needed */
  if ( !mesh->dim )  mesh->dim = 3;
  else if ( mesh->dim != 3 ) {
      fprintf(stdout,"  ** 3 DIMENSIONAL MESH NEEDED. Exit program.\n");
      return(0);
  }

  if ( !met->dim )  met->dim = 3;
  else if ( met->dim != 3 ) {
      fprintf(stdout,"  ** WRONG DIMENSION FOR METRIC. Exit program.\n");
      return(0);
  }
  if ( !mesh->ver )  mesh->ver = 2;
  if ( !met ->ver )  met ->ver = 2;

  /*  Check mesh data */
  if ( info.ddebug ) {
    if ( (!mesh->np) || (!mesh->point) ||
         (!mesh->ne) || (!mesh->tetra) ) {
      fprintf(stdout,"  ** MISSING DATA. Exit program.\n");
      return(0);
    }
  }
  mesh->base = mesh->mark = 0;

  mesh->npi   = mesh->np;
  mesh->nei   = mesh->ne;
  mesh->nai   = mesh->na;

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->ne + 1;
  for (k=mesh->npnil; k<mesh->npmax-1; k++) {
    mesh->point[k].tmp  = k+1;
  }
  for (k=mesh->nenil; k<mesh->nemax-1; k++) {
    mesh->tetra[k].v[3] = k+1;
  }

  /* tag points*/
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->tag  = MG_NUL;
    ppt->flag = 0;
  }
  for (k=1; k<=mesh->ne; k++) {
    for (i=0; i<4; i++) {
      ppt = &mesh->point[mesh->tetra[k].v[i]];
      ppt->tag &= ~MG_NUL;
    }
  }

  return(1);
}

static inline
int packMesh(pMesh mesh,pSol met) {
  pTetra	pt,ptnew;
  pPoint	ppt,pptnew;
  hgeom   *ph;
  int     np,nc,nr, k,ne,nbl,imet,imetnew,i;
  int     iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;

  /* compact vertices */
  np = nc = nr = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->tmp = ++np;
    if ( ppt->tag & MG_CRN )  nc++;
  }

  /* compact tetrahedra */
  ne  = 0;
  nbl = 1;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    pt->v[0] = mesh->point[pt->v[0]].tmp;
    pt->v[1] = mesh->point[pt->v[1]].tmp;
    pt->v[2] = mesh->point[pt->v[2]].tmp;
    pt->v[3] = mesh->point[pt->v[3]].tmp;
    ne++;
    if ( k!=nbl ) {
      ptnew = &mesh->tetra[nbl];
      memcpy(ptnew,pt,sizeof(Tetra));

      iadr = 4*(k-1) + 1;
      adja = &mesh->adja[iadr];
      iadrnew = 4*(nbl-1) + 1;
      adjanew = &mesh->adja[iadrnew];
      for(i=0 ; i<4 ; i++) {
        adjanew[i] = adja[i];
        if(!adja[i]) continue;
        iadrv = 4*(adja[i]/4-1) +1;
        adjav = &mesh->adja[iadrv];
        voy = i;
        adjav[adja[i]%4] = 4*nbl + voy;
      }
    }
    nbl++;
  }
  mesh->ne = ne;

  /* compact metric */
  nbl = 1;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    imet    = (k-1) * met->size + 1;
    imetnew = (nbl-1) * met->size + 1;

    for (i=0; i<met->size; i++)
      met->m[imetnew + i] = met->m[imet + i];
    ++nbl;
  }

  /*compact vertices*/
  np  = 0;
  nbl = 1;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    np++;
    if ( k!=nbl ) {
      pptnew = &mesh->point[nbl];
      memmove(pptnew,ppt,sizeof(Point));
      memset(ppt,0,sizeof(Point));
      ppt->tag    = MG_NUL;
    }
    nbl++;
  }
  mesh->np = np;
  met->np  = np;

  /* rebuild triangles*/
  mesh->nt = 0;
  if ( !bdryTria(mesh) ) {
    fprintf(stdout," ## Error: unable to rebuild triangles\n");
    return(0);
  }

  /* build hash table for edges */
  if ( mesh->htab.geom ) {
    free(mesh->htab.geom);
    mesh->htab.geom=NULL;
  }
  hNew(&mesh->htab,3*(mesh->xt),9*(mesh->xt));
  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    if ( MG_EOK(pt) &&  pt->xt ) {
      for (i=0; i<6; i++) {
        if ( mesh->xtetra[pt->xt].edg[i] ||
             ( MG_EDG(mesh->xtetra[pt->xt].tag[i] ) ||
               (mesh->xtetra[pt->xt].tag[i] & MG_REQ) ) )
          hEdge(&mesh->htab,pt->v[iare[i][0]],pt->v[iare[i][1]],
                mesh->xtetra[pt->xt].edg[i],mesh->xtetra[pt->xt].tag[i]);
      }
    }
  }

  /* edges + ridges + required edges */
  mesh->na = 0;
  for (k=0; k<=mesh->htab.max; k++) {
    ph = &mesh->htab.geom[k];
    if ( !(ph->a) )  continue;
    mesh->na++;
  }
  if ( mesh->na ) {
    mesh->edge = (pEdge)calloc(mesh->na+1,sizeof(Edge));
    if ( !mesh->edge ) {
      fprintf(stdout,"  ## Allocation problem (edge), not enough memory.\n");
      fprintf(stdout,"  ## Check the mesh size or increase the allocated memory");
      fprintf(stdout," with the -m option.\n");
      fprintf(stdout,"  Exit program.\n");
      return(0);
    }
    mesh->na = 0;
    for (k=0; k<=mesh->htab.max; k++) {
      ph = &mesh->htab.geom[k];
      if ( !ph->a )  continue;
      mesh->na++;
      mesh->edge[mesh->na ].a  = mesh->point[ph->a].tmp;
      mesh->edge[mesh->na ].b  = mesh->point[ph->b].tmp;
      mesh->edge[mesh->na].tag = ( ph->tag | MG_REF ) ;
      mesh->edge[mesh->na].ref = ph->ref;
      if ( MG_GEO & ph->tag ) nr++;
    }
  }

  for(k=1 ; k<=mesh->np ; k++)
    mesh->point[k].tmp = 0;

  mesh->npnil = mesh->np + 1;
  for (k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  mesh->nenil = mesh->ne + 1;
  for (k=mesh->nenil; k<mesh->nemax-1; k++)
    mesh->tetra[k].v[3] = k+1;

  /* to could save the mesh, the adjacency have to be correct */
  if ( info.ddebug && (!chkmsh(mesh,1,1) ) ) {
    fprintf(stdout,"  ##  Problem. Invalid mesh.\n");
    return(0);
  }
  if ( info.imprim ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",mesh->np,nc);
    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d   RIDGES  %8d\n",mesh->na,nr);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",mesh->nt);
    fprintf(stdout,"     NUMBER OF ELEMENTS   %8d\n",mesh->ne);
  }
  return(1);
}

/** Initialization of option tables with default values */
void mmg3dinit(int opt_i[10],double opt_d[6]) {

  /* default values for first tab (integer) */
  opt_i[   MMG5_imprim] = -99; /**< [-10..10],Tune level of imprim */
  opt_i[      MMG5_mem] = -1;  /**< [n/-1]   ,Set memory size to n Mbytes/keep the default value */
  opt_i[    MMG5_debug] =  0;  /**< [1/0]    ,Turn on/off debug mode */
  opt_i[    MMG5_angle] =  1;  /**< [1/0]    ,Turn on/off angle detection */
  opt_i[      MMG5_iso] =  0;  /**< [1/0]    ,Turn on/off levelset meshing */
  opt_i[ MMG5_noinsert] =  0;  /**< [1/0]    ,avoid/allow point insertion/deletion */
  opt_i[   MMG5_noswap] =  0;  /**< [1/0]    ,avoid/allow edge or face flipping */
  opt_i[   MMG5_nomove] =  0;  /**< [1/0]    ,avoid/allow point relocation */
#ifdef USE_SCOTCH
  opt_i[    MMG5_renum] = 1;   /**< [1/0]    , Turn on/off the renumbering using SCOTCH; */
#else
  opt_i[    MMG5_renum] = 0;   /**< [1/0]    , Turn on/off the renumbering using SCOTCH; */
#endif

  /* default values for second tab (double) */
  opt_d[  MMG5_dhd] = 45;       /**< angle detection; */
  opt_d[ MMG5_hmin] = 0.0;      /**< minimal mesh size; */
  opt_d[ MMG5_hmax] = FLT_MAX;  /**< maximal mesh size; */
  opt_d[MMG5_hausd] = 0.01;     /**< control Hausdorff */
  opt_d[MMG5_hgrad] = exp(0.1); /**< control gradation; */
  opt_d[   MMG5_ls] = 0.0;      /**< level set value */
}

/** Store user options in the info structure */
void stockOption(int opt_i[10],double opt_d[6],pMesh mesh){

  /* recovering of first option table (integers) */
  info.imprim   = opt_i[MMG5_imprim];
  info.mem      = opt_i[MMG5_mem];
  info.ddebug   = opt_i[MMG5_debug];
  if ( !opt_i[MMG5_angle] )
    info.dhd = -1.0;
  else {
    info.dhd = opt_d[MMG5_dhd];
    info.dhd = MG_MAX(0.0, MG_MIN(180.0,info.dhd));
    info.dhd = cos(info.dhd*M_PI/180.0);
  }

  info.iso      = opt_i[MMG5_iso];
  info.noinsert = opt_i[MMG5_noinsert];
  info.noswap   = opt_i[MMG5_noswap];
  info.nomove   = opt_i[MMG5_nomove];
#ifdef USE_SCOTCH
  info.renum    = opt_i[MMG5_renum];
#else
  info.renum    = 0;
#endif
#ifdef SINGUL
  info.sing     = opt_i[MMG5_sing];
#else
  info.sing     = 0;
#endif

  /* recovering of second option table (doubles) */
  info.hmin     = opt_d[MMG5_hmin];
  info.hmax     = opt_d[MMG5_hmax];
  info.hausd    = opt_d[MMG5_hausd];
  info.hgrad    = opt_d[MMG5_hgrad];
  if ( info.hgrad < 0.0 )
    info.hgrad = -1.0;
  else
    info.hgrad = log(info.hgrad);

  info.ls       = opt_d[MMG5_ls];

  /* other options */
  info.fem      = 0;

}

int mmg3dlib(int opt_i[10],double opt_d[6],pMesh mesh,pSol met
#ifdef SINGUL
             ,pSingul sing
#endif
             ) {

  char      stim[32];
  int       ier;
#ifndef SINGUL
  pSingul   sing;
#endif

  fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"    %s\n",COMPIL);

  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);

  tminit(info.ctim,TIMEMAX);
  chrono(ON,&(info.ctim[0]));

  stockOption(opt_i,opt_d,mesh);

  /* load data */
  fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");
  chrono(ON,&(info.ctim[1]));
  /* input data */
  if ( !inputdata(mesh,met) ) return(MMG5_STRONGFAILURE);

  met->npmax = mesh->npmax;
  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    free(met->m);
    met->m   = NULL;
    met->np = 0;
  }
  else if ( met->size!=1 ) {
    fprintf(stdout,"  ## ERROR: ANISOTROPIC METRIC NOT IMPLEMENTED.\n");
    return(MMG5_STRONGFAILURE);
  }
#ifdef SINGUL
  if ( info.sing ) {
    if ( !info.iso ) {
      if ( !sing->namein )
        fprintf(stdout,"  ## WARNING: NO SINGULARITIES PROVIDED.\n");
      else
        if ( !loadSingul(sing) ) return(MMG5_STRONGFAILURE);
    }
    else if ( sing->namein ) {
      fprintf(stdout,"  ## WARNING: SINGULARITIES MUST BE INSERTED IN");
      fprintf(stdout," A PRE-REMESHING PROCESS.\n");
      fprintf(stdout,"              FILE %s IGNORED\n",sing->namein);
    }
  }
#endif

  chrono(OFF,&(info.ctim[1]));
  printim(info.ctim[1].gdif,stim);
  fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(info.ctim[2]));
  setfunc(mesh,met);
  if ( abs(info.imprim) > 0 )  outqua(mesh,met);
  fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
  if ( info.imprim )  fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");

  if ( !scaleMesh(mesh,met,sing) ) return(MMG5_STRONGFAILURE);
  if ( info.iso ) {
    if ( !met->np ) {
      fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
      return(MMG5_STRONGFAILURE);
    }
    if ( !mmg3d2(mesh,met) ) return(MMG5_STRONGFAILURE);
  }

#ifdef SINGUL
  if ( info.sing ) {
    if ( !info.iso ) {
      if ( !met->np && !DoSol(mesh,met,&info) )
        return(MMG5_LOWFAILURE);
      if ( !( ier=inserSingul(mesh,met,sing) ) )
        return(MMG5_STRONGFAILURE);
      else if (ier > 0 ) {
        chrono(OFF,&info.ctim[2]);
        printim(info.ctim[2].gdif,stim);
        fprintf(stdout,"  -- INSERTION OF SINGULARITIES COMPLETED.     %s\n\n",stim);
        chrono(ON,&info.ctim[2]);
      }
    }
  }
#endif

#ifdef DEBUG
  if ( !met->np && !DoSol(mesh,met,&info) ) {
  if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
    return(mesh,met,MMG5_LOWFAILURE);
  }
#endif
  if ( !analys(mesh) ) {
  if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
    return(MMG5_LOWFAILURE);
  }

  if ( info.imprim > 4 && !info.iso && met->m ) prilen(mesh,met);

  chrono(OFF,&(info.ctim[2]));
  printim(info.ctim[2].gdif,stim);
  if ( info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(info.ctim[3]));
  if ( info.imprim )
    fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");

#ifdef SINGUL
  if ( info.sing && (!info.iso) ) {
    if ( colSing(mesh,met)<0 ) {
      fprintf(stdout,"  ## Collapse of singularities problem.\n");
      // return(MMG5_STRONGFAILURE);
    }
  }
#endif

  if ( !mmg3d1(mesh,met) ){
    if ( !(mesh->adja) && !hashTetra(mesh) ) {
      fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
      return(MMG5_STRONGFAILURE);
    }
    if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
    return(MMG5_LOWFAILURE);
  }

#ifdef SINGUL
  if ( info.sing && (!info.iso) ) {
    if ( !solveUnsignedTet(mesh,met) ) {
      fprintf(stdout,"  ## Solve of undetermined tetrahedra problem.\n");
      if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
      return(MMG5_LOWFAILURE);
    }
  }
#endif

  chrono(OFF,&(info.ctim[3]));
  printim(info.ctim[3].gdif,stim);
  if ( info.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  fprintf(stdout,"\n  %s\n   END OF MODULE MMG3d: IMB-LJLL \n  %s\n",MG_STR,MG_STR);

  /* save file */
  outqua(mesh,met);
  if ( info.imprim > 4 && !info.iso )
    prilen(mesh,met);

  chrono(ON,&(info.ctim[1]));
  if ( info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");
  if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
  if ( !packMesh(mesh,met) )          return(MMG5_STRONGFAILURE);
  met->np = mesh->np;
  chrono(OFF,&(info.ctim[1]));

  chrono(OFF,&info.ctim[0]);
  printim(info.ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3DLIB: ELAPSED TIME  %s\n",stim);
  return(MMG5_SUCCESS);
}
