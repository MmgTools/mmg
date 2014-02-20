 /**
 *
 * Written by Cecile Dobrzynski (IMB), Charles Dapogny,
 * Pascal Frey (LJLL) and Algiane Froehly
 * Copyright (c) 2004- IMB/LJLL.
 * All rights reserved.
 *
 * MMG5_mmg3dlib(MMG5_pMesh mesh,MMG5_pSol met [,MMG5_pSingul singul] ):
 *    to use mmg3d via a library
 *
 * Integers parameters:
 *    MMG5_IPARAM_verbose            = [-10..10] , Tune level of verbosity;
 *    MMG5_IPARAM_mem                = [n/-1]    , Set maximal memory size to n Mbytes/keep the default value;
 *    MMG5_IPARAM_debug              = [1/0]     , Turn on/off debug mode;
 *    MMG5_IPARAM_angle              = [1/0]     , Turn on/off angle detection;
 *    MMG5_IPARAM_iso                = [1/0]     , Turn on/off levelset meshing;
 *    MMG5_IPARAM_noinsert           = [1/0]     , avoid/allow point insertion/deletion;
 *    MMG5_IPARAM_noswap             = [1/0]     , avoid/allow edge or face flipping;
 *    MMG5_IPARAM_nomove             = [1/0]     , avoid/allow point relocation;
 *    MMG5_IPARAM_numberOflocalParam = [n]       , number of local parameters;
 *    MMG5_IPARAM_renum              = [1/0]     , Turn on/off the renumbering using SCOTCH;
 *    MMG5_IPARAM_sing               = [1/0]     , Turn on/off the insertion of singularities
 *                                        (need to compile with -DSINGUL flag);
 * Double parameters:
 *    MMG5_DPARAM_dhd   = [val]     , angle detection;
 *    MMG5_DPARAM_hmin  = [val]     , minimal mesh size;
 *    MMG5_DPARAM_hmax  = [val]     , maximal mesh size;
 *    MMG5_DPARAM_hausd = [val]     , control global Hausdorff distance
 *                                    (on all the boundary surfaces of the mesh);
 *    MMG5_DPARAM_hgrad = [val]     , control gradation;
 *    MMG5_DPARAM_ls    = [val]     , level set value;
 **/

#include "mmg3d.h"
#include "shared_func.h"

#define RETURN_AND_PACK(mesh,met,val)do    \
    {                                           \
      packMesh(mesh,met);                       \
      return(val);                              \
    }while(0)

/** Deallocations before return */
void Free_all(pMesh mesh,pSol met
#ifdef SINGUL
             ,pSingul singul
#endif
             ){

#ifdef SINGUL
  Free_structures(mesh,met,singul);
#else
  Free_structures(mesh,met);
#endif

#ifdef SINGUL
  SAFE_FREE(singul);
#endif
  SAFE_FREE(met);
  SAFE_FREE(mesh);
}

/** set pointer for MMG5_saveMesh function */
void Set_saveFunc(pMesh mesh) {
  MMG5_saveMesh = saveLibraryMesh;
}

/** Free adja, xtetra and xpoint tables */
static inline
void Free_topoTables(pMesh mesh) {
  int k;

  mesh->xp = 0;
  if ( mesh->adja )
    DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

  freeXTets(mesh);

  if ( mesh->xpoint )
    DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(xPoint));
  for(k=1; k <=mesh->np; k++) {
    mesh->point[k].xp = 0;
  }

  return;
}

/** pack the sparse mesh and create triangles and edges before getting
    out of library */
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
  if ( met->m ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      imet    = (k-1) * met->size + 1;
      imetnew = (nbl-1) * met->size + 1;

      for (i=0; i<met->size; i++)
        met->m[imetnew + i] = met->m[imet + i];
      ++nbl;
    }
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
  chkNumberOfTri(mesh);
  if ( !bdryTria(mesh) ) {
    fprintf(stdout," ## Error: unable to rebuild triangles\n");
    return(0);
  }

  /* build hash table for edges */
  if ( mesh->htab.geom )
    DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(hgeom));

  mesh->na = 0;
  /* in the wost case (all edges are marked), we will have around 1 edge per *
   * triangle (we count edges only one time) */
  mesh->memCur += (long long)((3*mesh->nt+2)*sizeof(hgeom));
  if ( (mesh->memCur) > (mesh->memMax) ) {
    mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(hgeom));
    fprintf(stdout,"  ## Warning:");
    fprintf(stdout," unable to allocate htab.\n");
  }
  else if ( hNew(&mesh->htab,mesh->nt,3*(mesh->nt),0) ) {
    for (k=1; k<=mesh->ne; k++) {
      pt   = &mesh->tetra[k];
      if ( MG_EOK(pt) &&  pt->xt ) {
        for (i=0; i<6; i++) {
          if ( mesh->xtetra[pt->xt].edg[i] ||
               ( MG_EDG(mesh->xtetra[pt->xt].tag[i] ) ||
                 (mesh->xtetra[pt->xt].tag[i] & MG_REQ) ) )
            hEdge(mesh,pt->v[iare[i][0]],pt->v[iare[i][1]],
                  mesh->xtetra[pt->xt].edg[i],mesh->xtetra[pt->xt].tag[i]);
        }
      }
    }

    /* edges + ridges + required edges */
    for (k=0; k<=mesh->htab.max; k++) {
      ph = &mesh->htab.geom[k];
      if ( !(ph->a) )  continue;
      mesh->na++;
    }
    if ( mesh->na ) {
      ADD_MEM(mesh,(mesh->na+1)*sizeof(Edge),"edges",);
      SAFE_CALLOC(mesh->edge,mesh->na+1,Edge);

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
    DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(hgeom));
  }
  else
    mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(hgeom));

  for(k=1 ; k<=mesh->np ; k++)
    mesh->point[k].tmp = 0;

  mesh->npnil = mesh->np + 1;
  for(k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  mesh->nenil = mesh->ne + 1;
  for(k=mesh->nenil; k<mesh->nemax-1; k++)
    mesh->tetra[k].v[3] = k+1;

  /* to could save the mesh, the adjacency have to be correct */
  if ( mesh->info.ddebug && (!chkmsh(mesh,1,1) ) ) {
    fprintf(stdout,"  ##  Problem. Invalid mesh.\n");
    return(0);
  }

  Free_topoTables(mesh);

  if ( mesh->info.imprim ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",mesh->np,nc);
    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d   RIDGES  %8d\n",mesh->na,nr);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",mesh->nt);
    fprintf(stdout,"     NUMBER OF ELEMENTS   %8d\n",mesh->ne);
  }
  return(1);
}

/** main programm */
int mmg3dlib(pMesh mesh,pSol met
#ifdef SINGUL
             ,pSingul sing
#endif
             ) {
  mytime    ctim[TIMEMAX];
  char      stim[32];
#ifdef SINGUL
  int       ier;
#else
  /* sing is not used but must be declared */
  pSingul   sing;
  Singul    singul;
  sing = &singul;
  memset(sing,0,sizeof(Singul));
#endif

  fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

#ifdef USE_SCOTCH
  warnScotch(mesh);
#endif

  fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  warnOrientation(mesh);

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
  setfunc(mesh,met);
  Set_saveFunc(mesh);
  if ( abs(mesh->info.imprim) > 0 )  outqua(mesh,met);
  fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");

  if ( !scaleMesh(mesh,met,sing) ) return(MMG5_STRONGFAILURE);
  if ( mesh->info.iso ) {
    if ( !met->np ) {
      fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
      return(MMG5_STRONGFAILURE);
    }
    if ( !mmg3d2(mesh,met) ) return(MMG5_STRONGFAILURE);
  }

#ifdef SINGUL
  if ( mesh->info.sing ) {
    if ( !mesh->info.iso ) {
      if ( !met->np && !DoSol(mesh,met) )
        RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
      if ( !( ier=inserSingul(mesh,met,sing) ) )
        return(MMG5_STRONGFAILURE);
      else if (ier > 0 ) {
        chrono(OFF,&ctim[2]);
        printim(ctim[2].gdif,stim);
        fprintf(stdout,"  -- INSERTION OF SINGULARITIES COMPLETED.     %s\n\n",stim);
        chrono(ON,&ctim[2]);
      }
    }
  }
#endif

#ifdef DEBUG
  if ( !met->np && !DoSol(mesh,met,&mesh->info) ) {
    if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
    RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
  }
#endif
  if ( !analys(mesh) ) {
  if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
  RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 4 && !mesh->info.iso && met->m ) prilen(mesh,met);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");

#ifdef USE_SCOTCH
    /*check enough vertex to renum*/
    if ( mesh->info.renum && (mesh->np/2. > BOXSIZE) && mesh->np>100000 ) {
      /* renumbering begin */
      if ( mesh->info.imprim > 5 )
        fprintf(stdout,"  -- RENUMBERING. \n");
      if ( !renumbering(BOXSIZE,mesh, met) ) {
        fprintf(stdout,"  ## Unable to renumbering mesh. \n");
        fprintf(stdout,"  ## Try to run without renumbering option (-rn 0)\n");
        return(0);
      }

      if ( mesh->info.imprim > 5) {
        fprintf(stdout,"  -- PHASE RENUMBERING COMPLETED. \n");
      }

      if ( mesh->info.ddebug )  chkmsh(mesh,1,0);
      /* renumbering end */
    }
#endif

#ifdef SINGUL
  if ( mesh->info.sing && (!mesh->info.iso) ) {
    if ( colSing(mesh,met)<0 ) {
      fprintf(stdout,"  ## Collapse of singularities problem.\n");
      // return(MMG5_STRONGFAILURE);
    }
  }
#endif


#ifdef PATTERN
  if ( !mmg3d1(mesh,met) ) {
    if ( !(mesh->adja) && !hashTetra(mesh,1) ) {
      fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
      return(MMG5_STRONGFAILURE);
    }
    if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
    RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
  }
#else
  /** Patterns in iso mode, delauney otherwise */
  if ( !mesh->info.iso ) {
    if ( !mmg3d1_delone(mesh,met) ) {
      if ( !(mesh->adja) && !hashTetra(mesh,1) ) {
        fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
        return(MMG5_STRONGFAILURE);
      }
      if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
      RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
    }
  }
  else {
    if ( !mmg3d1(mesh,met) ) {
      if ( !(mesh->adja) && !hashTetra(mesh,1) ) {
        fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
        return(MMG5_STRONGFAILURE);
      }
      if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
      RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
    }
  }

#endif

#ifdef SINGUL
  if ( mesh->info.sing && (!mesh->info.iso) ) {
    if ( !solveUnsignedTet(mesh,met) ) {
      fprintf(stdout,"  ## Solve of undetermined tetrahedra problem.\n");
      if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
      RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
    }
  }
#endif

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  fprintf(stdout,"\n  %s\n   END OF MODULE MMG3d: IMB-LJLL \n  %s\n",MG_STR,MG_STR);

  /* save file */
  outqua(mesh,met);
  if ( mesh->info.imprim > 4 && !mesh->info.iso )
    prilen(mesh,met);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");
  if ( !unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
  if ( !packMesh(mesh,met) )     return(MMG5_STRONGFAILURE);
  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3DLIB: ELAPSED TIME  %s\n",stim);
  return(MMG5_SUCCESS);
}
