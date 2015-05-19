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

static void endcod() {
  double   ttot,ttim[TIMEMAX];
  int      k,call[TIMEMAX];

  /* chrono(OFF,&ctim[0]); */
  /* for (k=0; k<TIMEMAX; k++) { */
  /*   call[k] = ctim[k].call; */
  /*   ttim[k] = ctim[k].call ? gttime(ctim[k]) : 0.0; */
  /* } */
  /* ttot    = ttim[1]+ttim[2]+ttim[3]+ttim[4]; */
  /* ttim[0] = M_MAX(ttim[0],ttot); */

  /* fprintf(stdout,"\n  -- CPU REQUIREMENTS\n"); */
  /* fprintf(stdout,"  in/out %8.2f %%    %3d. calls,   %7.2f sec/call\n", */
  /*         100.*ttim[1]/ttim[0],call[1],ttim[1]/(float)call[1]); */
  /* fprintf(stdout,"  analysis %8.2f %%    %3d. calls,   %7.2f sec/call\n", */
  /*         100.*ttim[2]/ttim[0],call[2],ttim[2]/(float)call[2]); */
  /* fprintf(stdout,"  remeshing  %8.2f %%    %3d. calls,   %7.2f sec/call\n", */
  /*         100.*ttim[3]/ttim[0],call[3],ttim[3]/(float)call[3]); */
  /* fprintf(stdout,"  total     %8.2f %%    %3d. calls,  %7.2f sec/call\n", */
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

int MMG2_inputdata(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pPoint   ppt;
  int   k;



  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nanil = mesh->na + 1;
  mesh->nenil = mesh->nt + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;
  for (k=mesh->nenil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;
  if ( mesh->na ) {
    mesh->nanil = mesh->na + 1;
    for (k=mesh->nanil; k<mesh->namax-1; k++)
      mesh->edge[k].b = k+1;
  }
  /*tag points*/
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->tag  = M_NUL;
  }
  return(1);
}

int MMG2_tassage(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pEdge         ped;
  MMG5_pTria    pt,ptnew;
  MMG5_pPoint ppt,pptnew;
  int     np,nt,k,nbl,isol,isolnew,i;
  int     iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;


  /* compact vertices */
  np=0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( ppt->tag & M_NUL )  continue;
    ppt->tmp = ++np;
  }


  /* compact edges */
  for (k=1; k<=mesh->na; k++) {
    ped  = &mesh->edge[k];
    ped->a = mesh->point[ped->a].tmp;
    ped->b = mesh->point[ped->b].tmp;
  }

  /* compact triangle */
  nt  = 0;
  nbl = 1;
  printf("on a %d triangle \n",mesh->nt);
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] )  {
      continue;
    }

    pt->v[0] = mesh->point[pt->v[0]].tmp;
    pt->v[1] = mesh->point[pt->v[1]].tmp;
    pt->v[2] = mesh->point[pt->v[2]].tmp;
    nt++;
    if(k!=nbl) {
      //printf("on voudrait.tmpser\n");
      ptnew = &mesh->tria[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tria));

      //and the adjacency
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      iadrnew = 3*(nbl-1) + 1;
      adjanew = &mesh->adja[iadrnew];
      for(i=0 ; i<2 ; i++) {
        adjanew[i] = adja[i];
        if(!adja[i]) continue;
        iadrv = 3*(adja[i]/3-1) +1;
        adjav = &mesh->adja[iadrv];
        voy = i;
        adjav[adja[i]%3] = 3*nbl + voy;
      }
    }
    nbl++;
  }
  mesh->nt = nt;
  printf("apres.tmpsage on a %d triangle \n",mesh->nt);

  /* compact metric */
  nbl = 1;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( ppt->tag & M_NUL )  continue;
    isol    = (k-1) * sol->size + 1;
    isolnew = (nbl-1) * sol->size + 1;

    for (i=0; i<sol->size; i++)
      sol->m[isolnew + i] = sol->m[isol + i];
    ++nbl;
  }


  /*compact vertices*/
  np  = 0;
  nbl = 1;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( ppt->tag & M_NUL )  continue;
    if(k!=nbl) {
      pptnew = &mesh->point[nbl];
      memcpy(pptnew,ppt,sizeof(MMG5_Point));
      ppt->tag   &= ~M_NUL;
      assert(ppt->tmp == nbl);
    }
    np++;
    if(k != nbl) {
      ppt = &mesh->point[k];
      memset(ppt,0,sizeof(MMG5_Point));
      ppt->tag    = M_NUL;
    }
    nbl++;
  }
  mesh->np = np;
  sol->np  = np;

  for(k=1 ; k<=mesh->np ; k++)
    mesh->point[k].tmp = 0;

  mesh->npnil = mesh->np + 1;
  for (k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  mesh->nanil = mesh->na + 1;
  for (k=mesh->nanil; k<mesh->namax-1; k++)
    mesh->edge[k].b = k+1;

  mesh->nenil = mesh->nt + 1;
  for (k=mesh->nenil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;

  return(1);
}

/*
  opt[0] = option
  opt[1] = ddebug
  opt[2] = noswap
  opt[3] = noinsert
  opt[4] = nomove
  opt[5] = imprim
  opt[6] = nr

  optdbl[0] = hgrad
  optdbl[1] =.dhd
*/
int MMG2_mmg2dlib(int opt[7],double optdbl[2],MMG5_pMesh mesh,MMG5_pSol sol,void (*titi)(int ,int ,int,int,int)) {
  double declic;
  int    nsiter;

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",M_VER,M_REL);
  fprintf(stdout,"     %s\n",M_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  MMG2_callbackinsert = titi;
  /* interrupts */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  atexit(endcod);

  // tminit(ctim,TIMEMAX);
  //chrono(ON,&ctim[0]);

  /* default values */
  mesh->info.imprim = opt[5];
  mesh->info.mem    = 0;
  mesh->info.ddebug = opt[1];
  mesh->info.iso = 0;
  mesh->info.lag = -1;
  mesh->info.hmin = -1;
  mesh->info.hmax = -1;
  mesh->info.hausd = 175.;
  switch(opt[0]) {
  case 0:
  case 1:
  case 2:
    break;
  case 6:
    mesh->info.iso = 1;
    break;
  case 9:
  case 99:
    mesh->info.lag = 0;
    break;
  default:
    fprintf(stdout,"option not recognized %d\n",opt[0]);
    exit(0);
  }
  mesh->info.noswap = opt[2];
  mesh->info.nomove = opt[4];
  mesh->info.noinsert = opt[3];
  mesh->info.hgrad  = optdbl[0];
  mesh->info.renum    = 0;
  if(opt[6])
    mesh->info.dhd  = -1;
  else
    mesh->info.dhd  = optdbl[1];
  /*this options are not used inside library version*/
  //qdegrad[0] = 10./ALPHA;
  //qdegrad[1] = 1.3;
  mesh->info.renum = 0;

  sol->type = 1;


  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  //chrono(ON,&ctim[1]);
  MMG2_inputdata(mesh,sol);
  if ( sol->np && sol->np != mesh->np ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    sol->np = 0;
  }

  setfunc(sol->size);//setfunc(sol.type);


  if(mesh->info.lag >= 0) {
    /*alloc Disp*/
    fprintf(stdout,"OPTION LAG NOT AVAILABLE WITH THE LIBRARY   \n");
    return(1);
  }
  //chrono(OFF,&ctim[1]);
  //fprintf(stdout,"  -- DATA READING COMPLETED.     %.2f sec.\n",
  //       gttime(ctim[1]));


  fprintf(stdout,"\n  %s\n   MODULE MMG2D-IMB/LJLL : %s (%s) %s\n  %s\n",
          M_STR,M_VER,M_REL,sol->size == 1 ? "ISO" : "ANISO",M_STR);
  fprintf(stdout,"  MAXIMUM NUMBER OF POINTS    (NPMAX) : %8d\n",mesh->npmax);
  fprintf(stdout,"  MAXIMUM NUMBER OF TRIANGLES (NTMAX) : %8d\n",mesh->ntmax);

  /* analysis */
  //chrono(ON,&ctim[2]);
  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : DATA ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");
  if ( !MMG2_scaleMesh(mesh,sol) )  return(1);
  MMG2_saveMesh(mesh,"tata.mesh");
  if ( mesh->nt && !MMG2_hashel(mesh) )  return(1);
  if ( !mesh->info.renum && !MMG2_chkmsh(mesh,1) )        return(1);
  /*geom : corner detection*/
  if ( mesh->info.dhd<0 )
    if( !MMG2_evalgeom(mesh) ) return(1);

  /*mesh gradation*/
  if( mesh->nt && mesh->info.hgrad > 0 ) {
    if ( mesh->info.imprim )   fprintf(stdout,"\n  -- GRADATION : %8f\n",mesh->info.hgrad);
    MMG2_lissmet(mesh,sol);
  }
  if ( mesh->nt && abs(mesh->info.imprim) > 3 )  MMG2_outqua(mesh,sol);

  if ( mesh->nt && abs(mesh->info.imprim) > 4 )  {
    MMG2_prilen(mesh,sol);
  }

  //chrono(OFF,&ctim[2]);
  //fprintf(stdout,"  -- PHASE 1 COMPLETED.     %.2f sec.\n",gttime(ctim[2]));

  /* remeshing */
  //chrono(ON,&ctim[3]);
  if ( mesh->info.iso ) {
    fprintf(stdout,"Fit an embedded mesh\n");
    MMG2_mmg2d6(mesh,sol);
    MMG2_saveMesh(mesh,mesh->nameout);
    return(0);
  } else if ( mesh->info.lag >= 0 ) {

#warning option 9
    printf("exit option 9 not implemented\n");
    exit(1);
  } else {

    if(!mesh->nt) {
      fprintf(stdout,"\n  -- PHASE 2 : MESH GENERATION\n");
      if ( !MMG2_mmg2d2(mesh,sol) )  return(1);
    } else {
      fprintf(stdout,"\n  -- PHASE 2 : MESH ADAPTATION\n");
      if ( (!mesh->info.noinsert) && !MMG2_mmg2d1(mesh,sol) )  return(1);
    }

    /* optimisation */
    //chrono(ON,&ctim[4]);
    fprintf(stdout,"\n  -- PHASE 3 : MESH OPTIMISATION\n");
    //if ( !optlap(&mesh,&sol) ) return(1);
    if ( !MMG2_mmg2d0(mesh,sol) )  return(1);
    // chrono(OFF,&ctim[4]);
    //fprintf(stdout,"  -- PHASE 3 COMPLETED.     %.2f sec.\n",gttime(ctim[4]));

  }
  //chrono(OFF,&ctim[3]);
  // fprintf(stdout,"  -- PHASE 2 COMPLETED.     %.2f sec.\n",gttime(ctim[3]));


  if ( !MMG2_unscaleMesh(mesh,sol) )  return(1);

  if ( abs(mesh->info.imprim) > 3 )  {
    MMG2_outqua(mesh,sol);
    MMG2_prilen(mesh,sol);
  }
  //chrono(ON,&ctim[1]);
  MMG2_saveMesh(mesh,"toto.mesh");
  MMG2_tassage(mesh,sol);
  //chrono(OFF,&ctim[1]);


  if ( mesh->info.imprim < -4 ) M_memDump();
  return(0);
}


