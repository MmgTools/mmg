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

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * Deallocations before return.
 *
 */
void MMG2D_Free_all(MMG5_pMesh mesh,MMG5_pSol met
  ){

  MMG2D_Free_structures(mesh,met);

}


int MMG2_tassage(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pEdge         ped;
  MMG5_pTria         pt,ptnew;
  MMG5_pPoint        ppt,pptnew;
  int                np,nt,k,nbl,isol,isolnew,i;
  int                iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;


  /* compact vertices */
  np=0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( ppt->tag & M_NUL )  continue;
    ppt->tmp = ++np;
  }


  /* compact edges */
  nbl = 0;
  for (k=1; k<=mesh->na; k++) {
    ped  = &mesh->edge[k];
    if(!ped->a) continue;
    ped->a = mesh->point[ped->a].tmp;
    ped->b = mesh->point[ped->b].tmp;
    /* impossible to do that without update triangle....*/
    /* if(k!=nbl) { */
    /*   pednew = &mesh->edge[nbl]; */
    /*   memcpy(pednew,ped,sizeof(MMG5_Edge)); */
    /*   memset(ped,0,sizeof(MMG5_Tria)); */
    /* } */
    /*nbl++;*/
  }
  /* mesh->na = nbl;*/

  /* compact triangle */
  nt  = 0;
  nbl = 1;
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
      ptnew = &mesh->tria[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tria));
      //and the adjacency
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      iadrnew = 3*(nbl-1) + 1;
      adjanew = &mesh->adja[iadrnew];
      for(i=0 ; i<3 ; i++) {
        adjanew[i] = adja[i];
        if(!adja[i]) continue;
        //if(adja[i]/3==2414 || nbl==2414) printf("adja of %d nbl %d : %d\n",k,nbl,adja[i]/3);
        iadrv = 3*(adja[i]/3-1) +1;
        adjav = &mesh->adja[iadrv];
        //if(k==2414 || nbl==2414) printf("on met %d pour %d\n",nbl,adja[i]/3);
        voy = i;
        adjav[adja[i]%3] = 3*nbl + voy;

        adja[i] = 0;
      }
      memset(pt,0,sizeof(MMG5_Tria));
    }
    nbl++;
  }
  mesh->nt = nt;

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

 if(mesh->np < mesh->npmax - 3) {
    mesh->npnil = mesh->np + 1;
    for (k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;
  } else {
    mesh->npnil = 0;
  }
 
 /*to do only if the edges are packed*/
  /* if(mesh->na < mesh->namax - 3) { */
  /*   mesh->nanil = mesh->na + 1; */
  /*   for (k=mesh->nanil; k<mesh->namax-1; k++) */
  /*     mesh->edge[k].b = k+1; */
  /* } else { */
  /*   mesh->nanil = 0; */
  /* } */

  if(mesh->nt < mesh->ntmax - 3) {
    mesh->nenil = mesh->nt + 1;
    for (k=mesh->nenil; k<mesh->ntmax-1; k++)
      mesh->tria[k].v[2] = k+1;
  } else {
    mesh->nenil = 0;
  }

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
  optdbl[1] =ar
*/
int MMG2D_mmg2dlib(MMG5_pMesh mesh,MMG5_pSol sol,void (*titi)(int ,int ,int,int,int)) {

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",M_VER,M_REL);
  fprintf(stdout,"     %s\n",M_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  MMG2D_callbackinsert = titi;
  /* interrupts */
  signal(SIGABRT,_MMG2_excfun);
  signal(SIGFPE,_MMG2_excfun);
  signal(SIGILL,_MMG2_excfun);
  signal(SIGSEGV,_MMG2_excfun);
  signal(SIGTERM,_MMG2_excfun);
  signal(SIGINT,_MMG2_excfun);

  sol->np = mesh->np;
  sol->size = 1;
  sol->ver  = mesh->ver;
  if ( !sol->m ) {
    /* mem alloc */
    _MMG5_SAFE_CALLOC(sol->m,sol->size*mesh->npmax,double);
    sol->np = 0;
  } else   if ( sol->np && (sol->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER : %d != %d\n",sol->np,mesh->np);
    //exit(1);
  }


  // tminit(ctim,TIMEMAX);
  //chrono(ON,&ctim[0]);

  /* /\* default values *\/ */
  /* mesh->info.imprim = opt[5]; */
  /* mesh->info.mem    = 0; */
  /* mesh->info.ddebug = opt[1]; */
  /* mesh->info.iso = 0; */
  /* mesh->info.lag = -1; */
  /* mesh->info.hmin = -1; */
  /* mesh->info.hmax = -1; */
  /* mesh->info.hausd = 0.01; */
  /* switch(opt[0]) { */
  /* case 0: */
  /* case 1: */
  /* case 2: */
  /*   break; */
  /* case 6: */
  /*   mesh->info.iso = 1; */
  /*   break; */
  /* case 9: */
  /* case 99: */
  /*   mesh->info.lag = 0; */
  /*   break; */
  /* default: */
  /*   fprintf(stdout,"option not recognized %d\n",opt[0]); */
  /*   exit(0); */
  /* } */
  /* mesh->info.noswap = opt[2]; */
  /* mesh->info.nomove = opt[4]; */
  /* mesh->info.noinsert = opt[3]; */
  /* mesh->info.hgrad  = optdbl[0]; */
  /* if(opt[6]) */
  /*   mesh->info.dhd  = -1; */
  /* else */
  /*   mesh->info.dhd  = 180.-optdbl[1]; */
  /* /\*this options are not used inside library version*\/ */
  /* //qdegrad[0] = 10./ALPHA; */
  /* //qdegrad[1] = 1.3; */
  /* mesh->info.renum = 0; */

  /* sol->type = 1; */

  MMG2D_setfunc(sol->size);
  _MMG2D_Set_commonFunc();
 
  if(mesh->info.lag >= 0) {
    /*alloc Disp*/
    fprintf(stdout,"OPTION LAG NOT AVAILABLE WITH THE LIBRARY   \n");
    return(MMG5_STRONGFAILURE);
  }
  //chrono(OFF,&ctim[1]);
  //fprintf(stdout,"  -- DATA READING COMPLETED.     %.2f sec.\n",
  //       gttime(ctim[1]));


  fprintf(stdout,"\n  %s\n   MODULE MMG2D-IMB/LJLL : %s (%s) %s\n  %s\n",
          M_STR,M_VER,M_REL,sol->size == 1 ? "ISO" : "ANISO",M_STR);
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM NUMBER OF POINTS    (NPMAX) : %8d\n",mesh->npmax);
    fprintf(stdout,"  MAXIMUM NUMBER OF TRIANGLES (NTMAX) : %8d\n",mesh->ntmax);
  }

  /* analysis */
  //chrono(ON,&ctim[2]);
  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : DATA ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");
  if ( !MMG2_scaleMesh(mesh,sol) )  return(MMG5_STRONGFAILURE);
  if ( !sol->np && !MMG2_doSol(mesh,sol) )  return(MMG5_STRONGFAILURE);

  if ( mesh->nt && !MMG2_hashel(mesh) )  return(MMG5_STRONGFAILURE);
  if ( mesh->info.ddebug && !_MMG5_chkmsh(mesh,1,0) )  return(MMG5_STRONGFAILURE);
  /*geom : corner detection*/
  if ( mesh->info.dhd>0 )
    if( !MMG2_evalgeom(mesh) ) return(MMG5_STRONGFAILURE);


  /*mesh gradation*/
  if( mesh->nt && mesh->info.hgrad > 0 ) {
    if ( mesh->info.imprim )   fprintf(stdout,"\n  -- GRADATION : %8f\n",mesh->info.hgrad);
    MMG2_lissmet(mesh,sol);
  }
  if ( mesh->nt )  MMG2_outqua(mesh,sol);

  if ( mesh->nt && abs(mesh->info.imprim) > 1 )  {
    MMG2_prilen(mesh,sol);
  }

  //chrono(OFF,&ctim[2]);
  //fprintf(stdout,"  -- PHASE 1 COMPLETED.     %.2f sec.\n",gttime(ctim[2]));

  /* remeshing */
  //chrono(ON,&ctim[3]);
  if ( mesh->info.iso ) {
    fprintf(stdout,"Fit an embedded mesh\n");
    MMG2_mmg2d6(mesh,sol);
    MMG2D_saveMesh(mesh,mesh->nameout);
    return(MMG5_SUCCESS);
  } else if ( mesh->info.lag >= 0 ) {

//#warning option 9
    printf("exit option 9 not implemented\n");
    exit(1);
  } else {

    if(!mesh->nt) {
      if ( mesh->info.imprim )
        fprintf(stdout,"\n  -- PHASE 2 : MESH GENERATION\n");
      if ( !MMG2_mmg2d2(mesh,sol) )  {
        if ( !MMG2_unscaleMesh(mesh,sol) )  return(MMG5_STRONGFAILURE);
        MMG2D_saveMesh(mesh,mesh->nameout);
        return(MMG5_LOWFAILURE);
      }
      /*geom : corner detection*/
      if ( mesh->info.dhd>0 )
        if( !MMG2_evalgeom(mesh) ) return(MMG5_STRONGFAILURE);


      /*mesh gradation*/
      if( mesh->nt && mesh->info.hgrad > 0 ) {
        if ( mesh->info.imprim )   fprintf(stdout,"\n  -- GRADATION : %8f\n",mesh->info.hgrad);
        MMG2_lissmet(mesh,sol);
      }
      if ( mesh->nt )  MMG2_outqua(mesh,sol);
      
      if ( mesh->nt && abs(mesh->info.imprim) > 1 )  {
        MMG2_prilen(mesh,sol);
      }
      if ( (!mesh->info.noinsert) && !MMG2_mmg2d1(mesh,sol) )  return(MMG5_LOWFAILURE);
      
    } else {
      if ( mesh->info.imprim )
        fprintf(stdout,"\n  -- PHASE 2 : MESH ADAPTATION\n");
      if ( (!mesh->info.noinsert) && !MMG2_mmg2d1(mesh,sol) )  return(MMG5_LOWFAILURE);
    }

    /* optimisation */
    //chrono(ON,&ctim[4]);
    if ( mesh->info.imprim )
      fprintf(stdout,"\n  -- PHASE 3 : MESH OPTIMISATION\n");
    //if ( !optlap(&mesh,&sol) ) return(1);
    if ( !MMG2_mmg2d0(mesh,sol) )  return(MMG5_LOWFAILURE);

    // chrono(OFF,&ctim[4]);
    //fprintf(stdout,"  -- PHASE 3 COMPLETED.     %.2f sec.\n",gttime(ctim[4]));

  }
  //chrono(OFF,&ctim[3]);
  // fprintf(stdout,"  -- PHASE 2 COMPLETED.     %.2f sec.\n",gttime(ctim[3]));


  if ( !MMG2_unscaleMesh(mesh,sol) )  return(MMG5_STRONGFAILURE);

  MMG2_outqua(mesh,sol);

  if ( abs(mesh->info.imprim) > 1 )  {
    MMG2_prilen(mesh,sol);
  }
  //chrono(ON,&ctim[1]);
  //MMG2D_saveMesh(mesh,"toto.mesh");
  MMG2_tassage(mesh,sol);
  //chrono(OFF,&ctim[1]);


  return(MMG5_SUCCESS);

}


