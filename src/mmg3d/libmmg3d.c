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
 * \file mmg3d/libmmg3d.c
 * \brief Private API functions for MMG3D library.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \todo documentation doxygen
 *
 * Private API functions for MMG3D library: incompatible functions
 * with the main binary.
 *
 */

#include "mmg3d.h"
#include "shared_func.h"

/**
 * Pack the mesh \a mesh and its associated metric \a met and return \a val.
 */
#define _MMG5_RETURN_AND_PACK(mesh,met,disp,val)do  \
  {                                                 \
    MMG5_packMesh(mesh,met,disp);                   \
    return(val);                                    \
  }while(0)

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure (metric or solution).
 * \param disp pointer toward a sol structure (displacement).
 *
 * Deallocations before return.
 *
 */
void MMG5_Free_all(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol disp
  ){

  MMG5_Free_structures(mesh,met,disp);

  _MMG5_SAFE_FREE(met);
  _MMG5_SAFE_FREE(mesh);
  _MMG5_SAFE_FREE(disp);
}

/** Free adja, xtetra and xpoint tables */
static inline
void MMG5_Free_topoTables(MMG5_pMesh mesh) {
  int k;

  mesh->xp = 0;
  if ( mesh->adja )
    _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

  _MMG5_freeXTets(mesh);

  if ( mesh->xpoint )
    _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
  for(k=1; k <=mesh->np; k++) {
    mesh->point[k].xp = 0;
  }

  return;
}

/**
 * \param mesh pointer toward the mesh structure (unused).
 * \param met pointer toward the solution (metric or level-set) structure.
 * \param disp pointer toward the solution (displacement) structure.
 *
 * Pack the sparse mesh and create triangles and edges before getting
 * out of library
 *
 */
static inline
int MMG5_packMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol disp) {
  MMG5_pTetra   pt,ptnew;
  MMG5_pPoint   ppt,pptnew;
  MMG5_hgeom   *ph;
  int     np,nc,nr, k,ne,nbl,imet,imetnew,i;
  int     iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;

  /* compact vertices */
  np = nc = nr = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->tmp = ++np;
    if ( ppt->tag & MG_CRN )  nc++;
    ppt->ref = abs(ppt->ref);
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
      memcpy(ptnew,pt,sizeof(MMG5_Tetra));

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
  if ( met && met->m ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      imet    = k   * met->size;
      imetnew = nbl * met->size;

      for (i=0; i<met->size; i++)
        met->m[imetnew + i] = met->m[imet + i];
      ++nbl;
    }
  }

  /* compact displacement */
  nbl = 1;
  if ( disp && disp->m ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      imet    = k   * disp->size;
      imetnew = nbl * disp->size;

      for (i=0; i<disp->size; i++)
        disp->m[imetnew + i] = disp->m[imet + i];
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
      memmove(pptnew,ppt,sizeof(MMG5_Point));
      memset(ppt,0,sizeof(MMG5_Point));
      ppt->tag    = MG_NUL;
    }
    nbl++;
  }
  mesh->np = np;
  if ( met && met->m )
    met->np  = np;
  if ( disp && disp->m )
    disp->np = np;

  /* rebuild triangles*/
  mesh->nt = 0;
  _MMG5_chkNumberOfTri(mesh);
  if ( !_MMG5_bdryTria(mesh) ) {
    fprintf(stdout," ## Error: unable to rebuild triangles\n");
    return(0);
  }

  /* build hash table for edges */
  if ( mesh->htab.geom )
    _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));

  mesh->na = 0;
  /* in the wost case (all edges are marked), we will have around 1 edge per *
   * triangle (we count edges only one time) */
  mesh->memCur += (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));
  if ( (mesh->memCur) > (mesh->memMax) ) {
    mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));
    fprintf(stdout,"  ## Warning:");
    fprintf(stdout," unable to allocate htab.\n");
  }
  else if ( _MMG5_hNew(&mesh->htab,mesh->nt,3*(mesh->nt),0) ) {
    for (k=1; k<=mesh->ne; k++) {
      pt   = &mesh->tetra[k];
      if ( MG_EOK(pt) &&  pt->xt ) {
        for (i=0; i<6; i++) {
          if ( mesh->xtetra[pt->xt].edg[i] ||
               ( MG_EDG(mesh->xtetra[pt->xt].tag[i] ) ||
                 (mesh->xtetra[pt->xt].tag[i] & MG_REQ) ) )
            _MMG5_hEdge(mesh,pt->v[_MMG5_iare[i][0]],pt->v[_MMG5_iare[i][1]],
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
      _MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"edges",);
      _MMG5_SAFE_CALLOC(mesh->edge,mesh->na+1,MMG5_Edge);

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
    _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));
  }
  else
    mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));

  for(k=1 ; k<=mesh->np ; k++)
    mesh->point[k].tmp = 0;

  mesh->npnil = mesh->np + 1;
  for(k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  mesh->nenil = mesh->ne + 1;
  for(k=mesh->nenil; k<mesh->nemax-1; k++)
    mesh->tetra[k].v[3] = k+1;

  /* to could save the mesh, the adjacency have to be correct */
  if ( mesh->info.ddebug && (!_MMG5_chkmsh(mesh,1,1) ) ) {
    fprintf(stdout,"  ##  Problem. Invalid mesh.\n");
    return(0);
  }

  MMG5_Free_topoTables(mesh);

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

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (metric or solution).
 * \param disp pointer toward a sol structure (displacement).
 * \return Return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 *
 * Main program for the library .
 *
 */
int MMG5_mmg3dlib(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol disp
  ) {
  mytime    ctim[TIMEMAX];
  char      stim[32];
  MMG5_pSol scaledSol;

  fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  _MMG5_Set_commonFunc();

  signal(SIGABRT,_MMG5_excfun);
  signal(SIGFPE,_MMG5_excfun);
  signal(SIGILL,_MMG5_excfun);
  signal(SIGSEGV,_MMG5_excfun);
  signal(SIGTERM,_MMG5_excfun);
  signal(SIGINT,_MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

#ifdef USE_SCOTCH
  _MMG5_warnScotch(mesh);
#endif

  fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  _MMG5_warnOrientation(mesh);

  if ( mesh->info.lag > -1 ) {
    if (disp->np && (disp->np != mesh->np) ) {
      fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
      _MMG5_DEL_MEM(mesh,disp->m,(disp->size*(disp->npmax+1))*sizeof(double));
      disp->np = 0;
    }
    else if (disp->size!=3) {
      fprintf(stdout,"  ## ERROR: LAGRANGIAN MOTION OPTION NEED A VECTORIAL DISPLACEMENT\n");
      return(MMG5_STRONGFAILURE);
    }
  }
  else {
    if ( met->np && (met->np != mesh->np) ) {
      fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
      _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));
      met->np = 0;
    }
    else if ( met->size!=1 && met->size!=6 ) {
      fprintf(stdout,"  ## ERROR: WRONG DATA TYPE.\n");
      return(MMG5_STRONGFAILURE);
    }
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  _MMG5_setfunc(mesh,met);
  if ( abs(mesh->info.imprim) > 0 )  _MMG5_outqua(mesh,met);
  fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");

 /* scaling mesh */
  if ( mesh->info.lag == -1 ) {
    if ( !_MMG5_scaleMesh(mesh,met) ) return(MMG5_STRONGFAILURE);
    scaledSol = met;
  }
  else {
    if ( !_MMG5_scaleMesh(mesh,disp) ) return(MMG5_STRONGFAILURE);
    scaledSol = disp;
  }

  /* specific meshing */
  if ( mesh->info.iso ) {
    if ( !met->np ) {
      fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
      return(MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_mmg3d2(mesh,met) ) return(MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.lag < 0 ) {
    if ( mesh->info.optim && !met->np && !_MMG5_DoSol(mesh,met) ) {
      if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
      _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
    }
  }

  /* mesh analysis */
  if ( !_MMG5_analys(mesh) ) {
    if ( !_MMG5_unscaleMesh(mesh,scaledSol) )  return(MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 4 && !mesh->info.iso && met->m ) _MMG5_prilen(mesh,met);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim ) {
    if ( mesh->info.lag < 0 )
      fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");
    else
      fprintf(stdout,"\n  -- PHASE 2 : LAGRANGIAN MOTION\n");
  }

  /* renumerotation if available */
  if ( !_MMG5_scotchCall(mesh,met) )
  {
    if ( !_MMG5_unscaleMesh(mesh,scaledSol) )  return(MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }

  /* Lagrangian mode */
#ifdef USE_SUSCELAS
  if ( mesh->info.lag >= 0 ) {
    if ( !_MMG5_mmg3d3(mesh,disp,met) ) {
      return(MMG5_STRONGFAILURE);
    }
    if ( !met->np && !_MMG5_DoSol(mesh,met) ) {
      if ( !_MMG5_unscaleMesh(mesh,disp) )  return(MMG5_STRONGFAILURE);
      _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
    }
  }
#endif


/* ******************* Part to skip in lag mode ? *************************** */
if ( mesh->info.lag == -1 ) {

#ifdef PATTERN
  if ( !_MMG5_mmg3d1_pattern(mesh,met) ) {
    if ( !(mesh->adja) && !_MMG5_hashTetra(mesh,1) ) {
      fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
      return(MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }
#else
  /** Patterns in iso mode, delauney otherwise */
  if ( !mesh->info.iso ) {
    if ( !_MMG5_mmg3d1_delone(mesh,met) ) {
      if ( !(mesh->adja) && !_MMG5_hashTetra(mesh,1) ) {
        fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
        return(MMG5_STRONGFAILURE);
      }
      if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
      _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
    }
  }
  else {
    if ( !_MMG5_mmg3d1_pattern(mesh,met) ) {
      if ( !(mesh->adja) && !_MMG5_hashTetra(mesh,1) ) {
        fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
        return(MMG5_STRONGFAILURE);
      }
      if ( !_MMG5_unscaleMesh(mesh,met) )  return(MMG5_STRONGFAILURE);
      _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
    }
  }
#endif

}

/* *************************************** End of part to skip in lag mode ? *************************** */

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  fprintf(stdout,"\n  %s\n   END OF MODULE MMG3d: IMB-LJLL \n  %s\n",MG_STR,MG_STR);

  /* save file */
  _MMG5_outqua(mesh,met);
  if ( mesh->info.imprim > 4 && !mesh->info.iso )
    _MMG5_prilen(mesh,met);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");
  if ( !_MMG5_unscaleMesh(mesh,scaledSol) )  return(MMG5_STRONGFAILURE);
  if ( !MMG5_packMesh(mesh,met,disp) )     return(MMG5_STRONGFAILURE);
  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3DLIB: ELAPSED TIME  %s\n",stim);
  return(MMG5_SUCCESS);
}
