/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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
 * \brief Most of the API functions of the MMG3D library.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning Use the MMG3D_ prefix: MMG5_ prefix will became obsolete soon...
 * \todo documentation doxygen
 *
 * Private API functions for MMG3D library: incompatible functions
 * with the main binary.
 *
 */

#include "mmg3d.h"

/**
 * Pack the mesh \a mesh and its associated metric \a met and return \a val.
 */
#define _MMG5_RETURN_AND_PACK(mesh,met,disp,val)do                      \
  {                                                                     \
    if ( !_MMG3D_packMesh(mesh,met,disp) )  {                           \
      mesh->npi = mesh->np;                                             \
      mesh->nti = mesh->nt;                                             \
      mesh->nai = mesh->na;                                             \
      mesh->nei = mesh->ne;                                             \
      met->npi  = met->np;                                              \
      return(MMG5_LOWFAILURE);                                          \
    }                                                                   \
    _LIBMMG5_RETURN(mesh,met,val);                                      \
  }while(0)

/** Free adja, xtetra and xpoint tables */
static inline
void _MMG3D_Free_topoTables(MMG5_pMesh mesh) {
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
 * \return 1 if success, 0 if chkmsh fail or if we are unable to build
 * triangles.
 *
 * Pack the sparse mesh and create triangles and edges before getting
 * out of library
 *
 */
static inline
int _MMG3D_packMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol disp) {
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
    if ( ppt->tag & MG_CRN ) {
      if ( mesh->info.nosurf && (ppt->tag & MG_REQ) ) {
        // Warning: in -nosurf option : we loose the corner and the required
        // points
        ppt->tag &= ~MG_CRN;
        ppt->tag &= ~MG_REQ;
        continue;
      }
      nc++;
    }
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
  /* in the worst case (all edges are marked), we will have around 1 edge per *
   * triangle (we count edges only one time) */
  mesh->memCur += (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));
  if ( (mesh->memCur) > (mesh->memMax) ) {
    mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));
    fprintf(stdout,"  ## Warning:");
    fprintf(stdout," unable to allocate an hash table to store the mesh edges."
      "\n  ## Warning: uncomplete mesh\n");
  }
  else if ( _MMG5_hNew(&mesh->htab,mesh->nt,3*(mesh->nt),0) ) {
    for (k=1; k<=mesh->ne; k++) {
      pt   = &mesh->tetra[k];
      if ( MG_EOK(pt) &&  pt->xt ) {
        for (i=0; i<6; i++) {
          if ( mesh->xtetra[pt->xt].edg[i] )
            _MMG5_hEdge(mesh,pt->v[_MMG5_iare[i][0]],pt->v[_MMG5_iare[i][1]],
                        mesh->xtetra[pt->xt].edg[i],mesh->xtetra[pt->xt].tag[i]);
          else if ( mesh->xtetra[pt->xt].tag[i] & MG_REQ ) {
            if ( mesh->info.nosurf ) {
              if ( mesh->xtetra[pt->xt].tag[i]==MG_REQ+MG_CRN+MG_GEO )
                continue;
              else if ( mesh->xtetra[pt->xt].tag[i] & MG_CRN ) {
                mesh->xtetra[pt->xt].tag[i] &= ~MG_CRN;
                mesh->xtetra[pt->xt].tag[i] &= ~MG_REQ;
              }
              _MMG5_hEdge(mesh,pt->v[_MMG5_iare[i][0]],pt->v[_MMG5_iare[i][1]],
                          mesh->xtetra[pt->xt].edg[i],mesh->xtetra[pt->xt].tag[i]);
            }
          }
          else if ( MG_EDG(mesh->xtetra[pt->xt].tag[i] ) ) {
            _MMG5_hEdge(mesh,pt->v[_MMG5_iare[i][0]],pt->v[_MMG5_iare[i][1]],
                        mesh->xtetra[pt->xt].edg[i],mesh->xtetra[pt->xt].tag[i]);
          }
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
      _MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"edges",
                    mesh->na = 0;
                    printf("  ## Warning: uncomplete mesh\n")
        );
    }

    if ( mesh->na ) {
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
 * \param met pointer toward a sol structure (metric).
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 *
 * Main program for the remesh library .
 *
 */
int MMG3D_mmg3dlib(MMG5_pMesh mesh,MMG5_pSol met) {
  mytime    ctim[TIMEMAX];
  char      stim[32];

  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",MG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
  }

  _MMG3D_Set_commonFunc();

  /** Free topologic tables (adja, xpoint, xtetra) resulting from a previous
   * run */
  _MMG3D_Free_topoTables(mesh);

  signal(SIGABRT,_MMG5_excfun);
  signal(SIGFPE,_MMG5_excfun);
  signal(SIGILL,_MMG5_excfun);
  signal(SIGSEGV,_MMG5_excfun);
  signal(SIGTERM,_MMG5_excfun);
  signal(SIGINT,_MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->info.lag > -1 ) {
    fprintf(stdout,"  ## Error: lagrangian mode unavailable (MMG3D_IPARAM_lag):\n"
            "            You must call the MMG3D_mmg3dmov function to move a rigidbody.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.iso ) {
    fprintf(stdout,"  ## Error: level-set discretisation unavailable"
            " (MMG3D_IPARAM_iso):\n"
            "          You must call the MMG3D_mmg3dmov function to use this option.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

#ifdef USE_SCOTCH
  _MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  _MMG5_warnOrientation(mesh);



  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));
    met->np = 0;
  }
  else if ( met->size!=1 && met->size!=6 ) {
    fprintf(stdout,"  ## ERROR: WRONG DATA TYPE.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMG3D_setfunc(mesh,met);

  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  /* scaling mesh */
  if ( !_MMG5_scaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  if ( abs(mesh->info.imprim) > 0 ) {
    if ( !_MMG3D_inqua(mesh,met) ) {
      if ( !_MMG5_unscaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
    }
  }

  /* specific meshing */
  if ( mesh->info.optim && !met->np && !MMG3D_DoSol(mesh,met) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
  }

  /* mesh analysis */
  if ( !_MMG3D_analys(mesh) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 1 && met->m ) _MMG3D_prilen(mesh,met,0);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");
  }

  /* renumerotation if available */
  if ( !_MMG5_scotchCall(mesh,met) )
  {
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }

#ifdef PATTERN
  if ( !_MMG5_mmg3d1_pattern(mesh,met) ) {
    if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }
#else
  if ( !_MMG5_mmg3d1_delone(mesh,met) ) {
    if ( (!mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }
#endif

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG3d: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  /* save file */
  if ( !_MMG3D_outqua(mesh,met) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 4 )
    _MMG3D_prilen(mesh,met,1);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");
  if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  if ( !_MMG3D_packMesh(mesh,met,NULL) )     _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMG3DLIB: ELAPSED TIME  %s\n",stim);
  _LIBMMG5_RETURN(mesh,met,MMG5_SUCCESS);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (level-set).
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 *
 * Main program for the levelset discretization library .
 *
 */
int MMG3D_mmg3dls(MMG5_pMesh mesh,MMG5_pSol met) {
  mytime    ctim[TIMEMAX];
  char      stim[32];

  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",MG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
  }

  _MMG3D_Set_commonFunc();

  signal(SIGABRT,_MMG5_excfun);
  signal(SIGFPE,_MMG5_excfun);
  signal(SIGILL,_MMG5_excfun);
  signal(SIGSEGV,_MMG5_excfun);
  signal(SIGTERM,_MMG5_excfun);
  signal(SIGINT,_MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->info.lag > -1 ) {
    fprintf(stdout,"  ## Error: lagrangian mode unavailable (MMG3D_IPARAM_lag):\n"
            "            You must call the MMG3D_mmg3dmov function to move a rigidbody.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

#ifdef USE_SCOTCH
  _MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG3DLS: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  _MMG5_warnOrientation(mesh);

  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));
    met->np = 0;
  }
  else if ( met->size!=1 ) {
    fprintf(stdout,"  ## ERROR: WRONG DATA TYPE.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMG3D_setfunc(mesh,met);

  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

 /* scaling mesh */
  if ( !_MMG5_scaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  if ( abs(mesh->info.imprim) > 0 ) {
    if ( !_MMG3D_inqua(mesh,met) ) {
      if ( !_MMG5_unscaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
    }
  }

  /* specific meshing */
  if ( !met->np ) {
    fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  if ( !_MMG5_mmg3d2(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  /* mesh analysis */
  if ( !_MMG3D_analys(mesh) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }


  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 2 : ISOTROPIC MESHING\n");
  }

  /* renumerotation if available */
  if ( !_MMG5_scotchCall(mesh,met) )
  {
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }

#ifdef PATTERN
  if ( !_MMG5_mmg3d1_pattern(mesh,met) ) {
    if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }
#else
  if ( !_MMG5_mmg3d1_pattern(mesh,met) ) {
    if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }
#endif

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG3d: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  /* save file */
  if ( !_MMG3D_outqua(mesh,met) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");
  if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  if ( !_MMG3D_packMesh(mesh,met,NULL) )     _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMG3DLS: ELAPSED TIME  %s\n",stim);
  _LIBMMG5_RETURN(mesh,met,MMG5_SUCCESS);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (for the output metric).
 * \param disp pointer toward a sol structure (displacement for the
 * lagrangian motion mode).
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 *
 * Main program for the rigidbody moving library.
 *
 * \warning set the disp->npi field at the good value at the good moment (needed
 * by API).
 */
int MMG3D_mmg3dmov(MMG5_pMesh mesh,MMG5_pSol met, MMG5_pSol disp) {
  mytime    ctim[TIMEMAX];
  char      stim[32];

  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",MG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
  }

  _MMG3D_Set_commonFunc();

  signal(SIGABRT,_MMG5_excfun);
  signal(SIGFPE,_MMG5_excfun);
  signal(SIGILL,_MMG5_excfun);
  signal(SIGSEGV,_MMG5_excfun);
  signal(SIGTERM,_MMG5_excfun);
  signal(SIGINT,_MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->info.iso ) {
    fprintf(stdout,"  ## Error: level-set discretisation unavailable"
            " (MMG3D_IPARAM_iso):\n"
            "          You must call the MMG3D_mmg3dmov function to use this option.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

#ifdef USE_SCOTCH
  _MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG3DMOV: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  _MMG5_warnOrientation(mesh);

  if ( mesh->info.lag == -1 ) {
    if ( mesh->info.imprim )
      fprintf(stdout,"  ## Warning: displacement mode for the rigidbody"
              " movement is not set.\n"
              "               Lagrangian displacement computed according"
              " to mode 1.\n");
    mesh->info.lag = 1;
  }

#ifndef USE_SUSCELAS
  fprintf(stdout,"  ## Error: you need to compile with the USE_SUSCELAS"
    " CMake's flag set to ON to use the rigidbody movement library.\n");
  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
#endif

  if ( !disp ) {
    fprintf(stdout,"  ## Error: in lagrangian mode, a structure of type"
            " \"MMG5_pSol\" is needed to store the displacement field.\n"
            "            This structure must be different from the one used"
            " to store the metric.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  if (disp->np && (disp->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,disp->m,(disp->size*(disp->npmax+1))*sizeof(double));
    disp->np = 0;
  }
  else if (disp->size!=3) {
    fprintf(stdout,"  ## ERROR: LAGRANGIAN MOTION OPTION NEED A VECTORIAL DISPLACEMENT\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMG3D_setfunc(mesh,met);

  if ( abs(mesh->info.imprim) > 0 ) {
    if ( !_MMG3D_inqua(mesh,met) ) {
      _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
    }
  }

  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

 /* scaling mesh */
  if ( !_MMG5_scaleMesh(mesh,disp) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);


  /* mesh analysis */
  if ( !_MMG3D_analys(mesh) ) {
    if ( !_MMG5_unscaleMesh(mesh,disp) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 4 && !mesh->info.iso && met->m ) _MMG3D_prilen(mesh,met,0);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 2 : LAGRANGIAN MOTION\n");
  }

  /* renumerotation if available */
  if ( !_MMG5_scotchCall(mesh,met) )
  {
    if ( !_MMG5_unscaleMesh(mesh,disp) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }

#ifdef USE_SUSCELAS
  /* Lagrangian mode */
  if ( !_MMG5_mmg3d3(mesh,disp,met) ) {
    disp->npi = disp->np;
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
#endif
  disp->npi = disp->np;

  if ( !met->np && !MMG3D_DoSol(mesh,met) ) {
    if ( !_MMG5_unscaleMesh(mesh,disp) ) {
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }

/*   /\* ******************* Part to skip in lag mode ? *************************** *\/ */
/* #ifdef PATTERN */
/*   if ( !_MMG5_mmg3d1_pattern(mesh,met) ) { */
/*     if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) { */
/*       fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n"); */
/*       disp->npi = disp->np; */
/*       _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE); */
/*     } */
/*     if ( !_MMG5_unscaleMesh(mesh,met) ) { */
/*       disp->npi = disp->np; */
/*       _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE); */
/*     } */
/*     _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE); */
/*   } */
/* #else */
/*   /\** Patterns in iso mode, delauney otherwise *\/ */
/*   if ( !mesh->info.iso ) { */
/*     if ( !_MMG5_mmg3d1_delone(mesh,met) ) { */
/*       if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) { */
/*         fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n"); */
/*         disp->npi = disp->np; */
/*         _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE); */
/*       } */
/*       if ( !_MMG5_unscaleMesh(mesh,met) ) { */
/*         disp->npi = disp->np; */
/*         _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE); */
/*       } */
/*       _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE); */
/*     } */
/*   } */
/*   else { */
/*     if ( !_MMG5_mmg3d1_pattern(mesh,met) ) { */
/*       if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) { */
/*         fprintf(stdout,"  ## Hashing problem. Invalid mesh.\n"); */
/*         disp->npi = disp->np; */
/*         _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE); */
/*       } */
/*       if ( !_MMG5_unscaleMesh(mesh,met) ) { */
/*         disp->npi = disp->np; */
/*         _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE); */
/*       } */
/*       _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE); */
/*     } */
/*   } */
/* #endif */

/* *************************************** End of part to skip in lag mode ? *************************** */

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG3D: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  /* save file */
  if ( !_MMG3D_outqua(mesh,met) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) ) {
      disp->npi = disp->np;
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 1 && !mesh->info.iso )
    _MMG3D_prilen(mesh,met,1);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");
  if ( !_MMG5_unscaleMesh(mesh,disp) ) {
    disp->npi = disp->np;
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  if ( !_MMG3D_packMesh(mesh,met,disp) ) {
    disp->npi = disp->np;
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMG3DMOV: ELAPSED TIME  %s\n",stim);
  disp->npi = disp->np;
  _LIBMMG5_RETURN(mesh,met,MMG5_SUCCESS);
}

/** Old API °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°*/

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (metric or solution).
 * lagrangian motion mode).
 * \return Return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 * \warning Do not use: will became obsolete
 *
 * Main program for the library .
 *
 */
int MMG5_mmg3dlib(MMG5_pMesh mesh,MMG5_pSol met)
{
  printf("  ## MMG5_mmg3dlib: "
         "MMG5_ API is deprecated (replaced by the MMG3D_ one) and will"
        " be removed soon\n." );
  return(MMG3D_mmg3dlib(mesh,met));
}
