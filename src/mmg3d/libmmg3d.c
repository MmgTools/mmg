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

#include "inlined_functions_3d.h"

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
inline
void _MMG3D_Free_topoTables(MMG5_pMesh mesh) {
  int k;

  mesh->xp = 0;
  if ( mesh->adja )
    _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

  _MMG5_freeXTets(mesh);

  if ( mesh->adjapr )
    _MMG5_DEL_MEM(mesh,mesh->adjapr,(5*mesh->nprism+6)*sizeof(int));

  _MMG5_freeXPrisms(mesh);

  if ( mesh->xpoint )
    _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));

  for(k=1; k <=mesh->np; k++) {
    mesh->point[k].xp = 0;
  }

  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the solution structure.
 *
 * Truncate the metric computed by the DoSol function by hmax and hmin values
 * (if setted by the user). Set hmin and hmax if they are not setted.
 *
 */
void _MMG3D_solTruncatureForOptim(MMG5_pMesh mesh, MMG5_pSol met) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  double      isqhmin, isqhmax;
  int         i,k,iadr,sethmin,sethmax;

  assert ( mesh->info.optim || mesh->info.hsiz > 0. );

  /* Detect the point used only by prisms */
  if ( mesh->nprism ) {
    for (k=1; k<=mesh->np; k++) {
      mesh->point[k].flag = 1;
    }
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;

      for (i=0; i<4; i++) {
        mesh->point[pt->v[i]].flag = 0;
      }
    }
  }


  /* If not provided by the user, compute hmin/hmax from the metric computed by
   * the DoSol function. */
  sethmin = sethmax = 1;
  if ( mesh->info.hmin < 0 ) {
    sethmin = 0;
    if ( met->size == 1 ) {
      mesh->info.hmin = FLT_MAX;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        mesh->info.hmin = MG_MIN(mesh->info.hmin,met->m[k]);
      }
    }
    else if ( met->size == 6 ){
      mesh->info.hmin = 0.;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        iadr = met->size*k;
        mesh->info.hmin = MG_MAX(mesh->info.hmin,met->m[iadr]);
        mesh->info.hmin = MG_MAX(mesh->info.hmin,met->m[iadr+3]);
        mesh->info.hmin = MG_MAX(mesh->info.hmin,met->m[iadr+5]);
      }
      mesh->info.hmin = 1./sqrt(mesh->info.hmin);
    }
  }

  if ( mesh->info.hmax < 0 ) {
    sethmax = 1;
    if ( met->size == 1 ) {
      mesh->info.hmax = 0.;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        mesh->info.hmax = MG_MAX(mesh->info.hmax,met->m[k]);
      }
    }
    else if ( met->size == 6 ){
      mesh->info.hmax = FLT_MAX;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        iadr = met->size*k;
        mesh->info.hmax = MG_MIN(mesh->info.hmax,met->m[iadr]);
        mesh->info.hmax = MG_MIN(mesh->info.hmax,met->m[iadr+3]);
        mesh->info.hmax = MG_MIN(mesh->info.hmax,met->m[iadr+5]);
      }
      mesh->info.hmax = 1./sqrt(mesh->info.hmax);
    }
  }


  if ( !sethmin ) {
    mesh->info.hmin *=.1;
    /* Check that user has not given a hmax value lower that the founded
     * hmin. */
    if ( mesh->info.hmin > mesh->info.hmax ) {
      mesh->info.hmin = 0.1*mesh->info.hmax;
    }
  }
  if ( !sethmax ) {
    mesh->info.hmax *=10.;
    /* Check that user has not given a hmin value bigger that the founded
     * hmax. */
    if ( mesh->info.hmax < mesh->info.hmin ) {
      mesh->info.hmax = 10.*mesh->info.hmin;
    }
  }

  /* vertex size */
  if ( met->size == 1 ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      met->m[k] = MG_MIN(mesh->info.hmax,MG_MAX(mesh->info.hmin,met->m[k]));
    }
  }
  else if ( met->size == 6 ) {
    isqhmin = 1./(mesh->info.hmin*mesh->info.hmin);
    isqhmax = 1./(mesh->info.hmax*mesh->info.hmax);

    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      iadr = 6*k;
      met->m[iadr]   = MG_MAX(isqhmax,MG_MIN(isqhmin,met->m[iadr]));
      met->m[iadr+3] = met->m[iadr];
      met->m[iadr+5] = met->m[iadr];
    }
  }

  return;
}

/**
 * \param mesh pointer toward the mesh structure (unused).
 * \return -1 if fail, the number of detected ridges otherwise
 *
 * Create the boundary entities of the mesh (triangles and edges).
 *
 * \warning mesh must be packed and hashed
 */
int _MMG3D_bdryBuild(MMG5_pMesh mesh) {
  MMG5_pTetra pt;
  MMG5_hgeom  *ph;
  int         k,i,nr;


  /* rebuild triangles*/
  mesh->nt = 0;
  if ( !_MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to rebuild triangles\n",__func__);
    return(-1);
  }

  /* build hash table for edges */
  if ( mesh->htab.geom )
    _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));

  mesh->na = nr = 0;
  /* in the worst case (all edges are marked), we will have around 1 edge per *
   * triangle (we count edges only one time) */
  if ( _MMG5_hNew(mesh,&mesh->htab,mesh->nt,3*(mesh->nt)) ) {
    for (k=1; k<=mesh->ne; k++) {
      pt   = &mesh->tetra[k];
      if ( MG_EOK(pt) &&  pt->xt ) {
        for (i=0; i<6; i++) {
          if ( mesh->xtetra[pt->xt].edg[i] ||
               ( mesh->xtetra[pt->xt].tag[i] & MG_REQ ||
                 MG_EDG(mesh->xtetra[pt->xt].tag[i])) )
            if ( !_MMG5_hEdge(mesh,&mesh->htab,pt->v[_MMG5_iare[i][0]],pt->v[_MMG5_iare[i][1]],
                              mesh->xtetra[pt->xt].edg[i],mesh->xtetra[pt->xt].tag[i]))
              return -1;
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
                    printf("  ## Warning: uncomplete mesh\n"));
    }

    if ( mesh->na ) {
      _MMG5_SAFE_CALLOC(mesh->edge,mesh->na+1,MMG5_Edge,-1);

      mesh->na = 0;
      for (k=0; k<=mesh->htab.max; k++) {
        ph = &mesh->htab.geom[k];
        if ( !ph->a )  continue;
        mesh->na++;
        mesh->edge[mesh->na ].a  = ph->a;
        mesh->edge[mesh->na ].b  = ph->b;
        mesh->edge[mesh->na].tag = ( ph->tag | MG_REF ) ;
        mesh->edge[mesh->na].ref = ph->ref;

        if ( MG_GEO & ph->tag ) nr++;
      }
    }
    _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));
  }
  else
    mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));

  if ( mesh->info.imprim ) {
    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d   RIDGES  %8d\n",mesh->na,nr);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",mesh->nt);
  }

  return nr;
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
int _MMG3D_packMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol disp) {
  MMG5_pTetra   pt,ptnew;
  MMG5_pPrism   pp;
  MMG5_pQuad    pq;
  MMG5_pPoint   ppt,pptnew;
  int     np,nc,nr, k,ne,nbl,imet,imetnew,i;
  int     iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;

  /* compact vertices */
  if ( !mesh->point ) {
    fprintf(stderr, "\n  ## Error: %s: points array not allocated.\n",
            __func__);
    return 0;
  }
  if ( !mesh->tetra ) {
    fprintf(stderr, "\n  ## Error: %s: tetra array not allocated.\n",
            __func__);
    return 0;
  }

  np = nc = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->tmp = ++np;

    if ( ppt->tag & MG_NOSURF ) {
      ppt->tag &= ~MG_NOSURF;
      ppt->tag &= ~MG_REQ;
    }

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

      if ( mesh->adja ) {
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
    }
    nbl++;
  }
  mesh->ne = ne;

  /* update prisms and quads vertex indices */
  for (k=1; k<=mesh->nprism; k++) {
    pp = &mesh->prism[k];
    if ( !MG_EOK(pp) )  continue;

    pp->v[0] = mesh->point[pp->v[0]].tmp;
    pp->v[1] = mesh->point[pp->v[1]].tmp;
    pp->v[2] = mesh->point[pp->v[2]].tmp;
    pp->v[3] = mesh->point[pp->v[3]].tmp;
    pp->v[4] = mesh->point[pp->v[4]].tmp;
    pp->v[5] = mesh->point[pp->v[5]].tmp;
  }
  for (k=1; k<=mesh->nquad; k++) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) )  continue;

    pq->v[0] = mesh->point[pq->v[0]].tmp;
    pq->v[1] = mesh->point[pq->v[1]].tmp;
    pq->v[2] = mesh->point[pq->v[2]].tmp;
    pq->v[3] = mesh->point[pq->v[3]].tmp;
  }

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
  mesh->nc1 = 0;


  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;

    if ( ppt->tag & MG_BDY &&
         !(ppt->tag & MG_CRN || ppt->tag & MG_NOM || MG_EDG(ppt->tag)) ) {
      assert ( ppt->xp );

      memcpy(ppt->n,mesh->xpoint[ppt->xp].n1,3*sizeof(double));
      ++mesh->nc1;
    }

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

  /* create prism adjacency */
  if ( !MMG3D_hashPrism(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: prism hashing problem. Exit program.\n",
            __func__);
    return(0);
  }

  /* Remove the MG_REQ tags added by the nosurf option */
  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    if ( MG_EOK(pt) &&  pt->xt ) {

      for (i=0; i<6; i++) {
        if ( mesh->xtetra[pt->xt].tag[i] & MG_NOSURF ) {
          mesh->xtetra[pt->xt].tag[i] &= ~MG_REQ;
          mesh->xtetra[pt->xt].tag[i] &= ~MG_NOSURF;
        }
      }
    }
  }

  if ( mesh->info.imprim ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",mesh->np,nc);
    fprintf(stdout,"     NUMBER OF TETRAHEDRA %8d\n",mesh->ne);
  }

  nr = _MMG3D_bdryBuild(mesh);
  if ( nr < 0 ) return 0;

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
    fprintf(stderr,"\n  ##  Warning: %s: invalid mesh.\n",__func__);
    return(0);
  }

  return(1);
}

int MMG3D_mmg3dlib(MMG5_pMesh mesh,MMG5_pSol met) {
  mytime    ctim[TIMEMAX];
  char      stim[32];

  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",MG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
  }

  _MMG3D_Set_commonFunc();


  _MMG5_warnOrientation(mesh);

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
    fprintf(stderr,"\n  ## ERROR: LAGRANGIAN MODE UNAVAILABLE (MMG3D_IPARAM_lag):\n"
            "            YOU MUST CALL THE MMG3D_MMG3DMOV FUNCTION TO MOVE A RIGIDBODY.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.iso ) {
    fprintf(stderr,"\n  ## ERROR: LEVEL-SET DISCRETISATION UNAVAILABLE"
            " (MMG3D_IPARAM_iso):\n"
            "          YOU MUST CALL THE MMG3D_MMG3DMOV FUNCTION TO USE THIS OPTION.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.optimLES && met->size==6 ) {
    fprintf(stdout,"\n  ## ERROR: STRONG MESH OPTIMIZATION FOR LES METHODS"
            " UNAVAILABLE (MMG3D_IPARAM_optimLES) WITH AN ANISOTROPIC METRIC.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");

  chrono(ON,&(ctim[1]));

  /* check input */
  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"\n  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));
    met->np = 0;
  }
  else if ( met->size!=1 && met->size!=6 ) {
    fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  /* specific meshing */
  if ( met->np ) {
    if ( mesh->info.optim ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: OPTIM OPTION CAN NOT BE USED"
             " WITH AN INPUT METRIC.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }

    if ( mesh->info.hsiz>0. ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ OPTION CAN NOT BE USED"
             " WITH AN INPUT METRIC.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
  }

  if ( mesh->info.optim &&  mesh->info.hsiz>0. ) {
    printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ AND OPTIM OPTIONS CAN NOT BE USED"
           " TOGETHER.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

#ifdef USE_SCOTCH
  _MMG5_warnScotch(mesh);
#endif

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  /* scaling mesh */
  if ( !_MMG5_scaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  /* specific meshing */
  if ( mesh->info.optim ) {
    if ( !MMG3D_doSol(mesh,met) ) {
      if ( !_MMG5_unscaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
    }
    _MMG3D_solTruncatureForOptim(mesh,met);
  }

  if ( mesh->info.hsiz > 0. ) {
    if ( !MMG3D_Set_constantSize(mesh,met) ) {
     if ( !_MMG5_unscaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
     _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
  }

  MMG3D_setfunc(mesh,met);

  if ( !_MMG3D_tetraQual(mesh,met,0) ) _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);

  if ( abs(mesh->info.imprim) > 0 ) {
    if ( !_MMG3D_inqua(mesh,met) ) {
      if ( !_MMG5_unscaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
    }
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
      fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }
#else
  if ( !_MMG5_mmg3d1_delone(mesh,met) ) {
    if ( (!mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
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
    fprintf(stderr,"\n  ## ERROR: LAGRANGIAN MODE UNAVAILABLE (MMG3D_IPARAM_lag):\n"
            "            YOU MUST CALL THE MMG3D_MMG3DMOV FUNCTION TO MOVE A RIGIDBODY.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.optimLES ) {
    fprintf(stdout,"\n  ## ERROR: STRONG MESH OPTIMIZATION FOR LES METHODS"
            " UNAVAILABLE (MMG3D_IPARAM_optimLES) IN ISOSURFACE"
            " DISCRETIZATION MODE.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  if ( mesh->info.optim ) {
    printf("\n  ## ERROR: OPTIM OPTION UNAVAILABLE IN ISOSURFACE"
           " DISCRETIZATION MODE\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  if ( mesh->info.hsiz>0. ) {
    printf("\n  ## ERROR: HSIZ OPTION UNAVAILABLE IN ISOSURFACE"
           " DISCRETIZATION MODE.\n");
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
    fprintf(stdout,"\n  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));
    met->np = 0;
  }
  else if ( met->size!=1 ) {
    fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  chrono(ON,&(ctim[2]));

  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ISOSURFACE DISCRETIZATION\n");
  }

  /* scaling mesh */
  if ( !_MMG5_scaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  MMG3D_setfunc(mesh,met);

  if ( !_MMG3D_tetraQual(mesh,met,0) ) _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);

  if ( abs(mesh->info.imprim) > 0 ) {
    if ( !_MMG3D_inqua(mesh,met) ) {
      if ( !_MMG5_unscaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
    }
  }

  /* specific meshing */
  if ( !met->np ) {
    fprintf(stderr,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  if ( !_MMG3D_mmg3d2(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 2 : ANALYSIS\n");
  }

  /* mesh analysis */
  if ( !_MMG3D_analys(mesh) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[4]));
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 3 : MESH IMPROVEMENT\n");
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
      fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }
#else
  if ( !_MMG5_mmg3d1_pattern(mesh,met) ) {
    if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMG5_RETURN_AND_PACK(mesh,met,NULL,MMG5_LOWFAILURE);
  }
#endif

  chrono(OFF,&(ctim[4]));
  printim(ctim[4].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 3 COMPLETED.     %s\n",stim);
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
    fprintf(stderr,"\n  ## ERROR: LEVEL-SET DISCRETISATION UNAVAILABLE"
            " (MMG3D_IPARAM_iso):\n"
            "          YOU MUST CALL THE MMG3D_mmg3dmov FUNCTION TO USE THIS OPTION.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.optimLES ) {
    fprintf(stdout,"\n  ## ERROR: STRONG MESH OPTIMIZATION FOR LES METHODS"
            " UNAVAILABLE (MMG3D_IPARAM_optimLES) IN LAGRANGIAN MODE.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  if ( mesh->info.optim ) {
    printf("\n  ## ERROR: OPTIM OPTION UNAVAILABLE IN LAGRANGIAN MODE\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  if ( mesh->info.hsiz>0. ) {
    printf("\n  ## ERROR: HSIZ OPTION UNAVAILABLE IN LAGRANGIAN MODE\n");
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
      fprintf(stdout,"\n  ## Warning: displacement mode for the rigidbody"
              " movement is not set.\n"
              "               Lagrangian displacement computed according"
              " to mode 1.\n");
    mesh->info.lag = 1;
  }

#ifndef USE_ELAS
  fprintf(stderr,"\n  ## ERROR: YOU NEED TO COMPILE WITH THE USE_ELAS"
    " CMake's FLAG SET TO ON TO USE THE RIGIDBODY MOVEMENT LIBRARY.\n");
  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
#endif

  if ( !disp ) {
    fprintf(stderr,"\n  ## ERROR: IN LAGRANGIAN MODE, A STRUCTURE OF TYPE"
            " \"MMG5_pSoL\" IS NEEDED TO STORE THE DISPLACEMENT FIELD.\n"
            "            THIS STRUCTURE MUST BE DIFFERENT FROM THE ONE USED"
            " TO STORE THE METRIC.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  if (disp->np && (disp->np != mesh->np) ) {
    fprintf(stdout,"\n  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,disp->m,(disp->size*(disp->npmax+1))*sizeof(double));
    disp->np = 0;
  }
  else if (disp->size!=3) {
    fprintf(stderr,"\n  ## ERROR: LAGRANGIAN MOTION OPTION NEED A VECTORIAL DISPLACEMENT\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));

  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  /* scaling mesh */
  if ( !_MMG5_scaleMesh(mesh,disp) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  MMG3D_setfunc(mesh,met);

  if ( !_MMG3D_tetraQual(mesh,met,0) ) _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);

  if ( abs(mesh->info.imprim) > 0 ) {
    if ( !_MMG3D_inqua(mesh,met) ) {
      _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
    }
  }

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

#ifdef USE_ELAS
  /* Lagrangian mode */
  if ( !_MMG5_mmg3d3(mesh,disp,met) ) {
    disp->npi = disp->np;
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
#endif
  disp->npi = disp->np;

  if ( mesh->info.optim ) {
    if ( !MMG3D_doSol(mesh,met) ) {
      if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
      _MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
    }
    _MMG3D_solTruncatureForOptim(mesh,met);
  }

/* ******************* Add mesh improvement ? *************************** */

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
