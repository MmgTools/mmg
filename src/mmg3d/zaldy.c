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
 * \file mmg3d/zaldy.c
 * \brief Memory management.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmg3d.h"

/** get new point address */
int _MMG5_newPt(MMG5_pMesh mesh,double c[3],char tag) {
  MMG5_pPoint  ppt;
  int     curpt;

  if ( !mesh->npnil )  return(0);
  curpt = mesh->npnil;
  if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
  ppt   = &mesh->point[curpt];
  memcpy(ppt->c,c,3*sizeof(double));
  mesh->npnil = ppt->tmp;
  ppt->tmp    = 0;

  ppt->ref = 0;
  ppt->xp = 0;
  ppt->flag = 0;
  /* point on geometry */
  if ( tag & MG_BDY ) {
    mesh->xp++;
    if(mesh->xp > mesh->xpmax){
      /* reallocation of xpoint table */
      _MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                         "larger xpoint table",
                         return(0));
    }
    ppt->xp  = mesh->xp;
  }
  assert(tag < 127);
  assert(tag >= 0);
  ppt->n   = 0;
  ppt->tag = tag;
  ppt->tagdel = 0;
  return(curpt);
}

void _MMG5_delPt(MMG5_pMesh mesh,int ip) {
  MMG5_pPoint   ppt;
  MMG5_xPoint  *pxp;

  ppt = &mesh->point[ip];
  if ( ppt->xp ) {
    pxp = &mesh->xpoint[ppt->xp];
    memset(pxp,0,sizeof(MMG5_xPoint));
  }
  memset(ppt,0,sizeof(MMG5_Point));
  ppt->tag    = MG_NUL;
  ppt->tmp    = mesh->npnil;
  mesh->npnil = ip;
  if ( ip == mesh->np ) {
    while ( !MG_VOK((&mesh->point[mesh->np])) )  mesh->np--;
  }
}

/** get new elt address */
int _MMG5_newElt(MMG5_pMesh mesh) {
  int     curiel;

  if ( !mesh->nenil )  return(0);
  curiel = mesh->nenil;

  if ( mesh->nenil > mesh->ne )  mesh->ne = mesh->nenil;
  mesh->nenil = mesh->tetra[curiel].v[3];
  mesh->tetra[curiel].v[3] = 0;
  mesh->tetra[curiel].mark=0;

  return(curiel);
}


void _MMG5_delElt(MMG5_pMesh mesh,int iel) {
  MMG5_pTetra   pt;
  int      iadr;

  pt = &mesh->tetra[iel];
  if ( !MG_EOK(pt) ) {
    fprintf(stdout,"  ## INVALID ELEMENT %d.\n",iel);
    exit(EXIT_FAILURE);
  }
  memset(pt,0,sizeof(MMG5_Tetra));
  pt->v[3] = mesh->nenil;
  iadr = 4*(iel-1) + 1;
  if ( mesh->adja )
    memset(&mesh->adja[iadr],0,4*sizeof(int));
  mesh->nenil = iel;
  if ( iel == mesh->ne ) {
    while ( !MG_EOK((&mesh->tetra[mesh->ne])) )  mesh->ne--;
  }
}

long long _MMG5_memSize (void) {
  long long mem;

#if (defined(__APPLE__) && defined(__MACH__))
  size_t size;

  size = sizeof(mem);
  if ( sysctlbyname("hw.memsize",&mem,&size,NULL,0) == -1)
    return(0);

#elif defined(__unix__) || defined(__unix) || defined(unix)
  mem = ((long long)sysconf(_SC_PHYS_PAGES))*
    ((long long)sysconf(_SC_PAGE_SIZE));
#else
  printf("  ## WARNING: UNKNOWN SYSTEM, RECOVER OF MAXIMAL MEMORY NOT AVAILABLE.\n");
  return(0);
#endif

  return(mem);
}

/** memory repartition for the -m option */
void _MMG5_memOption(MMG5_pMesh mesh) {
  long long  million = 1048576L;
  int        ctri,npask,bytes,memtmp;

  mesh->memMax = _MMG5_memSize();

  mesh->npmax = MG_MAX(1.5*mesh->np,_MMG5_NPMAX);
  mesh->nemax = MG_MAX(1.5*mesh->ne,_MMG5_NEMAX);
  mesh->ntmax = MG_MAX(1.5*mesh->nt,_MMG5_NTMAX);

  if ( mesh->info.mem <= 0 ) {
    if ( mesh->memMax )
      /* maximal memory = 50% of total physical memory */
      mesh->memMax = mesh->memMax*50/100;
    else {
      /* default value = 800 Mo */
      printf("  Maximum memory set to default value: %d Mo.\n",_MMG5_MEMMAX);
      mesh->memMax = _MMG5_MEMMAX*million;
    }
  }
  else {
    /* memory asked by user if possible, otherwise total physical memory */
    if ( (long long)mesh->info.mem*million > mesh->memMax && mesh->memMax ) {
      fprintf(stdout,"  ## Warning: asking for %d Mo of memory ",mesh->info.mem);
      fprintf(stdout,"when only %lld available.\n",mesh->memMax/million);
    }
    else {
      mesh->memMax= (long long)(mesh->info.mem)*million;
    }

    /* if asked memory is lower than default _MMG5_NPMAX/_MMG5_NEMAX/_MMG5_NTMAX we take lower values */
    ctri = 2;

    /* Euler-poincare: ne = 6*np; nt = 2*np; na = np/5 *
     * point+tria+tets+adja+adjt+sol+item *
     * warning: we exceed memory in saveMesh when we call _MMG5_hNew */
    bytes = sizeof(MMG5_Point) + sizeof(MMG5_xPoint) +
      6*sizeof(MMG5_Tetra) + ctri*sizeof(MMG5_xTetra) +
      4*6*sizeof(int) + ctri*3*sizeof(int) +
      sizeof(MMG5_Sol)+4*sizeof(_MMG5_hedge);
#ifdef USE_SCOTCH
    /* bytes = bytes + vertTab + edgeTab + PermVrtTab *
     * + vertOldTab + sortPartTab - adja */
    bytes = bytes + 3*6*sizeof(int);
#endif

    /*init allocation need 38Mo*/
    npask = (double)(mesh->info.mem-38) / bytes * (int)million;
    mesh->npmax = MG_MIN(npask,mesh->npmax);
    mesh->ntmax = MG_MIN(ctri*npask,mesh->ntmax);
    mesh->nemax = MG_MIN(6*npask,mesh->nemax);
    /*check if the memory asked is enough to load the mesh*/
    if(mesh->np &&
       (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne)) {
      memtmp = mesh->np * bytes /(int)million + 38;
      memtmp = MG_MAX(memtmp, mesh->nt * bytes /(ctri* (int)million) + 38);
      memtmp = MG_MAX(memtmp, mesh->ne * bytes /(6*(int)million) + 38);
      mesh->memMax = (long long) memtmp+1;
      fprintf(stdout,"  ## ERROR: asking for %d Mo of memory ",mesh->info.mem);
      fprintf(stdout,"is not enough to load mesh. You need to ask %d Mo minimum\n",
              memtmp+1);
    }
    if(mesh->info.mem < 39) {
      mesh->memMax = (long long) 39;
      fprintf(stdout,"  ## ERROR: asking for %d Mo of memory ",mesh->info.mem);
      fprintf(stdout,"is not enough to load mesh. You need to ask %d Mo minimum\n",
              39);
    }
  }

  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED (Mo)    %lld\n",
            mesh->memMax/million);

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  _MMG5_NPMAX    %d\n",mesh->npmax);
    fprintf(stdout,"  _MMG5_NTMAX    %d\n",mesh->ntmax);
    fprintf(stdout,"  _MMG5_NEMAX    %d\n",mesh->nemax);
  }

  return;
}

/** allocate main structure */
int _MMG5_zaldy(MMG5_pMesh mesh) {
  int     k;

  _MMG5_memOption(mesh);

  _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point);

  _MMG5_ADD_MEM(mesh,(mesh->nemax+1)*sizeof(MMG5_Tetra),"initial tetrahedra",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->tetra,mesh->nemax+1,MMG5_Tetra);

  if ( mesh->nt ) {
    _MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(MMG5_Tria),"initial triangles",return(0));
    _MMG5_SAFE_CALLOC(mesh->tria,mesh->nt+1,MMG5_Tria);
    memset(&mesh->tria[0],0,sizeof(MMG5_Tria));
  }
  if ( mesh->na ) {
    _MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"initial edges",return(0));
    _MMG5_SAFE_CALLOC(mesh->edge,(mesh->na+1),MMG5_Edge);
  }

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->ne + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++) {
    /* Set tangent field of point to 0 */
    mesh->point[k].n[0] = 0;
    mesh->point[k].n[1] = 0;
    mesh->point[k].n[2] = 0;
    /* link */
    mesh->point[k].tmp  = k+1;
  }

  for (k=mesh->nenil; k<mesh->nemax-1; k++)
    mesh->tetra[k].v[3] = k+1;

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Free xtetra structure.
 *
 */
void _MMG5_freeXTets(MMG5_pMesh mesh) {
  MMG5_pTetra pt;
  int    k;

  for (k=1; k<=mesh->ne; k++) {
    pt     = &mesh->tetra[k];
    pt->xt = 0;
  }
  if ( mesh->xtetra )
    _MMG5_DEL_MEM(mesh,mesh->xtetra,(mesh->xtmax+1)*sizeof(MMG5_xTetra));
  mesh->xt = 0;
}
