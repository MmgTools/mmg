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
 * \file mmgs/zaldy_s.c
 * \brief Memory management.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmgs.h"

/* get new point address */
int _MMGS_newPt(MMG5_pMesh mesh,double c[3],double n[3]) {
  MMG5_pPoint  ppt;
  int     curpt;

  if ( !mesh->npnil )  return(0);

  curpt = mesh->npnil;
  if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
  ppt   = &mesh->point[curpt];
  memcpy(ppt->c,c,3*sizeof(double));
  if ( n )
    memcpy(ppt->n,n,3*sizeof(double));
  ppt->tag   &= ~MG_NUL;
  mesh->npnil = ppt->tmp;
  ppt->tmp    = 0;

  return(curpt);
}

void _MMGS_delPt(MMG5_pMesh mesh,int ip) {
  MMG5_pPoint   ppt;

  ppt = &mesh->point[ip];
  memset(ppt,0,sizeof(MMG5_Point));
  ppt->tag    = MG_NUL;
  ppt->tmp    = mesh->npnil;
  mesh->npnil = ip;
  if ( ip == mesh->np ) {
    while ( !MG_VOK((&mesh->point[mesh->np])) )  mesh->np--;
  }
}

int _MMGS_newElt(MMG5_pMesh mesh) {
  int     curiel;

  if ( !mesh->nenil )  return(0);
  curiel = mesh->nenil;

  if ( mesh->nenil > mesh->nt )  mesh->nt = mesh->nenil;
  mesh->nenil = mesh->tria[curiel].v[2];
  mesh->tria[curiel].v[2] = 0;

  return(curiel);
}

/**
 * \param mesh pointer toward the mesh
 * \param iel index of the element to delete
 *
 * \return 1 if success, 0 if fail
 *
 * Delete the element \a iel
 *
 */
int _MMGS_delElt(MMG5_pMesh mesh,int iel) {
  MMG5_pTria    pt;

  pt = &mesh->tria[iel];
  if ( !MG_EOK(pt) ) {
    fprintf(stderr,"\n  ## INVALID ELEMENT %d.\n",iel);
    return 0;
  }
  memset(pt,0,sizeof(MMG5_Tria));
  pt->v[2] = mesh->nenil;
  if ( mesh->adja )
    memset(&mesh->adja[3*(iel-1)+1],0,3*sizeof(int));
  mesh->nenil = iel;
  if ( iel == mesh->nt ) {
    while ( !MG_EOK((&mesh->tria[mesh->nt])) )  mesh->nt--;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * memory repartition for the -m option with initial values of memMax,npmax,
 * and ntmax.
 *
 */
static inline
int _MMGS_memOption_memRepartition(MMG5_pMesh mesh) {
  long long  million = 1048576L;
  long long  usedMem,avMem,reservedMem;
  long       castedVal;
  int        npadd,bytes;

  if ( mesh->info.mem <= 0 ) {
    if ( mesh->memMax )
      /* maximal memory = 50% of total physical memory */
      mesh->memMax = (long long)(mesh->memMax*_MMG5_MEMPERCENT);
    else {
      /* default value = 800 Mo */
      printf("  Maximum memory set to default value: %d Mo.\n",_MMG5_MEMMAX);
      mesh->memMax = _MMG5_MEMMAX*million;
    }
  }
  else {
    /* memory asked by user if possible, otherwise total physical memory */
    if ( (long long)(mesh->info.mem)*million > mesh->memMax && mesh->memMax ) {
      fprintf(stderr,"\n  ## Warning: %s: asking for %d Mo of memory ",
              __func__,mesh->info.mem);
      castedVal = _MMG5_SAFELL2LCAST(mesh->memMax/million);
      fprintf(stderr,"when only %ld available.\n",castedVal);
    }
    else {
      mesh->memMax= (long long)(mesh->info.mem)*million;
    }
  }

  /* init allocation need 38 Mo */
  reservedMem = 38*million;

  /* Compute the needed initial memory */
  usedMem = reservedMem + (mesh->np+1)*sizeof(MMG5_Point)
    + (mesh->nt+1)*sizeof(MMG5_Tria) + (3*mesh->nt+1)*sizeof(int)
    + (mesh->np+1)*sizeof(double);

  if ( usedMem > mesh->memMax  ) {
    fprintf(stderr,"\n  ## Error: %s: %lld Mo of memory ",__func__,mesh->memMax/million);
    castedVal =  _MMG5_SAFELL2LCAST(usedMem/million+1);
    fprintf(stderr,"is not enough to load mesh. You need to ask %ld Mo minimum\n",
            castedVal);
    return 0;
  }

  /* point+tria+adja */
  bytes = sizeof(MMG5_Point) + sizeof(MMG5_xPoint) +
    2*sizeof(MMG5_Tria) + 3*sizeof(int) + sizeof(MMG5_Sol);

  avMem = mesh->memMax-usedMem;

  npadd = (int) ( (double)avMem/bytes );
  mesh->npmax = MG_MIN(mesh->npmax,mesh->np+npadd);
  mesh->ntmax = MG_MIN(mesh->ntmax,2*npadd+mesh->nt);

  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug ) {
    castedVal = _MMG5_SAFELL2LCAST(mesh->memMax/million);
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED (Mo)    %ld\n",castedVal);
  }

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  _MMG2D_NPMAX    %d\n",mesh->npmax);
    fprintf(stdout,"  _MMG2D_NTMAX    %d\n",mesh->ntmax);
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * memory repartition for the -m option
 *
 */
int _MMGS_memOption(MMG5_pMesh mesh) {

  mesh->memMax = _MMG5_memSize();

  mesh->npmax = MG_MAX(1.5*mesh->np,_MMGS_NPMAX);
  mesh->ntmax = MG_MAX(1.5*mesh->nt,_MMGS_NTMAX);

  return ( _MMGS_memOption_memRepartition(mesh) );
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Allocation of the array fields of the mesh.
 *
 */
int MMGS_setMeshSize_alloc( MMG5_pMesh mesh ) {
  int k;

  _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point,0);
  _MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),"initial triangles",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  _MMG5_SAFE_CALLOC(mesh->tria,mesh->ntmax+1,MMG5_Tria,0);


  mesh->namax = mesh->na;
  if ( mesh->na ) {
    _MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"initial edges",return(0));
    _MMG5_SAFE_CALLOC(mesh->edge,(mesh->na+1),MMG5_Edge,0);
  }

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->nt + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  for (k=mesh->nenil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 *
 * \return 1 if success, 0 if fail
 *
 * allocate main structure
 *
 */
int _MMGS_zaldy(MMG5_pMesh mesh) {

  if ( !_MMGS_memOption(mesh) )  return 0;

  return ( MMGS_setMeshSize_alloc(mesh) );
}
