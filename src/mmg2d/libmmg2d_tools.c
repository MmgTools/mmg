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
 * \file mmg2d/libmmg2d_tools.c
 * \brief Tools functions for the mmg2d library.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/

#include "mmg2d.h"

void MMG2D_setfunc(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( met->size == 3 ) {
    MMG2D_lencurv  = _MMG2_lencurv_ani;
    MMG2D_defsiz   = _MMG2_defsiz_ani;
    MMG2D_gradsiz  = _MMG2_gradsiz_ani;
    MMG2D_caltri   = _MMG2_caltri_ani;
    MMG2D_intmet   = _MMG2_intmet_ani;
    //    MMG2_optlen    = optlen_ani;
  }
  else {
    MMG2D_lencurv   = _MMG2_lencurv_iso;
    MMG2D_defsiz    = _MMG2_defsiz_iso;
    MMG2D_gradsiz   = _MMG2_gradsiz_iso;
    MMG2D_caltri    = _MMG2_caltri_iso;
    MMG2D_intmet    = _MMG2_intmet_iso;
  }
  return;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return 1 if success, 0 if fail
 *
 * Read parameter file DEFAULT.mmg2d
 *
 */
int MMG2_parsop(MMG5_pMesh mesh,MMG5_pSol met) {
  int       ret,i;
  char     *ptr,data[256];
  FILE     *in;
  MMG5_pMat pm;
  fpos_t    position;
  
  /* Check for parameter file */
  strcpy(data,mesh->namein);
  ptr = strstr(data,".mesh");
  if ( ptr ) *ptr = '\0';
  strcat(data,".mmg2d");
  in = fopen(data,"rb");
  
  if ( !in ) {
    sprintf(data,"%s","DEFAULT.mmg2d");
    in = fopen(data,"rb");
    if ( !in )
      return(1);
  }
  fprintf(stdout,"\n  %%%% %s OPENED\n",data);
  
  /* Read parameters */
  while ( !feof(in) ) {
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) ) break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);
    
    /* Read user defined references for the LS mode */
    if ( !strcmp(data,"lsreferences") ) {
      ret = fscanf(in,"%d",&mesh->info.nmat);
      
      if ( !ret ) {
        fprintf(stderr,"  %%%% Wrong format: %d\n",mesh->info.nmat);
        return 0;
      }

      if ( mesh->info.nmat ) {
        _MMG5_SAFE_CALLOC(mesh->info.mat,mesh->info.nmat,MMG5_Mat,0);
        for (i=0; i<mesh->info.nmat; i++) {
          pm = &mesh->info.mat[i];
          fscanf(in,"%d",&pm->ref);
          fgetpos(in,&position);
          fscanf(in,"%s",data);
          if ( !strcmp(data,"nosplit") ) {
            pm->dospl = 0;
            pm->rin = pm->ref;
            pm->rex = pm->ref;
          }
          else {
            fsetpos(in,&position);
            fscanf(in,"%d",&pm->rin);
            fscanf(in,"%d",&pm->rex);
            pm->dospl = 1;
          }
        }
      }
    }
    else {
      fprintf(stderr,"  %%%% Wrong format: %s\n",data);
      return 0;
    }
  }
  
  fclose(in);
  return(1);
}

int MMG2D_Get_adjaTri(MMG5_pMesh mesh, int kel, int listri[3]) {

  if ( ! mesh->adja ) {
    if (! MMG2_hashTria(mesh))
      return(0);
  }

  listri[0] = mesh->adja[3*(kel-1)+1]/3;
  listri[1] = mesh->adja[3*(kel-1)+2]/3;
  listri[2] = mesh->adja[3*(kel-1)+3]/3;

  return(1);
}

inline
int MMG2D_Get_adjaVertices(MMG5_pMesh mesh, int ip, int lispoi[MMG2D_LMAX])
{
  int start;

  if ( !mesh->tria ) return 0;

  start=MMG2_findTria(mesh,ip);
  if ( !start ) return 0;

  return MMG2D_Get_adjaVerticesFast(mesh,ip,start,lispoi);
}

inline
int MMG2D_Get_adjaVerticesFast(MMG5_pMesh mesh, int ip,int start, int lispoi[MMG2D_LMAX])
{
  MMG5_pTria pt;
  int k,prevk,nbpoi,iploc,i,i1,i2,*adja;

  pt   = &mesh->tria[start];

  for ( iploc=0; iploc<3; ++iploc ) {
    if ( pt->v[iploc] == ip ) break;
  }

  assert(iploc!=3);

  k = start;
  i = iploc;
  nbpoi = 0;
  do {
    if ( nbpoi == MMG2D_LMAX ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to compute adjacent"
              " vertices of the vertex %d:\nthe ball of point contain too many"
              " elements.\n",__func__,ip);
      return(0);
    }
    i1 = _MMG5_inxt2[i];
    lispoi[nbpoi] = mesh->tria[k].v[i1];
    ++nbpoi;

    adja = &mesh->adja[3*(k-1)+1];
    prevk = k;
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = _MMG5_inxt2[i];
  }
  while ( k && k != start );

  if ( k > 0 ) return(nbpoi);

  /* store the last point of the boundary triangle */
  if ( nbpoi == MMG2D_LMAX ) {
    fprintf(stderr,"\n  ## Warning: %s: unable to compute adjacent vertices of the"
            " vertex %d:\nthe ball of point contain too many elements.\n",
            __func__,ip);
    return(0);
  }
  i1 = _MMG5_inxt2[i1];
  lispoi[nbpoi] = mesh->tria[prevk].v[i1];
  ++nbpoi;

  /* check if boundary hit */
  k = start;
  i = iploc;
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i2 = _MMG5_iprv2[i];
    k  = adja[i2] / 3;
    if ( k == 0 )  break;

    if ( nbpoi == MMG2D_LMAX ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to compute adjacent vertices of the"
              " vertex %d:\nthe ball of point contain too many elements.\n",
              __func__,ip);
      return(0);
    }
    i  = adja[i2] % 3;
    lispoi[nbpoi] = mesh->tria[k].v[i];
    ++nbpoi;

    i  = _MMG5_iprv2[i];
  }
  while ( k );

  return nbpoi;
}

int MMG2D_Get_triFromEdge(MMG5_pMesh mesh, int ked, int *ktri, int *ied)
{
  int val;

  val = mesh->edge[ked].base;

  if ( !val ) return(0);

  *ktri = val/3;

  *ied = val%3;

  return 1;


}

int MMG2D_Set_constantSize(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pPoint ppt;
  double      hsiz;
  int         k,iadr;

  met->np = mesh->np;

  if ( !MMG5_Compute_constantSize(mesh,met,&hsiz) )
    return 0;

  if ( met->size == 1 ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      met->m[k] = hsiz;
    }
  }
  else {
    hsiz    = 1./(hsiz*hsiz);

    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;

      iadr           = met->size*k;
      met->m[iadr]   = hsiz;
      met->m[iadr+2] = hsiz;
    }
  }
  return 1;
}



void MMG2D_Reset_verticestags(MMG5_pMesh mesh) {
  int k;

  for ( k=1; k<=mesh->np;  ++k ) {
    mesh->point[k].tag = 0;
  }

}

void MMG2D_Free_triangles(MMG5_pMesh mesh) {

  if ( mesh->adja )
    _MMG5_DEL_MEM(mesh,mesh->adja,(3*mesh->ntmax+5)*sizeof(int));

  if ( mesh->tria )
    _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->ntmax+1)*sizeof(MMG5_Tria));

  mesh->nt = 0;
  mesh->nti = 0;
  mesh->nenil = 0;

  return;
}

void MMG2D_Free_edges(MMG5_pMesh mesh) {

  if ( mesh->edge )
    _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->namax+1)*sizeof(MMG5_Edge));

  if ( mesh->xpoint )
    _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));

  mesh->na = 0;
  mesh->nai = 0;
  mesh->nanil = 0;

  mesh->xp = 0;

  return;
}

void MMG2D_Free_solutions(MMG5_pMesh mesh,MMG5_pSol sol) {

  /* sol */
  if ( sol && sol->m )
    _MMG5_DEL_MEM(mesh,sol->m,(sol->size*(sol->npmax+1))*sizeof(double));

  return;
}
