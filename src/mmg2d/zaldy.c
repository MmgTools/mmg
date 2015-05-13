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


/* get new point address */
int MMG2_newPt(MMG5_pMesh mesh,double c[2]) {
  MMG5_pPoint  ppt;
  int     curpt;

  if ( !mesh->npnil )  return(0);

  curpt = mesh->npnil;
  if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
  ppt   = &mesh->point[curpt];
  memcpy(ppt->c,c,2*sizeof(double));
  ppt->tag   &= ~M_NUL;
  mesh->npnil = ppt->tmp;
  ppt->tmp    = 0;
  ppt->xp     = 0;
  //ppt->fla   = mesh->flag;

  return(curpt);
}


void MMG2_delPt(MMG5_pMesh mesh,int ip) {
  MMG5_pPoint   ppt;
  MMG5_pxPoint  pxp;

  ppt = &mesh->point[ip];
  if ( ppt->xp ) {
    pxp = &mesh->xpoint[ppt->xp];
    memset(pxp,0,sizeof(MMG5_xPoint));
  }

  memset(ppt,0,sizeof(MMG5_Point));
  ppt->tag    = M_NUL;
  ppt->tmp    = mesh->npnil; 
  
  mesh->npnil = ip;
  if ( ip == mesh->np )  mesh->np--;
}

/* get new elt address */
int MMG2_newEdge(MMG5_pMesh mesh) {
  int     curiel;

  if ( !mesh->nanil ) {
    fprintf(stdout,"  ## UNABLE TO ALLOCATE NEW ELEMENT.\n");
    return(0);
  }
  curiel = mesh->nanil;
  if ( mesh->nanil > mesh->na )  mesh->na = mesh->nanil;
  mesh->nanil = mesh->edge[curiel].b;
  mesh->edge[curiel].b = 0;

  return(curiel);
}


void MMG2_delEdge(MMG5_pMesh mesh,int iel) {
  MMG5_pEdge    pt;

  pt = &mesh->edge[iel];
  if ( !pt->a ) {
    fprintf(stdout,"  ## INVALID EDGE.\n");
    return;
  }
  memset(pt,0,sizeof(MMG5_Edge));
  pt->b = mesh->nanil;
  mesh->nanil = iel;
  if ( iel == mesh->na )  mesh->na--;
}

/* get new elt address */
int MMG2_newElt(MMG5_pMesh mesh) {
  int     curiel;

  if ( !mesh->nenil ) {
    fprintf(stdout,"  ## UNABLE TO ALLOCATE NEW ELEMENT.\n");
    return(0);
  }
  curiel = mesh->nenil;
  if ( mesh->nenil > mesh->nt )  mesh->nt = mesh->nenil;
  mesh->nenil = mesh->tria[curiel].v[2];
  mesh->tria[curiel].v[2] = 0;
  mesh->tria[curiel].ref = 0;
  mesh->tria[curiel].base = 0;
  

  return(curiel);
}


void MMG2_delElt(MMG5_pMesh mesh,int iel) {
  MMG5_pTria    pt;
  int      iadr;

  pt = &mesh->tria[iel];
  if ( !M_EOK(pt) ) {
    fprintf(stdout,"  ## INVALID ELEMENT.\n");
    return;
  }
  memset(pt,0,sizeof(MMG5_Tria));
  pt->v[2] = mesh->nenil;
  pt->qual = 0.0;
  iadr = (iel-1)*3 + 1;
  memset(&mesh->adja[iadr],0,3*sizeof(int));

  mesh->nenil = iel;
  if ( iel == mesh->nt )  mesh->nt--;
}


/* check if n elets available */
int MMG2_getnElt(MMG5_pMesh mesh,int n) {
  int     curiel;

  if ( !mesh->nenil )  return(0);
  curiel = mesh->nenil;
  do {
    curiel = mesh->tria[curiel].v[2];
  }
  while (--n);

  return(n == 0);
}


/* allocate main structure */
int MMG2_zaldy(MMG5_pMesh mesh) {
  int     million = 1048576L;
  int     k,npask;

  if ( mesh->info.mem < 0 ) {
    mesh->npmax  = M_MAX(1.5*mesh->np,NPMAX);
    mesh->xpmax  = M_MAX(0.1*mesh->xp,0.1*NPMAX);
    mesh->namax = M_MAX(1.5*mesh->na,NEDMAX);
    mesh->ntmax  = M_MAX(1.5*mesh->nt,NEMAX);
  }
  else {
    /* point+tria+adja+sol+bucket+queue */
    int bytes = sizeof(MMG5_Point) +  0.1*sizeof(MMG5_xPoint) + 2*sizeof(MMG5_Tria) + 3*sizeof(int)
      + sizeof(MMG5_Sol) /*+ sizeof(Displ)*/
                + sizeof(int) + 5*sizeof(int);

    npask = (double)mesh->info.mem / bytes * million;
    mesh->npmax = M_MAX(1.5*mesh->np,npask);
    mesh->xpmax = M_MAX(0.1*mesh->xp,0.1*npask);
    mesh->namax = M_MAX(1.5*mesh->na,2*npask);
    mesh->ntmax = M_MAX(1.5*mesh->nt,2*npask);
  }
  mesh->point = (MMG5_pPoint)M_calloc(mesh->npmax+1,sizeof(MMG5_Point),"zaldy.point");
  assert(mesh->point);
  mesh->xpoint = (MMG5_pxPoint)M_calloc(mesh->npmax+1,sizeof(MMG5_xPoint),"zaldy.xpoint");
  assert(mesh->point);
  mesh->tria  = (MMG5_pTria)M_calloc(mesh->ntmax+1,sizeof(MMG5_Tria),"zaldy.tria");
  assert(mesh->tria);
  mesh->edge  = (MMG5_pEdge)M_calloc(mesh->namax+1,sizeof(MMG5_Edge),"zaldy.edge");
  assert(mesh->edge);
  mesh->adja = (int*)calloc(3*mesh->ntmax+5,sizeof(int));
  assert(mesh->adja);

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nanil = mesh->na + 1;
  mesh->nenil = mesh->nt + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  for (k=mesh->nanil; k<mesh->namax-1; k++)
    mesh->edge[k].b = k+1;

  for (k=mesh->nenil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;


  return(1);
}

