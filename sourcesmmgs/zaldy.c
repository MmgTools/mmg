#include "mmgs.h"

extern Info  info;


/* get new point address */
int newPt(pMesh mesh,double c[3],double n[3]) {
  pPoint  ppt;
  int     curpt;

  if ( !mesh->npnil )  return(0);

  curpt = mesh->npnil;
  if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
  ppt   = &mesh->point[curpt];
  memcpy(ppt->c,c,3*sizeof(double));
	memcpy(ppt->n,n,3*sizeof(double));
  ppt->tag   &= ~MS_NUL;
  mesh->npnil = ppt->tmp;
  ppt->tmp    = 0;

  return(curpt);
}

void delPt(pMesh mesh,int ip) {
  pPoint   ppt;

  ppt = &mesh->point[ip];
  memset(ppt,0,sizeof(Point));
  ppt->tag    = MS_NUL;
  ppt->tmp    = mesh->npnil;
  mesh->npnil = ip;
  if ( ip == mesh->np ) {
	  while ( !MS_VOK((&mesh->point[mesh->np])) )  mesh->np--;
  }
}

int newElt(pMesh mesh) {
  int     curiel;

  if ( !mesh->ntnil )  return(0);
  curiel = mesh->ntnil;

  if ( mesh->ntnil > mesh->nt )  mesh->nt = mesh->ntnil;
  mesh->ntnil = mesh->tria[curiel].v[2];
  mesh->tria[curiel].v[2] = 0;

  return(curiel);
}

void delElt(pMesh mesh,int iel) {
  pTria    pt;

  pt = &mesh->tria[iel];
  if ( !MS_EOK(pt) ) {
    fprintf(stdout,"  ## INVALID ELEMENT: %d.\n",iel);
    return;
  }
  memset(pt,0,sizeof(Tria));
  pt->v[2] = mesh->ntnil;
  if ( mesh->adja )
    memset(&mesh->adja[3*(iel-1)+1],0,3*sizeof(int));
  mesh->ntnil = iel;
	if ( iel == mesh->nt ) {
		while ( !MS_EOK((&mesh->tria[mesh->nt])) )  mesh->nt--;       
	}
}


int zaldy(pMesh mesh) {
  int     million = 1048576L;
  int     k,npask,bytes;

  if ( info.mem < 0 ) {
    mesh->npmax = MS_MAX(1.5*mesh->np,NPMAX);
    mesh->ntmax = MS_MAX(1.5*mesh->nt,NTMAX);
  }
  else {
    /* point+tria+adja */
    bytes = sizeof(Point) + 2*sizeof(Tria) + 3*sizeof(int);

    npask = (int)((double)info.mem / bytes * million);
    mesh->npmax = MS_MAX(1.5*mesh->np,npask);
    mesh->ntmax = MS_MAX(1.5*mesh->nt,2*npask);
  }

  mesh->point = (pPoint)calloc(mesh->npmax+1,sizeof(Point));
  assert(mesh->point);
  mesh->tria  = (pTria)calloc(mesh->ntmax+1,sizeof(Tria));
  assert(mesh->tria);

  /* store empty links */
  mesh->npnil = mesh->np + 1;
  mesh->ntnil = mesh->nt + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  for (k=mesh->ntnil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;

  return(1);
}
