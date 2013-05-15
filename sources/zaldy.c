#include "mmg3d.h"

extern Info  info;

/** get new point address */
int newPt(pMesh mesh,double c[3],char tag) {
  pPoint  ppt;
  int     curpt;

  if ( !mesh->npnil )  return(0);
  curpt = mesh->npnil;
  if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
  ppt   = &mesh->point[curpt];
  memcpy(ppt->c,c,3*sizeof(double));
  ppt->tag   &= ~MG_NUL;
  mesh->npnil = ppt->tmp;
  ppt->tmp    = 0;

  /* point on geometry */
  if ( tag & MG_BDY ) {
    mesh->xp++;
    if(mesh->xp >= mesh->xpmax){
      printf("%s:%d: Error: unable to allocate a new xpoint 0\n",
             __FILE__,__LINE__);
      exit(EXIT_FAILURE);
    }
    ppt->tag = MG_BDY;
    ppt->xp  = mesh->xp;
  }
  return(curpt);
}

void delPt(pMesh mesh,int ip) {
  pPoint   ppt;
  xPoint  *pxp;

  ppt = &mesh->point[ip];
  if ( ppt->xp ) {
    pxp = &mesh->xpoint[ppt->xp];
    memset(pxp,0,sizeof(xPoint));
  }
  memset(ppt,0,sizeof(Point));
  ppt->tag    = MG_NUL;
  ppt->tmp    = mesh->npnil;
  mesh->npnil = ip;
  if ( ip == mesh->np ) {
    while ( !MG_VOK((&mesh->point[mesh->np])) )  mesh->np--;
  }
}

/** get new elt address */
int newElt(pMesh mesh) {
  int     curiel;

  if ( !mesh->nenil )  return(0);
  curiel = mesh->nenil;

  if ( mesh->nenil > mesh->ne )  mesh->ne = mesh->nenil;
  mesh->nenil = mesh->tetra[curiel].v[3];
  mesh->tetra[curiel].v[3] = 0;

  return(curiel);
}


void delElt(pMesh mesh,int iel) {
  pTetra   pt;
  int      iadr;

  pt = &mesh->tetra[iel];
  if ( !MG_EOK(pt) ) {
    fprintf(stdout,"  ## INVALID ELEMENT %d.\n",iel);
    exit(EXIT_FAILURE);
  }
  memset(pt,0,sizeof(Tetra));
  pt->v[3] = mesh->nenil;
  iadr = 4*(iel-1) + 1;
  if ( mesh->adja )
    memset(&mesh->adja[iadr],0,4*sizeof(int));
  mesh->nenil = iel;
  if ( iel == mesh->ne ) {
    while ( !MG_EOK((&mesh->tetra[mesh->ne])) )  mesh->ne--;
  }
}


/** allocate main structure */
int zaldy(pMesh mesh) {
  int     million = 1048576L;
  int     k,npask,bytes;

  if ( info.mem < 0 ) {
    mesh->npmax = MG_MAX(1.5*mesh->np,NPMAX);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,NTMAX);
    mesh->nemax = MG_MAX(1.5*mesh->ne,NEMAX);
  }
  else {
    /* point+tria+tets+adja+sol */
    bytes = 2*sizeof(Point) + 6*sizeof(Tetra) + 4*sizeof(int);

    npask = (double)info.mem / bytes * million;
    mesh->npmax = MG_MAX(1.5*mesh->np,npask);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,2*npask);
    mesh->nemax = MG_MAX(1.5*mesh->ne,6*npask);
  }

  mesh->point = (pPoint)calloc(mesh->npmax+1,sizeof(Point));
  assert(mesh->point);
  mesh->tetra = (pTetra)calloc(mesh->nemax+1,sizeof(Tetra));
  assert(mesh->tetra);
  if ( mesh->nt ) {
    mesh->tria = (pTria)calloc(mesh->ntmax+1,sizeof(Tria));
    assert(mesh->tria);
  }
  if ( mesh->na ) {
    mesh->edge = (pEdge)calloc(mesh->na+1,sizeof(Edge));
    assert(mesh->edge);
  }
  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->ne + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  for (k=mesh->nenil; k<mesh->nemax-1; k++)
    mesh->tetra[k].v[3] = k+1;

  return(1);
}
