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
      fprintf(stdout,"  ## Allocation problem (xpoint), not enough memory.\n");
      fprintf(stdout,"  ## Check the mesh size or ");
      fprintf(stdout,"increase the allocated memory with the -m option.\n");
      return(0);
    }
    ppt->xp  = mesh->xp;
  }
  ppt->tag = tag;
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
  mesh->tetra[curiel].mark=0;

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
  int     k,npask,bytes,ctri;

  if ( info.mem < 0 ) {
    mesh->npmax = MG_MAX(1.5*mesh->np,NPMAX);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,NTMAX);
    mesh->nemax = MG_MAX(1.5*mesh->ne,NEMAX);
  }
  else {
#ifdef SINGUL
    if ( info.sing )
      ctri = 4;
    else
      ctri = 2;
#else
      ctri = 2;
#endif

    /* Euler-poincare: ne = 6*np; nt = 2*np; na = np/5 *
     * point+tria+tets+adja+adjt+sol+item *
     * warning: we exceed memory in saveMesh when we call hNew */
    bytes = sizeof(Point) + sizeof(xPoint) +
      6*sizeof(Tetra) + ctri*sizeof(xTetra) +
      4*6*sizeof(int) + ctri*3*sizeof(int) +
      sizeof(Sol)+4*sizeof(hedge);

    npask = (double)info.mem / bytes * million;
    mesh->npmax = npask;
    mesh->ntmax = ctri*npask;
    mesh->nemax = 6*npask;
  }

  mesh->point = (pPoint)calloc(mesh->npmax+1,sizeof(Point));
  if ( !mesh->point ){
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }

  mesh->tetra = (pTetra)calloc(mesh->nemax+1,sizeof(Tetra));
  if ( !mesh->tetra ){
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  if ( mesh->nt ) {
    mesh->tria = (pTria)calloc(mesh->ntmax+1,sizeof(Tria));
    if ( !mesh->tria ){
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
  }
  if ( mesh->na ) {
    mesh->edge = (pEdge)calloc(mesh->na+1,sizeof(Edge));
    if ( !mesh->edge ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
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

/** free xtetra */
void freeXTets(pMesh mesh) {
  pTetra pt;
  int    k;

  for (k=1; k<=mesh->ne; k++) {
    pt     = &mesh->tetra[k];
    pt->xt = 0;
  }
  free(mesh->xtetra);
  mesh->xtetra = 0;
  mesh->xt = 0;
}
