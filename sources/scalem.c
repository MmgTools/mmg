#include "mmg3d.h"

extern Info  info;


int scaleMesh(pMesh mesh,pSol met,pSingul sing) {
  pPoint    ppt;
#ifdef SINGUL
  psPoint   ppts;
  double    deltb,delta[3];
#endif
  double    dd;
  int       i,k;
  pPar      par;

  /* compute bounding box */
  for (i=0; i<3; i++) {
    info.min[i] =  DBL_MAX;
    info.max[i] = -DBL_MAX;
  }
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    for (i=0; i<3; i++) {
      if ( ppt->c[i] > info.max[i] )  info.max[i] = ppt->c[i];
      if ( ppt->c[i] < info.min[i] )  info.min[i] = ppt->c[i];
    }
  }
  info.delta = 0.0;
  for (i=0; i<3; i++) {
    dd = info.max[i] - info.min[i];
    if ( dd > info.delta )  info.delta = dd;
  }
  if ( info.delta < EPSD ) {
    fprintf(stdout,"  ## Unable to scale mesh.\n");
    return(0);
  }

  /* normalize coordinates */
  dd = 1.0 / info.delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->c[0] = dd * (ppt->c[0] - info.min[0]);
    ppt->c[1] = dd * (ppt->c[1] - info.min[1]);
    ppt->c[2] = dd * (ppt->c[2] - info.min[2]);
  }

  /* normalize values */
  info.hmin  *= dd;
  info.hmax  *= dd;
  info.hausd *= dd;

  /* normalize sizes */
  if ( met->size == 1 && met->m ) {
    for (k=1; k<=mesh->np; k++)
      met->m[k] *= dd;
  }

#ifdef SINGUL
  /* 2nd mesh (sing) is quarter sized */
  /* Get the size of sing in every direction */
  if ( info.sing && sing->ns ) {
    deltb = 0.0;

    for (i=0; i<mesh->dim; i++) {
      delta[i] = fabs(sing->max[i]-sing->min[i]);
      if ( delta[i] > deltb )  deltb = delta[i];   // deltb = max dimension
    }
    if ( deltb < EPSD ) {
      fprintf(stdout,"  ## Unable to scale mesh\n");
      return(0);
    }

    /* centering */
    dd = 1.0 / deltb;
    for (k=1; k<=sing->ns; k++) {
      ppts = &sing->point[k];
      for (i=0; i<mesh->dim; i++) {
        ppts->c[i] = SIZE * (dd * (ppts->c[i]-sing->min[i])) +
          0.5 * (1.0 - SIZE*dd*delta[i]);
      }
    }
  }
#endif

  /* normalize local parameters */
	for (k=0; k<info.npar; k++) {
		par = &info.par[k];
		par->hmin  *= dd;
		par->hmax  *= dd;
		par->hausd *= dd;
	}

  return(1);
}

int unscaleMesh(pMesh mesh,pSol met) {
  pPoint     ppt;
  double     dd;
  int        k;

  /* de-normalize coordinates */
  dd = info.delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->c[0] = ppt->c[0] * dd + info.min[0];
    ppt->c[1] = ppt->c[1] * dd + info.min[1];
    ppt->c[2] = ppt->c[2] * dd + info.min[2];
  }

  /* unscale sizes */
  if(met->m){
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) )	met->m[k] *= dd;
    }
  }
  return(1);
}
