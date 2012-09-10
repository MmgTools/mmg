#include "mmg3d.h"

extern Info  info;


int scaleMesh(pMesh mesh,pSol met) {
  pPoint    ppt;
  double    dd;
  int       i,k;

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

	return(1);
}