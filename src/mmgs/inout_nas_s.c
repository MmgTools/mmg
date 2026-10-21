/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file inout_nas_s.c
 *  \brief Low-order Nastran bulk-data input/output for MMGS.
 */

#include "libmmgs.h"
#include "libmmgs_private.h"
#include "inout_nas.h"

static int MMGS_nasValidPoint(MMG5_pMesh mesh,MMG5_int index) {
  return index > 0 && index <= mesh->np && mesh->point[index].tmp >= 0;
}

int MMGS_loadNastranMesh(MMG5_pMesh mesh,const char *filename) {
  MMG5_NastranMesh    nas;
  MMG5_NastranElement *element;
  MMG5_int            i,j,nt=0,vertex[8];
  int                 nvertex,ier;

  ier = MMG5_loadNastranMeshData(filename,&nas);
  if ( ier < 1 ) return ier;
  for ( i=0; i<nas.ne; ++i ) {
    element = &nas.elements[i];
    if ( element->type == MMG5_NAS_CTRIA3 ) {
      if ( nt == MMG5_INTMAX ) goto error;
      ++nt;
    }
    else if ( element->type == MMG5_NAS_CQUAD4 ) {
      if ( nt > MMG5_INTMAX-2 ) goto error;
      nt += 2;
    }
    else {
      fprintf(stderr,"  ## Error: solid Nastran card in an MMGS mesh.\n");
      goto error;
    }
    nvertex = element->type == MMG5_NAS_CTRIA3 ? 3 : 4;
    for ( j=0; j<nvertex; ++j ) {
      if ( !MMG5_nastranPointIndex(&nas,element->vertex[j]) ) goto error;
    }
  }
  if ( !MMGS_Set_meshSize(mesh,nas.np,nt,0) ) goto error;
  for ( i=0; i<nas.np; ++i ) {
    if ( !MMGS_Set_vertex(mesh,nas.points[i].c[0],nas.points[i].c[1],
                         nas.points[i].c[2],0,i+1) ) goto error;
  }
  nt = 0;
  for ( i=0; i<nas.ne; ++i ) {
    element = &nas.elements[i];
    nvertex = element->type == MMG5_NAS_CTRIA3 ? 3 : 4;
    for ( j=0; j<nvertex; ++j ) {
      vertex[j] = MMG5_nastranPointIndex(&nas,element->vertex[j]);
    }
    if ( element->type == MMG5_NAS_CTRIA3 &&
         !MMGS_Set_triangle(mesh,vertex[0],vertex[1],vertex[2],element->ref,
                           ++nt) ) goto error;
    if ( element->type == MMG5_NAS_CQUAD4 &&
         (!MMGS_Set_triangle(mesh,vertex[0],vertex[1],vertex[2],element->ref,
                            ++nt) ||
          !MMGS_Set_triangle(mesh,vertex[0],vertex[2],vertex[3],element->ref,
                            ++nt)) ) goto error;
  }
  MMG5_freeNastranMeshData(&nas);
  MMG5_check_readedMesh(mesh,0);
  return 1;

error:
  fprintf(stderr,"  ## Error: unable to import Nastran mesh %s.\n",filename);
  MMG5_freeNastranMeshData(&nas);
  return -1;
}

static int MMGS_nasAddRef(MMG5_int *refs,MMG5_int *nref,MMG5_int ref) {
  MMG5_int i;

  for ( i=0; i<*nref; ++i ) {
    if ( refs[i] == ref ) return 1;
  }
  refs[(*nref)++] = ref;
  return 1;
}

static int MMGS_nasAssignPids(const MMG5_int *refs,MMG5_int *pids,
                              MMG5_int nref) {
  MMG5_int i,j,next=1;

  for ( i=0; i<nref; ++i ) pids[i] = refs[i] > 0 ? refs[i] : 0;
  for ( i=0; i<nref; ++i ) {
    if ( pids[i] ) continue;
    while ( 1 ) {
      for ( j=0; j<nref && pids[j] != next; ++j ) {}
      if ( j == nref ) break;
      if ( next == MMG5_INTMAX ) return 0;
      ++next;
    }
    pids[i] = next;
  }
  return 1;
}

static MMG5_int MMGS_nasPid(const MMG5_int *refs,const MMG5_int *pids,
                            MMG5_int nref,MMG5_int ref) {
  MMG5_int i;

  for ( i=0; i<nref; ++i ) {
    if ( refs[i] == ref ) return pids[i];
  }
  return 0;
}

int MMGS_saveNastranMesh(MMG5_pMesh mesh,const char *filename) {
  MMG5_int   *refs=NULL,*pids=NULL;
  MMG5_int   i,np=0,nref=0,eid=1,pid;
  MMG5_pPoint point;
  MMG5_pTria triangle;
  FILE       *out;
  int        failed;

  if ( !filename || !(out=fopen(filename,"w")) ) return 0;
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    point->tmp = MG_VOK(point) ? ++np : -1;
  }
  if ( mesh->nt ) {
    MMG5_SAFE_MALLOC(refs,mesh->nt,MMG5_int,goto error);
    MMG5_SAFE_MALLOC(pids,mesh->nt,MMG5_int,goto error);
  }
  for ( i=1; i<=mesh->nt; ++i ) {
    triangle = &mesh->tria[i];
    if ( MG_EOK(triangle) && MMGS_nasValidPoint(mesh,triangle->v[0]) &&
         MMGS_nasValidPoint(mesh,triangle->v[1]) &&
         MMGS_nasValidPoint(mesh,triangle->v[2]) ) {
      MMGS_nasAddRef(refs,&nref,triangle->ref);
    }
  }
  if ( !MMGS_nasAssignPids(refs,pids,nref) ) goto error;

  fprintf(out,"$ Nastran bulk-data mesh written by MMGS\n");
  for ( i=0; i<nref; ++i ) {
    /* PIDs must be positive. This comment preserves arbitrary signed Mmg
     * references while remaining harmless to other Nastran readers. */
    if ( pids[i] != refs[i] ) {
      fprintf(out,"$MMG_REF,%" MMG5_PRId ",%" MMG5_PRId "\n",
              pids[i],refs[i]);
    }
  }
  fprintf(out,"BEGIN BULK\n");
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    if ( !MG_VOK(point) ) continue;
    fprintf(out,"GRID,%" MMG5_PRId ",,%.17g,%.17g,%.17g\n",point->tmp,
            point->c[0],point->c[1],point->c[2]);
  }
  for ( i=1; i<=mesh->nt; ++i ) {
    triangle = &mesh->tria[i];
    if ( !MG_EOK(triangle) || !MMGS_nasValidPoint(mesh,triangle->v[0]) ||
         !MMGS_nasValidPoint(mesh,triangle->v[1]) ||
         !MMGS_nasValidPoint(mesh,triangle->v[2]) ) continue;
    pid = MMGS_nasPid(refs,pids,nref,triangle->ref);
    fprintf(out,"CTRIA3,%" MMG5_PRId ",%" MMG5_PRId ",%" MMG5_PRId
            ",%" MMG5_PRId ",%" MMG5_PRId "\n",eid++,pid,
            mesh->point[triangle->v[0]].tmp,
            mesh->point[triangle->v[1]].tmp,
            mesh->point[triangle->v[2]].tmp);
  }
  fprintf(out,"ENDDATA\n");
  failed = ferror(out);
  MMG5_SAFE_FREE(refs);
  MMG5_SAFE_FREE(pids);
  return !fclose(out) && !failed;

error:
  MMG5_SAFE_FREE(refs);
  MMG5_SAFE_FREE(pids);
  fclose(out);
  return 0;
}
