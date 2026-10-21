/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file inout_nas_3d.c
 *  \brief Low-order Nastran bulk-data input/output for MMG3D.
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"
#include "inout_nas.h"

static int MMG3D_nasVertexCount(enum MMG5_NastranElementType type) {
  switch ( type ) {
  case MMG5_NAS_CTRIA3: return 3;
  case MMG5_NAS_CQUAD4: return 4;
  case MMG5_NAS_CTETRA: return 4;
  case MMG5_NAS_CPENTA: return 6;
  case MMG5_NAS_CPYRAM: return 5;
  case MMG5_NAS_CHEXA: return 8;
  }
  return 0;
}

static MMG5_int MMG3D_nasConvertedTetraCount(
  enum MMG5_NastranElementType type) {
  switch ( type ) {
  case MMG5_NAS_CTETRA: return 1;
  case MMG5_NAS_CPENTA: return 8;
  case MMG5_NAS_CPYRAM: return 6;
  case MMG5_NAS_CHEXA: return 12;
  default: return 0;
  }
}

static int MMG3D_nasSetFace(MMG5_pMesh mesh,MMG5_int center,
                            const MMG5_int vertex[4],int size,MMG5_int ref,
                            MMG5_int *tetraIndex) {
  MMG5_int a=vertex[0],b=vertex[1],c=vertex[2],d;
  MMG5_int ac0,ac1,bd0,bd1;

  if ( size == 3 ) {
    return MMG3D_Set_tetrahedron(mesh,center,a,b,c,ref,++(*tetraIndex));
  }
  d = vertex[3];
  ac0 = MG_MIN(a,c); ac1 = MG_MAX(a,c);
  bd0 = MG_MIN(b,d); bd1 = MG_MAX(b,d);
  /* A diagonal chosen from global point IDs is independent of local face
   * orientation, so adjacent converted cells triangulate their shared quad
   * identically. */
  if ( ac0 < bd0 || (ac0 == bd0 && ac1 < bd1) ) {
    return MMG3D_Set_tetrahedron(mesh,center,a,b,c,ref,++(*tetraIndex)) &&
           MMG3D_Set_tetrahedron(mesh,center,a,c,d,ref,++(*tetraIndex));
  }
  return MMG3D_Set_tetrahedron(mesh,center,b,c,d,ref,++(*tetraIndex)) &&
         MMG3D_Set_tetrahedron(mesh,center,b,d,a,ref,++(*tetraIndex));
}

static int MMG3D_nasSetConvertedCell(MMG5_pMesh mesh,
                                     enum MMG5_NastranElementType type,
                                     const MMG5_int vertex[8],MMG5_int center,
                                     MMG5_int ref,MMG5_int *tetraIndex) {
  static const int prismFaces[5][4] = {
    {0,2,1,-1},{3,4,5,-1},{0,1,4,3},{1,2,5,4},{2,0,3,5}
  };
  static const int pyramidFaces[5][4] = {
    {0,3,2,1},{0,1,4,-1},{1,2,4,-1},{2,3,4,-1},{3,0,4,-1}
  };
  static const int hexFaces[6][4] = {
    {0,3,2,1},{4,5,6,7},{0,1,5,4},
    {1,2,6,5},{2,3,7,6},{3,0,4,7}
  };
  const int (*faces)[4];
  MMG5_int face[4];
  int      i,j,faceCount,faceSize;

  if ( type == MMG5_NAS_CPENTA ) {
    faces = prismFaces; faceCount = 5;
  }
  else if ( type == MMG5_NAS_CPYRAM ) {
    faces = pyramidFaces; faceCount = 5;
  }
  else if ( type == MMG5_NAS_CHEXA ) {
    faces = hexFaces; faceCount = 6;
  }
  else return 0;

  /* Coning a consistently triangulated cell boundary to its center yields a
   * conforming tetrahedralization for convex low-order Nastran solids. */
  for ( i=0; i<faceCount; ++i ) {
    faceSize = faces[i][3] < 0 ? 3 : 4;
    for ( j=0; j<faceSize; ++j ) face[j] = vertex[faces[i][j]];
    if ( !MMG3D_nasSetFace(mesh,center,face,faceSize,ref,tetraIndex) ) return 0;
  }
  return 1;
}

int MMG3D_loadNastranMesh(MMG5_pMesh mesh,const char *filename) {
  MMG5_NastranMesh     nas;
  MMG5_NastranElement *element;
  MMG5_int             i,j,npOut,neOut=0,ne=0,nprism=0,nt=0,nquad=0;
  MMG5_int             ncenters=0,vertex[8];
  int                  convertHybrid=0,nvertex,k,ier;

  ier = MMG5_loadNastranMeshData(filename,&nas);
  if ( ier < 1 ) return ier;
  for ( i=0; i<nas.ne; ++i ) {
    element = &nas.elements[i];
    nvertex = MMG3D_nasVertexCount(element->type);
    for ( j=0; j<nvertex; ++j ) {
      if ( !MMG5_nastranPointIndex(&nas,element->vertex[j]) ) goto error;
    }
    switch ( element->type ) {
    case MMG5_NAS_CTRIA3: ++nt; break;
    case MMG5_NAS_CQUAD4: ++nquad; break;
    case MMG5_NAS_CTETRA: ++ne; break;
    case MMG5_NAS_CPENTA: ++nprism; break;
    case MMG5_NAS_CPYRAM: case MMG5_NAS_CHEXA:
      convertHybrid = 1; ++ncenters; break;
    }
  }
  /* A split hex or pyramid quad cannot meet an unsplit prism quad. Convert
   * every prism whenever either unsupported cell type is present. */
  if ( convertHybrid ) {
    if ( nprism > MMG5_INTMAX-ncenters ) goto error;
    ncenters += nprism;
  }
  if ( ncenters > MMG5_INTMAX-nas.np ) goto error;
  npOut = nas.np+ncenters;
  if ( convertHybrid ) {
    for ( i=0; i<nas.ne; ++i ) {
      MMG5_int increment = MMG3D_nasConvertedTetraCount(nas.elements[i].type);

      if ( increment > MMG5_INTMAX-neOut ) goto error;
      neOut += increment;
    }
    nprism = 0;
  }
  else neOut = ne;
  if ( !MMG3D_Set_meshSize(mesh,npOut,neOut,nprism,nt,nquad,0) ) goto error;
  for ( i=0; i<nas.np; ++i ) {
    if ( !MMG3D_Set_vertex(mesh,nas.points[i].c[0],nas.points[i].c[1],
                          nas.points[i].c[2],0,i+1) ) goto error;
  }

  ne = nprism = nt = nquad = 0;
  npOut = nas.np;
  for ( i=0; i<nas.ne; ++i ) {
    double center[3] = {0.,0.,0.};

    element = &nas.elements[i];
    nvertex = MMG3D_nasVertexCount(element->type);
    for ( j=0; j<nvertex; ++j ) {
      vertex[j] = MMG5_nastranPointIndex(&nas,element->vertex[j]);
    }
    if ( element->type == MMG5_NAS_CTRIA3 ) {
      if ( !MMG3D_Set_triangle(mesh,vertex[0],vertex[1],vertex[2],element->ref,
                              ++nt) ) goto error;
    }
    else if ( element->type == MMG5_NAS_CQUAD4 ) {
      if ( !MMG3D_Set_quadrilateral(mesh,vertex[0],vertex[1],vertex[2],
                                   vertex[3],element->ref,++nquad) ) goto error;
    }
    else if ( element->type == MMG5_NAS_CTETRA ) {
      if ( !MMG3D_Set_tetrahedron(mesh,vertex[0],vertex[1],vertex[2],vertex[3],
                                 element->ref,++ne) ) goto error;
    }
    else if ( !convertHybrid ) {
      /* Nastran/VTK and Mmg orient the two prism triangles oppositely. */
      if ( !MMG3D_Set_prism(mesh,vertex[0],vertex[2],vertex[1],vertex[3],
                           vertex[5],vertex[4],element->ref,++nprism) ) {
        goto error;
      }
    }
    else {
      MMG5_int centerIndex = ++npOut;

      for ( j=0; j<nvertex; ++j ) {
        for ( k=0; k<3; ++k ) center[k] += mesh->point[vertex[j]].c[k];
      }
      for ( k=0; k<3; ++k ) center[k] /= nvertex;
      if ( !MMG3D_Set_vertex(mesh,center[0],center[1],center[2],0,
                            centerIndex) ||
           !MMG3D_nasSetConvertedCell(mesh,element->type,vertex,centerIndex,
                                     element->ref,&ne) ) goto error;
    }
  }
  if ( convertHybrid && mesh->info.imprim >= 0 ) {
    fprintf(stdout,"  ## Warning: non-tetrahedral Nastran cells were converted "
            "to conforming tetrahedra.\n");
  }
  MMG5_freeNastranMeshData(&nas);
  MMG5_check_readedMesh(mesh,0);
  return 1;

error:
  fprintf(stderr,"  ## Error: unable to import Nastran mesh %s.\n",filename);
  MMG5_freeNastranMeshData(&nas);
  return -1;
}

static int MMG3D_nasValidPoint(MMG5_pMesh mesh,MMG5_int index) {
  return index > 0 && index <= mesh->np && mesh->point[index].tmp >= 0;
}

static int MMG3D_nasValidVertices(MMG5_pMesh mesh,const MMG5_int *vertex,
                                  int count) {
  int i;

  for ( i=0; i<count; ++i ) {
    if ( !MMG3D_nasValidPoint(mesh,vertex[i]) ) return 0;
  }
  return 1;
}

static int MMG3D_nasAddRef(MMG5_int **refs,MMG5_int *nref,
                           MMG5_int *capacity,MMG5_int ref) {
  MMG5_int i;

  for ( i=0; i<*nref; ++i ) {
    if ( (*refs)[i] == ref ) return 1;
  }
  if ( *nref == *capacity ) {
    MMG5_int old = *capacity;

    if ( old > MMG5_INTMAX/2 ) return 0;
    *capacity = old ? 2*old : 8;
    MMG5_SAFE_REALLOC(*refs,old,*capacity,MMG5_int,"Nastran references",
                      return 0);
  }
  (*refs)[(*nref)++] = ref;
  return 1;
}

static int MMG3D_nasAssignPids(const MMG5_int *refs,MMG5_int *pids,
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

static MMG5_int MMG3D_nasPid(const MMG5_int *refs,const MMG5_int *pids,
                             MMG5_int nref,MMG5_int ref) {
  MMG5_int i;

  for ( i=0; i<nref; ++i ) {
    if ( refs[i] == ref ) return pids[i];
  }
  return 0;
}

int MMG3D_saveNastranMesh(MMG5_pMesh mesh,const char *filename) {
  MMG5_int    *refs=NULL,*pids=NULL;
  MMG5_int    i,np=0,nref=0,capacity=0,eid=1,pid;
  MMG5_pPoint point;
  MMG5_pTetra tetra;
  MMG5_pPrism prism;
  MMG5_pTria  triangle;
  MMG5_pQuad  quad;
  FILE        *out;
  int         failed;

  if ( !filename || !(out=fopen(filename,"w")) ) return 0;
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    point->tmp = MG_VOK(point) ? ++np : -1;
  }
  for ( i=1; i<=mesh->ne; ++i ) {
    tetra = &mesh->tetra[i];
    if ( MG_EOK(tetra) && MMG3D_nasValidVertices(mesh,tetra->v,4) &&
         !MMG3D_nasAddRef(&refs,&nref,&capacity,tetra->ref) ) goto error;
  }
  for ( i=1; i<=mesh->nprism; ++i ) {
    prism = &mesh->prism[i];
    if ( MG_EOK(prism) && MMG3D_nasValidVertices(mesh,prism->v,6) &&
         !MMG3D_nasAddRef(&refs,&nref,&capacity,prism->ref) ) goto error;
  }
  for ( i=1; i<=mesh->nt; ++i ) {
    triangle = &mesh->tria[i];
    if ( MG_EOK(triangle) && MMG3D_nasValidVertices(mesh,triangle->v,3) &&
         !MMG3D_nasAddRef(&refs,&nref,&capacity,triangle->ref) ) goto error;
  }
  for ( i=1; i<=mesh->nquad; ++i ) {
    quad = &mesh->quadra[i];
    if ( MG_EOK(quad) && MMG3D_nasValidVertices(mesh,quad->v,4) &&
         !MMG3D_nasAddRef(&refs,&nref,&capacity,quad->ref) ) goto error;
  }
  if ( nref ) MMG5_SAFE_MALLOC(pids,nref,MMG5_int,goto error);
  if ( !MMG3D_nasAssignPids(refs,pids,nref) ) goto error;

  fprintf(out,"$ Nastran bulk-data mesh written by MMG3D\n");
  for ( i=0; i<nref; ++i ) {
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
  for ( i=1; i<=mesh->ne; ++i ) {
    tetra = &mesh->tetra[i];
    if ( !MG_EOK(tetra) || !MMG3D_nasValidVertices(mesh,tetra->v,4) ) continue;
    pid = MMG3D_nasPid(refs,pids,nref,tetra->ref);
    fprintf(out,"CTETRA,%" MMG5_PRId ",%" MMG5_PRId ",%" MMG5_PRId
            ",%" MMG5_PRId ",%" MMG5_PRId ",%" MMG5_PRId "\n",eid++,pid,
            mesh->point[tetra->v[0]].tmp,mesh->point[tetra->v[1]].tmp,
            mesh->point[tetra->v[2]].tmp,mesh->point[tetra->v[3]].tmp);
  }
  for ( i=1; i<=mesh->nprism; ++i ) {
    prism = &mesh->prism[i];
    if ( !MG_EOK(prism) || !MMG3D_nasValidVertices(mesh,prism->v,6) ) continue;
    pid = MMG3D_nasPid(refs,pids,nref,prism->ref);
    /* Convert Mmg's prism convention back to Nastran/VTK ordering. */
    fprintf(out,"CPENTA,%" MMG5_PRId ",%" MMG5_PRId ",%" MMG5_PRId
            ",%" MMG5_PRId ",%" MMG5_PRId ",%" MMG5_PRId
            ",%" MMG5_PRId ",%" MMG5_PRId "\n",eid++,pid,
            mesh->point[prism->v[0]].tmp,mesh->point[prism->v[2]].tmp,
            mesh->point[prism->v[1]].tmp,mesh->point[prism->v[3]].tmp,
            mesh->point[prism->v[5]].tmp,mesh->point[prism->v[4]].tmp);
  }
  for ( i=1; i<=mesh->nt; ++i ) {
    triangle = &mesh->tria[i];
    if ( !MG_EOK(triangle) ||
         !MMG3D_nasValidVertices(mesh,triangle->v,3) ) continue;
    pid = MMG3D_nasPid(refs,pids,nref,triangle->ref);
    fprintf(out,"CTRIA3,%" MMG5_PRId ",%" MMG5_PRId ",%" MMG5_PRId
            ",%" MMG5_PRId ",%" MMG5_PRId "\n",eid++,pid,
            mesh->point[triangle->v[0]].tmp,
            mesh->point[triangle->v[1]].tmp,
            mesh->point[triangle->v[2]].tmp);
  }
  for ( i=1; i<=mesh->nquad; ++i ) {
    quad = &mesh->quadra[i];
    if ( !MG_EOK(quad) || !MMG3D_nasValidVertices(mesh,quad->v,4) ) continue;
    pid = MMG3D_nasPid(refs,pids,nref,quad->ref);
    fprintf(out,"CQUAD4,%" MMG5_PRId ",%" MMG5_PRId ",%" MMG5_PRId
            ",%" MMG5_PRId ",%" MMG5_PRId ",%" MMG5_PRId "\n",eid++,pid,
            mesh->point[quad->v[0]].tmp,mesh->point[quad->v[1]].tmp,
            mesh->point[quad->v[2]].tmp,mesh->point[quad->v[3]].tmp);
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
