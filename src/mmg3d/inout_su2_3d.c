/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file inout_su2_3d.c
 *  \brief SU2 mesh input/output for MMG3D.
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>

#define MMG3D_SU2_LINE_LENGTH 4096

typedef struct {
  int      type;
  MMG5_int vertex[8];
} MMG3D_Su2Cell;

typedef struct {
  int      type;
  MMG5_int vertex[4];
  MMG5_int ref;
} MMG3D_Su2Boundary;

static char *MMG3D_su2Trim(char *line) {
  char *end;

  while ( isspace((unsigned char)*line) ) ++line;
  end = line+strlen(line);
  while ( end > line && isspace((unsigned char)end[-1]) ) --end;
  *end = '\0';
  return line;
}

static int MMG3D_su2Line(FILE *inm,char line[MMG3D_SU2_LINE_LENGTH]) {
  char *trimmed;

  while ( fgets(line,MMG3D_SU2_LINE_LENGTH,inm) ) {
    if ( !strchr(line,'\n') && !feof(inm) ) return -1;
    trimmed = MMG3D_su2Trim(line);
    if ( !*trimmed || *trimmed == '%' ) continue;
    if ( trimmed != line ) memmove(line,trimmed,strlen(trimmed)+1);
    return 1;
  }
  return ferror(inm) ? -1 : 0;
}

static int MMG3D_su2Integer(const char *line,const char *key,MMG5_int *value) {
  const char *equal;
  char       *end;
  long long  parsed;
  size_t     keyLength = strlen(key);

  if ( strncmp(line,key,keyLength) ||
       (line[keyLength] != '=' && !isspace((unsigned char)line[keyLength])) ) {
    return 0;
  }
  equal = strchr(line+keyLength,'=');
  if ( !equal ) return 0;
  errno = 0;
  parsed = strtoll(equal+1,&end,10);
  while ( isspace((unsigned char)*end) ) ++end;
  if ( errno || end == equal+1 || (*end && *end != '%') || parsed < 0 ||
       (sizeof(MMG5_int) == 4 && parsed > INT32_MAX) ) return -1;
  *value = (MMG5_int)parsed;
  return 1;
}

static int MMG3D_su2Vertex(long long value,MMG5_int *vertex) {
  if ( value < 0 || value == LLONG_MAX ||
       (sizeof(MMG5_int) == 4 && value >= INT32_MAX) ) {
    return 0;
  }
  *vertex = (MMG5_int)value+1;
  return 1;
}

static int MMG3D_su2MarkerRef(const char *name,MMG5_int fallback,
                              MMG5_int *reference) {
  const char *number = name;
  char       *end;
  long long  parsed;

  if ( !strncmp(name,"mmg_ref_",8) ) number = name+8;
  if ( number == name ) {
    *reference = fallback;
    return 1;
  }
  errno = 0;
  parsed = strtoll(number,&end,10);
  while ( isspace((unsigned char)*end) ) ++end;
  if ( errno || end == number || *end ||
       (sizeof(MMG5_int) == 4 &&
        (parsed < INT32_MIN || parsed > INT32_MAX)) ) return 0;
  *reference = (MMG5_int)parsed;
  return 1;
}

static int MMG3D_su2CellSize(int type) {
  switch ( type ) {
  case 10: /* tetrahedron */
    return 4;
  case 12: /* hexahedron */
    return 8;
  case 13: /* prism */
    return 6;
  case 14: /* pyramid */
    return 5;
  default:
    return 0;
  }
}

static MMG5_int MMG3D_su2ConvertedTetraCount(int type) {
  switch ( type ) {
  case 10:
    return 1;
  case 12:
    return 12;
  case 13:
    return 8;
  case 14:
    return 6;
  default:
    return 0;
  }
}

static int MMG3D_su2SetFace(MMG5_pMesh mesh,MMG5_int center,
                            const MMG5_int vertex[4],int size,
                            MMG5_int *tetraIndex) {
  MMG5_int a,b,c,d;
  MMG5_int ac0,ac1,bd0,bd1;

  a = vertex[0];
  b = vertex[1];
  c = vertex[2];
  if ( size == 3 ) {
    return MMG3D_Set_tetrahedron(mesh,center,a,b,c,0,++(*tetraIndex));
  }

  d = vertex[3];
  ac0 = MG_MIN(a,c);
  ac1 = MG_MAX(a,c);
  bd0 = MG_MIN(b,d);
  bd1 = MG_MAX(b,d);
  /* Select the diagonal from global vertex IDs, not from the local face
   * orientation. Thus two cells sharing this quadrilateral always generate
   * the same pair of triangles, even if their face order is reversed. */
  if ( ac0 < bd0 || (ac0 == bd0 && ac1 < bd1) ) {
    return MMG3D_Set_tetrahedron(mesh,center,a,b,c,0,++(*tetraIndex)) &&
           MMG3D_Set_tetrahedron(mesh,center,a,c,d,0,++(*tetraIndex));
  }
  return MMG3D_Set_tetrahedron(mesh,center,b,c,d,0,++(*tetraIndex)) &&
         MMG3D_Set_tetrahedron(mesh,center,b,d,a,0,++(*tetraIndex));
}

static int MMG3D_su2SetConvertedCell(MMG5_pMesh mesh,
                                     const MMG3D_Su2Cell *cell,
                                     MMG5_int center,MMG5_int *tetraIndex) {
  static const int prismFaces[5][4] = {
    {0,2,1,-1},{3,4,5,-1},{0,1,4,3},{1,2,5,4},{2,0,3,5}
  };
  static const int hexFaces[6][4] = {
    {0,3,2,1},{4,5,6,7},{0,1,5,4},
    {1,2,6,5},{2,3,7,6},{3,0,4,7}
  };
  static const int pyramidFaces[5][4] = {
    {0,3,2,1},{0,1,4,-1},{1,2,4,-1},{2,3,4,-1},{3,0,4,-1}
  };
  const int (*faces)[4];
  MMG5_int face[4];
  int      i,j,faceCount,faceSize;

  switch ( cell->type ) {
  case 12:
    faces = hexFaces;
    faceCount = 6;
    break;
  case 13:
    faces = prismFaces;
    faceCount = 5;
    break;
  case 14:
    faces = pyramidFaces;
    faceCount = 5;
    break;
  default:
    return 0;
  }

  /* Cone a consistently triangulated boundary of the cell to its center.
   * This produces a conforming tetrahedralization for convex SU2 cells. */
  for ( i=0; i<faceCount; ++i ) {
    faceSize = faces[i][3] < 0 ? 3 : 4;
    for ( j=0; j<faceSize; ++j ) face[j] = cell->vertex[faces[i][j]];
    if ( !MMG3D_su2SetFace(mesh,center,face,faceSize,tetraIndex) ) return 0;
  }
  return 1;
}

int MMG3D_loadSu2Mesh(MMG5_pMesh mesh,const char *filename) {
  MMG3D_Su2Cell     *cells = NULL;
  MMG3D_Su2Boundary *boundary = NULL;
  double             *points = NULL;
  MMG5_int           ndim=0,nelem=0,np=0,nmark=0,ne=0,nprism=0;
  MMG5_int           nt=0,nquad=0,i,j,ncenters=0,boundarySize=0;
  MMG5_int           boundaryCapacity=0,npOut,neOut;
  char               line[MMG3D_SU2_LINE_LENGTH];
  char               tag[MMG3D_SU2_LINE_LENGTH];
  FILE               *inm;
  int                status,convertHybrid=0,haveCells=0,havePoints=0;

  if ( !filename || !(inm=fopen(filename,"r")) ) return 0;
  if ( MMG3D_su2Line(inm,line) < 1 ||
       MMG3D_su2Integer(line,"NDIME",&ndim) != 1 || ndim != 3 ) {
    goto parse_error;
  }
  /* SU2 tools conventionally put NELEM first, while some third-party writers
   * put NPOIN first. Both sections are self-describing, so accept either. */
  while ( !haveCells || !havePoints ) {
    if ( MMG3D_su2Line(inm,line) < 1 ) goto parse_error;
    if ( !haveCells && MMG3D_su2Integer(line,"NELEM",&nelem) == 1 ) {
      MMG5_SAFE_MALLOC(cells,nelem,MMG3D_Su2Cell,goto memory_error);
      for ( i=0; i<nelem; ++i ) {
        long long v[8];
        int       count;

        if ( MMG3D_su2Line(inm,line) < 1 ||
             sscanf(line,"%d",&cells[i].type) != 1 ) goto parse_error;
        count = MMG3D_su2CellSize(cells[i].type);
        switch ( count ) {
        case 4:
          if ( sscanf(line,"%*d %lld %lld %lld %lld",&v[0],&v[1],&v[2],
                      &v[3]) != 4 ) goto parse_error;
          ++ne;
          break;
        case 5:
          if ( sscanf(line,"%*d %lld %lld %lld %lld %lld",&v[0],&v[1],
                      &v[2],&v[3],&v[4]) != 5 ) goto parse_error;
          convertHybrid = 1;
          ++ncenters;
          break;
        case 6:
          if ( sscanf(line,"%*d %lld %lld %lld %lld %lld %lld",&v[0],&v[1],
                      &v[2],&v[3],&v[4],&v[5]) != 6 ) goto parse_error;
          ++nprism;
          break;
        case 8:
          if ( sscanf(line,"%*d %lld %lld %lld %lld %lld %lld %lld %lld",
                      &v[0],&v[1],&v[2],&v[3],&v[4],&v[5],&v[6],
                      &v[7]) != 8 ) goto parse_error;
          convertHybrid = 1;
          ++ncenters;
          break;
        default:
          fprintf(stderr,"  ## Error: unsupported SU2 3D element code %d.\n",
                  cells[i].type);
          goto parse_error;
        }
        for ( j=0; j<count; ++j ) {
          if ( !MMG3D_su2Vertex(v[j],&cells[i].vertex[j]) ) goto parse_error;
        }
      }
      haveCells = 1;
    }
    else if ( !havePoints && MMG3D_su2Integer(line,"NPOIN",&np) == 1 && np &&
              np <= MMG5_INTMAX/3 ) {
      MMG5_SAFE_MALLOC(points,3*np,double,goto memory_error);
      for ( i=0; i<np; ++i ) {
        if ( MMG3D_su2Line(inm,line) < 1 ||
             sscanf(line,"%lf %lf %lf",&points[3*i],&points[3*i+1],
                    &points[3*i+2]) != 3 || !isfinite(points[3*i]) ||
             !isfinite(points[3*i+1]) || !isfinite(points[3*i+2]) ) {
          goto parse_error;
        }
      }
      havePoints = 1;
    }
    else {
      goto parse_error;
    }
  }

  /* A converted hex face cannot meet an unsplit prism quad conformingly.
   * Therefore, as soon as a hex or pyramid is present, convert every prism
   * too and obtain an entirely tetrahedral volume mesh. Pure tet/prism files
   * retain their native cells. */
  if ( convertHybrid ) {
    if ( nprism > MMG5_INTMAX-ncenters ) goto parse_error;
    ncenters += nprism;
  }
  if ( ncenters > MMG5_INTMAX-np || np+ncenters > MMG5_INTMAX/3 ) {
    goto parse_error;
  }
  npOut = np+ncenters;
  if ( ncenters ) {
    MMG5_SAFE_REALLOC(points,3*np,3*npOut,double,"SU2 points",
                      goto memory_error);
  }

  status = MMG3D_su2Line(inm,line);
  if ( status < 0 ) goto parse_error;
  if ( status > 0 ) {
    if ( MMG3D_su2Integer(line,"NMARK",&nmark) != 1 ) goto parse_error;
    for ( i=0; i<nmark; ++i ) {
      MMG5_int markerElements,reference;
      char     *equal;

      if ( MMG3D_su2Line(inm,line) < 1 || strncmp(line,"MARKER_TAG",10) ||
           !(equal=strchr(line,'=')) ) goto parse_error;
      strncpy(tag,MMG3D_su2Trim(equal+1),sizeof(tag)-1);
      tag[sizeof(tag)-1] = '\0';
      if ( !MMG3D_su2MarkerRef(tag,i+1,&reference) ||
           MMG3D_su2Line(inm,line) < 1 ||
           MMG3D_su2Integer(line,"MARKER_ELEMS",&markerElements) != 1 ) {
        goto parse_error;
      }
      if ( markerElements > MMG5_INTMAX-boundarySize ) goto parse_error;
      if ( boundarySize+markerElements > boundaryCapacity ) {
        MMG5_int old = boundaryCapacity;

        boundaryCapacity = boundarySize+markerElements;
        MMG5_SAFE_REALLOC(boundary,old,boundaryCapacity,MMG3D_Su2Boundary,
                          "SU2 boundary",goto memory_error);
      }
      for ( j=0; j<markerElements; ++j ) {
        MMG3D_Su2Boundary *entity = &boundary[boundarySize];
        long long          v[4];
        int                k,count;

        if ( MMG3D_su2Line(inm,line) < 1 ||
             sscanf(line,"%d",&entity->type) != 1 ) goto parse_error;
        if ( entity->type == 5 ) {
          count = 3;
          ++nt;
          if ( sscanf(line,"%*d %lld %lld %lld",&v[0],&v[1],&v[2]) != 3 ) {
            goto parse_error;
          }
        }
        else if ( entity->type == 9 ) {
          count = 4;
          ++nquad;
          if ( sscanf(line,"%*d %lld %lld %lld %lld",&v[0],&v[1],&v[2],
                      &v[3]) != 4 ) goto parse_error;
        }
        else {
          fprintf(stderr,"  ## Error: unsupported SU2 3D boundary element "
                  "code %d.\n",entity->type);
          goto parse_error;
        }
        entity->ref = reference;
        for ( k=0; k<count; ++k ) {
          if ( !MMG3D_su2Vertex(v[k],&entity->vertex[k]) ||
               entity->vertex[k] > np ) goto parse_error;
        }
        ++boundarySize;
      }
    }
  }

  for ( i=0; i<nelem; ++i ) {
    int count = MMG3D_su2CellSize(cells[i].type);

    for ( j=0; j<count; ++j ) {
      if ( cells[i].vertex[j] > np ) goto parse_error;
    }
  }

  if ( convertHybrid ) {
    neOut = 0;
    for ( i=0; i<nelem; ++i ) {
      MMG5_int increment = MMG3D_su2ConvertedTetraCount(cells[i].type);

      if ( increment > MMG5_INTMAX-neOut ) goto parse_error;
      neOut += increment;
    }
    nprism = 0;
  }
  else {
    neOut = ne;
  }
  if ( !MMG3D_Set_meshSize(mesh,npOut,neOut,nprism,nt,nquad,0) ) {
    goto memory_error;
  }
  for ( i=0; i<np; ++i ) {
    if ( !MMG3D_Set_vertex(mesh,points[3*i],points[3*i+1],points[3*i+2],
                          0,i+1) ) goto parse_error;
  }

  ne = nprism = 0;
  npOut = np;
  for ( i=0; i<nelem; ++i ) {
    if ( cells[i].type == 10 ) {
      if ( !MMG3D_Set_tetrahedron(mesh,cells[i].vertex[0],cells[i].vertex[1],
                                 cells[i].vertex[2],cells[i].vertex[3],0,
                                 ++ne) ) goto parse_error;
    }
    else if ( !convertHybrid ) {
      /* SU2 and Mmg orient the two triangular prism faces oppositely. */
      if ( !MMG3D_Set_prism(mesh,cells[i].vertex[0],cells[i].vertex[2],
                           cells[i].vertex[1],cells[i].vertex[3],
                           cells[i].vertex[5],cells[i].vertex[4],0,
                           ++nprism) ) goto parse_error;
    }
    else {
      double   center[3] = {0.,0.,0.};
      MMG5_int centerIndex = ++npOut;
      int      count = MMG3D_su2CellSize(cells[i].type);
      int      k;

      for ( j=0; j<count; ++j ) {
        for ( k=0; k<3; ++k ) {
          center[k] += points[3*(cells[i].vertex[j]-1)+k];
        }
      }
      for ( k=0; k<3; ++k ) center[k] /= count;
      if ( !MMG3D_Set_vertex(mesh,center[0],center[1],center[2],0,
                            centerIndex) ||
           !MMG3D_su2SetConvertedCell(mesh,&cells[i],centerIndex,&ne) ) {
        goto parse_error;
      }
    }
  }

  nt = nquad = 0;
  for ( i=0; i<boundarySize; ++i ) {
    if ( boundary[i].type == 5 &&
         !MMG3D_Set_triangle(mesh,boundary[i].vertex[0],boundary[i].vertex[1],
                            boundary[i].vertex[2],boundary[i].ref,
                            ++nt) ) goto parse_error;
    if ( boundary[i].type == 9 &&
         !MMG3D_Set_quadrilateral(mesh,boundary[i].vertex[0],
                                 boundary[i].vertex[1],boundary[i].vertex[2],
                                 boundary[i].vertex[3],boundary[i].ref,
                                 ++nquad) ) goto parse_error;
  }

  if ( convertHybrid && mesh->info.imprim >= 0 ) {
    fprintf(stdout,"  ## Warning: non-tetrahedral SU2 cells were converted "
            "to conforming tetrahedra.\n");
  }
  fclose(inm);
  MMG5_SAFE_FREE(cells);
  MMG5_SAFE_FREE(boundary);
  MMG5_SAFE_FREE(points);
  MMG5_check_readedMesh(mesh,0);
  return 1;

parse_error:
  fprintf(stderr,"  ## Error: unable to parse SU2 mesh %s.\n",filename);
memory_error:
  fclose(inm);
  MMG5_SAFE_FREE(cells);
  MMG5_SAFE_FREE(boundary);
  MMG5_SAFE_FREE(points);
  return -1;
}

static int MMG3D_su2ValidPoint(MMG5_pMesh mesh,MMG5_int index) {
  return index > 0 && index <= mesh->np && mesh->point[index].tmp >= 0;
}

static int MMG3D_su2ValidTetra(MMG5_pMesh mesh,MMG5_pTetra tetra) {
  int i;

  if ( !MG_EOK(tetra) ) return 0;
  for ( i=0; i<4; ++i ) {
    if ( !MMG3D_su2ValidPoint(mesh,tetra->v[i]) ) return 0;
  }
  return 1;
}

static int MMG3D_su2ValidPrism(MMG5_pMesh mesh,MMG5_pPrism prism) {
  int i;

  if ( !MG_EOK(prism) ) return 0;
  for ( i=0; i<6; ++i ) {
    if ( !MMG3D_su2ValidPoint(mesh,prism->v[i]) ) return 0;
  }
  return 1;
}

static int MMG3D_su2ValidTriangle(MMG5_pMesh mesh,MMG5_pTria triangle) {
  int i;

  if ( !MG_EOK(triangle) ) return 0;
  for ( i=0; i<3; ++i ) {
    if ( !MMG3D_su2ValidPoint(mesh,triangle->v[i]) ) return 0;
  }
  return 1;
}

static int MMG3D_su2ValidQuad(MMG5_pMesh mesh,MMG5_pQuad quad) {
  int i;

  if ( !MG_EOK(quad) ) return 0;
  for ( i=0; i<4; ++i ) {
    if ( !MMG3D_su2ValidPoint(mesh,quad->v[i]) ) return 0;
  }
  return 1;
}

static int MMG3D_su2AddRef(MMG5_int **refs,MMG5_int *size,
                           MMG5_int *capacity,MMG5_int ref) {
  MMG5_int i;

  for ( i=0; i<*size; ++i ) {
    if ( (*refs)[i] == ref ) return 1;
  }
  if ( *size == *capacity ) {
    MMG5_int old = *capacity;

    if ( old > MMG5_INTMAX/2 ) return 0;
    *capacity = old ? 2*old : 8;
    MMG5_SAFE_REALLOC(*refs,old,*capacity,MMG5_int,"SU2 marker refs",
                      return 0);
  }
  (*refs)[(*size)++] = ref;
  return 1;
}

int MMG3D_saveSu2Mesh(MMG5_pMesh mesh,const char *filename) {
  MMG5_int    *refs = NULL;
  MMG5_int    i,j,np=0,nelem=0,nmark=0,refCapacity=0,elementId=0;
  MMG5_pPoint point;
  MMG5_pTetra tetra;
  MMG5_pPrism prism;
  MMG5_pTria  triangle;
  MMG5_pQuad  quad;
  FILE        *out;

  if ( !filename || !(out=fopen(filename,"w")) ) return 0;
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    point->tmp = MG_VOK(point) ? np++ : -1;
  }
  for ( i=1; i<=mesh->ne; ++i ) {
    if ( MMG3D_su2ValidTetra(mesh,&mesh->tetra[i]) ) ++nelem;
  }
  for ( i=1; i<=mesh->nprism; ++i ) {
    if ( MMG3D_su2ValidPrism(mesh,&mesh->prism[i]) ) ++nelem;
  }
  for ( i=1; i<=mesh->nt; ++i ) {
    triangle = &mesh->tria[i];
    if ( MMG3D_su2ValidTriangle(mesh,triangle) &&
         !MMG3D_su2AddRef(&refs,&nmark,&refCapacity,triangle->ref) ) goto error;
  }
  for ( i=1; i<=mesh->nquad; ++i ) {
    quad = &mesh->quadra[i];
    if ( MMG3D_su2ValidQuad(mesh,quad) &&
         !MMG3D_su2AddRef(&refs,&nmark,&refCapacity,quad->ref) ) goto error;
  }

  fprintf(out,"NDIME= 3\nNELEM= %" MMG5_PRId "\n",nelem);
  for ( i=1; i<=mesh->ne; ++i ) {
    tetra = &mesh->tetra[i];
    if ( !MMG3D_su2ValidTetra(mesh,tetra) ) continue;
    fprintf(out,"10 %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId
            " %" MMG5_PRId " %" MMG5_PRId "\n",
            mesh->point[tetra->v[0]].tmp,mesh->point[tetra->v[1]].tmp,
            mesh->point[tetra->v[2]].tmp,mesh->point[tetra->v[3]].tmp,
            elementId++);
  }
  for ( i=1; i<=mesh->nprism; ++i ) {
    prism = &mesh->prism[i];
    if ( !MMG3D_su2ValidPrism(mesh,prism) ) continue;
    /* SU2 and Mmg orient the two triangular prism faces oppositely. */
    fprintf(out,"13 %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId
            " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId
            " %" MMG5_PRId "\n",mesh->point[prism->v[0]].tmp,
            mesh->point[prism->v[2]].tmp,mesh->point[prism->v[1]].tmp,
            mesh->point[prism->v[3]].tmp,mesh->point[prism->v[5]].tmp,
            mesh->point[prism->v[4]].tmp,elementId++);
  }
  fprintf(out,"NPOIN= %" MMG5_PRId "\n",np);
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    if ( !MG_VOK(point) ) continue;
    fprintf(out,"%.17g %.17g %.17g %" MMG5_PRId "\n",point->c[0],
            point->c[1],point->c[2],point->tmp);
  }
  fprintf(out,"NMARK= %" MMG5_PRId "\n",nmark);
  for ( j=0; j<nmark; ++j ) {
    MMG5_int count = 0;

    for ( i=1; i<=mesh->nt; ++i ) {
      triangle = &mesh->tria[i];
      if ( triangle->ref == refs[j] &&
           MMG3D_su2ValidTriangle(mesh,triangle) ) ++count;
    }
    for ( i=1; i<=mesh->nquad; ++i ) {
      quad = &mesh->quadra[i];
      if ( quad->ref == refs[j] && MMG3D_su2ValidQuad(mesh,quad) ) ++count;
    }
    fprintf(out,"MARKER_TAG= mmg_ref_%" MMG5_PRId "\n"
            "MARKER_ELEMS= %" MMG5_PRId "\n",refs[j],count);
    for ( i=1; i<=mesh->nt; ++i ) {
      triangle = &mesh->tria[i];
      if ( triangle->ref != refs[j] ||
           !MMG3D_su2ValidTriangle(mesh,triangle) ) continue;
      fprintf(out,"5 %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
              mesh->point[triangle->v[0]].tmp,
              mesh->point[triangle->v[1]].tmp,
              mesh->point[triangle->v[2]].tmp);
    }
    for ( i=1; i<=mesh->nquad; ++i ) {
      quad = &mesh->quadra[i];
      if ( quad->ref != refs[j] || !MMG3D_su2ValidQuad(mesh,quad) ) continue;
      fprintf(out,"9 %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId
              " %" MMG5_PRId "\n",mesh->point[quad->v[0]].tmp,
              mesh->point[quad->v[1]].tmp,mesh->point[quad->v[2]].tmp,
              mesh->point[quad->v[3]].tmp);
    }
  }
  i = ferror(out);
  MMG5_SAFE_FREE(refs);
  return !fclose(out) && !i;

error:
  MMG5_SAFE_FREE(refs);
  fclose(out);
  return 0;
}
