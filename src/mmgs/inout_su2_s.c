/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file inout_su2_s.c
 *  \brief SU2 surface mesh input/output for MMGS.
 */

#include "libmmgs.h"
#include "libmmgs_private.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>

#define MMGS_SU2_LINE_LENGTH 4096

typedef struct {
  int      type;
  MMG5_int vertex[4];
  MMG5_int ref;
} MMGS_Su2Face;

static char *MMGS_su2Trim(char *line) {
  char *end;

  while ( isspace((unsigned char)*line) ) ++line;
  end = line+strlen(line);
  while ( end > line && isspace((unsigned char)end[-1]) ) --end;
  *end = '\0';
  return line;
}

static int MMGS_su2Line(FILE *inm,char line[MMGS_SU2_LINE_LENGTH]) {
  char *trimmed;

  while ( fgets(line,MMGS_SU2_LINE_LENGTH,inm) ) {
    if ( !strchr(line,'\n') && !feof(inm) ) return -1;
    trimmed = MMGS_su2Trim(line);
    if ( !*trimmed || *trimmed == '%' ) continue;
    if ( trimmed != line ) memmove(line,trimmed,strlen(trimmed)+1);
    return 1;
  }
  return ferror(inm) ? -1 : 0;
}

static int MMGS_su2Integer(const char *line,const char *key,MMG5_int *value) {
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

static int MMGS_su2Vertex(long long value,MMG5_int *vertex) {
  if ( value < 0 || value == LLONG_MAX ||
       (sizeof(MMG5_int) == 4 && value >= INT32_MAX) ) return 0;
  *vertex = (MMG5_int)value+1;
  return 1;
}

static int MMGS_su2CellSize(int type) {
  switch ( type ) {
  case 3: return 2;
  case 5: return 3;
  case 9: case 10: return 4;
  case 14: return 5;
  case 13: return 6;
  case 12: return 8;
  default: return 0;
  }
}

static int MMGS_su2MarkerRef(const char *name,MMG5_int *reference,
                             unsigned char *isExplicit) {
  const char *number = name;
  char       *end;
  long long  parsed;

  if ( !strncmp(name,"mmg_ref_",8) ) number = name+8;
  if ( number == name ) {
    *reference = 0;
    *isExplicit = 0;
    return 1;
  }
  errno = 0;
  parsed = strtoll(number,&end,10);
  while ( isspace((unsigned char)*end) ) ++end;
  if ( errno || end == number || *end ||
       (sizeof(MMG5_int) == 4 &&
        (parsed < INT32_MIN || parsed > INT32_MAX)) ) return 0;
  *reference = (MMG5_int)parsed;
  *isExplicit = 1;
  return 1;
}

static int MMGS_su2CompareRef(const void *left,const void *right) {
  const MMG5_int a = *(const MMG5_int *)left;
  const MMG5_int b = *(const MMG5_int *)right;

  return (a > b)-(a < b);
}

/** Assign fallback marker references after every explicit mmg_ref_N value is
 * known, so an ordinary marker cannot accidentally consume a reserved value. */
static int MMGS_su2ResolveMarkerRefs(MMG5_int *refs,
                                    const unsigned char *isExplicit,
                                    MMG5_int *reserved,MMG5_int nmark) {
  MMG5_int candidate=1,i,nreserved=0,reservedIndex=0;

  for ( i=0; i<nmark; ++i ) {
    if ( isExplicit[i] && refs[i] > 0 ) reserved[nreserved++] = refs[i];
  }
  qsort(reserved,(size_t)nreserved,sizeof(MMG5_int),MMGS_su2CompareRef);
  for ( i=0; i<nmark; ++i ) {
    if ( isExplicit[i] ) continue;
    while ( reservedIndex < nreserved && reserved[reservedIndex] < candidate ) {
      ++reservedIndex;
    }
    while ( reservedIndex < nreserved && reserved[reservedIndex] == candidate ) {
      while ( reservedIndex < nreserved && reserved[reservedIndex] == candidate ) {
        ++reservedIndex;
      }
      if ( candidate == MMG5_INTMAX ) return 0;
      ++candidate;
    }
    refs[i] = candidate;
    if ( candidate == MMG5_INTMAX ) {
      for ( ++i; i<nmark; ++i ) if ( !isExplicit[i] ) return 0;
      break;
    }
    ++candidate;
  }
  return 1;
}

static int MMGS_su2ReadFace(const char *line,MMGS_Su2Face *face) {
  long long v[4];
  int       count,i;

  if ( sscanf(line,"%d",&face->type) != 1 ) return 0;
  count = face->type == 5 ? 3 : (face->type == 9 ? 4 : 0);
  if ( count == 3 ) {
    if ( sscanf(line,"%*d %lld %lld %lld",&v[0],&v[1],&v[2]) != 3 ) return 0;
  }
  else if ( count == 4 ) {
    if ( sscanf(line,"%*d %lld %lld %lld %lld",&v[0],&v[1],&v[2],&v[3]) != 4 ) {
      return 0;
    }
  }
  else return 0;
  for ( i=0; i<count; ++i ) {
    if ( !MMGS_su2Vertex(v[i],&face->vertex[i]) ) return 0;
  }
  return 1;
}

static int MMGS_su2SetFace(MMG5_pMesh mesh,const MMGS_Su2Face *face,
                           MMG5_int *triangleIndex) {
  MMG5_int a,b,c,d;
  MMG5_int ac0,ac1,bd0,bd1;

  a = face->vertex[0];
  b = face->vertex[1];
  c = face->vertex[2];
  if ( face->type == 5 ) {
    return MMGS_Set_triangle(mesh,a,b,c,face->ref,++(*triangleIndex));
  }
  d = face->vertex[3];
  ac0 = MG_MIN(a,c); ac1 = MG_MAX(a,c);
  bd0 = MG_MIN(b,d); bd1 = MG_MAX(b,d);
  /* Pick a diagonal from global IDs so the split is deterministic even when
   * equivalent files use a different local orientation. */
  if ( ac0 < bd0 || (ac0 == bd0 && ac1 < bd1) ) {
    return MMGS_Set_triangle(mesh,a,b,c,face->ref,++(*triangleIndex)) &&
           MMGS_Set_triangle(mesh,a,c,d,face->ref,++(*triangleIndex));
  }
  return MMGS_Set_triangle(mesh,b,c,d,face->ref,++(*triangleIndex)) &&
         MMGS_Set_triangle(mesh,b,d,a,face->ref,++(*triangleIndex));
}

int MMGS_loadSu2Mesh(MMG5_pMesh mesh,const char *filename) {
  MMGS_Su2Face *domainFaces=NULL,*markerFaces=NULL,*selected;
  MMG5_int      *markerRefs=NULL,*reservedRefs=NULL;
  unsigned char *explicitRefs=NULL;
  double        *points=NULL;
  MMG5_int      ndim=0,nelem=0,np=0,nmark=0,ndomain=0,nmarker=0;
  MMG5_int      markerCapacity=0,i,j,nt=0;
  char          line[MMGS_SU2_LINE_LENGTH],tag[MMGS_SU2_LINE_LENGTH];
  FILE          *inm;
  int           status,haveCells=0,havePoints=0,hasVolume=0;

  if ( !filename || !(inm=fopen(filename,"r")) ) return 0;
  if ( MMGS_su2Line(inm,line) < 1 ||
       MMGS_su2Integer(line,"NDIME",&ndim) != 1 || ndim != 3 ) goto parse_error;

  /* Accept both the SU2 section order and the NPOIN-first order emitted by
   * some third-party writers. */
  while ( !haveCells || !havePoints ) {
    if ( MMGS_su2Line(inm,line) < 1 ) goto parse_error;
    if ( !haveCells && MMGS_su2Integer(line,"NELEM",&nelem) == 1 ) {
      if ( nelem ) MMG5_SAFE_MALLOC(domainFaces,nelem,MMGS_Su2Face,
                                    goto memory_error);
      for ( i=0; i<nelem; ++i ) {
        long long v[8];
        int       type,count,nread,k;

        if ( MMGS_su2Line(inm,line) < 1 || sscanf(line,"%d",&type) != 1 ||
             !(count=MMGS_su2CellSize(type)) ) goto parse_error;
        nread = sscanf(line,"%*d %lld %lld %lld %lld %lld %lld %lld %lld",
                       &v[0],&v[1],&v[2],&v[3],&v[4],&v[5],&v[6],&v[7]);
        if ( nread < count ) goto parse_error;
        if ( type == 5 || type == 9 ) {
          domainFaces[ndomain].type = type;
          domainFaces[ndomain].ref = 0;
          for ( k=0; k<count; ++k ) {
            if ( !MMGS_su2Vertex(v[k],&domainFaces[ndomain].vertex[k]) ) {
              goto parse_error;
            }
          }
          ++ndomain;
        }
        else if ( type == 10 || type == 12 || type == 13 || type == 14 ) {
          hasVolume = 1;
        }
      }
      haveCells = 1;
    }
    else if ( !havePoints && MMGS_su2Integer(line,"NPOIN",&np) == 1 && np &&
              np <= MMG5_INTMAX/3 ) {
      MMG5_SAFE_MALLOC(points,3*np,double,goto memory_error);
      for ( i=0; i<np; ++i ) {
        if ( MMGS_su2Line(inm,line) < 1 ||
             sscanf(line,"%lf %lf %lf",&points[3*i],&points[3*i+1],
                    &points[3*i+2]) != 3 || !isfinite(points[3*i]) ||
             !isfinite(points[3*i+1]) || !isfinite(points[3*i+2]) ) {
          goto parse_error;
        }
      }
      havePoints = 1;
    }
    else goto parse_error;
  }

  status = MMGS_su2Line(inm,line);
  if ( status < 0 ) goto parse_error;
  if ( status > 0 ) {
    if ( MMGS_su2Integer(line,"NMARK",&nmark) != 1 ) goto parse_error;
    if ( nmark ) {
      MMG5_SAFE_MALLOC(markerRefs,nmark,MMG5_int,goto memory_error);
      MMG5_SAFE_MALLOC(reservedRefs,nmark,MMG5_int,goto memory_error);
      MMG5_SAFE_MALLOC(explicitRefs,nmark,unsigned char,goto memory_error);
    }
    for ( i=0; i<nmark; ++i ) {
      MMG5_int markerElements;
      char     *equal;

      if ( MMGS_su2Line(inm,line) < 1 || strncmp(line,"MARKER_TAG",10) ||
           !(equal=strchr(line,'=')) ) goto parse_error;
      strncpy(tag,MMGS_su2Trim(equal+1),sizeof(tag)-1);
      tag[sizeof(tag)-1] = '\0';
      if ( !MMGS_su2MarkerRef(tag,&markerRefs[i],&explicitRefs[i]) ||
           MMGS_su2Line(inm,line) < 1 ||
           MMGS_su2Integer(line,"MARKER_ELEMS",&markerElements) != 1 ) {
        goto parse_error;
      }
      if ( markerElements > MMG5_INTMAX-nmarker ) goto parse_error;
      if ( nmarker+markerElements > markerCapacity ) {
        MMG5_int old=markerCapacity;
        markerCapacity=nmarker+markerElements;
        MMG5_SAFE_REALLOC(markerFaces,old,markerCapacity,MMGS_Su2Face,
                          "SU2 surface faces",goto memory_error);
      }
      for ( j=0; j<markerElements; ++j ) {
        if ( MMGS_su2Line(inm,line) < 1 ||
             !MMGS_su2ReadFace(line,&markerFaces[nmarker]) ) goto parse_error;
        /* Keep the ordinal until explicit references from all later markers
         * are known and collision-free fallbacks can be assigned. */
        markerFaces[nmarker++].ref=i;
      }
    }
    if ( !MMGS_su2ResolveMarkerRefs(markerRefs,explicitRefs,reservedRefs,nmark) ) {
      goto parse_error;
    }
    for ( i=0; i<nmarker; ++i ) markerFaces[i].ref=markerRefs[markerFaces[i].ref];
  }

  /* Standard 3D SU2 files store their boundary surface in NMARK. For
   * interoperability, also accept standalone surface files whose producer
   * placed triangles/quads in NELEM, but never mix those with volume cells. */
  if ( nmarker ) {
    selected = markerFaces;
    i = nmarker;
  }
  else if ( ndomain && !hasVolume ) {
    selected = domainFaces;
    i = ndomain;
  }
  else goto parse_error;
  for ( j=0; j<i; ++j ) {
    MMG5_int increment=selected[j].type == 5 ? 1 : 2;
    int      count=selected[j].type == 5 ? 3 : 4,k;
    if ( nt > MMG5_INTMAX-increment ) goto parse_error;
    nt += increment;
    for ( k=0; k<count; ++k ) {
      if ( selected[j].vertex[k] > np ) goto parse_error;
    }
  }
  if ( !MMGS_Set_meshSize(mesh,np,nt,0) ) goto memory_error;
  for ( j=0; j<np; ++j ) {
    if ( !MMGS_Set_vertex(mesh,points[3*j],points[3*j+1],points[3*j+2],0,j+1) ) {
      goto parse_error;
    }
  }
  nt = 0;
  for ( j=0; j<i; ++j ) {
    if ( !MMGS_su2SetFace(mesh,&selected[j],&nt) ) goto parse_error;
  }

  fclose(inm);
  MMG5_SAFE_FREE(domainFaces);
  MMG5_SAFE_FREE(markerFaces);
  MMG5_SAFE_FREE(markerRefs);
  MMG5_SAFE_FREE(reservedRefs);
  MMG5_SAFE_FREE(explicitRefs);
  MMG5_SAFE_FREE(points);
  MMG5_check_readedMesh(mesh,0);
  return 1;

parse_error:
  fprintf(stderr,"  ## Error: unable to parse SU2 surface mesh %s.\n",filename);
memory_error:
  fclose(inm);
  MMG5_SAFE_FREE(domainFaces);
  MMG5_SAFE_FREE(markerFaces);
  MMG5_SAFE_FREE(markerRefs);
  MMG5_SAFE_FREE(reservedRefs);
  MMG5_SAFE_FREE(explicitRefs);
  MMG5_SAFE_FREE(points);
  return -1;
}

static int MMGS_su2ValidTriangle(MMG5_pMesh mesh,MMG5_pTria triangle) {
  int i;

  if ( !MG_EOK(triangle) ) return 0;
  for ( i=0; i<3; ++i ) {
    MMG5_int vertex = triangle->v[i];
    if ( vertex <= 0 || vertex > mesh->np || mesh->point[vertex].tmp < 0 ) {
      return 0;
    }
  }
  return 1;
}

static int MMGS_su2AddRef(MMG5_int **refs,MMG5_int *size,
                          MMG5_int *capacity,MMG5_int ref) {
  MMG5_int i;

  for ( i=0; i<*size; ++i ) if ( (*refs)[i] == ref ) return 1;
  if ( *size == *capacity ) {
    MMG5_int old = *capacity;
    if ( old > MMG5_INTMAX/2 ) return 0;
    *capacity = old ? 2*old : 8;
    MMG5_SAFE_REALLOC(*refs,old,*capacity,MMG5_int,"SU2 marker refs",return 0);
  }
  (*refs)[(*size)++] = ref;
  return 1;
}

int MMGS_saveSu2Mesh(MMG5_pMesh mesh,const char *filename) {
  MMG5_int    *refs=NULL;
  MMG5_int    i,j,np=0,nmark=0,refCapacity=0;
  MMG5_pPoint point;
  MMG5_pTria  triangle;
  FILE        *out;

  if ( !filename || !(out=fopen(filename,"w")) ) return 0;
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    point->tmp = MG_VOK(point) ? np++ : -1;
  }
  for ( i=1; i<=mesh->nt; ++i ) {
    triangle = &mesh->tria[i];
    if ( MMGS_su2ValidTriangle(mesh,triangle) &&
         !MMGS_su2AddRef(&refs,&nmark,&refCapacity,triangle->ref) ) goto error;
  }

  /* A 3D SU2 boundary surface belongs in marker sections. NELEM is therefore
   * empty, matching the representation written by common SU2 tooling. */
  fprintf(out,"NDIME= 3\nNELEM= 0\nNPOIN= %" MMG5_PRId "\n",np);
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    if ( !MG_VOK(point) ) continue;
    fprintf(out,"%.17g %.17g %.17g %" MMG5_PRId "\n",point->c[0],
            point->c[1],point->c[2],point->tmp);
  }
  fprintf(out,"NMARK= %" MMG5_PRId "\n",nmark);
  for ( j=0; j<nmark; ++j ) {
    MMG5_int count=0;
    for ( i=1; i<=mesh->nt; ++i ) {
      triangle = &mesh->tria[i];
      if ( triangle->ref == refs[j] && MMGS_su2ValidTriangle(mesh,triangle) ) {
        ++count;
      }
    }
    fprintf(out,"MARKER_TAG= mmg_ref_%" MMG5_PRId "\n"
            "MARKER_ELEMS= %" MMG5_PRId "\n",refs[j],count);
    for ( i=1; i<=mesh->nt; ++i ) {
      triangle = &mesh->tria[i];
      if ( triangle->ref != refs[j] || !MMGS_su2ValidTriangle(mesh,triangle) ) {
        continue;
      }
      fprintf(out,"5 %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
              mesh->point[triangle->v[0]].tmp,
              mesh->point[triangle->v[1]].tmp,
              mesh->point[triangle->v[2]].tmp);
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
