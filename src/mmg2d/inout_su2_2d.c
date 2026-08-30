/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file inout_su2_2d.c
 *  \brief SU2 mesh input/output for MMG2D.
 */

#include "libmmg2d.h"
#include "libmmg2d_private.h"

#include <ctype.h>
#include <errno.h>

#define MMG2D_SU2_LINE_LENGTH 4096

typedef struct {
  int      type;
  MMG5_int vertex[4];
} MMG2D_Su2Cell;

static char *MMG2D_su2Trim(char *line) {
  char *end;

  while ( isspace((unsigned char)*line) ) ++line;
  end = line+strlen(line);
  while ( end > line && isspace((unsigned char)end[-1]) ) --end;
  *end = '\0';
  return line;
}

static int MMG2D_su2Line(FILE *inm,char line[MMG2D_SU2_LINE_LENGTH]) {
  char *trimmed;

  while ( fgets(line,MMG2D_SU2_LINE_LENGTH,inm) ) {
    if ( !strchr(line,'\n') && !feof(inm) ) return -1;
    trimmed = MMG2D_su2Trim(line);
    if ( !*trimmed || *trimmed == '%' ) continue;
    if ( trimmed != line ) memmove(line,trimmed,strlen(trimmed)+1);
    return 1;
  }
  return ferror(inm) ? -1 : 0;
}

static int MMG2D_su2Integer(const char *line,const char *key,MMG5_int *value) {
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

static int MMG2D_su2Vertex(long long value,MMG5_int *vertex) {
  if ( value < 0 || value == LLONG_MAX ||
       (sizeof(MMG5_int) == 4 && value >= INT32_MAX) ) {
    return 0;
  }

  *vertex = (MMG5_int)value+1;
  return 1;
}

static int MMG2D_su2MarkerRef(const char *name,MMG5_int fallback,
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
       (sizeof(MMG5_int) == 4 && (parsed < INT32_MIN || parsed > INT32_MAX)) ) {
    return 0;
  }
  *reference = (MMG5_int)parsed;
  return 1;
}

int MMG2D_loadSu2Mesh(MMG5_pMesh mesh,const char *filename) {
  MMG2D_Su2Cell *cells = NULL;
  MMG5_int      *edges = NULL,*edgeRefs = NULL;
  double        *points = NULL;
  MMG5_int      ndim=0,nelem=0,np=0,nmark=0,nt=0,nquad=0,na=0;
  MMG5_int      i,j,edgeCapacity=0;
  char          line[MMG2D_SU2_LINE_LENGTH],tag[MMG2D_SU2_LINE_LENGTH];
  FILE          *inm;
  int           status,haveCells=0,havePoints=0;

  if ( !filename || !(inm=fopen(filename,"r")) ) return 0;
  if ( MMG2D_su2Line(inm,line) < 1 ||
       MMG2D_su2Integer(line,"NDIME",&ndim) != 1 || ndim != 2 ) goto parse_error;
  /* SU2 tools conventionally put NELEM first, while some third-party writers
   * put NPOIN first. Both sections are self-describing, so accept either. */
  while ( !haveCells || !havePoints ) {
    if ( MMG2D_su2Line(inm,line) < 1 ) goto parse_error;
    if ( !haveCells && MMG2D_su2Integer(line,"NELEM",&nelem) == 1 ) {
      MMG5_SAFE_MALLOC(cells,nelem,MMG2D_Su2Cell,goto memory_error);
      for ( i=0; i<nelem; ++i ) {
        long long v[4];

        if ( MMG2D_su2Line(inm,line) < 1 ||
             sscanf(line,"%d",&cells[i].type) != 1 ) goto parse_error;
        if ( cells[i].type == 5 ) {
          if ( sscanf(line,"%*d %lld %lld %lld",&v[0],&v[1],&v[2]) != 3 ) {
            goto parse_error;
          }
          ++nt;
          for ( j=0; j<3; ++j ) {
            if ( !MMG2D_su2Vertex(v[j],&cells[i].vertex[j]) ) goto parse_error;
          }
        }
        else if ( cells[i].type == 9 ) {
          if ( sscanf(line,"%*d %lld %lld %lld %lld",&v[0],&v[1],
                      &v[2],&v[3]) != 4 ) goto parse_error;
          ++nquad;
          for ( j=0; j<4; ++j ) {
            if ( !MMG2D_su2Vertex(v[j],&cells[i].vertex[j]) ) {
              goto parse_error;
            }
          }
        }
        else {
          fprintf(stderr,"  ## Error: unsupported SU2 2D element code %d.\n",
                  cells[i].type);
          goto parse_error;
        }
      }
      haveCells = 1;
    }
    else if ( !havePoints && MMG2D_su2Integer(line,"NPOIN",&np) == 1 && np &&
              np <= MMG5_INTMAX/2 ) {
      MMG5_SAFE_MALLOC(points,2*np,double,goto memory_error);
      for ( i=0; i<np; ++i ) {
        if ( MMG2D_su2Line(inm,line) < 1 ||
             sscanf(line,"%lf %lf",&points[2*i],&points[2*i+1]) != 2 ) {
          goto parse_error;
        }
      }
      havePoints = 1;
    }
    else {
      goto parse_error;
    }
  }
  status = MMG2D_su2Line(inm,line);
  if ( status < 0 ) goto parse_error;
  if ( status > 0 ) {
    if ( MMG2D_su2Integer(line,"NMARK",&nmark) != 1 ) goto parse_error;
    for ( i=0; i<nmark; ++i ) {
      MMG5_int markerElements,reference;
      char *equal;
      if ( MMG2D_su2Line(inm,line) < 1 || strncmp(line,"MARKER_TAG",10) ||
           !(equal=strchr(line,'=')) ) goto parse_error;
      strncpy(tag,MMG2D_su2Trim(equal+1),sizeof(tag)-1);
      tag[sizeof(tag)-1] = '\0';
      if ( !MMG2D_su2MarkerRef(tag,i+1,&reference) ||
           MMG2D_su2Line(inm,line) < 1 ||
           MMG2D_su2Integer(line,"MARKER_ELEMS",&markerElements) != 1 ) goto parse_error;
      if ( markerElements > MMG5_INTMAX-na ||
           na+markerElements > MMG5_INTMAX/2 ) goto parse_error;
      if ( na+markerElements > edgeCapacity ) {
        MMG5_int old = edgeCapacity;

        edgeCapacity = na+markerElements;
        MMG5_SAFE_REALLOC(edges,2*old,2*edgeCapacity,MMG5_int,"SU2 edges",
                          goto memory_error);
        MMG5_SAFE_REALLOC(edgeRefs,old,edgeCapacity,MMG5_int,"SU2 edge refs",
                          goto memory_error);
      }
      for ( j=0; j<markerElements; ++j ) {
        int       type;
        long long a,b;

        if ( MMG2D_su2Line(inm,line) < 1 ||
             sscanf(line,"%d %lld %lld",&type,&a,&b) != 3 || type != 3 ) goto parse_error;
        if ( !MMG2D_su2Vertex(a,&edges[2*na]) ||
             !MMG2D_su2Vertex(b,&edges[2*na+1]) ||
             edges[2*na] > np || edges[2*na+1] > np ) goto parse_error;
        edgeRefs[na++] = reference;
      }
    }
  }
  for ( i=0; i<nelem; ++i ) {
    for ( j=0; j<(cells[i].type == 5 ? 3 : 4); ++j ) {
      if ( cells[i].vertex[j] > np ) goto parse_error;
    }
  }
  if ( !MMG2D_Set_meshSize(mesh,np,nt,nquad,na) ) goto memory_error;
  for ( i=0; i<np; ++i ) {
    if ( !MMG2D_Set_vertex(mesh,points[2*i],points[2*i+1],0,i+1) ) {
      goto parse_error;
    }
  }
  nt = nquad = 0;
  for ( i=0; i<nelem; ++i ) {
    if ( cells[i].type == 5 &&
         !MMG2D_Set_triangle(mesh,cells[i].vertex[0],cells[i].vertex[1],
                            cells[i].vertex[2],0,++nt) ) goto parse_error;
    if ( cells[i].type == 9 &&
         !MMG2D_Set_quadrilateral(mesh,cells[i].vertex[0],cells[i].vertex[1],
                                 cells[i].vertex[2],cells[i].vertex[3],0,
                                 ++nquad) ) goto parse_error;
  }
  for ( i=0; i<na; ++i ) {
    if ( !MMG2D_Set_edge(mesh,edges[2*i],edges[2*i+1],edgeRefs[i],i+1) ) {
      goto parse_error;
    }
  }
  fclose(inm);
  MMG5_SAFE_FREE(cells);
  MMG5_SAFE_FREE(points);
  MMG5_SAFE_FREE(edges);
  MMG5_SAFE_FREE(edgeRefs);
  MMG5_check_readedMesh(mesh,0);
  return 1;
parse_error:
  fprintf(stderr,"  ## Error: unable to parse SU2 mesh %s.\n",filename);
memory_error:
  fclose(inm);
  MMG5_SAFE_FREE(cells);
  MMG5_SAFE_FREE(points);
  MMG5_SAFE_FREE(edges);
  MMG5_SAFE_FREE(edgeRefs);
  return -1;
}

static int MMG2D_su2ValidPoint(MMG5_pMesh mesh,MMG5_int index) {
  return index > 0 && index <= mesh->np && mesh->point[index].tmp >= 0;
}

int MMG2D_saveSu2Mesh(MMG5_pMesh mesh,const char *filename) {
  MMG5_pPoint point;
  MMG5_pTria  triangle;
  MMG5_pQuad  quadrilateral;
  MMG5_pEdge  edge;
  MMG5_int    *references = NULL;
  MMG5_int    i,j,np=0,nelem=0,nmark=0,elementId=0;
  FILE        *out;

  if ( !filename || !(out=fopen(filename,"w")) ) return 0;
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    point->tmp = MG_VOK(point) ? np++ : -1;
  }
  for ( i=1; i<=mesh->nt; ++i ) {
    triangle=&mesh->tria[i];
    if ( MG_EOK(triangle) && MMG2D_su2ValidPoint(mesh,triangle->v[0]) &&
         MMG2D_su2ValidPoint(mesh,triangle->v[1]) &&
         MMG2D_su2ValidPoint(mesh,triangle->v[2]) ) ++nelem;
  }
  for ( i=1; i<=mesh->nquad; ++i ) {
    quadrilateral=&mesh->quadra[i];
    if ( MG_EOK(quadrilateral) &&
         MMG2D_su2ValidPoint(mesh,quadrilateral->v[0]) &&
         MMG2D_su2ValidPoint(mesh,quadrilateral->v[1]) &&
         MMG2D_su2ValidPoint(mesh,quadrilateral->v[2]) &&
         MMG2D_su2ValidPoint(mesh,quadrilateral->v[3]) ) ++nelem;
  }
  if ( mesh->na ) MMG5_SAFE_MALLOC(references,mesh->na,MMG5_int,goto error);
  for ( i=1; i<=mesh->na; ++i ) {
    edge=&mesh->edge[i];
    if ( !MMG2D_su2ValidPoint(mesh,edge->a) ||
         !MMG2D_su2ValidPoint(mesh,edge->b) ) continue;
    for ( j=0; j<nmark && references[j] != edge->ref; ++j ) {}
    if ( j == nmark ) references[nmark++] = edge->ref;
  }
  fprintf(out,"NDIME= 2\nNELEM= %" MMG5_PRId "\n",nelem);
  for ( i=1; i<=mesh->nt; ++i ) {
    triangle=&mesh->tria[i];
    if ( !MG_EOK(triangle) || !MMG2D_su2ValidPoint(mesh,triangle->v[0]) ||
         !MMG2D_su2ValidPoint(mesh,triangle->v[1]) ||
         !MMG2D_su2ValidPoint(mesh,triangle->v[2]) ) continue;
    fprintf(out,"5 %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
            mesh->point[triangle->v[0]].tmp,mesh->point[triangle->v[1]].tmp,
            mesh->point[triangle->v[2]].tmp,elementId++);
  }
  for ( i=1; i<=mesh->nquad; ++i ) {
    quadrilateral=&mesh->quadra[i];
    if ( !MG_EOK(quadrilateral) ||
         !MMG2D_su2ValidPoint(mesh,quadrilateral->v[0]) ||
         !MMG2D_su2ValidPoint(mesh,quadrilateral->v[1]) ||
         !MMG2D_su2ValidPoint(mesh,quadrilateral->v[2]) ||
         !MMG2D_su2ValidPoint(mesh,quadrilateral->v[3]) ) continue;
    fprintf(out,"9 %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
            mesh->point[quadrilateral->v[0]].tmp,mesh->point[quadrilateral->v[1]].tmp,
            mesh->point[quadrilateral->v[2]].tmp,mesh->point[quadrilateral->v[3]].tmp,elementId++);
  }
  fprintf(out,"NPOIN= %" MMG5_PRId "\n",np);
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    if ( !MG_VOK(point) ) continue;
    fprintf(out,"%.17g %.17g %" MMG5_PRId "\n",point->c[0],point->c[1],point->tmp);
  }
  fprintf(out,"NMARK= %" MMG5_PRId "\n",nmark);
  for ( j=0; j<nmark; ++j ) {
    MMG5_int count=0;
    for ( i=1; i<=mesh->na; ++i ) {
      edge=&mesh->edge[i];
      if ( edge->ref == references[j] &&
           MMG2D_su2ValidPoint(mesh,edge->a) &&
           MMG2D_su2ValidPoint(mesh,edge->b) ) ++count;
    }
    fprintf(out,"MARKER_TAG= mmg_ref_%" MMG5_PRId "\n"
            "MARKER_ELEMS= %" MMG5_PRId "\n",references[j],count);
    for ( i=1; i<=mesh->na; ++i ) {
      edge = &mesh->edge[i];
      if ( edge->ref != references[j] ||
           !MMG2D_su2ValidPoint(mesh,edge->a) ||
           !MMG2D_su2ValidPoint(mesh,edge->b) ) continue;
      fprintf(out,"3 %" MMG5_PRId " %" MMG5_PRId "\n",
              mesh->point[edge->a].tmp,mesh->point[edge->b].tmp);
    }
  }
  i = ferror(out);
  MMG5_SAFE_FREE(references);
  return !fclose(out) && !i;
error:
  MMG5_SAFE_FREE(references);
  fclose(out);
  return 0;
}
