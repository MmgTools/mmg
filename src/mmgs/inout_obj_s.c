/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
**
**  Mmg is free software: you can redistribute it and/or modify it under the
**  terms of the GNU Lesser General Public License as published by the Free
**  Software Foundation, either version 3 of the License, or any later version.
** =============================================================================
*/

/**
 * \file mmgs/inout_obj_s.c
 * \brief Wavefront OBJ input/output for MMGS.
 */

#include "libmmgs.h"
#include "libmmgs_private.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>

typedef struct {
  char     **names;
  MMG5_int *refs;
  size_t   size;
  size_t   capacity;
  MMG5_int nextRef;
} MMGS_ObjGroups;

/** Read an arbitrarily long text line. */
static int MMGS_objReadLine(FILE *inm,char **line,size_t *capacity) {
  size_t length;
  int    c;

  if ( !*line ) {
    *capacity = 256;
    MMG5_SAFE_MALLOC(*line,*capacity,char,return -1);
  }

  length = 0;
  while ( (c = fgetc(inm)) != EOF ) {
    if ( length+1 >= *capacity ) {
      size_t oldCapacity = *capacity;

      if ( oldCapacity > SIZE_MAX/2 ) return -1;
      *capacity *= 2;
      MMG5_SAFE_REALLOC(*line,oldCapacity,*capacity,char,"OBJ line",
                        return -1);
    }
    (*line)[length++] = (char)c;
    if ( c == '\n' ) break;
  }

  if ( c == EOF && !length ) return 0;
  (*line)[length] = '\0';
  return 1;
}

/** Return the payload of a line beginning with keyword, or NULL. */
static char *MMGS_objPayload(char *line,const char *keyword) {
  size_t length = strlen(keyword);

  while ( isspace((unsigned char)*line) ) ++line;
  if ( strncmp(line,keyword,length) ) return NULL;
  if ( line[length] && !isspace((unsigned char)line[length]) ) return NULL;

  line += length;
  while ( isspace((unsigned char)*line) ) ++line;
  return line;
}

/** Extract the next whitespace-delimited token, stopping at an OBJ comment. */
static char *MMGS_objNextToken(char **cursor) {
  char *begin = *cursor;

  while ( isspace((unsigned char)*begin) ) ++begin;
  if ( !*begin || *begin == '#' ) {
    *cursor = begin;
    return NULL;
  }

  *cursor = begin;
  while ( **cursor && !isspace((unsigned char)**cursor) && **cursor != '#' ) {
    ++*cursor;
  }
  if ( **cursor ) *(*cursor)++ = '\0';
  return begin;
}

/** Convert the vertex part of an OBJ face token to a one-based index. */
static int MMGS_objVertexIndex(char *token,MMG5_int npoint,
                               MMG5_int *index) {
  long long value;
  char      *end;

  errno = 0;
  value = strtoll(token,&end,10);
  if ( errno == ERANGE || end == token || (*end && *end != '/') || !value ) {
    return 0;
  }

  if ( value > 0 ) {
    if ( value > (long long)npoint ) return 0;
    *index = (MMG5_int)value;
  }
  else {
    if ( value < -(long long)npoint ) return 0;
    *index = npoint+(MMG5_int)value+1;
  }
  return 1;
}

static void MMGS_objFreeGroups(MMGS_ObjGroups *groups) {
  size_t i;

  for ( i=0; i<groups->size; ++i ) MMG5_SAFE_FREE(groups->names[i]);
  MMG5_SAFE_FREE(groups->names);
  MMG5_SAFE_FREE(groups->refs);
}

static int MMGS_objRefIsUsed(const MMGS_ObjGroups *groups,MMG5_int ref) {
  size_t i;

  for ( i=0; i<groups->size; ++i ) {
    if ( groups->refs[i] == ref ) return 1;
  }
  return 0;
}

/** Decode an explicit mmg_ref_<integer> group name. */
static int MMGS_objExplicitRef(const char *name,MMG5_int *ref) {
  char      *end;
  long long value;

  if ( strncmp(name,"mmg_ref_",8) ) return 0;
  errno = 0;
  value = strtoll(name+8,&end,10);
  if ( errno || end == name+8 || *end ||
       (sizeof(MMG5_int) == 4 &&
        (value < INT32_MIN || value > INT32_MAX)) ) return 0;
  *ref = (MMG5_int)value;
  return 1;
}

/** Get a stable reference for a group name. */
static int MMGS_objGroupRef(MMGS_ObjGroups *groups,const char *name,
                            MMG5_int *ref) {
  size_t     i,length;

  for ( i=0; i<groups->size; ++i ) {
    if ( !strcmp(groups->names[i],name) ) {
      *ref = groups->refs[i];
      return 1;
    }
  }

  if ( !MMGS_objExplicitRef(name,ref) ) {
    while ( MMGS_objRefIsUsed(groups,groups->nextRef) ) {
      if ( groups->nextRef == MMG5_INTMAX ) return 0;
      ++groups->nextRef;
    }
    *ref = groups->nextRef;
    if ( groups->nextRef < MMG5_INTMAX ) ++groups->nextRef;
  }

  if ( groups->size == groups->capacity ) {
    size_t oldCapacity = groups->capacity;
    size_t newCapacity;

    if ( oldCapacity > SIZE_MAX/2 ) return 0;
    newCapacity = oldCapacity ? 2*oldCapacity : 8;

    MMG5_SAFE_REALLOC(groups->names,oldCapacity,newCapacity,char *,
                      "OBJ groups",return 0);
    MMG5_SAFE_REALLOC(groups->refs,oldCapacity,newCapacity,MMG5_int,
                      "OBJ group references",return 0);
    groups->capacity = newCapacity;
  }

  length = strlen(name)+1;
  MMG5_SAFE_MALLOC(groups->names[groups->size],length,char,return 0);
  memcpy(groups->names[groups->size],name,length);
  groups->refs[groups->size] = *ref;
  ++groups->size;
  return 1;
}

/** Trim a group payload in place. */
static char *MMGS_objGroupName(char *payload) {
  char *end = strchr(payload,'#');

  if ( end ) *end = '\0';
  end = payload+strlen(payload);
  while ( end > payload && isspace((unsigned char)end[-1]) ) --end;
  *end = '\0';
  return payload;
}

int MMGS_loadObjMesh(MMG5_pMesh mesh,const char *filename) {
  MMGS_ObjGroups groups = { NULL,NULL,0,0,1 };
  MMG5_int       np = 0,nt = 0,ip = 0,it = 0,currentRef = 0;
  MMG5_int       groupRef = 0,materialRef = 0,nvertex;
  size_t         capacity = 0;
  char           *line = NULL,*payload,*cursor,*token;
  FILE           *inm;
  int            status,hasGroup = 0;

  if ( !filename || !(inm = fopen(filename,"r")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",filename ? filename : "(null)");
    return 0;
  }

  /* Count vertices and the triangles produced by polygon fan triangulation. */
  while ( (status = MMGS_objReadLine(inm,&line,&capacity)) > 0 ) {
    if ( MMGS_objPayload(line,"v") ) {
      if ( np == MMG5_INTMAX ) {
        status = -1;
        break;
      }
      ++np;
    }
    else if ( (payload = MMGS_objPayload(line,"f")) ) {
      cursor = payload;
      nvertex = 0;
      while ( MMGS_objNextToken(&cursor) ) {
        if ( nvertex == MMG5_INTMAX ) {
          status = -1;
          break;
        }
        ++nvertex;
      }
      if ( status < 0 ) break;
      if ( nvertex < 3 ) {
        fprintf(stderr,"  ## Error: OBJ face with fewer than 3 vertices.\n");
        status = -1;
        break;
      }
      if ( nvertex-2 > MMG5_INTMAX-nt ) {
        status = -1;
        break;
      }
      nt += nvertex-2;
    }
    else if ( (payload = MMGS_objPayload(line,"g")) ||
              (payload = MMGS_objPayload(line,"usemtl")) ) {
      MMG5_int explicitRef;

      payload = MMGS_objGroupName(payload);
      /* Reserve every explicit reference before assigning fallbacks.  This
       * prevents an earlier arbitrary group from silently taking mmg_ref_N. */
      if ( MMGS_objExplicitRef(payload,&explicitRef) &&
           !MMGS_objGroupRef(&groups,payload,&explicitRef) ) {
        status = -1;
        break;
      }
    }
  }
  if ( status < 0 || ferror(inm) || !np || !nt ) {
    fprintf(stderr,"  ## Error: invalid or empty OBJ mesh in %s.\n",filename);
    fclose(inm);
    MMG5_SAFE_FREE(line);
    MMGS_objFreeGroups(&groups);
    return -1;
  }

  rewind(inm);
  if ( !MMGS_Set_meshSize(mesh,np,nt,0) ) {
    fclose(inm);
    MMG5_SAFE_FREE(line);
    MMGS_objFreeGroups(&groups);
    return -1;
  }

  while ( (status = MMGS_objReadLine(inm,&line,&capacity)) > 0 ) {
    if ( (payload = MMGS_objPayload(line,"v")) ) {
      double x,y,z;
      char   *end;

      errno = 0;
      x = strtod(payload,&end);
      if ( end == payload ) goto parse_error;
      payload = end;
      y = strtod(payload,&end);
      if ( end == payload ) goto parse_error;
      payload = end;
      z = strtod(payload,&end);
      if ( end == payload || errno == ERANGE || !isfinite(x) ||
           !isfinite(y) || !isfinite(z) ||
           !MMGS_Set_vertex(mesh,x,y,z,0,++ip) ) goto parse_error;
    }
    else if ( (payload = MMGS_objPayload(line,"g")) ) {
      payload = MMGS_objGroupName(payload);
      if ( *payload && strcmp(payload,"off") &&
           !MMGS_objGroupRef(&groups,payload,&groupRef) ) {
        goto memory_error;
      }
      hasGroup = *payload && strcmp(payload,"off");
      currentRef = hasGroup ? groupRef : materialRef;
    }
    else if ( (payload = MMGS_objPayload(line,"usemtl")) ) {
      payload = MMGS_objGroupName(payload);
      if ( *payload && strcmp(payload,"off") &&
           !MMGS_objGroupRef(&groups,payload,&materialRef) ) {
        goto memory_error;
      }
      if ( !*payload || !strcmp(payload,"off") ) materialRef = 0;
      /* OBJ groups remain authoritative. Materials provide a reference only
       * when no group is active, which preserves existing group semantics. */
      if ( !hasGroup ) currentRef = materialRef;
    }
    else if ( (payload = MMGS_objPayload(line,"f")) ) {
      MMG5_int first,previous,current;

      cursor = payload;
      token = MMGS_objNextToken(&cursor);
      if ( !token || !MMGS_objVertexIndex(token,ip,&first) ) goto parse_error;
      token = MMGS_objNextToken(&cursor);
      if ( !token || !MMGS_objVertexIndex(token,ip,&previous) ) {
        goto parse_error;
      }
      while ( (token = MMGS_objNextToken(&cursor)) ) {
        if ( !MMGS_objVertexIndex(token,ip,&current) ||
             !MMGS_Set_triangle(mesh,first,previous,current,currentRef,++it) ) {
          goto parse_error;
        }
        previous = current;
      }
    }
  }

  if ( status < 0 || ferror(inm) || ip != np || it != nt ) goto parse_error;
  fclose(inm);
  MMG5_SAFE_FREE(line);
  MMGS_objFreeGroups(&groups);
  MMG5_check_readedMesh(mesh,0);
  return 1;

parse_error:
  fprintf(stderr,"  ## Error: unable to parse OBJ mesh %s.\n",filename);
memory_error:
  fclose(inm);
  MMG5_SAFE_FREE(line);
  MMGS_objFreeGroups(&groups);
  return -1;
}

static int MMGS_objAddRef(MMG5_int *refs,MMG5_int *nref,MMG5_int ref) {
  MMG5_int i;

  for ( i=0; i<*nref; ++i ) {
    if ( refs[i] == ref ) return 1;
  }
  refs[(*nref)++] = ref;
  return 1;
}

static int MMGS_objMtlPath(const char *filename,char **path,
                           const char **basename) {
  const char *slash,*backslash,*separator;
  char       *extension;
  size_t     length = strlen(filename);

  if ( length > SIZE_MAX-5 ) return 0;
  MMG5_SAFE_MALLOC(*path,length+5,char,return 0);
  strcpy(*path,filename);
  slash = strrchr(*path,'/');
  backslash = strrchr(*path,'\\');
  separator = slash;
  if ( backslash && (!separator || backslash > separator) ) {
    separator = backslash;
  }
  extension = strrchr(*path,'.');
  if ( extension && (!separator || extension > separator) ) {
    strcpy(extension,".mtl");
  }
  else strcat(*path,".mtl");
  slash = strrchr(*path,'/');
  backslash = strrchr(*path,'\\');
  separator = slash;
  if ( backslash && (!separator || backslash > separator) ) {
    separator = backslash;
  }
  *basename = separator ? separator+1 : *path;
  return 1;
}

/** Generate stable display colors without attaching physical meaning. */
static void MMGS_objMaterialColor(MMG5_int ref,double color[3]) {
  uint64_t hash = (uint64_t)(int64_t)ref+UINT64_C(0x9e3779b97f4a7c15);
  int      i;

  hash = (hash^(hash >> 30))*UINT64_C(0xbf58476d1ce4e5b9);
  hash = (hash^(hash >> 27))*UINT64_C(0x94d049bb133111eb);
  hash ^= hash >> 31;
  for ( i=0; i<3; ++i ) {
    color[i] = .25+.65*(double)((hash >> (8*i)) & UINT64_C(255))/255.;
  }
}

static int MMGS_objWriteMtl(const char *filename,const MMG5_int *refs,
                            MMG5_int nref) {
  double color[3];
  FILE   *out;
  int    failed;
  MMG5_int i;

  if ( !(out=fopen(filename,"w")) ) return 0;
  fprintf(out,"# Wavefront material library written by MMGS\n");
  for ( i=0; i<nref; ++i ) {
    MMGS_objMaterialColor(refs[i],color);
    fprintf(out,"\nnewmtl mmg_ref_%" MMG5_PRId "\n",refs[i]);
    fprintf(out,"Ka %.8g %.8g %.8g\n",.2*color[0],.2*color[1],.2*color[2]);
    fprintf(out,"Kd %.8g %.8g %.8g\n",color[0],color[1],color[2]);
    fprintf(out,"Ks 0.1 0.1 0.1\nNs 10\nd 1\nillum 2\n");
  }
  failed = ferror(out);
  return !fclose(out) && !failed;
}

int MMGS_saveObjMesh(MMG5_pMesh mesh,const char *filename) {
  MMG5_pPoint ppt;
  MMG5_pTria  ptt;
  MMG5_int    *refs = NULL;
  MMG5_int    k,np = 0,nt = 0,nref = 0,lastRef = 0;
  char        *mtlPath = NULL;
  const char  *mtlName;
  FILE        *inm;
  int         hasLastRef = 0,failed;

  if ( !filename ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",filename ? filename : "(null)");
    return 0;
  }

  for ( k=1; k<=mesh->np; ++k ) {
    ppt = &mesh->point[k];
    ppt->tmp = 0;
    if ( MG_VOK(ppt) ) ppt->tmp = ++np;
  }
  for ( k=1; k<=mesh->nt; ++k ) {
    ptt = &mesh->tria[k];
    if ( MG_EOK(ptt) && mesh->point[ptt->v[0]].tmp &&
         mesh->point[ptt->v[1]].tmp && mesh->point[ptt->v[2]].tmp ) ++nt;
  }
  if ( nt ) MMG5_SAFE_MALLOC(refs,nt,MMG5_int,goto error);
  for ( k=1; k<=mesh->nt; ++k ) {
    ptt = &mesh->tria[k];
    if ( MG_EOK(ptt) && mesh->point[ptt->v[0]].tmp &&
         mesh->point[ptt->v[1]].tmp && mesh->point[ptt->v[2]].tmp ) {
      MMGS_objAddRef(refs,&nref,ptt->ref);
    }
  }
  if ( !MMGS_objMtlPath(filename,&mtlPath,&mtlName) ||
       !MMGS_objWriteMtl(mtlPath,refs,nref) || !(inm=fopen(filename,"w")) ) {
    fprintf(stderr,"  ** UNABLE TO WRITE OBJ/MTL OUTPUT %s.\n",filename);
    goto error;
  }

  fprintf(inm,"# Wavefront OBJ written by MMGS\n");
  fprintf(inm,"# vertices: %" MMG5_PRId ", triangles: %" MMG5_PRId "\n",np,nt);
  fprintf(inm,"mtllib %s\n",mtlName);
  for ( k=1; k<=mesh->np; ++k ) {
    ppt = &mesh->point[k];
    if ( ppt->tmp ) {
      fprintf(inm,"v %.17g %.17g %.17g\n",ppt->c[0],ppt->c[1],ppt->c[2]);
    }
  }

  for ( k=1; k<=mesh->nt; ++k ) {
    ptt = &mesh->tria[k];
    if ( !MG_EOK(ptt) || !mesh->point[ptt->v[0]].tmp ||
         !mesh->point[ptt->v[1]].tmp || !mesh->point[ptt->v[2]].tmp ) {
      continue;
    }
    if ( !hasLastRef || ptt->ref != lastRef ) {
      fprintf(inm,"g mmg_ref_%" MMG5_PRId "\n",ptt->ref);
      fprintf(inm,"usemtl mmg_ref_%" MMG5_PRId "\n",ptt->ref);
      lastRef = ptt->ref;
      hasLastRef = 1;
    }
    fprintf(inm,"f %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
            mesh->point[ptt->v[0]].tmp,mesh->point[ptt->v[1]].tmp,
            mesh->point[ptt->v[2]].tmp);
  }

  failed = ferror(inm);
  MMG5_SAFE_FREE(refs);
  MMG5_SAFE_FREE(mtlPath);
  return !fclose(inm) && !failed;

error:
  MMG5_SAFE_FREE(refs);
  MMG5_SAFE_FREE(mtlPath);
  return 0;
}
