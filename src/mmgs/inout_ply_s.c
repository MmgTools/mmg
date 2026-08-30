/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file inout_ply_s.c
 *  \brief ASCII and binary little-endian PLY input/output for MMGS.
 */

#include "libmmgs.h"
#include "libmmgs_private.h"

#include <errno.h>
#include <math.h>

#define MMGS_PLY_MAX_PROPERTIES 64

enum MMGS_PlyType {
  MMGS_PLY_I8,MMGS_PLY_U8,MMGS_PLY_I16,MMGS_PLY_U16,
  MMGS_PLY_I32,MMGS_PLY_U32,MMGS_PLY_I64,MMGS_PLY_U64,
  MMGS_PLY_F32,MMGS_PLY_F64,MMGS_PLY_INVALID
};

typedef struct {
  enum MMGS_PlyType type,countType;
  int isList;
  char name[64];
} MMGS_PlyProperty;

typedef struct {
  int binary;
  MMG5_int np,nface;
  MMGS_PlyProperty vertex[MMGS_PLY_MAX_PROPERTIES];
  MMGS_PlyProperty face[MMGS_PLY_MAX_PROPERTIES];
  int nVertexProperties,nFaceProperties;
} MMGS_PlyHeader;

typedef struct {
  MMG5_int *vertices,*refs;
  size_t size,capacity;
} MMGS_PlyTriangles;

static enum MMGS_PlyType MMGS_plyType(const char *name) {
  if ( !strcmp(name,"char") || !strcmp(name,"int8") ) return MMGS_PLY_I8;
  if ( !strcmp(name,"uchar") || !strcmp(name,"uint8") ) return MMGS_PLY_U8;
  if ( !strcmp(name,"short") || !strcmp(name,"int16") ) return MMGS_PLY_I16;
  if ( !strcmp(name,"ushort") || !strcmp(name,"uint16") ) return MMGS_PLY_U16;
  if ( !strcmp(name,"int") || !strcmp(name,"int32") ) return MMGS_PLY_I32;
  if ( !strcmp(name,"uint") || !strcmp(name,"uint32") ) return MMGS_PLY_U32;
  if ( !strcmp(name,"int64") ) return MMGS_PLY_I64;
  if ( !strcmp(name,"uint64") ) return MMGS_PLY_U64;
  if ( !strcmp(name,"float") || !strcmp(name,"float32") ) return MMGS_PLY_F32;
  if ( !strcmp(name,"double") || !strcmp(name,"float64") ) return MMGS_PLY_F64;
  return MMGS_PLY_INVALID;
}

static int MMGS_plyIntegerType(enum MMGS_PlyType type) {
  return type >= MMGS_PLY_I8 && type <= MMGS_PLY_U64;
}

static int MMGS_plyHeader(FILE *inm,MMGS_PlyHeader *header) {
  char line[512],word[64],type[64],countType[64],name[64];
  int  element = 0,otherCount = 0,formatFound = 0;
  int  vertexFound = 0,faceFound = 0;

  if ( !fgets(line,sizeof(line),inm) || strncmp(line,"ply",3) ||
       (line[3] && line[3] != '\n' && line[3] != '\r') ) return 0;
  while ( fgets(line,sizeof(line),inm) ) {
    if ( sscanf(line,"%63s",word) != 1 ) continue;
    if ( !strcmp(word,"comment") || !strcmp(word,"obj_info") ) continue;
    if ( !strcmp(word,"format") ) {
      if ( sscanf(line,"format %63s",type) != 1 ) return 0;
      if ( !strcmp(type,"ascii") ) header->binary = 0;
      else if ( !strcmp(type,"binary_little_endian") ) header->binary = 1;
      else return 0;
      formatFound = 1;
    }
    else if ( !strcmp(word,"element") ) {
      long long count;
      if ( sscanf(line,"element %63s %lld",name,&count) != 2 || count < 0 ||
           count > MMG5_INTMAX ) return 0;
      if ( !strcmp(name,"vertex") ) {
        if ( vertexFound || faceFound ) return 0;
        element = 1; vertexFound = 1; header->np = (MMG5_int)count;
      }
      else if ( !strcmp(name,"face") ) {
        if ( !vertexFound || faceFound ) return 0;
        element = 2; faceFound = 1; header->nface = (MMG5_int)count;
      }
      else { element = 3; otherCount += count != 0; }
    }
    else if ( !strcmp(word,"property") ) {
      MMGS_PlyProperty *property;
      int *number;

      if ( element == 1 ) { property = header->vertex; number = &header->nVertexProperties; }
      else if ( element == 2 ) { property = header->face; number = &header->nFaceProperties; }
      else continue;
      if ( *number >= MMGS_PLY_MAX_PROPERTIES ) return 0;
      property += (*number)++;
      if ( sscanf(line,"property list %63s %63s %63s",countType,type,name) == 3 ) {
        property->isList = 1;
        property->countType = MMGS_plyType(countType);
        property->type = MMGS_plyType(type);
      }
      else if ( sscanf(line,"property %63s %63s",type,name) == 2 ) {
        property->isList = 0;
        property->type = MMGS_plyType(type);
        property->countType = MMGS_PLY_INVALID;
      }
      else return 0;
      if ( property->type == MMGS_PLY_INVALID ||
           (property->isList && property->countType == MMGS_PLY_INVALID) ) return 0;
      if ( element == 1 && property->isList ) return 0;
      if ( property->isList && !MMGS_plyIntegerType(property->countType) ) return 0;
      strncpy(property->name,name,sizeof(property->name)-1);
      property->name[sizeof(property->name)-1] = '\0';
    }
    else if ( !strcmp(word,"end_header") ) {
      return formatFound && !otherCount && header->np > 0 && header->nface > 0;
    }
  }
  return 0;
}

static int MMGS_plyReadBytes(FILE *inm,int n,uint64_t *bits) {
  unsigned char bytes[8];
  int i;
  *bits = 0;
  if ( fread(bytes,n,1,inm) != 1 ) return 0;
  for ( i=0; i<n; ++i ) *bits |= (uint64_t)bytes[i] << (8*i);
  return 1;
}

static int MMGS_plyReadValue(FILE *inm,int binary,enum MMGS_PlyType type,
                             double *real,int64_t *integer) {
  uint64_t bits;
  char token[256],*end;
  int nbyte;

  if ( !binary ) {
    if ( fscanf(inm,"%255s",token) != 1 ) return 0;
    errno = 0;
    if ( type == MMGS_PLY_F32 || type == MMGS_PLY_F64 ) {
      *real = strtod(token,&end);
      *integer = 0;
      return !errno && end != token && !*end && isfinite(*real);
    }
    else if ( type == MMGS_PLY_U8 || type == MMGS_PLY_U16 ||
              type == MMGS_PLY_U32 || type == MMGS_PLY_U64 ) {
      unsigned long long parsed;

      if ( token[0] == '-' ) return 0;
      parsed = strtoull(token,&end,10);
      if ( errno || end == token || *end || parsed > INT64_MAX ) return 0;
      *integer = (int64_t)parsed;
    }
    else {
      *integer = strtoll(token,&end,10);
    }
    if ( errno || end == token || *end ) return 0;
    if ( (type == MMGS_PLY_I8 && (*integer < INT8_MIN || *integer > INT8_MAX)) ||
         (type == MMGS_PLY_U8 && *integer > UINT8_MAX) ||
         (type == MMGS_PLY_I16 &&
          (*integer < INT16_MIN || *integer > INT16_MAX)) ||
         (type == MMGS_PLY_U16 && *integer > UINT16_MAX) ||
         (type == MMGS_PLY_I32 &&
          (*integer < INT32_MIN || *integer > INT32_MAX)) ||
         (type == MMGS_PLY_U32 && *integer > UINT32_MAX) ) return 0;
    *real = (double)*integer;
    return 1;
  }
  nbyte = (type == MMGS_PLY_I8 || type == MMGS_PLY_U8) ? 1 :
    (type == MMGS_PLY_I16 || type == MMGS_PLY_U16) ? 2 :
    (type == MMGS_PLY_I64 || type == MMGS_PLY_U64 ||
     type == MMGS_PLY_F64) ? 8 : 4;
  if ( !MMGS_plyReadBytes(inm,nbyte,&bits) ) return 0;
  switch ( type ) {
  case MMGS_PLY_I8:  *integer = (int8_t)bits; break;
  case MMGS_PLY_U8:  *integer = (uint8_t)bits; break;
  case MMGS_PLY_I16: *integer = (int16_t)bits; break;
  case MMGS_PLY_U16: *integer = (uint16_t)bits; break;
  case MMGS_PLY_I32: *integer = (int32_t)bits; break;
  case MMGS_PLY_U32: *integer = (uint32_t)bits; break;
  case MMGS_PLY_I64: *integer = (int64_t)bits; break;
  case MMGS_PLY_U64:
    if ( bits > INT64_MAX ) return 0;
    *integer = (int64_t)bits; break;
  case MMGS_PLY_F32: {
    uint32_t b = (uint32_t)bits;
    float    f;

    memcpy(&f,&b,sizeof(f));
    *real = f;
    *integer = 0;
    return isfinite(*real);
  }
  case MMGS_PLY_F64: {
    double d;

    memcpy(&d,&bits,sizeof(d));
    *real = d;
    *integer = 0;
    return isfinite(*real);
  }
  default: return 0;
  }
  *real = (double)*integer;
  return 1;
}

static int MMGS_plyAppendTriangle(MMGS_PlyTriangles *triangles,
                                  MMG5_int a,MMG5_int b,MMG5_int c,MMG5_int ref) {
  if ( triangles->size >= (size_t)MMG5_INTMAX ) return 0;
  if ( triangles->size == triangles->capacity ) {
    size_t old = triangles->capacity;
    size_t newCapacity;

    if ( old > SIZE_MAX/2 ) return 0;
    newCapacity = old ? 2*old : 256;
    if ( newCapacity > SIZE_MAX/(3*sizeof(MMG5_int)) ||
         newCapacity > SIZE_MAX/sizeof(MMG5_int) ) return 0;
    MMG5_SAFE_REALLOC(triangles->vertices,3*old,3*newCapacity,MMG5_int,
                      "PLY triangles",return 0);
    MMG5_SAFE_REALLOC(triangles->refs,old,newCapacity,MMG5_int,
                      "PLY references",return 0);
    triangles->capacity = newCapacity;
  }
  triangles->vertices[3*triangles->size] = a;
  triangles->vertices[3*triangles->size+1] = b;
  triangles->vertices[3*triangles->size+2] = c;
  triangles->refs[triangles->size++] = ref;
  return 1;
}

int MMGS_loadPlyMesh(MMG5_pMesh mesh,const char *filename) {
  MMGS_PlyHeader header = {0};
  MMGS_PlyTriangles triangles = {NULL,NULL,0,0};
  double *points = NULL,real;
  MMG5_int i,j,*polygon = NULL;
  size_t polygonCapacity = 0;
  int64_t integer;
  FILE *inm;
  int xProperty=-1,yProperty=-1,zProperty=-1,indexProperty=-1,refProperty=-1;

  if ( !filename || !(inm=fopen(filename,"rb")) ) return 0;
  if ( !MMGS_plyHeader(inm,&header) ) goto parse_error;
  for ( j=0; j<header.nVertexProperties; ++j ) {
    if ( !strcmp(header.vertex[j].name,"x") ) xProperty=j;
    else if ( !strcmp(header.vertex[j].name,"y") ) yProperty=j;
    else if ( !strcmp(header.vertex[j].name,"z") ) zProperty=j;
  }
  for ( j=0; j<header.nFaceProperties; ++j ) {
    if ( header.face[j].isList &&
         (!strcmp(header.face[j].name,"vertex_indices") ||
          !strcmp(header.face[j].name,"vertex_index")) ) indexProperty=j;
    else if ( !header.face[j].isList &&
              (!strcmp(header.face[j].name,"ref") ||
               !strcmp(header.face[j].name,"reference") ||
               !strcmp(header.face[j].name,"material_index")) ) refProperty=j;
  }
  if ( xProperty<0 || yProperty<0 || zProperty<0 || indexProperty<0 ) goto parse_error;
  if ( !MMGS_plyIntegerType(header.face[indexProperty].type) ||
       (refProperty >= 0 && !MMGS_plyIntegerType(header.face[refProperty].type)) ) {
    goto parse_error;
  }
  if ( (size_t)header.np > SIZE_MAX/(3*sizeof(double)) ) goto parse_error;
  MMG5_SAFE_MALLOC(points,3*(size_t)header.np,double,goto memory_error);
  for ( i=0; i<header.np; ++i ) {
    for ( j=0; j<header.nVertexProperties; ++j ) {
      if ( header.vertex[j].isList || !MMGS_plyReadValue(inm,header.binary,
           header.vertex[j].type,&real,&integer) ) goto parse_error;
      if ( j==xProperty ) points[3*i]=real;
      else if ( j==yProperty ) points[3*i+1]=real;
      else if ( j==zProperty ) points[3*i+2]=real;
    }
  }
  for ( i=0; i<header.nface; ++i ) {
    MMG5_int polygonSize=0,ref=0;
    for ( j=0; j<header.nFaceProperties; ++j ) {
      if ( header.face[j].isList ) {
        MMG5_int k,count;
        if ( !MMGS_plyReadValue(inm,header.binary,header.face[j].countType,
                                 &real,&integer) || integer<0 ) goto parse_error;
        if ( integer > MMG5_INTMAX ||
             (uint64_t)integer > SIZE_MAX/sizeof(MMG5_int) ) goto parse_error;
        count=(MMG5_int)integer;
        if ( j==indexProperty && (size_t)count>polygonCapacity ) {
          size_t old=polygonCapacity;
          MMG5_SAFE_REALLOC(polygon,old,(size_t)count,MMG5_int,
                            "PLY polygon",goto memory_error);
          polygonCapacity = (size_t)count;
        }
        for ( k=0; k<count; ++k ) {
          if ( !MMGS_plyReadValue(inm,header.binary,header.face[j].type,
                                   &real,&integer) ) goto parse_error;
          if ( j==indexProperty ) {
            if ( integer<0 || integer>=header.np ) goto parse_error;
            polygon[k]=(MMG5_int)integer+1;
          }
        }
        if ( j==indexProperty ) polygonSize=count;
      }
      else {
        if ( !MMGS_plyReadValue(inm,header.binary,header.face[j].type,
                                 &real,&integer) ) goto parse_error;
        if ( j==refProperty ) {
          if ( integer < -MMG5_INTMAX-1 || integer > MMG5_INTMAX ) {
            goto parse_error;
          }
          ref=(MMG5_int)integer;
        }
      }
    }
    if ( polygonSize<3 ) goto parse_error;
    for ( j=1; j<polygonSize-1; ++j ) {
      if ( !MMGS_plyAppendTriangle(&triangles,polygon[0],polygon[j],
                                    polygon[j+1],ref) ) goto memory_error;
    }
  }
  if ( triangles.size > (size_t)MMG5_INTMAX ||
       !MMGS_Set_meshSize(mesh,header.np,(MMG5_int)triangles.size,0) ) {
    goto memory_error;
  }
  for ( i=0; i<header.np; ++i ) {
    if ( !MMGS_Set_vertex(mesh,points[3*i],points[3*i+1],points[3*i+2],
                          0,i+1) ) goto parse_error;
  }
  for ( i=0; i<(MMG5_int)triangles.size; ++i ) {
    if ( !MMGS_Set_triangle(mesh,triangles.vertices[3*i],
                            triangles.vertices[3*i+1],
                            triangles.vertices[3*i+2],triangles.refs[i],
                            i+1) ) goto parse_error;
  }
  fclose(inm); MMG5_SAFE_FREE(points); MMG5_SAFE_FREE(polygon);
  MMG5_SAFE_FREE(triangles.vertices); MMG5_SAFE_FREE(triangles.refs);
  MMG5_check_readedMesh(mesh,0); return 1;

parse_error:
  fprintf(stderr,"  ## Error: unable to parse PLY mesh %s.\n",filename);
memory_error:
  fclose(inm); MMG5_SAFE_FREE(points); MMG5_SAFE_FREE(polygon);
  MMG5_SAFE_FREE(triangles.vertices); MMG5_SAFE_FREE(triangles.refs); return -1;
}

static int MMGS_plyWriteBytes(FILE *out,uint64_t value,int n) {
  unsigned char bytes[8];
  int           i;

  for ( i=0; i<n; ++i ) bytes[i] = (unsigned char)(value >> (8*i));
  return fwrite(bytes,n,1,out) == 1;
}

static int MMGS_plyWriteDouble(FILE *out,double value) {
  uint64_t bits;

  memcpy(&bits,&value,sizeof(bits));
  return MMGS_plyWriteBytes(out,bits,8);
}

int MMGS_savePlyMesh(MMG5_pMesh mesh,const char *filename) {
  MMG5_pPoint point;
  MMG5_pTria  triangle;
  MMG5_int    i,np=0,nt=0;
  const char  *extension;
  FILE        *out;
  int         ascii,wideIndex = 0,wideRef = 0;

  if ( !filename ) return 0;
  extension = strrchr(filename,'.');
  ascii = extension && !strcmp(extension,".plya");
  if ( !(out=fopen(filename,ascii ? "w" : "wb")) ) return 0;
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    point->tmp = MG_VOK(point) ? ++np : 0;
  }
  for ( i=1; i<=mesh->nt; ++i ) {
    triangle = &mesh->tria[i];
    if ( MG_EOK(triangle) && mesh->point[triangle->v[0]].tmp &&
         mesh->point[triangle->v[1]].tmp &&
         mesh->point[triangle->v[2]].tmp ) {
      ++nt;
      if ( sizeof(MMG5_int) > 4 &&
           (triangle->ref < INT32_MIN || triangle->ref > INT32_MAX) ) {
        wideRef = 1;
      }
    }
  }
  /* PLY's `int` type is exactly 32 bits.  Advertise and write 64-bit values
   * whenever an MMG5_int cannot be represented by that property type. */
  wideIndex = sizeof(MMG5_int) > 4 && np > INT32_MAX;
  fprintf(out,"ply\nformat %s 1.0\ncomment written by MMGS\n",
          ascii ? "ascii" : "binary_little_endian");
  fprintf(out,"element vertex %" MMG5_PRId "\nproperty double x\n"
          "property double y\nproperty double z\n",np);
  fprintf(out,"element face %" MMG5_PRId "\nproperty list uchar %s "
          "vertex_indices\nproperty %s ref\nend_header\n",nt,
          wideIndex ? "int64" : "int",wideRef ? "int64" : "int");
  for ( i=1; i<=mesh->np; ++i ) {
    point = &mesh->point[i];
    if ( !point->tmp ) continue;
    if ( ascii ) {
      fprintf(out,"%.17g %.17g %.17g\n",point->c[0],point->c[1],point->c[2]);
    }
    else if ( !MMGS_plyWriteDouble(out,point->c[0]) ||
              !MMGS_plyWriteDouble(out,point->c[1]) ||
              !MMGS_plyWriteDouble(out,point->c[2]) ) goto error;
  }
  for ( i=1; i<=mesh->nt; ++i ) {
    MMG5_int indices[3];
    int      j;

    triangle = &mesh->tria[i];
    if ( !MG_EOK(triangle) || !mesh->point[triangle->v[0]].tmp ||
         !mesh->point[triangle->v[1]].tmp ||
         !mesh->point[triangle->v[2]].tmp ) continue;
    for ( j=0; j<3; ++j ) {
      indices[j] = mesh->point[triangle->v[j]].tmp-1;
    }
    if ( ascii ) {
      fprintf(out,"3 %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId
              " %" MMG5_PRId "\n",indices[0],indices[1],indices[2],
              triangle->ref);
      continue;
    }
    if ( !MMGS_plyWriteBytes(out,3,1) ) goto error;
    for ( j=0; j<3; ++j ) {
      if ( !MMGS_plyWriteBytes(out,(uint64_t)indices[j],wideIndex ? 8 : 4) ) {
        goto error;
      }
    }
    if ( !MMGS_plyWriteBytes(out,(uint64_t)(int64_t)triangle->ref,
                             wideRef ? 8 : 4) ) goto error;
  }
  return !fclose(out);
error:
  fclose(out);
  return 0;
}
