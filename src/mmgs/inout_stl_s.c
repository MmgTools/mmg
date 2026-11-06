/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file inout_stl_s.c
 *  \brief ASCII/binary STL input and binary STL output for MMGS.
 */

#include "libmmgs.h"
#include "libmmgs_private.h"

#include <float.h>
#include <math.h>

typedef struct {
  double *coordinates;
  size_t size;
  size_t capacity;
} MMGS_StlFacets;

static uint32_t MMGS_stlUint32(const unsigned char *bytes) {
  return (uint32_t)bytes[0] | ((uint32_t)bytes[1] << 8) |
    ((uint32_t)bytes[2] << 16) | ((uint32_t)bytes[3] << 24);
}

static float MMGS_stlFloat(const unsigned char *bytes) {
  uint32_t bits = MMGS_stlUint32(bytes);
  float value;

  memcpy(&value,&bits,sizeof(value));
  return value;
}

static int MMGS_stlAppend(MMGS_StlFacets *facets,const double coordinates[9]) {
  if ( facets->size == facets->capacity ) {
    size_t oldCapacity = facets->capacity;
    size_t newCapacity;

    if ( oldCapacity > SIZE_MAX/2 ) return 0;
    newCapacity = oldCapacity ? 2*oldCapacity : 256;
    if ( newCapacity > SIZE_MAX/(9*sizeof(double)) ) return 0;
    MMG5_SAFE_REALLOC(facets->coordinates,9*oldCapacity,9*newCapacity,
                      double,"STL facets",return 0);
    facets->capacity = newCapacity;
  }
  memcpy(facets->coordinates+9*facets->size,coordinates,9*sizeof(double));
  ++facets->size;
  return 1;
}

static int MMGS_stlReadBinary(FILE *inm,uint32_t nfacet,
                              MMGS_StlFacets *facets) {
  unsigned char record[50];
  double        coordinates[9];
  uint32_t      k;
  int           i;

  /* An STL with no facets cannot describe an MMGS mesh. */
  if ( !nfacet || fseek(inm,84,SEEK_SET) ) return 0;
  for ( k=0; k<nfacet; ++k ) {
    if ( fread(record,sizeof(record),1,inm) != 1 ) return 0;
    for ( i=0; i<9; ++i ) {
      coordinates[i] = (double)MMGS_stlFloat(record+12+4*i);
      if ( !isfinite(coordinates[i]) ) return 0;
    }
    if ( !MMGS_stlAppend(facets,coordinates) ) return 0;
  }
  return 1;
}

static int MMGS_stlReadAscii(FILE *inm,MMGS_StlFacets *facets) {
  char   word[256];
  double coordinates[9];
  int    nvertex = 0;

  rewind(inm);
  while ( fscanf(inm,"%255s",word) == 1 ) {
    if ( strcmp(word,"vertex") ) continue;
    if ( fscanf(inm,"%lf %lf %lf",&coordinates[3*nvertex],
                &coordinates[3*nvertex+1],&coordinates[3*nvertex+2]) != 3 ) {
      return 0;
    }
    if ( !isfinite(coordinates[3*nvertex]) ||
         !isfinite(coordinates[3*nvertex+1]) ||
         !isfinite(coordinates[3*nvertex+2]) ) return 0;
    if ( ++nvertex == 3 ) {
      if ( !MMGS_stlAppend(facets,coordinates) ) return 0;
      nvertex = 0;
    }
  }
  return !ferror(inm) && !nvertex && facets->size;
}

static uint64_t MMGS_stlHash(int64_t x,int64_t y,int64_t z) {
  uint64_t hash = (uint64_t)x*UINT64_C(0x9e3779b185ebca87);

  hash ^= (uint64_t)y*UINT64_C(0xc2b2ae3d27d4eb4f);
  hash ^= (uint64_t)z*UINT64_C(0x165667b19e3779f9);
  hash ^= hash >> 33;
  return hash;
}

/** Weld facet vertices and transfer the resulting indexed mesh to MMGS. */
static int MMGS_stlBuildMesh(MMG5_pMesh mesh,const MMGS_StlFacets *facets) {
  double      *points = NULL,tolerance,min[3],max[3],scale = 0.0;
  MMG5_int    *triangles = NULL,*heads = NULL,*next = NULL,*map = NULL;
  int64_t     *cells = NULL,q[3];
  unsigned char *used = NULL;
  size_t      nraw,nhash,i,j,npoint = 0;
  MMG5_int    ntriangle = 0,nused = 0;
  int         axis,dx,dy,dz;

  if ( !facets->size || !facets->coordinates ||
       facets->size > SIZE_MAX/3 || facets->size > (size_t)MMG5_INTMAX/3 ) {
    return -1;
  }
  nraw = 3*facets->size;
  if ( nraw > (SIZE_MAX-1)/2 ) return -1;
  nhash = 2*nraw+1;

  for ( axis=0; axis<3; ++axis ) {
    min[axis] = max[axis] = facets->coordinates[axis];
  }
  for ( i=0; i<nraw; ++i ) {
    for ( axis=0; axis<3; ++axis ) {
      double value = facets->coordinates[3*i+axis];
      min[axis] = MG_MIN(min[axis],value);
      max[axis] = MG_MAX(max[axis],value);
    }
  }
  /* Base welding on the mesh extent, not its distance from the origin.  An
   * absolute-coordinate scale would merge unrelated vertices after a large
   * translation of an otherwise unchanged mesh. */
  for ( axis=0; axis<3; ++axis ) {
    double extent = max[axis]-min[axis];

    if ( !isfinite(extent) ) return -1;
    scale = MG_MAX(scale,extent);
  }
  /* A unit-sized lower bound would collapse otherwise valid meshes whose
   * coordinates and extent are both very small.  Degenerate, zero-extent
   * facet sets cannot produce a surface mesh and would also make the spatial
   * hash quantization divide by zero. */
  if ( scale <= 0.0 ) return -1;
  tolerance = 64.0*DBL_EPSILON*scale;
  if ( tolerance <= 0.0 ) return -1;

  MMG5_SAFE_MALLOC(points,3*nraw,double,goto memory_error);
  MMG5_SAFE_MALLOC(triangles,nraw,MMG5_int,goto memory_error);
  MMG5_SAFE_MALLOC(heads,nhash,MMG5_int,goto memory_error);
  MMG5_SAFE_MALLOC(next,nraw,MMG5_int,goto memory_error);
  MMG5_SAFE_MALLOC(cells,3*nraw,int64_t,goto memory_error);
  MMG5_SAFE_CALLOC(used,nraw,unsigned char,goto memory_error);
  MMG5_SAFE_MALLOC(map,nraw,MMG5_int,goto memory_error);
  for ( i=0; i<nhash; ++i ) heads[i] = -1;

  for ( i=0; i<nraw; ++i ) {
    MMG5_int found = -1;

    for ( axis=0; axis<3; ++axis ) {
      q[axis] = (int64_t)floor((facets->coordinates[3*i+axis]-min[axis]) /
                                tolerance);
    }
    for ( dx=-1; dx<=1 && found < 0; ++dx ) {
      for ( dy=-1; dy<=1 && found < 0; ++dy ) {
        for ( dz=-1; dz<=1 && found < 0; ++dz ) {
          size_t bucket = MMGS_stlHash(q[0]+dx,q[1]+dy,q[2]+dz)%nhash;
          MMG5_int candidate;

          for ( candidate=heads[bucket]; candidate>=0; candidate=next[candidate] ) {
            if ( cells[3*candidate] != q[0]+dx ||
                 cells[3*candidate+1] != q[1]+dy ||
                 cells[3*candidate+2] != q[2]+dz ) continue;
            if ( fabs(points[3*candidate]-facets->coordinates[3*i]) <= tolerance &&
                 fabs(points[3*candidate+1]-facets->coordinates[3*i+1]) <= tolerance &&
                 fabs(points[3*candidate+2]-facets->coordinates[3*i+2]) <= tolerance ) {
              found = candidate;
              break;
            }
          }
        }
      }
    }
    if ( found < 0 ) {
      size_t bucket = MMGS_stlHash(q[0],q[1],q[2])%nhash;

      found = (MMG5_int)npoint++;
      memcpy(points+3*found,facets->coordinates+3*i,3*sizeof(double));
      memcpy(cells+3*found,q,3*sizeof(int64_t));
      next[found] = heads[bucket];
      heads[bucket] = found;
    }
    triangles[i] = found;
  }

  for ( i=0; i<facets->size; ++i ) {
    MMG5_int *triangle = triangles+3*i;
    if ( triangle[0] == triangle[1] || triangle[1] == triangle[2] ||
         triangle[2] == triangle[0] ) continue;
    ++ntriangle;
    used[triangle[0]] = used[triangle[1]] = used[triangle[2]] = 1;
  }
  for ( i=0; i<npoint; ++i ) map[i] = used[i] ? ++nused : 0;
  if ( !nused || !ntriangle || !MMGS_Set_meshSize(mesh,nused,ntriangle,0) ) {
    goto data_error;
  }
  for ( i=0; i<npoint; ++i ) {
    if ( map[i] && !MMGS_Set_vertex(mesh,points[3*i],points[3*i+1],
                                     points[3*i+2],0,map[i]) ) goto data_error;
  }
  j = 0;
  for ( i=0; i<facets->size; ++i ) {
    MMG5_int *triangle = triangles+3*i;
    if ( triangle[0] == triangle[1] || triangle[1] == triangle[2] ||
         triangle[2] == triangle[0] ) continue;
    if ( !MMGS_Set_triangle(mesh,map[triangle[0]],map[triangle[1]],
                            map[triangle[2]],0,++j) ) goto data_error;
  }

  MMG5_SAFE_FREE(points); MMG5_SAFE_FREE(triangles); MMG5_SAFE_FREE(heads);
  MMG5_SAFE_FREE(next); MMG5_SAFE_FREE(cells); MMG5_SAFE_FREE(used);
  MMG5_SAFE_FREE(map);
  MMG5_check_readedMesh(mesh,0);
  return 1;

memory_error:
  fprintf(stderr,"  ## Error: unable to allocate the STL mesh.\n");
data_error:
  MMG5_SAFE_FREE(points); MMG5_SAFE_FREE(triangles); MMG5_SAFE_FREE(heads);
  MMG5_SAFE_FREE(next); MMG5_SAFE_FREE(cells); MMG5_SAFE_FREE(used);
  MMG5_SAFE_FREE(map);
  return -1;
}

int MMGS_loadStlMesh(MMG5_pMesh mesh,const char *filename) {
  MMGS_StlFacets facets = { NULL,0,0 };
  unsigned char header[84];
  uint32_t      nfacet = 0;
  long          fileSize;
  FILE          *inm;
  int           binary = 0,ier;

  if ( !filename || !(inm = fopen(filename,"rb")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",filename ? filename : "(null)");
    return 0;
  }
  if ( !fseek(inm,0,SEEK_END) ) fileSize = ftell(inm);
  else fileSize = -1;
  rewind(inm);
  if ( fileSize >= 84 && fread(header,sizeof(header),1,inm) == 1 ) {
    nfacet = MMGS_stlUint32(header+80);
    binary = (uint64_t)fileSize == UINT64_C(84)+UINT64_C(50)*nfacet;
  }

  ier = binary ? MMGS_stlReadBinary(inm,nfacet,&facets) :
    MMGS_stlReadAscii(inm,&facets);
  fclose(inm);
  if ( !ier ) {
    fprintf(stderr,"  ## Error: unable to parse STL mesh %s.\n",filename);
    MMG5_SAFE_FREE(facets.coordinates);
    return -1;
  }
  ier = MMGS_stlBuildMesh(mesh,&facets);
  MMG5_SAFE_FREE(facets.coordinates);
  return ier;
}

static void MMGS_stlNormal(MMG5_pMesh mesh,MMG5_pTria triangle,double normal[3]) {
  const double *a = mesh->point[triangle->v[0]].c;
  const double *b = mesh->point[triangle->v[1]].c;
  const double *c = mesh->point[triangle->v[2]].c;
  double u[3] = {b[0]-a[0],b[1]-a[1],b[2]-a[2]};
  double v[3] = {c[0]-a[0],c[1]-a[1],c[2]-a[2]};
  double length;

  normal[0] = u[1]*v[2]-u[2]*v[1];
  normal[1] = u[2]*v[0]-u[0]*v[2];
  normal[2] = u[0]*v[1]-u[1]*v[0];
  length = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  if ( length > 0.0 ) {
    normal[0] /= length; normal[1] /= length; normal[2] /= length;
  }
}

static int MMGS_stlWriteUint32(FILE *out,uint32_t value) {
  unsigned char bytes[4] = {(unsigned char)value,(unsigned char)(value>>8),
    (unsigned char)(value>>16),(unsigned char)(value>>24)};
  return fwrite(bytes,sizeof(bytes),1,out) == 1;
}

static int MMGS_stlWriteFloat(FILE *out,double value) {
  float         converted = (float)value;
  uint32_t      bits;
  unsigned char bytes[4];

  if ( !isfinite(converted) ) return 0;
  memcpy(&bits,&converted,sizeof(bits));
  bytes[0] = (unsigned char)bits; bytes[1] = (unsigned char)(bits>>8);
  bytes[2] = (unsigned char)(bits>>16); bytes[3] = (unsigned char)(bits>>24);
  return fwrite(bytes,sizeof(bytes),1,out) == 1;
}

int MMGS_saveStlMesh(MMG5_pMesh mesh,const char *filename) {
  unsigned char header[80] = {0},attribute[2] = {0,0};
  MMG5_pTria   triangle;
  MMG5_int     k,ntriangle = 0;
  FILE         *out;

  if ( !filename ) return 0;
  if ( !(out = fopen(filename,"wb")) ) return 0;
  for ( k=1; k<=mesh->nt; ++k ) if ( MG_EOK(&mesh->tria[k]) ) ++ntriangle;

  if ( ntriangle > UINT32_MAX ) { fclose(out); return 0; }
  snprintf((char *)header,sizeof(header),"Binary STL written by MMGS");
  if ( fwrite(header,sizeof(header),1,out) != 1 ||
       !MMGS_stlWriteUint32(out,(uint32_t)ntriangle) ) {
    fclose(out); return 0;
  }
  for ( k=1; k<=mesh->nt; ++k ) {
    double normal[3];
    int i,j;

    triangle = &mesh->tria[k];
    if ( !MG_EOK(triangle) ) continue;
    MMGS_stlNormal(mesh,triangle,normal);
    for ( i=0; i<3; ++i ) if ( !MMGS_stlWriteFloat(out,normal[i]) ) goto error;
    for ( i=0; i<3; ++i ) {
      const double *point = mesh->point[triangle->v[i]].c;
      for ( j=0; j<3; ++j ) if ( !MMGS_stlWriteFloat(out,point[j]) ) goto error;
    }
    if ( fwrite(attribute,sizeof(attribute),1,out) != 1 ) goto error;
  }
  return !fclose(out);

error:
  fclose(out);
  return 0;
}
