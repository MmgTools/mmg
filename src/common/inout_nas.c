/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file inout_nas.c
 *  \brief Shared parser for low-order Nastran bulk-data meshes.
 */

#include "inout_nas.h"
#include "mmgcommon_private.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>

#define MMG5_NAS_MAX_FIELDS 24
#define MMG5_NAS_FIELD_SIZE 64

typedef struct {
  char name[9];
  char field[MMG5_NAS_MAX_FIELDS][MMG5_NAS_FIELD_SIZE];
  int  nfield;
  int  active;
} MMG5_NastranCard;

typedef struct {
  MMG5_int pid;
  MMG5_int ref;
} MMG5_NastranRef;

static char *MMG5_nasTrim(char *text) {
  char *end;

  while ( isspace((unsigned char)*text) ) ++text;
  end = text+strlen(text);
  while ( end > text && isspace((unsigned char)end[-1]) ) --end;
  *end = '\0';
  return text;
}

static int MMG5_nasReadLine(FILE *inm,char **line,size_t *capacity) {
  size_t length = 0;
  int    c;

  if ( !*line ) {
    *capacity = 256;
    MMG5_SAFE_MALLOC(*line,*capacity,char,return -1);
  }
  while ( (c=fgetc(inm)) != EOF ) {
    if ( length+1 >= *capacity ) {
      size_t old = *capacity;

      if ( old > SIZE_MAX/2 ) return -1;
      *capacity = 2*old;
      MMG5_SAFE_REALLOC(*line,old,*capacity,char,"Nastran line",return -1);
    }
    (*line)[length++] = (char)c;
    if ( c == '\n' ) break;
  }
  if ( c == EOF && !length ) return 0;
  (*line)[length] = '\0';
  return 1;
}

static void MMG5_nasUpper(char *text) {
  while ( *text ) {
    *text = (char)toupper((unsigned char)*text);
    ++text;
  }
}

/** Portable case-insensitive comparison for Nastran keywords. */
static int MMG5_nasEqualNoCase(const char *left,const char *right) {
  while ( *left && *right ) {
    if ( toupper((unsigned char)*left) != toupper((unsigned char)*right) ) {
      return 0;
    }
    ++left;
    ++right;
  }
  return !*left && !*right;
}

static int MMG5_nasKnownCard(const char *name) {
  return !strcmp(name,"GRID") || !strcmp(name,"CTRIA3") ||
         !strcmp(name,"CTRIAR") || !strcmp(name,"CQUAD4") ||
         !strcmp(name,"CQUADR") || !strcmp(name,"CTETRA") ||
         !strcmp(name,"CPENTA") || !strcmp(name,"CPYRAM") ||
         !strcmp(name,"CPYRA") || !strcmp(name,"CHEXA");
}

static int MMG5_nasUnsupportedCard(const char *name) {
  return !strcmp(name,"CTRIA6") || !strcmp(name,"CQUAD8") ||
         !strcmp(name,"CQUAD9");
}

static int MMG5_nasAppendField(MMG5_NastranCard *card,const char *begin,
                               size_t length) {
  char *trimmed;

  if ( card->nfield == MMG5_NAS_MAX_FIELDS ||
       length >= MMG5_NAS_FIELD_SIZE ) return 0;
  memcpy(card->field[card->nfield],begin,length);
  card->field[card->nfield][length] = '\0';
  trimmed = MMG5_nasTrim(card->field[card->nfield]);
  if ( trimmed != card->field[card->nfield] ) {
    memmove(card->field[card->nfield],trimmed,strlen(trimmed)+1);
  }
  ++card->nfield;
  return 1;
}

static int MMG5_nasInteger(const char *field,MMG5_int *value,int blank) {
  char      *end;
  long long parsed;

  if ( !*field ) {
    if ( !blank ) return 0;
    *value = 0;
    return 1;
  }
  errno = 0;
  parsed = strtoll(field,&end,10);
  while ( isspace((unsigned char)*end) ) ++end;
  if ( errno || end == field || *end ||
       (sizeof(MMG5_int) == 4 &&
        (parsed < INT32_MIN || parsed > INT32_MAX)) ) return 0;
  *value = (MMG5_int)parsed;
  return 1;
}

/** Parse standard, D-exponent, and exponent-without-E Nastran reals. */
static int MMG5_nasReal(const char *field,double *value) {
  char normalized[MMG5_NAS_FIELD_SIZE+2];
  char *end;
  int  i,length,exponent = -1;

  if ( !*field ) {
    *value = 0.;
    return 1;
  }
  length = (int)strlen(field);
  if ( length >= MMG5_NAS_FIELD_SIZE ) return 0;
  memcpy(normalized,field,(size_t)length+1);
  for ( i=0; i<length; ++i ) {
    if ( normalized[i] == 'd' || normalized[i] == 'D' ) normalized[i] = 'E';
  }
  errno = 0;
  *value = strtod(normalized,&end);
  if ( !errno && end != normalized && !*end ) return isfinite(*value);

  /* Nastran permits forms such as 1.25-3 and -.7+2. */
  for ( i=1; i<length; ++i ) {
    if ( (normalized[i] == '+' || normalized[i] == '-') &&
         normalized[i-1] != 'E' && normalized[i-1] != 'e' ) exponent = i;
  }
  if ( exponent < 0 || length+1 >= (int)sizeof(normalized) ) return 0;
  memmove(normalized+exponent+1,normalized+exponent,
          (size_t)(length-exponent)+1);
  normalized[exponent] = 'E';
  errno = 0;
  *value = strtod(normalized,&end);
  return !errno && end != normalized && !*end && isfinite(*value);
}

static int MMG5_nasElementInfo(const char *name,
                               enum MMG5_NastranElementType *type,
                               int *nvertex) {
  if ( !strcmp(name,"CTRIA3") || !strcmp(name,"CTRIAR") ) {
    *type = MMG5_NAS_CTRIA3; *nvertex = 3;
  }
  else if ( !strcmp(name,"CQUAD4") || !strcmp(name,"CQUADR") ) {
    *type = MMG5_NAS_CQUAD4; *nvertex = 4;
  }
  else if ( !strcmp(name,"CTETRA") ) {
    *type = MMG5_NAS_CTETRA; *nvertex = 4;
  }
  else if ( !strcmp(name,"CPENTA") ) {
    *type = MMG5_NAS_CPENTA; *nvertex = 6;
  }
  else if ( !strcmp(name,"CPYRAM") || !strcmp(name,"CPYRA") ) {
    *type = MMG5_NAS_CPYRAM; *nvertex = 5;
  }
  else if ( !strcmp(name,"CHEXA") ) {
    *type = MMG5_NAS_CHEXA; *nvertex = 8;
  }
  else return 0;
  return 1;
}

static int MMG5_nasAddPoint(MMG5_NastranMesh *nas,MMG5_int *capacity,
                            const MMG5_NastranCard *card) {
  MMG5_NastranPoint *point;
  MMG5_int          cp;

  if ( card->nfield < 5 ) return 0;
  if ( nas->np == *capacity ) {
    MMG5_int old = *capacity;

    if ( old > MMG5_INTMAX/2 ) return 0;
    *capacity = old ? 2*old : 128;
    MMG5_SAFE_REALLOC(nas->points,old,*capacity,MMG5_NastranPoint,
                      "Nastran points",return 0);
  }
  point = &nas->points[nas->np];
  if ( !MMG5_nasInteger(card->field[0],&point->id,0) || point->id <= 0 ||
       !MMG5_nasInteger(card->field[1],&cp,1) ||
       !MMG5_nasReal(card->field[2],&point->c[0]) ||
       !MMG5_nasReal(card->field[3],&point->c[1]) ||
       !MMG5_nasReal(card->field[4],&point->c[2]) ) return 0;
  if ( cp ) {
    fprintf(stderr,"  ## Error: non-basic Nastran GRID coordinates are "
            "unsupported.\n");
    return 0;
  }
  ++nas->np;
  return 1;
}

static int MMG5_nasAddElement(MMG5_NastranMesh *nas,MMG5_int *capacity,
                              const MMG5_NastranCard *card) {
  MMG5_NastranElement *element;
  int                 i,j,nvertex;

  if ( nas->ne == *capacity ) {
    MMG5_int old = *capacity;

    if ( old > MMG5_INTMAX/2 ) return 0;
    *capacity = old ? 2*old : 128;
    MMG5_SAFE_REALLOC(nas->elements,old,*capacity,MMG5_NastranElement,
                      "Nastran elements",return 0);
  }
  element = &nas->elements[nas->ne];
  if ( !MMG5_nasElementInfo(card->name,&element->type,&nvertex) ||
       card->nfield < nvertex+2 ||
       !MMG5_nasInteger(card->field[0],&element->id,0) || element->id <= 0 ||
       !MMG5_nasInteger(card->field[1],&element->ref,1) ) return 0;
  for ( i=0; i<nvertex; ++i ) {
    if ( !MMG5_nasInteger(card->field[i+2],&element->vertex[i],0) ||
         element->vertex[i] <= 0 ) return 0;
    for ( j=0; j<i; ++j ) {
      if ( element->vertex[j] == element->vertex[i] ) return 0;
    }
  }
  /* Solid fields after the corner nodes are midside nodes, which are outside
   * the deliberately low-order scope of this reader. */
  if ( element->type >= MMG5_NAS_CTETRA ) {
    for ( i=nvertex+2; i<card->nfield; ++i ) {
      if ( card->field[i][0] ) {
        fprintf(stderr,"  ## Error: quadratic Nastran solids are unsupported.\n");
        return 0;
      }
    }
  }
  ++nas->ne;
  return 1;
}

static int MMG5_nasFinishCard(MMG5_NastranMesh *nas,
                              MMG5_NastranCard *card,
                              MMG5_int *pointCapacity,
                              MMG5_int *elementCapacity) {
  int success = 1;

  if ( !card->active ) return 1;
  if ( !strcmp(card->name,"GRID") ) {
    success = MMG5_nasAddPoint(nas,pointCapacity,card);
  }
  else {
    success = MMG5_nasAddElement(nas,elementCapacity,card);
  }
  card->active = 0;
  card->nfield = 0;
  return success;
}

static int MMG5_nasRef(const MMG5_NastranRef *refs,MMG5_int nref,
                       MMG5_int pid,MMG5_int *ref) {
  MMG5_int i;

  for ( i=0; i<nref; ++i ) {
    if ( refs[i].pid == pid ) {
      *ref = refs[i].ref;
      return 1;
    }
  }
  *ref = pid;
  return 1;
}

static int MMG5_nasPointCompare(const void *left,const void *right) {
  const MMG5_NastranPoint *a = (const MMG5_NastranPoint *)left;
  const MMG5_NastranPoint *b = (const MMG5_NastranPoint *)right;

  return (a->id > b->id)-(a->id < b->id);
}

static int MMG5_nasMetadata(char *line,MMG5_NastranRef **refs,
                            MMG5_int *nref,MMG5_int *capacity) {
  MMG5_int pid,ref,i;
  long long rawPid,rawRef;
  char      extra;

  if ( strncmp(line,"$MMG_REF,",9) ) return 1;
  if ( sscanf(line,"$MMG_REF,%lld,%lld %c",&rawPid,&rawRef,&extra) != 2 ||
       rawPid <= 0 ||
       (sizeof(MMG5_int) == 4 &&
        (rawPid > INT32_MAX || rawRef < INT32_MIN || rawRef > INT32_MAX)) ) {
    return 0;
  }
  pid = (MMG5_int)rawPid;
  ref = (MMG5_int)rawRef;
  for ( i=0; i<*nref; ++i ) {
    if ( (*refs)[i].pid == pid ) return (*refs)[i].ref == ref;
  }
  if ( *nref == *capacity ) {
    MMG5_int old = *capacity;

    if ( old > MMG5_INTMAX/2 ) return 0;
    *capacity = old ? 2*old : 8;
    MMG5_SAFE_REALLOC(*refs,old,*capacity,MMG5_NastranRef,
                      "Nastran reference metadata",return 0);
  }
  (*refs)[*nref].pid = pid;
  (*refs)[*nref].ref = ref;
  ++(*nref);
  return 1;
}

static int MMG5_nasPhysicalCard(char *line,MMG5_NastranCard *card,
                                MMG5_NastranMesh *nas,
                                MMG5_int *pointCapacity,
                                MMG5_int *elementCapacity) {
  char   first[MMG5_NAS_FIELD_SIZE],*cursor,*comma,*comment,*trimmed;
  size_t length,lineLength,offset,width,count,i;
  int    continuation,large;

  comment = strchr(line,'$');
  if ( comment ) *comment = '\0';
  lineLength = strlen(line);
  while ( lineLength && (line[lineLength-1] == '\n' ||
                         line[lineLength-1] == '\r') ) {
    line[--lineLength] = '\0';
  }
  trimmed = MMG5_nasTrim(line);
  if ( !*trimmed ) return 1;
  comma = strchr(line,',');
  if ( comma ) {
    cursor = line;
    comma = strchr(cursor,',');
    length = comma ? (size_t)(comma-cursor) : strlen(cursor);
    if ( length >= sizeof(first) ) return 0;
    memcpy(first,cursor,length); first[length] = '\0';
    trimmed = MMG5_nasTrim(first);
    if ( trimmed != first ) memmove(first,trimmed,strlen(trimmed)+1);
    continuation = !first[0] || first[0] == '+' || first[0] == '*';
    if ( !continuation ) {
      if ( !MMG5_nasFinishCard(nas,card,pointCapacity,elementCapacity) ) {
        return 0;
      }
      MMG5_nasUpper(first);
      large = strlen(first) && first[strlen(first)-1] == '*';
      if ( large ) first[strlen(first)-1] = '\0';
      if ( MMG5_nasUnsupportedCard(first) ) {
        fprintf(stderr,"  ## Error: high-order Nastran shell card %s is "
                "unsupported.\n",first);
        return 0;
      }
      card->active = MMG5_nasKnownCard(first);
      card->nfield = 0;
      if ( card->active ) strcpy(card->name,first);
    }
    if ( !card->active ) return 1;
    cursor = comma ? comma+1 : cursor+length;
    while ( 1 ) {
      comma = strchr(cursor,',');
      length = comma ? (size_t)(comma-cursor) : strlen(cursor);
      if ( !MMG5_nasAppendField(card,cursor,length) ) return 0;
      if ( !comma ) break;
      cursor = comma+1;
    }
    return 1;
  }

  length = MG_MIN((size_t)8,lineLength);
  memcpy(first,line,length); first[length] = '\0';
  trimmed = MMG5_nasTrim(first);
  if ( trimmed != first ) memmove(first,trimmed,strlen(trimmed)+1);
  continuation = !first[0] || first[0] == '+' || first[0] == '*';
  large = first[0] == '*' ||
          (!continuation && strlen(first) && first[strlen(first)-1] == '*');
  if ( !continuation ) {
    if ( !MMG5_nasFinishCard(nas,card,pointCapacity,elementCapacity) ) return 0;
    MMG5_nasUpper(first);
    if ( large ) first[strlen(first)-1] = '\0';
    if ( MMG5_nasUnsupportedCard(first) ) {
      fprintf(stderr,"  ## Error: high-order Nastran shell card %s is "
              "unsupported.\n",first);
      return 0;
    }
    card->active = MMG5_nasKnownCard(first);
    card->nfield = 0;
    if ( card->active ) strcpy(card->name,first);
  }
  if ( !card->active ) return 1;
  width = large ? 16 : 8;
  count = large ? 4 : 8;
  for ( i=0; i<count; ++i ) {
    offset = 8+i*width;
    length = offset < lineLength ? MG_MIN(width,lineLength-offset) : 0;
    if ( !MMG5_nasAppendField(card,line+MG_MIN(offset,lineLength),length) ) {
      return 0;
    }
  }
  return 1;
}

int MMG5_loadNastranMeshData(const char *filename,MMG5_NastranMesh *nas) {
  MMG5_NastranCard card = {{0},{{0}},0,0};
  MMG5_NastranRef  *refs = NULL;
  MMG5_int         pointCapacity=0,elementCapacity=0,nref=0,refCapacity=0;
  MMG5_int         i;
  char             *line=NULL,*trimmed;
  size_t           capacity=0;
  FILE             *inm;
  int              status=1;

  memset(nas,0,sizeof(*nas));
  if ( !filename || !(inm=fopen(filename,"r")) ) return 0;
  while ( (status=MMG5_nasReadLine(inm,&line,&capacity)) > 0 ) {
    trimmed = MMG5_nasTrim(line);
    if ( !strncmp(trimmed,"$MMG_REF,",9) ) {
      if ( !MMG5_nasMetadata(trimmed,&refs,&nref,&refCapacity) ) {
        status = -1; break;
      }
      continue;
    }
    if ( *trimmed == '$' ) continue;
    if ( MMG5_nasEqualNoCase(trimmed,"ENDDATA") ) break;
    if ( !MMG5_nasPhysicalCard(line,&card,nas,&pointCapacity,
                               &elementCapacity) ) {
      status = -1; break;
    }
  }
  if ( status >= 0 &&
       !MMG5_nasFinishCard(nas,&card,&pointCapacity,&elementCapacity) ) {
    status = -1;
  }
  if ( status >= 0 && nas->np ) {
    qsort(nas->points,(size_t)nas->np,sizeof(MMG5_NastranPoint),
          MMG5_nasPointCompare);
    for ( i=1; i<nas->np; ++i ) {
      if ( nas->points[i-1].id == nas->points[i].id ) status = -1;
    }
    for ( i=0; i<nas->ne; ++i ) {
      MMG5_nasRef(refs,nref,nas->elements[i].ref,&nas->elements[i].ref);
    }
  }
  fclose(inm);
  MMG5_SAFE_FREE(line);
  MMG5_SAFE_FREE(refs);
  if ( status < 0 || !nas->np ) {
    fprintf(stderr,"  ## Error: unable to parse Nastran mesh %s.\n",filename);
    MMG5_freeNastranMeshData(nas);
    return -1;
  }
  return 1;
}

void MMG5_freeNastranMeshData(MMG5_NastranMesh *nas) {
  MMG5_SAFE_FREE(nas->points);
  MMG5_SAFE_FREE(nas->elements);
  nas->np = nas->ne = 0;
}

MMG5_int MMG5_nastranPointIndex(const MMG5_NastranMesh *nas,MMG5_int id) {
  MMG5_int begin=0,end=nas->np;

  while ( begin < end ) {
    MMG5_int middle = begin+(end-begin)/2;

    if ( nas->points[middle].id < id ) begin = middle+1;
    else end = middle;
  }
  return begin < nas->np && nas->points[begin].id == id ? begin+1 : 0;
}
