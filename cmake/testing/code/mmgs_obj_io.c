/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file mmgs_obj_io.c
 *  \brief Test native OBJ input/output through the public MMGS API.
 */

#include "mmg/mmgs/libmmgs.h"

#include <stdio.h>
#include <string.h>

static int writeInput(const char *filename) {
  FILE *out = fopen(filename,"w");

  if ( !out ) return 0;
  fprintf(out,"# polygon triangulation and OBJ index syntax\n");
  fprintf(out,"v 0 0 0\n");
  fprintf(out,"v 1 0 0\n");
  fprintf(out,"v 1 1 0\n");
  fprintf(out,"v 0 1 0\n");
  fprintf(out,"v 0 0 1\n");
  fprintf(out,"mtllib ignored-input.mtl\n");
  fprintf(out,"usemtl skin\n");
  fprintf(out,"f 1/1 2/2 3/3 4/4\n");
  /* This explicit reference occurs after `skin`; it must still be reserved. */
  fprintf(out,"g mmg_ref_1\n");
  /* An active OBJ group is authoritative over a different material. */
  fprintf(out,"usemtl paint\n");
  fprintf(out,"f -5//1 -2//1 -1//1\n");
  fprintf(out,"g off\nusemtl mmg_ref_-7\n");
  fprintf(out,"f 2 3 5\n");
  return !fclose(out);
}

static int checkMesh(MMG5_pMesh mesh) {
  MMG5_int np,nt,na,v0,v1,v2,ref;
  int      required,nref1 = 0,nref2 = 0,nrefMinus7 = 0;

  if ( !MMGS_Get_meshSize(mesh,&np,&nt,&na) || np != 5 || nt != 4 || na ) {
    return 0;
  }
  while ( nt-- ) {
    if ( !MMGS_Get_triangle(mesh,&v0,&v1,&v2,&ref,&required) ) return 0;
    if ( ref == 1 ) ++nref1;
    else if ( ref == 2 ) ++nref2;
    else if ( ref == -7 ) ++nrefMinus7;
    else return 0;
  }
  return nref1 == 1 && nref2 == 2 && nrefMinus7 == 1;
}

static int checkOutputFiles(const char *objName,const char *mtlName) {
  const char *basename = strrchr(mtlName,'/');
  char       line[256],expected[256];
  FILE       *in;
  long long  ref;
  int        hasLibrary=0,hasUseMtl=0,nref1=0,nref2=0,nrefMinus7=0;

  basename = basename ? basename+1 : mtlName;
  if ( snprintf(expected,sizeof(expected),"mtllib %s",basename) < 0 ) return 0;
  if ( !(in=fopen(objName,"r")) ) return 0;
  while ( fgets(line,sizeof(line),in) ) {
    line[strcspn(line,"\r\n")] = '\0';
    if ( !strcmp(line,expected) ) hasLibrary = 1;
    if ( !strncmp(line,"usemtl mmg_ref_",15) ) hasUseMtl = 1;
  }
  if ( fclose(in) || !hasLibrary || !hasUseMtl || !(in=fopen(mtlName,"r")) ) {
    return 0;
  }
  while ( fgets(line,sizeof(line),in) ) {
    if ( sscanf(line,"newmtl mmg_ref_%lld",&ref) != 1 ) continue;
    if ( ref == 1 ) ++nref1;
    else if ( ref == 2 ) ++nref2;
    else if ( ref == -7 ) ++nrefMinus7;
    else {
      fclose(in);
      return 0;
    }
  }
  return !fclose(in) && nref1 == 1 && nref2 == 1 && nrefMinus7 == 1;
}

int main(int argc,char **argv) {
  MMG5_pMesh mesh = NULL;
  int        ier = 1;

  if ( argc != 4 || !writeInput(argv[1]) ) return 1;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);

  if ( MMGS_loadGenericMesh(mesh,NULL,NULL,argv[1]) != 1 ) goto cleanup;
  if ( !checkMesh(mesh) ) goto cleanup;
  if ( MMGS_saveGenericMesh(mesh,NULL,argv[2]) != 1 ||
       !checkOutputFiles(argv[2],argv[3]) ) goto cleanup;

  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  mesh = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  if ( MMGS_loadObjMesh(mesh,argv[2]) != 1 || !checkMesh(mesh) ) goto cleanup;
  ier = 0;

cleanup:
  MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_end);
  return ier;
}
