/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file mmgs/API_functions.c
 * \brief C API functions definitions for MMGS library.
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * \note This file contains some internal functions for the API, see
 * the \ref mmgs/libmmgs.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * C API for MMGS library.
 *
 */

#include "mmgs.h"

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 * \todo fill function
 */
void MMG5_MMGS_Init_parameters(MMG5_pMesh mesh) {
  printf("WARNING %s:%d: empty function, to fill.\n",__FILE__,__LINE__);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * Set the name of output mesh file.
 *
 */
int MMG5_MMGS_Set_outputMeshName(MMG5_pMesh mesh, char* meshout) {
  char *ptr;

  if ( mesh->nameout )
    _MMG5_DEL_MEM(mesh,mesh->nameout,(strlen(mesh->nameout)+1)*sizeof(char));

  if ( strlen(meshout) ) {
    _MMG5_ADD_MEM(mesh,(strlen(meshout)+1)*sizeof(char),"output mesh name",
                  printf("  Exit program.\n");
                  exit(EXIT_FAILURE));
    _MMG5_SAFE_CALLOC(mesh->nameout,strlen(meshout)+1,char);
    strcpy(mesh->nameout,meshout);
  }
  else {
    if ( strlen(mesh->namein) ) {
      _MMG5_ADD_MEM(mesh,(strlen(mesh->namein)+3)*sizeof(char),"output mesh name",
                    printf("  Exit program.\n");
                    exit(EXIT_FAILURE));
      _MMG5_SAFE_CALLOC(mesh->nameout,strlen(mesh->namein)+3,char);
      strcpy(mesh->nameout,mesh->namein);
      ptr = strstr(mesh->nameout,".mesh");
      if ( !ptr ) {
        /* filename without extension */
        strcat(mesh->nameout,".d");
      }
      else {
        *ptr = '\0';
        strcat(mesh->nameout,".d.mesh");
      }
      ptr = strstr(mesh->namein,".meshb");
      if ( ptr ) {
        /* filename with .meshb extention */
        strcat(mesh->nameout,"b");
      }
    }
    else {
      _MMG5_ADD_MEM(mesh,7*sizeof(char),"output mesh name",
                    printf("  Exit program.\n");
                    exit(EXIT_FAILURE));
      _MMG5_SAFE_CALLOC(mesh->nameout,7,char);
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
        fprintf(stdout,"  ## Warning: no name given for output mesh.\n");
        fprintf(stdout,"     Use of default value \"mesh.o.mesh\".\n");
      }
      strcpy(mesh->nameout,"mesh.o.mesh");
    }
  }
  return(1);
}
