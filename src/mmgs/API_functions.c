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
 * \todo try to remive paramters that do not coincide with mmg3d.
 */
void _MMG5_Init_parameters(MMG5_pMesh mesh) {

  memset(&mesh->info,0, sizeof(MMG5_Info));

  /* default values for integers */
  /** MMG5_IPARAM_verbose = -99 */
  mesh->info.imprim   =  -99;  /* [-10..10],Tune level of imprim */
  /** MMG5_IPARAM_mem = -1 */
  mesh->info.mem      = -1;  /* [n/-1]   ,Set memory size to n Mbytes/keep the default value */
  /** MMG5_IPARAM_debug = 0 */
  mesh->info.ddebug   =  0;  /* [0/1]    ,Turn on/off debug mode */
  /** MMG5_IPARAM_npar = 0 */
  mesh->info.npar     =  0;  /* [n]      ,number of local parameters */

  /* default values for doubles */
  /** MMG5_DPARAM_angleDetection = \ref _MMG5_ANGEDG */
  mesh->info.dhd      = _MMG5_ANGEDG;   /* angle detection; */
  /** MMG5_DPARAM_hmin = 0.0 */
  mesh->info.hmin     = 0.0;      /* minimal mesh size; */
  /** MMG5_DPARAM_hmax = \f$ \infty \f$ */
  mesh->info.hmax     = FLT_MAX;  /* maximal mesh size; */
  /** MMG5_DPARAM_hausd = 0.01 */
  mesh->info.hausd    = 0.01;     /* control Hausdorff */
  /** MMG5_DPARAM_hgrad = 0.1 */
  mesh->info.hgrad    = 0.1;      /* control gradation; */

  /* initial value for memMax and gap */
  mesh->gap = 0.2;
  mesh->memMax = _MMG5_memSize();
  if ( mesh->memMax )
    /* maximal memory = 50% of total physical memory */
    mesh->memMax = mesh->memMax*50/100;
  else {
    /* default value = 800 Mo */
    printf("  Maximum memory set to default value: %d Mo.\n",_MMG5_MEMMAX);
    mesh->memMax = _MMG5_MEMMAX << 20;
  }

}

/**
 * \param mesh pointer toward the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * Set the name of output mesh file.
 *
 */
int _MMG5_Set_outputMeshName(MMG5_pMesh mesh, char* meshout) {
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
