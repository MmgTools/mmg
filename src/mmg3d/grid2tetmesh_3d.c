/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
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
 * \file grid2tetmesh_3d.c
 * \brief Conversion of a structured grid into a tetrahedral mesh.
 * \author Juliette Busquet (Enseirb-Matmeca)
 * \author Antoine Huc (Enseirb-Matmeca)
 * \author Fanny Kuhn (Enseirb-Matmeca)
 * \author Algiane Froehly (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Conversion of a structured grid into a tetrahedral mesh.
 *
 */

#include "mmg3d.h"

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Create the initial octree from the strcuture grid implicitely loaded.
 *
 * Input data :
 * - mesh->freeint[i] is the number of cells in the i-direction
 * - mesh->info.min[i] stores the origin of the grid;
 * - mesh->info.max[i] stores the size of a cell in the i-direction
 *
 * \remark for now, we suppose that the grid is aligned with the canonical
 * directions
 */

static inline
int MMG3D_convert_grid2smallOctree(MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_MOctree_s *po;
  double         length[3];
  int            i,ip,depth_max;

  /** Step 1: Allocation and initialization of the octree root */
  /* Creation of the bottom-left-front corner of the root cell (grid origin) and
   * computation of the octree length */

  //  /!!!!!!!!!! IL FAUT LIRE LE DEPTH MAX!!!!!!!!!!!!!!!!
  depth_max=6;
  ip = 1;
  for ( i=0; i<3; ++i ) {
    mesh->point[ip].c[i] = mesh->info.min[i];
    // Remark : the level-set value associated to th point ip is stored in sol->m[ip].
    // It implies that you don't have the choice for the point numbering,
    // it will be in the same order than in the VTK file

    length[i] = mesh->info.max[i] * (double)mesh->freeint[i];
  }

  /* Computation of the octree length */
  /* Octree cell initialization */
  if ( !MMG3D_init_MOctree(mesh,mesh->octree,ip,length) ) return 0;
  po = mesh->octree->root;
  po->leaf=0;

  /** Step 2: Octree subdivision until reaching the grid size */
  MMG3D_split_MOctree_s (mesh, po, depth_max);

  /** Step 3: Giving the right split_ls to all the cells*/


  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Create a coarse octree from the initial octree.
 *
 */
static inline
int MMG3D_coarsen_octree(MMG5_pMesh mesh, MMG5_pSol sol) {

  printf ( " %s:%s: TO IMPLEMENT\n",__FILE__,__func__ ); return 0;

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure that contains the solution at
 * grid centroids at the beginning and at mesh nodes at the end.
 *
 * \return 1 if success, 0 if fail.
 *
 * Balance an unbalanced octree in order that 2 adjacent cells have at most 1
 * level of depth of difference (2:1 balancing).
 *
 */
static inline
int MMG3D_balance_octree(MMG5_pMesh mesh, MMG5_pSol sol) {

  printf ( " %s:%s: TO IMPLEMENT\n",__FILE__,__func__ ); return 0;

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Convert a balanced octree (2:1) into a tetrahedral mesh.
 *
 */
static inline
int MMG3D_convert_octree2tetmesh(MMG5_pMesh mesh, MMG5_pSol sol) {

  printf ( " %s:%s: TO IMPLEMENT\n",__FILE__,__func__ ); return 0;

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure that contains the solution at
 * grid centroids at the beginning and at mesh nodes at the end.
 *
 * \return 1 if success, 0 if fail.
 *
 * Convert a strcutured grid into an octree, then transform this octree into a
 * tetrahedral mesh.
 *
 */
int MMG3D_convert_grid2tetmesh(MMG5_pMesh mesh, MMG5_pSol sol) {

  /**--- stage 1: Octree computation */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** OCTREE INITIALIZATION\n");

  /* Conversion of the grid into an octree */
  if ( !MMG3D_convert_grid2smallOctree(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree initialization problem. Exit program.\n");
    return 0;
  }

  /* Creation of the coarse octree */
  if ( !MMG3D_coarsen_octree(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree coarsening problem. Exit program.\n");
    return 0;
  }

  /* Octree balancing */
  if ( !MMG3D_balance_octree(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree balancing problem. Exit program.\n");
    return 0;
  }

  /**--- stage 2: Tetrahedralization */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** OCTREE TETRAHEDRALIZATION\n");

  if ( !MMG3D_convert_octree2tetmesh(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree tetrahedralization problem. Exit program.\n");
    return 0;
  }

  return 1;
}
