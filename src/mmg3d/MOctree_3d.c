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
 * \file MOctree_3d.c
 * \brief Tools to manage an octree mesh.
 * \author Juliette Busquet (Enseirb-Matmeca)
 * \author Antoine Huc (Enseirb-Matmeca)
 * \author Fanny Kuhn (Enseirb-Matmeca)
 * \author Algiane Froehly (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Functions related to the creation of an octree from a grid and the associated
 * signed distance function. The octree cells are allowed to merge if they do
 * not cross the 0 levels-set of the signed distance function. At the end, the
 * octree must be 2:1 balanced (only 1 depth level of difference betwwen to
 * adjacent cells).
 *
 */

#include "mmg3d.h"

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MOctree cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Allocate and init an MOctree structure.
 *
 */
int MMG3D_init_MOctree  ( MMG5_pMesh mesh, MMG5_pMOctree* q, int depth_max ) {

  printf ( " %s:%s: TO IMPLEMENT\n",__FILE__,__func__ ); return 0;

  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MOctree cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Allocate and init an MOctree Cell
 *
 */
int MMG3D_init_MOctree_s( MMG5_pMesh mesh, MMG5_MOctree_s* q ) {

  printf ( " %s:%s: TO IMPLEMENT\n",__FILE__,__func__ ); return 0;

  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MOctree cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Split an MOctree cell \ref q into \ref MMG3D_SIZE_OCTREESONS MOctree Cells.
 * \ref q must be a leaf.
 *
 */
int  MMG3D_split_MOctree_s ( MMG5_pMesh mesh, MMG5_MOctree_s* q ) {

  printf ( " %s:%s: TO IMPLEMENT\n",__FILE__,__func__ ); return 0;

  return 1;
}

/**
 * \param q pointer toward the MOctree
 *
 * \return 1 if success, 0 if fail.
 *
 * Free a MOctree.
 *
 */
int MMG3D_free_MOctree  ( MMG5_pMOctree** q ) {

  printf ( " %s:%s: TO IMPLEMENT\n",__FILE__,__func__ ); return 0;

  return 1;
}

/**
 * \param q pointer toward the MOctree cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Free a MOctree cell.
 *
 */
int MMG3D_free_MOctree_s( MMG5_MOctree_s* q ) {
  int i;

  printf ( " %s:%s: TO IMPLEMENT\n",__FILE__,__func__ ); return 0;

  return 1;
}

/**
 * \param q pointer toward the MOctree cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Merge the \ref q MOctree cell (remove the \ref MMG3D_SIZE OCTREESONS sons of
 * \ref q, all sons being leafs).
 *
 */
int  MMG3D_merge_MOctree_s ( MMG5_MOctree_s* q ) {

  printf ( " %s:%s: TO IMPLEMENT\n",__FILE__,__func__ ); return 0;

  return 1;
}
