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
 * \param ip index of the bottom-left-front corner of the root cell
 * \param length length of the octree in each direction
 *
 * \return 1 if success, 0 if fail.
 *
 * Allocate and init an MOctree structure.
 *
 */
int MMG3D_init_MOctree  ( MMG5_pMesh mesh, MMG5_pMOctree q, int ip, double length[3] ) {

  q->length[0] = length[0];
  q->length[1] = length[1];
  q->length[2] = length[2];


  /** Check that we have enough memory to allocate a new cell */
  MMG5_ADD_MEM(mesh,sizeof(MMG5_MOctree_s),"initial MOctree cell",
               return 0);

  /** New cell allocation */
  MMG5_SAFE_MALLOC( q->root,1, MMG5_MOctree_s, return 0);
  MMG3D_init_MOctree_s( mesh,q->root,ip,1,0);
  q->root->nsons = 0;

  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MOctree cell
 * \param ip index of the bottom-left-front corner of the cell
 * \param depth cell's depth
 * \param split_ls 1 if the cell is intersected by the level-set
 *
 * \return 1 if success, 0 if fail.
 *
 * Allocate and init an MOctree Cell
 *
 */
int MMG3D_init_MOctree_s( MMG5_pMesh mesh, MMG5_MOctree_s* q,int ip, int depth,int8_t split_ls ) {

  q->father = NULL;
  q->sons   = NULL;
  q->nsons  = 0;

  q->blf_ip = ip;
  q->depth  = depth;

  q->split_ls  = split_ls;
  q->leaf = 0;
  q->coordoct[0]=0;
  q->coordoct[1]=0;
  q->coordoct[2]=0;
  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MOctree cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Set the parameter split_ls of a cell to 1 if the level-set intersect the cell for every cell of the octree max.
 *
 */
int  MMG3D_set_splitls_MOctree ( MMG5_pMesh mesh, MMG5_MOctree_s* q, MMG5_pSol sol, double max_distance) {

  if(sol->m[q->blf_ip] <= max_distance)
  {
    q->split_ls=1;
    // if(q->depth != 0)
    // {
    //   q->father->split_ls=1;
    // }
  }
  else
  {
    q->split_ls=0;
  }

  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MOctree cell
 * \param depth_max the depth maximum of the octree.

 *
 * \return 1 if success, 0 if fail.
 *
 * Split an MOctree cell \ref q into \ref MMG3D_SIZE_OCTREESONS MOctree Cells.
 * \ref q must be a leaf.
 *
 */
int  MMG3D_split_MOctree_s ( MMG5_pMesh mesh, MMG5_MOctree_s* q, int depth_max, MMG5_pSol sol, double max_distance) {

  // double dx = mesh->info.max[0];
  // double dy = mesh->info.max[1];
  // double dz = mesh->info.max[2];

  int ip;
  int i;
  MMG5_SAFE_MALLOC(q->sons,q->nsons, MMG5_MOctree_s, return 0);
  for(i=0; i<q->nsons; i++)
  {
    MMG3D_init_MOctree_s(mesh, &q->sons[i], 0, q->depth + 1, 0);
    q->sons[i].father = q;
    //calculus of octree coordinates and ip

    q->sons[i].coordoct[0]=q->coordoct[0];
    q->sons[i].coordoct[1]=q->coordoct[1];
    q->sons[i].coordoct[2]=q->coordoct[2];
    int power = (2^(depth_max)/2^(q->sons[i].depth));
    if(i==1)
    {
      q->sons[i].coordoct[0]+=power;
    }
    else if(i==2)
    {
      q->sons[i].coordoct[2]+=power;
    }
    else if(i==3)
    {
      q->sons[i].coordoct[0]+=power;
      q->sons[i].coordoct[2]+=power;
    }
    else if(i==4)
    {
      q->sons[i].coordoct[1]+=power;
    }
    else if(i==5)
    {
      q->sons[i].coordoct[0]+=power;
      q->sons[i].coordoct[1]+=power;
    }
    else if(i==6)
    {
      q->sons[i].coordoct[1]+=power;
      q->sons[i].coordoct[2]+=power;
    }
    else if(i==7)
    {
      q->sons[i].coordoct[0]+=power;
      q->sons[i].coordoct[1]+=power;
      q->sons[i].coordoct[2]+=power;
    }

    if(q->sons[i].depth < depth_max)
    {
      q->sons[i].nsons = 8;
      MMG3D_split_MOctree_s(mesh, &q->sons[i], depth_max, sol, max_distance);
    }
    else{
      /*calculus of ip for leaves*/
      q->sons[i].blf_ip=q->sons[i].coordoct[2]*pow(2,2*(depth_max-1))+q->sons[i].coordoct[1]*pow(2,depth_max-1)+q->sons[i].coordoct[0]+1;
      //moins couteux ?
      //q->sons[i].blf_ip=q->sons[i].coordoct[2]*dx*dy+q->sons[i].coordoct[1]*dx+q->sons[i].coordoct[0]+1;
      q->sons[i].leaf=1;
      MMG3D_set_splitls_MOctree (mesh, &q->sons[i], sol, max_distance);
    }
  }
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
int MMG3D_free_MOctree  ( MMG5_pMOctree** q, MMG5_pMesh mesh) {
  MMG5_DEL_MEM(mesh,*q);
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
int MMG3D_free_MOctree_s( MMG5_MOctree_s* q, MMG5_pMesh mesh) {
  MMG5_DEL_MEM(mesh,q);
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
int  MMG3D_merge_MOctree_s ( MMG5_MOctree_s* q, MMG5_pMesh mesh) {

  q->nsons = 0;
  q->leaf = 1;
  q->blf_ip=q->sons[0].blf_ip;//récupérer la ls à l'origine de la cellule
  MMG3D_free_MOctree_s(q->sons, mesh);

  return 1;
}
