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
  int            i,ip,depth_max,max_dim;

  /** Step 1: Allocation and initialization of the octree root */
  /* Creation of the bottom-left-front corner of the root cell (grid origin) and
   * computation of the octree length */

  max_dim=0;
  ip = 1;
  for ( i=0; i<3; ++i ) {
    length[i] = mesh->info.max[i] * (double)mesh->freeint[i];
    if(max_dim < mesh->freeint[i])
    {
      max_dim = mesh->freeint[i];
    }
  }
  /* Begin to work on the dual grid => we will have one cellule less in each
   * direction */
  max_dim--;

  /* set max dim to the next power of 2 */
  max_dim--;
  max_dim |= max_dim >> 1;
  max_dim |= max_dim >> 2;
  max_dim |= max_dim >> 4;
  max_dim |= max_dim >> 8;
  max_dim |= max_dim >> 16;
  max_dim++;

  depth_max=log(max_dim)/log(2);

  /* Computation of the octree length */
  /* Octree cell initialization */
  if ( !MMG3D_init_MOctree(mesh,&mesh->octree,ip,length,depth_max) ) return 0;
  po = mesh->octree->root;

  if(po->depth != depth_max)
  {
    po->nsons = 8;
  }
  else {
    po->nsons = 0;
    po->leaf=1;
  }

  /** Step 2: Octree subdivision until reaching the grid size */
  MMG3D_split_MOctree_s (mesh, po, sol);

  return 1;
}

/**
 *\param mesh toward the mesh structure
 *\param q pointer toward the MOctree cell
 *\param depth_max the depth maximum of the octree.
 *
 * \return 1 if success, 0 if fail.
 *
 * Merge once the small octree if ls=0 is not in the small octree.
 *
 */
static inline
int MMG3D_build_coarsen_octree_first_time(MMG5_pMesh mesh, MMG5_MOctree_s* q, int depth_max) {
  int i;
  if (q->depth < depth_max-1)
  {
    for (i=0; i< q->nsons; i++)
    {
      MMG3D_build_coarsen_octree_first_time(mesh, &q->sons[i],depth_max);
    }
  }
  else if (q->depth == depth_max-1)
  {
    int sum_ls;
    sum_ls=0;
    for(i=0; i<q->nsons; i++)
    {
      sum_ls += q->sons[i].split_ls;
    }
    if (sum_ls==0)
    {
      MMG3D_merge_MOctree_s (q, mesh);
    }
  }
  return 1;
}

/**
 *\param q pointer toward the MOctree cell
 *\param depth_max the depth maximum of the octree.
 *\param depth  the depth of the parent that we want to see if we can merge his sons.
 * \return 1 if success, 0 if fail.
 *
 * Merge the sons of the parent whose depth=depth and balance to have only one depth at maximum beteween them.
 *
 */
static inline
int MMG3D_build_coarsen_octree(MMG5_pMesh mesh, MMG5_MOctree_s* q, int depth_max, int depth) {
  int i;
  if (q->depth < depth)
  {
    for (i=0; i< q->nsons; i++)
    {
      MMG3D_build_coarsen_octree(mesh, &q->sons[i], depth_max, depth);
    }
  }
  else
  {
    MMG5_MOctree_s* Neighbour;
    MMG5_ADD_MEM(mesh,sizeof(MMG5_MOctree_s),"MOctree neighbour",
                 return 0);
    MMG5_SAFE_MALLOC(Neighbour,1, MMG5_MOctree_s, return 0);
    MMG3D_init_MOctree_s(mesh, Neighbour, 0, 0, 0 );

    int sum_leaf;
    sum_leaf=0;
    for(i=0; i<q->nsons; i++)
    {
      sum_leaf += q->sons[i].leaf;
    }

    if (sum_leaf==8)
    {
      for(i=0; i<q->nsons; i++)
      {
        int dir;
        for (dir=0; dir<18; dir++)
        {
          MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh , &q->sons[i], dir, Neighbour);
          if(Neighbour->leaf==0 && Neighbour->depth!=0)
          {
            MMG5_DEL_MEM (mesh,Neighbour);
            return 1;
          }
        }
      }
      MMG3D_merge_MOctree_s (q, mesh);
      MMG5_DEL_MEM (mesh,Neighbour);
      return 1;
    }
    else
    {
      MMG5_DEL_MEM (mesh,Neighbour);
    }
  }
  return 0;
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
  MMG5_MOctree_s *po;
  po=mesh->octree->root;

  int i, depth_max;
  int max_dim=0;
  for ( i=0; i<3; ++i ) {
    if(max_dim < mesh->freeint[i])
    {
      max_dim = mesh->freeint[i];
    }
  }

  max_dim--;

  /* set max dim to the next power of 2 */
  max_dim--;
  max_dim |= max_dim >> 1;
  max_dim |= max_dim >> 2;
  max_dim |= max_dim >> 4;
  max_dim |= max_dim >> 8;
  max_dim |= max_dim >> 16;
  max_dim++;

  depth_max=log(max_dim)/log(2);
  MMG3D_build_coarsen_octree_first_time(mesh, po, depth_max);
  int depth;
  for (depth=depth_max-2; depth>=0; depth --)
  {
    MMG3D_build_coarsen_octree(mesh, po, depth_max, depth);
  }

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Delete the bounding box created for the delaunay mesh generation
 *
 */
static inline
int MMG3D_delete_bounding_box ( MMG5_pMesh mesh ) {
  MMG5_pTetra pt;
  int ip[8];
  int k,i,j,iel,ifac,base;

  /** Step 1: if a tetra has at least 1 vertex of the bounding box, it must be
   * removed */
  for ( i=0; i<8; ++i ) {
    ip[i] = (mesh->np-i);
  }

  base = ++mesh->base;
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /* Mark the tetra with 1 bb vertex */
    for ( i=0; i<8; ++i ) {
      for ( j=0; j<4; ++j ) {
        if ( pt->v[j] == ip[i] ) {
          pt->flag = -base;
          break;
        }
      }
      if ( j != 4 ) break;
    }
  }

  /** Step 2: delete the BB tetra and update the adja array */
  for ( k=1; k<=mesh->ne; ++k ) {
    pt = &mesh->tetra[k];

    if ( !MG_EOK(pt) ) continue;

    if ( pt->flag != -base ) continue;

    for ( i=0; i<4; ++i ) {
      iel  = mesh->adja[ 4*(k-1) + 1 + i ];

      if ( !iel ) continue;

      ifac = iel % 4;
      iel /= 4;

      mesh->adja[ 4*(iel-1) + 1 + ifac ] = 0;
      mesh->adja[ 4*(k-1) + 1 + i ] = 0;

    }
    MMG3D_delElt ( mesh, k );
  }


  /** Step 3: delete the BB points */
  for ( i=0; i<8; ++i ) {
    MMG3D_delPt ( mesh,ip[i] );
  }

  return 1;
}


/**
 * \param mesh pointer toward a mesh structure.
 * \param q pointer toward to the MOctree cell
 * \param depth_max maximum depth of the octree
 * \param depth depth of the octree to reach
 *
 *
 * \return 1 if success, 0 if fail.
 *
 * Delete the octree cells at depth
 *
 */
static inline
int MMG3D_delete_MOctree(MMG5_pMesh mesh, MMG5_MOctree_s* q) {
  int i;

#ifndef NDEBUG
  /* Security check */
  if ( q->nsons ) {
    assert ( q->sons );
  }
  if ( q->sons ) {
    assert ( q->nsons );
  }
#endif

  /* Free the subcell */
  for ( i=0; i<q->nsons; i++ ) {
    if ( q->sons[i].nsons )  {
      MMG3D_delete_MOctree(mesh, &q->sons[i]);
    }
  }

  /* Free our sons */
  MMG3D_free_MOctree_s ( q, mesh );

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Delete the octree and free memory
 *
 */
static inline
int MMG3D_delete_octree ( MMG5_pMesh mesh ) {

  MMG3D_delete_MOctree(mesh,mesh->octree->root);

  MMG3D_free_MOctree( &mesh->octree,mesh);

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
  int i, depth_max;
  int max_dim=0;
  int ip_bb_pt_list[8],ip_bb_elt_list[5];
  int ier;

  /* Mark all the points as unused */
  for ( i=1; i<=mesh->np; ++i ) {
    mesh->point[i].tag = MG_NUL;
  }

  /* Get the maximal dimension */
  for ( i=0; i<3; ++i ) {
    if(max_dim < mesh->freeint[i])
    {
      max_dim = mesh->freeint[i];
    }
  }

  max_dim--;

  /* set max dim to the next power of 2 */
  max_dim--;
  max_dim |= max_dim >> 1;
  max_dim |= max_dim >> 2;
  max_dim |= max_dim >> 4;
  max_dim |= max_dim >> 8;
  max_dim |= max_dim >> 16;
  max_dim++;

  depth_max=log(max_dim)/log(2);

  if ( mesh->info.imprim > 4 ) {
    printf("  ** BOUNDING BOX GENERATION \n");
  }

  ier = MMG3D_build_bounding_box (mesh,sol, ip_bb_pt_list, ip_bb_elt_list);
  if ( !ier ) {
    fprintf (stderr,"\n  ## Warning: %s: unable to create the mesh bounding box.\n",__func__);
    return 0;
  }

  mesh->ntmax = MMG3D_NTMAX;

  if ( !MMG3D_analys(mesh) ) return 0;

  if ( mesh->info.imprim > 4 ) {
    printf("  ** MESH GENERATION \n");
  }

  MMG3D_add_Boundary (mesh, sol, depth_max);

  /* Delete the mesh bounding box */
  if ( !MMG3D_delete_bounding_box (mesh) ) {
    fprintf (stderr,"\n  ## Warning: %s: unable to delete the mesh bounding box.\n",__func__);
    return 0;
  }

  /* Reset mesh informations */
  MMG3D_delete_octree (mesh);

  mesh->mark = 0;
  MMG5_freeXTets(mesh);

  sol->np = mesh->np;

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Convert a balanced octree (2:1) into a tetrahedral mesh using tetgen.
 *
 */
static inline
int MMG3D_convert_octree2tetmesh_with_tetgen(MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_MOctree_s   *q;
  MMG5_pPoint      ppt;
  int              i,ier,span,np,nc;
  FILE             *inm;
  char             *filename="tmp_tetgen";
  char             tetgenfile[256],command[256];

  /* Mark all the points as unused */
  for ( i=1; i<=mesh->np; ++i ) {
    mesh->point[i].tag = MG_NUL;
  }

  /* Mark the octree points as used and pack the mesh vertices */
  q = mesh->octree->root;

  span = mesh->octree->nspan_at_root;
  np = nc = 0;
  ier = MMG3D_mark_MOctreeCellCorners(mesh,q,span,&np,&nc);
  if ( !ier ) {
    fprintf(stderr,"\n  ## Error: %s: unable to mark the octree cell corners as"
            " used.\n",__func__);
    return 0;
  }
  MMG3D_delete_octree (mesh);

  /* compact metric */
  if ( sol && sol->m ) {
    if ( !MMG3D_pack_sol(mesh,sol) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to pack the solution.\n",__func__);
      return 0;
    }
  }

  nc = MMG3D_pack_points(mesh);
  if ( nc<0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to pack the mesh vertices.\n",__func__);
    return 0;
  }

  /* MMG3D_saveMesh(mesh,"tt.mesh"); */
  /* MMG3D_saveSol(mesh,sol,"tt.sol"); */

  /* Save the points in a .node file */
  sprintf(tetgenfile, "%s%s",filename,".node");
  if( !(inm = fopen(tetgenfile,"w")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",tetgenfile);
    return 0;
  }

  /* .node file : np points, in dim mesh->dim, without attribute and without markers */
  fprintf(inm,"%d %d %d %d\n",np,mesh->dim, 0, 0);
  np = 0;
  for ( i=1; i<=mesh->np; ++i ) {
    ppt = &mesh->point[i];
    if ( MG_VOK(ppt) ) {
      ppt->tmp = np++;
      fprintf(inm," %d %.15lg %.15lg %.15lg\n",ppt->tmp, ppt->c[0],ppt->c[1],ppt->c[2]);
    }
  }
  fclose(inm);

  /* run tetgen on the .node file */
  printf("  --> Tetgen executable: %s\n",TETGEN);
  sprintf(command, "%s -BANEF -Q -g %s", TETGEN, tetgenfile);
  ier = system(command);
  if ( ier != 0 ) {
    printf("  ## Error:%s: Tetgen error.\n",__func__);
    return 0;
  }
  MMG3D_saveSol(mesh,sol,"tmp_tetgen.1.sol");

  /* Reset mesh informations */
  mesh->mark = 0;

  if ( mesh->tetra ) {
    MMG5_DEL_MEM(mesh,mesh->tetra);
    mesh->ne = 0;
  }
  if ( mesh->point ) {
    MMG5_DEL_MEM(mesh,mesh->point);
    mesh->np = 0;
  }
  assert ( !mesh->xtetra );
  assert ( !mesh->xpoint );
  assert ( !mesh->tria );

  /* Read the tetgen mesh */
  sprintf(tetgenfile, "%s%s",filename,".1.mesh");
  ier = MMG3D_loadMesh ( mesh,tetgenfile );
  if ( ier<= 0 ) return 0;


  /* Clean the tetgen mesh */
  /* Remove spurious triangles */
  if ( mesh->tria ) {
    MMG5_DEL_MEM(mesh,mesh->tria);
    mesh->nt = 0;
  }
  /* Remove spurious tags on mesh vertices */
  for ( i=1; i<=mesh->np; ++i ) {
    ppt = &mesh->point[i];
    if ( !MG_VOK(ppt) ) continue;

    ppt->tag = MG_NOTAG;
  }

  return 1;
}

static inline int MMG3D_check_octreeSons ( MMG5_MOctree_s *q ) {
  int i;

  /* Security check */
  if ( q->nsons ) {
    assert ( q->sons );
  }
  if ( q->sons ) {
    assert ( q->nsons );
  }

  for ( i=0; i<q->nsons; ++i ) {
    MMG3D_check_octreeSons (  &q->sons[i] );
  }

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure that contains the solution at
 * grid centroids at the beginning and at mesh nodes at the end.
 *
 * \return 1 if success, 0 if fail.
 *
 * Convert a structured grid into an octree, then transform this octree into a
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

#ifndef NDEBUG
  MMG3D_check_octreeSons ( mesh->octree->root );
#endif

  /* Creation of the coarse octree if an isovalue is provided */
  if ( sol && sol->m && !MMG3D_coarsen_octree(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree coarsening problem. Exit program.\n");
    return 0;
  }

#ifndef NDEBUG
  MMG3D_check_octreeSons ( mesh->octree->root );
  if ( !MMG3D_saveVTKOctree(mesh,sol,mesh->nameout) ) {
    fprintf(stderr,"\n  ## Warning: unable to save the coarsen octree\n");
  }
#endif

  /**--- stage 2: Tetrahedralization */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** OCTREE TETRAHEDRALIZATION\n");

  if ( !strcmp(TETGEN,"TETGEN_EXEC-NOTFOUND" ) ) {
    if ( !MMG3D_convert_octree2tetmesh(mesh,sol) ) {
      fprintf(stderr,"\n  ## Octree tetrahedralization problem. Exit program.\n");
      return 0;
    }
  }
  else {
    if ( !MMG3D_convert_octree2tetmesh_with_tetgen(mesh,sol) ) {
      fprintf(stderr,"\n  ## Octree tetrahedralization using tetgen problem. Exit program.\n");
      return 0;
    }
  }
  return 1;
}
