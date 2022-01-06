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
 * \file mmg/mmg2.c
 * \author Algiane Froehly (Bx INP/Inria/UBordeaux)
 * \author Luca Cirrottola (Inria)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * Common functions for ls discretization.
 *
 */

#include "mmgcommon.h"

/**
 * \param pim    multimaterials inverse data table.
 * \param ref    material reference.
 * \return the key of the material in the lookup/hash table.
 *
 * Compute key for the material in the hash table.
 */
static int MMG5_InvMat_key(MMG5_pInvMat pim,int ref) {
  return (ref - pim->offset);
}

/**
 * \param k      index of the material.
 * \param attr   the attribute of the material (nosplit/split/plus/minus).
 * \return the code to be stored in the inverse data table entry.
 *
 * Compute the material entry in the inverse data table, storing the index of
 * the parent material and the attribute of the child material.
 */
static int MMG5_InvMat_code(int k,int attr) {
  return 4*(k+1)+attr;
}

/**
 * \param pim multimaterials inverse data table.
 * \param ref    material reference.
 * \return Index of the material.
 *
 * Get index of the parent material from lookup table.
 */
static int MMG5_InvMat_getIndex(MMG5_pInvMat pim,int ref) {
  int key = MMG5_InvMat_key(pim,ref);
  /* The parent index is stored as 4*(k+1) */
  return (pim->lookup[key] / 4 - 1);
}

/**
 * \param mesh   pointer toward the mesh structure.
 * \param pim    multimaterials inverse data table.
 * \param ref    material reference.
 * \return the nosplit/split/plus/minus attribute of the material.
 *
 * Get attribute of the child material (nosplit/split/plus/minus) from lookup
 * table.
 */
static int MMG5_InvMat_getAttrib(MMG5_pInvMat pim,int ref) {
  int key = MMG5_InvMat_key(pim,ref);
  /* The nosplit/split/plus/minus attribute is stored as the rest of the
   * integer division. */
  return (pim->lookup[key] % 4);
}

/**
 * \param pim    multimaterials inverse data table.
 * \param key    material key in the inverse data table.
 * \return 0 if an entry for the material already exists, 1 otherwise.
 *
 * Check if a material reference already exists in the material lookup table.
 */
static int MMG5_InvMat_check(MMG5_pInvMat pim,int key) {
  return pim->lookup[key] ? 0 : 1;
}

/**
 * \param pim    multimaterials inverse data table.
 * \param ref    reference of the material..
 * \param k      index of the material in the direct table.
 * \return 0 if an entry for the material already exists, 1 otherwise.
 *
 * Raise an error if trying to overwrite a reference entry in the material
 * lookup table.
 */
static void MMG5_InvMat_error(MMG5_pInvMat pim,int ref,int k) {
  fprintf(stderr,"\n   ## Error: Trying to overwrite material reference %d"
    " (from LSReferences line %d) with another entry from LSReferences line %d."
    ,ref,MMG5_InvMat_getIndex(pim,ref)+1,k+1);
  fprintf(stderr,"\n             Check your LSReferences table: each material"
    " reference should be unique!\n");
}

/**
 * \param mesh   pointer toward the mesh structure.
 * \param pim    multimaterials inverse data table.
 * \param k      index of the material in the input table.
 *
 * Set materials lookup table entry.
 */
static int MMG5_InvMat_set(MMG5_pMesh mesh,MMG5_pInvMat pim,int k) {
  MMG5_pMat pm;
  int       key;

  /* Get material */
  pm = &mesh->info.mat[k];

  /** Store the dosplit attribute of the parent material */
  key = MMG5_InvMat_key(pim,pm->ref);
  if( MMG5_InvMat_check(pim,key) ) {
    pim->lookup[key] = MMG5_InvMat_code(k,pm->dospl);
  } else {
    MMG5_InvMat_error(pim,pm->ref,k);
    return 0;
  }

  /** Store the child material sign with the parent material index (in the
   *  lookup table).
   *  1) 0 is a legit material index, so store the parent as 4*(k+1).
   *  2) If a child material has the same reference as the parent, this
   *     effectively overwrites the result of the previous instruction.
   *  3) No different child materials are allowed to have the same reference,
   *     and this must have already been checked. */
  if( pm->dospl ) {
    key = MMG5_InvMat_key(pim,pm->rin);
    if( MMG5_InvMat_check(pim,key) ) {
      pim->lookup[key] = MMG5_InvMat_code(k,MG_MINUS);
    } else {
      MMG5_InvMat_error(pim,pm->rin,k);
      return 0;
    }
    key = MMG5_InvMat_key(pim,pm->rex);
    if( MMG5_InvMat_check(pim,key) ) {
      pim->lookup[key] = MMG5_InvMat_code(k,MG_PLUS);
    } else {
      MMG5_InvMat_error(pim,pm->rex,k);
      return 0;
    }
  }

  return 1;
}

/**
 * \param mesh   pointer toward the mesh structure.
 * \param pim    multimaterials inverse data table.
 * \param ref    material reference.
 * \param pref   pointer to the parent material reference.
 * \return 1 if found, 0 if not found.
 *
 * Get reference of the parent material in multimaterial mode.
 *.
 */
static int MMG5_InvMat_getParent(MMG5_pMesh mesh,MMG5_pInvMat pim,int ref,int *pref) {
  MMG5_pMat pm;
  int       k;

  /* The parent index */
  k = MMG5_InvMat_getIndex(pim,ref);

  /* Material not found in the table */
  if( k == -1 ) return 0;

  /* Get the material in the lookup table and return the parent reference */
  pm = &mesh->info.mat[k];
  *pref = pm->ref;
  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param ref  final reference for which we are searching the initial one
 * \param pref pointer to the reference of the parent material.
 * \return 1 if found, 0 otherwise.
 *
 * Retrieve the starting domain reference (parent material) associated to the
 * split reference ref. Allow the call in non-multimaterial mode.
 *
 */
int MMG5_getStartRef(MMG5_pMesh mesh,int ref,int *pref) {
  MMG5_pInvMat pim;

  /* No multi-materials nor single material reference preservation */
  if( !mesh->info.nmat ) {
    *pref = 0;
    return 1;
  }

  /* Get parent of material */
  pim = &mesh->info.invmat;

  /* Return 0 if the material does not exist, 1 otherwise */
  if( !MMG5_InvMat_getParent(mesh,pim,ref,pref) )
    return 0;
  else
    return 1;
}

/**
 * \param mesh   pointer toward the mesh structure.
 * \param pim    multimaterials inverse data table.
 *
 * Print materials lookup table.
 */
static void MMG5_InvMat_print(MMG5_pMesh mesh,MMG5_pInvMat pim) {
  int ref,pref;

  /* Scan all references in the table limits, some may not exist */
  for( ref = pim->offset; ref < pim->offset + pim->size; ref++ ) {
    if( !MMG5_InvMat_getParent(mesh,pim,ref,&pref) ) continue;
    printf("%d (%d): %d %d\n",ref,MMG5_InvMat_key(pim,ref),pref,
        MMG5_InvMat_getAttrib(pim,ref));
  }
}

/**
 * \param mesh   pointer toward the mesh structure.
 * \return 1 if success, 0 if fail.
 *
 *
 * Initialize handling of multimaterial mode.
 *
 * An indexed table of materials has been provided by the MMG5_Mat datatype in
 * the form:
 *
 *   index | dospl       | ref       | rin       | rex
 *   -------------------------------------------------------
 *   0     | dospl_0     | ref_0     | rin_0     | rex_0
 *   ...   | ...         | ...       | ...       |
 *   k     | dospl_k     | ref_k     | rin_k     | rex_k
 *   ...   | ...         | ...       | ...       |
 *   n-1   | dospl_{n-1} | ref_{n-1} | rin_{n-1} | rex_{n-1}
 *
 * where dospl is the split/preserve attribute of the material, and rin,rex
 * are its child materials (if dospl). Viceversa, ref is the parent material
 * for rin,rex.
 *
 * Here a lookup table for material references is built through trivial hashing
 * for all references (both parent and child materials) with the key:
 *
 *   key = ref - ref_min,    ref_min = min_{k = 0,n-1} (ref_k, rin_k, rex_k)
 *
 * For all references, it is important to store
 * - the index k of the row in the original table,
 * - the characteristic attribute (parent-split, parent-preserve, child-interior,
 *   child-exterior) of the material.
 *
 * Since
 *   dospl = 0 or 1,  MG_PLUS = 2, and MG_MINUS = 3
 *
 * we can store the attribute as dospl (for the parent material) or MG_MINUS/
 * MG_PLUS (for the child materials), with its value ranging between 0 and 3.
 *
 * A convenient entry to store both the index and the attribute in the lookup
 * table is thus:
 *
 *   entry = 4*(index+1) + attribute
 *
 * leading to a lookup table in the form (the key ordering is only an example):
 *
 *   key   | entry
 *   ------------------------
 *   ...   | ...
 *   ref_k | 4*(k+1)+dospl_k
 *   ...   | ...
 *   rin_k | 4*(k+1)+MG_MINUS
 *   ...   | ...
 *   rex_k | 4*(k+1)+MG_PLUS
 *   ...   | ...
 *
 * where the index and the attribute of the material can be retrieved as
 *
 *   index     = entry / 4 -1
 *   attribute = entry % 4
 *
 *
 * What if two materials have the same reference?
 *  - child references should be distinct (the entry in the lookup table would
 *    be overwritten),
 *  - a parent material can have itself as child (a positive attribute would say
 *    it should be split, the attribute value would say if it is interior or
 *    exterior).
 * Thus, each of the maps parent->child_in and parent->child_ext is injective,
 * but a fixed point is allowed.
 *
 * Why child materials should be different?
 *  - because on failure it is necessary to recover the parent references and
 *    apply them on the tetrahedra to restore the starting materials
 *    distribution.
 *
 */
int MMG5_MultiMat_init(MMG5_pMesh mesh) {
  MMG5_pMat    pm;
  MMG5_pInvMat pim;
  int          k;
  int          refmax,refmin;

  /* Nothing to do if no multi-material option */
  if( !mesh->info.nmat ) return 1;

  /* Error if all the materials have not been set */
  if( mesh->info.nmati < mesh->info.nmat ) {
    fprintf(stderr,"\n ## Error: %s: Only %d materials out of %d have been set.\n",
        __func__,mesh->info.nmati,mesh->info.nmat);
    return 0;
  }

  /* Get pointer to the structure */
  pim = &mesh->info.invmat;

  /* Initialize the max and min reference */
  refmax = 0;
  refmin = INT_MAX;

  /* Look for the max/min reference */
  for( k = 0; k < mesh->info.nmat; k++ ) {
    pm = &mesh->info.mat[k];
    /* Update max and min val for original ref */
    if( pm->ref > refmax ) refmax = pm->ref;
    if( pm->ref < refmin ) refmin = pm->ref;
    if( !pm->dospl ) continue;
    /* Update max and min val with interior ref */
    if( pm->rin > refmax ) refmax = pm->rin;
    if( pm->rin < refmin ) refmin = pm->rin;
    /* Update max and min val with exterior ref */
    if( pm->rex > refmax ) refmax = pm->rex;
    if( pm->rex < refmin ) refmin = pm->rex;
  }

  /* Get span of the lookup table */
  pim->offset = refmin;
  pim->size   = refmax - refmin + 1;
  assert( pim->size > 0 );

  /* Allocate lookup table */
  MMG5_ADD_MEM(mesh,pim->size*sizeof(int),"materials lookup table",return 0);
  MMG5_SAFE_CALLOC(pim->lookup,pim->size,int,return 0);

  /* Fill lookup table */
  for( k = 0; k < mesh->info.nmat; k++ ) {
    if( !MMG5_InvMat_set(mesh,pim,k) )
      return 0;
  }

//  MMG5_InvMat_print(mesh,pim);
  return 1;
}

/**
 * \param mesh   pointer toward the mesh structure.
 * \param ref    initial reference.
 * \param refint internal reference after ls discretization.
 * \param refint internal reference after ls discretization.
 * \return 1 if entity can be splitted, 0 if cannot be splitted.
 *
 * Identify whether an entity with reference ref should be split, and the
 * labels of the resulting entities.
 *
 */
int MMG5_isSplit(MMG5_pMesh mesh,int ref,int *refint,int *refext) {
  MMG5_pInvMat pim;
  MMG5_pMat    pm;
  int          k;

  /* Default case: split with references MG_MINUS, MG_PLUS */
  if( !mesh->info.nmat ) {
    *refint = MG_MINUS;
    *refext = MG_PLUS;
    return 1;
  }

  /* Check in the info->mat table if reference ref is supplied by the user */
  pim = &mesh->info.invmat;
  k = MMG5_InvMat_getIndex(pim,ref);

  assert( k != -1 );
  pm = &mesh->info.mat[k];

  if ( !pm->dospl ) {
    return 0;
  } else {
    *refint = pm->rin;
    *refext = pm->rex;
    return 1;
  }
}

/**
 * \param mesh   pointer toward the mesh structure.
 * \param ref    initial reference.
 * \return 1 if entity cannot be split, 0 if can be split.
 *
 * Identify whether an entity with reference ref should not be split.
 *
 */
int MMG5_isNotSplit(MMG5_pMesh mesh,int ref) {
  MMG5_pInvMat pim;

  /* Split material by default if not in multi-material mode */
  if( !mesh->info.nmat ) return 0;

  /* Look in the table otherwise */
  pim = &mesh->info.invmat;
  if( !MMG5_InvMat_getAttrib(pim,ref) )
    return 0;
  else
    return 1;

}

/**
 * \param mesh   pointer toward the mesh structure.
 * \param ref0   reference of the first tetrahedron sharing the face.
 * \param ref1   reference of the second tetrahedron sharing the face..
 * \return 1 if face is on the discrete level set, 0 if not.
 *
 * Identify whether a face is on the discrete level set or not.
 *
 */
int MMG5_isLevelSet(MMG5_pMesh mesh,int ref0,int ref1) {
  MMG5_pInvMat pim;
  int8_t       found0,found1;

  /* Check whether multimaterial case or not */
  if( mesh->info.nmat ) {
    /* Retrieve levelset information from the lookup table */
    pim = &mesh->info.invmat;
    found0 = MMG5_InvMat_getAttrib(pim,ref0);
    found1 = MMG5_InvMat_getAttrib(pim,ref1);

    if( (found0+found1) == (MG_MINUS+MG_PLUS) ) return 1;
    else return 0;

  } else {
    /* Single material, check references directly */
    if( ( ref0 == MG_MINUS && ref1 == MG_PLUS ) ||
        ( ref1 == MG_MINUS && ref0 == MG_PLUS ) ) return 1;
    else return 0;
  }
}
