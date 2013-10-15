/**
 * Logiciel initial: MMG3D Version 4.0
 * Co-auteurs : Cecile Dobrzynski et Pascal Frey.
 * Propriétaires :IPB - UPMC -INRIA.
 *
 * Copyright © 2004-2005-2006-2007-2008-2009-2010-2011,
 * diffusé sous les termes et conditions de la licence publique générale de GNU
 * Version 3 ou toute version ultérieure.
 *
 * Ce fichier est une partie de MMG3D.
 * MMG3D est un logiciel libre ; vous pouvez le redistribuer et/ou le modifier
 * suivant les termes de la licence publique générale de GNU
 * Version 3 ou toute version ultérieure.
 * MMG3D est distribué dans l'espoir qu'il sera utile, mais SANS
 * AUCUNE GARANTIE ; sans même garantie de valeur marchande.
 * Voir la licence publique générale de GNU pour plus de détails.
 * MMG3D est diffusé en espérant qu’il sera utile,
 * mais SANS AUCUNE GARANTIE, ni explicite ni implicite,
 * y compris les garanties de commercialisation ou
 * d’adaptation dans un but spécifique.
 * Reportez-vous à la licence publique générale de GNU pour plus de détails.
 * Vous devez avoir reçu une copie de la licence publique générale de GNU
 * en même temps que ce document.
 * Si ce n’est pas le cas, aller voir <http://www.gnu.org/licenses/>.**/
/**
 * Initial software: MMG3D Version 4.0
 * Co-authors: Cecile Dobrzynski et Pascal Frey.
 * Owners: IPB - UPMC -INRIA.
 *
 * Copyright © 2004-2005-2006-2007-2008-2009-2010-2011,
 * spread under the terms and conditions of the license GNU General Public License
 * as published Version 3, or (at your option) any later version.
 *
 * This file is part of MMG3D
 * MMG3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * MMG3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with MMG3D. If not, see <http://www.gnu.org/licenses/>.
 **/
/** librnbg
 *
 * Written by Cedric Lachat
 **/
#include "mmg3d.h"

#ifdef USE_SCOTCH

#include <time.h>
#include "librnbg.h"

#define CHECK_SCOTCH(t,m,e) if(0!=t){perror(m);exit(e);}




/** Internal function : biPartBoxCompute
 * it computes a new numbering of graph vertices, using a bipartitioning.
 *
 *  - graf : the input graph
 *  - vertNbr : the number of vertices
 *  - boxVertNbr : the number of vertices of each box
 *  - permVrtTab : the new numbering
 *
 *  returning 0 if OK, 1 else
 **/
int biPartBoxCompute(SCOTCH_Graph graf, int vertNbr, int boxVertNbr, SCOTCH_Num *permVrtTab) {
  int boxNbr, vertIdx, boxIdx;
  SCOTCH_Num tmp, tmp2, *partTab, *partNumTab, *partPrmTab;
  SCOTCH_Strat strat ;

  /* Computing the number of boxes */
  boxNbr = vertNbr / boxVertNbr;
  if (boxNbr * boxVertNbr != vertNbr) {
    boxNbr = boxNbr + 1;
  }


  /* Initializing SCOTCH functions */
  CHECK_SCOTCH(SCOTCH_stratInit(&strat), "scotch_stratInit", 0) ;
  CHECK_SCOTCH(SCOTCH_stratGraphMap(&strat, "r{job=t,map=t,poli=S,sep=m{,vert=80,low=h{pass=10}f{bal=0.005,move=0},asc=b{bnd=f{bal=0.05,move=0},org=f{bal=0.05,move=0}}}|m{,vert=80,low=h{pass=10}f{bal=0.005,move=0},asc=b{bnd=f{bal=0.05,move=0},org=f{bal=0.05,move=0}}}}"), "scotch_stratGraphMap", 0) ;

  partTab = (SCOTCH_Num *)calloc(vertNbr, sizeof(SCOTCH_Num));


  /* Partionning the graph */
  CHECK_SCOTCH(SCOTCH_graphPart(&graf, boxNbr, &strat, partTab), "scotch_graphPart", 0);

  partNumTab = (SCOTCH_Num *)calloc(boxNbr, sizeof(SCOTCH_Num));

  if (!memset(partNumTab, 0, boxNbr*sizeof(SCOTCH_Num))) {
    perror("memset");
    return 0;
  }

  /* Computing the number of elements of each box */
  for( vertIdx = 0 ; vertIdx< vertNbr ;vertIdx++)
    partNumTab[partTab[vertIdx]] += 1;


  /* partition permutation tabular */
  partPrmTab = (SCOTCH_Num *)calloc(vertNbr + 1, sizeof(SCOTCH_Num));


  /* Copying the previous tabular in order to have the index of the first
   * element of each box
   * */
  tmp = partNumTab[0];
  partNumTab[0] = 0;
  for(boxIdx = 1; boxIdx < boxNbr ; boxIdx++) {
    tmp2 = partNumTab[boxIdx];
    partNumTab[boxIdx] = partNumTab[boxIdx-1] + tmp;
    tmp = tmp2;
  }

  /* partPrmTab is built such as each vertex belongs to his box */
  for( vertIdx = 0;vertIdx< vertNbr;vertIdx++)
    partPrmTab[partNumTab[partTab[vertIdx]]++] = vertIdx;


  /* Infering the new numbering */
  for (vertIdx = 0; vertIdx < vertNbr ; vertIdx++)
    permVrtTab[partPrmTab[vertIdx] + 1] = vertIdx + 1;

  free(partTab);
  free(partNumTab);
  free(partPrmTab);

  SCOTCH_stratExit(&strat) ;
  return 0;
}




/* Internal function : kPartBoxCompute
 * it computes a new numbering of graph vertices, using a k-partitioning.
 * Assuming that baseval of the graph is 1
 *
 *  - graf : the input graph
 *  - vertNbr : the number of vertices
 *  - boxVertNbr : the number of vertices of each box
 *  - permVrtTab : the new numbering
 *
 *  returning 0 if OK, 1 else
 */
int kPartBoxCompute(SCOTCH_Graph graf, int vertNbr, int boxVertNbr, SCOTCH_Num *permVrtTab) {
  int boxNbr, vertIdx;
  SCOTCH_Num logMaxVal, SupMaxVal, InfMaxVal, maxVal;
  char s[200];
  SCOTCH_Num *sortPartTb;
  SCOTCH_Strat strat ;
  SCOTCH_Arch arch;

  /* Computing the number of boxes */
  boxNbr = vertNbr / boxVertNbr;
  if (boxNbr * boxVertNbr != vertNbr) {
    boxNbr = boxNbr + 1;
  }


  /* Initializing SCOTCH functions */
  CHECK_SCOTCH(SCOTCH_stratInit(&strat), "scotch_stratInit", 0) ;
  CHECK_SCOTCH(SCOTCH_archCmplt(&arch, boxNbr), "scotch_archVcmplt", 0) ;

  sprintf(s, "m{vert=%d,low=r{job=t,map=t,poli=S,sep=m{vert=80,low=h{pass=10}f{bal=0.0005,move=80},asc=f{bal=0.005,move=80}}}}", vertNbr / boxVertNbr);
  CHECK_SCOTCH(SCOTCH_stratGraphMap(&strat, s), "scotch_stratGraphMap", 0) ;


  sortPartTb= (SCOTCH_Num *)calloc(2*vertNbr, sizeof(SCOTCH_Num));


  /* Partionning the graph */
  CHECK_SCOTCH(SCOTCH_graphMap(&graf, &arch, &strat, sortPartTb), "scotch_graphMap", 0);


  for (vertIdx = vertNbr - 1 ; vertIdx >= 0 ; vertIdx--) {
    sortPartTb[2*vertIdx] = sortPartTb[vertIdx];
    sortPartTb[2*vertIdx+1] = vertIdx + 1;
  }

  // Sorting the tabular, which contains box values and vertex numbers
  _SCOTCHintSort2asc1(sortPartTb, vertNbr);


  /* Infering the new numbering */
  for (vertIdx = 0; vertIdx < vertNbr ; vertIdx++) {
    permVrtTab[sortPartTb[2*vertIdx + 1]] = vertIdx + 1;
  }

  SCOTCH_stratExit(&strat) ;
  SCOTCH_archExit(&arch) ;

  free(sortPartTb);

  return 0;
}



/* Function : renumbering
 *  it modifies the numbering of each node to prevent from cache missing.
 *
 *  - boxVertNbr : number of vertices by box
 *  - mesh : the input mesh which is modified
 *
 *  returning 0 if OK, 1 else
 */
int renumbering(int boxVertNbr, pMesh mesh, pSol sol) {
  pPoint ppt;
  pPoint points;
  pTetra ptet, tetras;
  SCOTCH_Num edgeNbr;
  SCOTCH_Num *vertTab, *vendTab, *edgeTab, *permVrtTab;
  SCOTCH_Graph graf ;
  int    vertNbr, nodeGlbIdx, tetraIdx, ballTetIdx;
  int    i, j, k, addrNew, addrOld;
  int    edgeSiz;
  int    *vertOldTab, *permNodTab, nereal, npreal;
  int    *adja,iadr;
  double *metNew;


  /* Computing the number of vertices and a contiguous tabular of vertices */
  vertNbr = 0;
  vertOldTab = (int *)calloc(mesh->ne + 1, sizeof(int));

  if (!memset(vertOldTab, 0, sizeof(int)*(mesh->ne+1))) {
    perror("memset");
    return 1;
  }

  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {

    /* Testing if the tetra exists */
    if (!mesh->tetra[tetraIdx].v[0]) continue;
    vertOldTab[tetraIdx] = ++vertNbr;
  }

  if ( vertNbr/2 < BOXSIZE ) {
    /* not enough tetra to renum */
    free(vertOldTab);
    return(1);
  }


  /* Allocating memory to compute adjacency lists */
  vertTab = (SCOTCH_Num *)calloc(vertNbr + 1, sizeof(SCOTCH_Num));

  if (!memset(vertTab, ~0, sizeof(SCOTCH_Num)*(vertNbr + 1))) {
    perror("memset");
    return 1;
  }

  vendTab = (SCOTCH_Num *)calloc(vertNbr + 1, sizeof(SCOTCH_Num));

  edgeNbr = 1;
  edgeSiz = vertNbr*2;
  edgeTab = (SCOTCH_Num *)calloc(edgeSiz, sizeof(SCOTCH_Num));



  /* Computing the adjacency list for each vertex */
  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {

    /* Testing if the tetra exists */
    if (!mesh->tetra[tetraIdx].v[0]) continue;

    iadr = 4*(tetraIdx-1) + 1;
    adja = &mesh->adja[iadr];
    for (i=0; i<4; i++) {
      ballTetIdx = adja[i] / 4;

      if (!ballTetIdx) continue;

      /* Testing if one neighbour of tetraIdx has already been added */
      if (vertTab[vertOldTab[tetraIdx]] < 0)
        vertTab[vertOldTab[tetraIdx]] = edgeNbr;
      vendTab[vertOldTab[tetraIdx]] = edgeNbr+1;

      /* Testing if edgeTab memory is enough */
      if (edgeNbr >= edgeSiz) {
        edgeSiz += EDGEGAP;
        edgeTab = (SCOTCH_Num *)realloc(edgeTab, edgeSiz * sizeof(SCOTCH_Num));
      }

      edgeTab[edgeNbr++] = vertOldTab[ballTetIdx];
    }
  }

  edgeNbr--;


  /* Building the graph by calling Scotch functions */

  SCOTCH_graphInit(&graf) ;
  CHECK_SCOTCH(SCOTCH_graphBuild(&graf, (SCOTCH_Num) 1, vertNbr, vertTab+1,
                                 vendTab+1, NULL, NULL, edgeNbr, edgeTab+1, NULL),
               "scotch_graphbuild", 0) ;
  CHECK_SCOTCH(SCOTCH_graphCheck(&graf), "scotch_graphcheck", 0) ;

  permVrtTab = (SCOTCH_Num *)calloc(vertNbr + 1, sizeof(SCOTCH_Num));

  CHECK_SCOTCH(kPartBoxCompute(graf, vertNbr, boxVertNbr, permVrtTab), "boxCompute", 0);

  SCOTCH_graphExit(&graf) ;

  free(vertTab);
  free(vendTab);
  free(edgeTab);


  permNodTab = (int *)calloc(mesh->np + 1, sizeof(int));

  /* Computing the new point list and modifying the sol structures*/
  tetras = (pTetra)calloc(mesh->nemax+1,sizeof(Tetra));

  points = (pPoint)calloc(mesh->npmax+1,sizeof(Point));

  metNew = (double*)calloc(sol->npmax+1,sol->size*sizeof(double));

  adja   = (int*)calloc(4*mesh->nemax+5,sizeof(int));
  memset(adja, 0, (4*mesh->nemax+5)*sizeof(int) );

  nereal = 0;
  npreal = 0;

  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {
    ptet = &mesh->tetra[tetraIdx];

    /* Testing if the tetra exists */
    if (!ptet->v[0]) continue;

    /* Building the new point list */
    tetras[permVrtTab[vertOldTab[tetraIdx]]] = *ptet;

    adja[4*(permVrtTab[vertOldTab[tetraIdx]]-1)+1] =
      4*(permVrtTab[vertOldTab[mesh->adja[4*(tetraIdx-1)+1]/4]])+mesh->adja[4*(tetraIdx-1)+1]%4;
    adja[4*(permVrtTab[vertOldTab[tetraIdx]]-1)+2] =
      4*(permVrtTab[vertOldTab[mesh->adja[4*(tetraIdx-1)+2]/4]])+mesh->adja[4*(tetraIdx-1)+2]%4;
    adja[4*(permVrtTab[vertOldTab[tetraIdx]]-1)+3] =
      4*(permVrtTab[vertOldTab[mesh->adja[4*(tetraIdx-1)+3]/4]])+mesh->adja[4*(tetraIdx-1)+3]%4;
    adja[4*(permVrtTab[vertOldTab[tetraIdx]]-1)+4] =
      4*(permVrtTab[vertOldTab[mesh->adja[4*(tetraIdx-1)+4]/4]])+mesh->adja[4*(tetraIdx-1)+4]%4;

    nereal++;

    for(j = 0 ; j <= 3 ; j++) {

      nodeGlbIdx = mesh->tetra[tetraIdx].v[j];

      if ( permNodTab[nodeGlbIdx] ) continue;

      ppt = &mesh->point[nodeGlbIdx];

      if (!(ppt->tag & MG_NUL)) {
        /* Building the new point list */
        permNodTab[nodeGlbIdx] = ++npreal;

        points[permNodTab[nodeGlbIdx]] = *ppt;

        /* Building the new sol met */
        addrOld = (nodeGlbIdx-1)*sol->size + 1;
        addrNew = (permNodTab[nodeGlbIdx]-1)*sol->size + 1;
        memcpy(&metNew[addrNew], &sol->m[addrOld], sol->size*sizeof(double));
      }
    }
  }

  free(mesh->tetra);
  mesh->tetra = tetras;
  mesh->ne = nereal;

  free(mesh->point);
  mesh->point = points;
  mesh->np    = npreal;

  free(sol->m);
  sol->m = metNew;

  free(mesh->adja);
  mesh->adja = adja;

  if ( mesh->np == mesh->npmax )
    mesh->npnil = 0;
  else
    mesh->npnil = mesh->np + 1;

  if ( mesh->ne == mesh->nemax )
    mesh->nenil = 0;
  else
    mesh->nenil = mesh->ne + 1;

  if ( mesh->npnil )
    for (k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;

  if ( mesh->nenil )
    for (k=mesh->nenil; k<mesh->nemax-1; k++)
      mesh->tetra[k].v[3] = k+1;

  /* Modifying the numbering of the nodes of each tetra */
  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {
    if (!mesh->tetra[tetraIdx].v[0]) continue;
    for(j = 0 ; j <= 3 ; j++) {
      mesh->tetra[tetraIdx].v[j] = permNodTab[mesh->tetra[tetraIdx].v[j]];
    }
  }

  free(permVrtTab);
  free(permNodTab);
  free(vertOldTab);

  return 1;
}
#endif
