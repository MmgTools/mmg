/**
 * Logiciel initial: MMG3D Version 4.0
 * Co-auteurs : Cecile Dobrzynski et Pascal Frey.
 * Proprietaires :IPB - UPMC -INRIA.
 *
 * Copyright \copyright 2004-2005-2006-2007-2008-2009-2010-2011,
 * diffuse sous les termes et conditions de la licence publique generale de GNU
 * Version 3 ou toute version ulterieure.
 *
 * Ce fichier est une partie de MMG3D.
 * MMG3D est un logiciel libre ; vous pouvez le redistribuer et/ou le modifier
 * suivant les termes de la licence publique generale de GNU
 * Version 3 ou toute version ulterieure.
 * MMG3D est distribue dans l espoir qu il sera utile, mais SANS
 * AUCUNE GARANTIE ; sans meme garantie de valeur marchande.
 * Voir la licence publique generale de GNU pour plus de details.
 * MMG3D est diffuse en esperant qu il sera utile,
 * mais SANS AUCUNE GARANTIE, ni explicite ni implicite,
 * y compris les garanties de commercialisation ou
 * d adaptation dans un but specifique.
 * Reportez-vous a la licence publique generale de GNU pour plus de details.
 * Vous devez avoir re\c{c}u une copie de la licence publique generale de GNU
 * en meme temps que ce document.
 * Si ce n est pas le cas, aller voir <http://www.gnu.org/licenses/>.**/
/**
 * Initial software: MMG3D Version 4.0
 * Co-authors: Cecile Dobrzynski et Pascal Frey.
 * Owners: IPB - UPMC -INRIA.
 *
 * Copyright \copyright 2004-2005-2006-2007-2008-2009-2010-2011,
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
 * Written by Cedric Lachat and Algiane Froehly
 **/
#include "mmg3d.h"

#ifdef USE_SCOTCH

#include <time.h>
#include "librnbg.h"

#define CHECK_SCOTCH(t,m,e) if(0!=t){perror(m);exit(e);}


/** Not used */
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
int biPartBoxCompute(SCOTCH_Graph graf, int vertNbr, int boxVertNbr, SCOTCH_Num *permVrtTab,MMG5_pMesh mesh) {
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
#if SCOTCH_VERSION>=6
  CHECK_SCOTCH(SCOTCH_stratGraphMap(&strat, "r{job=t,map=t,poli=S,sep=m{,vert=80,low=h{pass=10}f{bal=0.005,move=0},asc=b{bnd=f{bal=0.05,move=0},org=f{bal=0.05,move=0}}}|m{,vert=80,low=h{pass=10}f{bal=0.005,move=0},asc=b{bnd=f{bal=0.05,move=0},org=f{bal=0.05,move=0}}}}"), "scotch_stratGraphMap", 0) ;
#else
  CHECK_SCOTCH(SCOTCH_stratGraphMap(&strat, "r{job=t,map=t,poli=S,sep=m{type=h,vert=80,low=h{pass=10}f{bal=0.005,move=0},asc=b{bnd=f{bal=0.05,move=0},org=f{bal=0.05,move=0}}}|m{,vert=80,low=h{pass=10}f{bal=0.005,move=0},asc=b{bnd=f{bal=0.05,move=0},org=f{bal=0.05,move=0}}}}"), "scotch_stratGraphMap", 0) ;
#endif

  ADD_MEM(mesh,vertNbr*sizeof(SCOTCH_Num),"partTab",return(1));
  SAFE_CALLOC(partTab,vertNbr,SCOTCH_Num);

  /* Partionning the graph */
  CHECK_SCOTCH(SCOTCH_graphPart(&graf, boxNbr, &strat, partTab), "scotch_graphPart", 0);

  ADD_MEM(mesh,boxNbr*sizeof(SCOTCH_Num),"partNumTab",return(1));
  SAFE_CALLOC(partNumTab,boxNbr,SCOTCH_Num);

  if (!memset(partNumTab, 0, boxNbr*sizeof(SCOTCH_Num))) {
    perror("  ## Memory problem: memset");
    return 1;
  }

  /* Computing the number of elements of each box */
  for( vertIdx = 0 ; vertIdx< vertNbr ;vertIdx++)
    partNumTab[partTab[vertIdx]] += 1;


  /* partition permutation tabular */
  ADD_MEM(mesh,(vertNbr+1)*sizeof(SCOTCH_Num),"partPrmTab",return(1));
  SAFE_CALLOC(partPrmTab,vertNbr+1,SCOTCH_Num);

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

  DEL_MEM(mesh,partTab,vertNbr*sizeof(SCOTCH_Num));
  DEL_MEM(mesh,partNumTab,boxNbr*sizeof(SCOTCH_Num));
  DEL_MEM(mesh,partPrmTab,(vertNbr+1)*sizeof(SCOTCH_Num));

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
int kPartBoxCompute(SCOTCH_Graph graf, int vertNbr, int boxVertNbr,
                    SCOTCH_Num *permVrtTab,MMG5_pMesh mesh) {
  int boxNbr, vertIdx;
#if SCOTCH_VERSION<6
  SCOTCH_Num logMaxVal, SupMaxVal, InfMaxVal, maxVal;
#endif
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
#if SCOTCH_VERSION>=6
  CHECK_SCOTCH(SCOTCH_archCmplt(&arch, boxNbr), "scotch_archCmplt", 0) ;
#else
  CHECK_SCOTCH(SCOTCH_archVcmplt(&arch), "scotch_archVcmplt", 0) ;
#endif
  sprintf(s, "m{vert=%d,low=r{job=t,map=t,poli=S,sep=m{vert=80,low=h{pass=10}f{bal=0.0005,move=80},asc=f{bal=0.005,move=80}}}}", vertNbr / boxVertNbr);
  CHECK_SCOTCH(SCOTCH_stratGraphMap(&strat, s), "scotch_stratGraphMap", 0) ;

  ADD_MEM(mesh,2*vertNbr*sizeof(SCOTCH_Num),"sortPartTb",return(1));
  SAFE_CALLOC(sortPartTb,2*vertNbr,SCOTCH_Num);

  /* Partionning the graph */
  CHECK_SCOTCH(SCOTCH_graphMap(&graf, &arch, &strat, sortPartTb), "scotch_graphMap", 0);


#if SCOTCH_VERSION<6
  // Looking for the max value in sortPartTb and computing sortPartTb as
  // followed :
  //  - sortPartTb[2i] is the box value
  //  - sortPartTb[2i+1] is the vertex number
  maxVal = sortPartTb[0];
#endif
  for (vertIdx = vertNbr - 1 ; vertIdx >= 0 ; vertIdx--) {
    sortPartTb[2*vertIdx] = sortPartTb[vertIdx];
    sortPartTb[2*vertIdx+1] = vertIdx + 1;
#if SCOTCH_VERSION<6
    if (sortPartTb[vertIdx] > maxVal)
      maxVal = sortPartTb[vertIdx];
#endif
  }

#if SCOTCH_VERSION<6
  // Determining the log of MaxVal
  logMaxVal = 0;
  while ( maxVal > 0) {
    logMaxVal++;
    maxVal >>= 1;
  }

  // Infering the interval in which box values will be
  InfMaxVal = logMaxVal << logMaxVal;
  SupMaxVal = (logMaxVal << (logMaxVal + 1)) - 1;

  // Increasing box values until they are in the previous interval
  for (vertIdx = 0 ; vertIdx < vertNbr ; vertIdx++) {
    while (!(sortPartTb[2*vertIdx] >= InfMaxVal && sortPartTb[2*vertIdx] <= SupMaxVal)) {
      sortPartTb[2*vertIdx] <<= 1;
    }
  }
#endif

  // Sorting the tabular, which contains box values and vertex numbers
  _SCOTCHintSort2asc1(sortPartTb, vertNbr);


  /* Infering the new numbering */
  for (vertIdx = 0; vertIdx < vertNbr ; vertIdx++) {
    permVrtTab[sortPartTb[2*vertIdx + 1]] = vertIdx + 1;
  }

  SCOTCH_stratExit(&strat) ;
  SCOTCH_archExit(&arch) ;

  DEL_MEM(mesh,sortPartTb,2*vertNbr*sizeof(SCOTCH_Num));

  return 0;
}

/** Swap two tetras in the table of tetrahedras */
static inline
void swapTet(pTetra tetras/*, int* adja*/, int* perm, int ind1, int ind2) {
  Tetra pttmp;
  int   tmp;

/* Commentated part: swap for adja table if we don't free it in renumbering *
 * function (faster but need of 4*mesh->nemax*sizeof(int) extra bytes ) */
  /* int   adjatmp,kadj,ifadj,j; */

  /* /\* 1-- swap the adja table *\/ */
  /* /\* First: replace ind2 by ind1 in adjacents tetras of ind2 *\/ */
  /* for ( j=1; j<5; j++ ) { */
  /*   if ( adja[4*(ind2-1)+j]/4 ) { */
  /*     kadj  = adja[4*(ind2-1)+j]/4; */
  /*     ifadj = adja[4*(ind2-1)+j]%4; */
  /*     adja[4*(kadj-1)+1+ifadj] = 4*ind1 + adja[4*(kadj-1)+1+ifadj]%4; */
  /*   } */
  /* } */

  /* /\* Second: replace ind1 by ind2 in adjacents tetras of ind1*\/ */
  /* for ( j=1; j<5; j++ ) { */
  /*   if ( adja[4*(ind1-1)+j]/4 ) { */
  /*     kadj  = adja[4*(ind1-1)+j]/4; */
  /*     ifadj = adja[4*(ind1-1)+j]%4; */
  /*     if ( kadj == ind1 ) */
	/*       adja[4*(ind2-1)+1+ifadj] = 4*ind2 + adja[4*(ind2-1)+1+ifadj]%4; */
  /*     else */
	/*       adja[4*(kadj-1)+1+ifadj] = 4*ind2 + adja[4*(kadj-1)+1+ifadj]%4; */
  /*   } */
  /* } */

  /* /\* Third: swap adjacents for ind1 and ind2 *\/ */
  /* for ( j=1; j<5; j++ ) { */
  /*   adjatmp = adja[4*(ind2-1)+j]; */
  /*   adja[4*(ind2-1)+j] = adja[4*(ind1-1)+j]; */
  /*   adja[4*(ind1-1)+j] = adjatmp; */
  /* } */

  /* 2-- swap the tetrahedras */
  memcpy(&pttmp       ,&tetras[ind2],sizeof(Tetra));
  memcpy(&tetras[ind2],&tetras[ind1],sizeof(Tetra));
  memcpy(&tetras[ind1],&pttmp       ,sizeof(Tetra));

  /* 3-- swap the permutation table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}

/** Swap two points in the table of points and in the sol table */
static inline
void swapNod(pPoint points, double* sols, int* perm,
             int ind1, int ind2, int solsiz) {
  Point ptttmp;
  Sol   soltmp;
  int   tmp,addr2,addr1;

  /* swap the points */
  memcpy(&ptttmp      ,&points[ind2],sizeof(Point));
  memcpy(&points[ind2],&points[ind1],sizeof(Point));
  memcpy(&points[ind1],&ptttmp      ,sizeof(Point));

  /* swap the sols */
  addr1 = (ind1-1)*solsiz + 1;
  addr2 = (ind2-1)*solsiz + 1;
  memcpy(&soltmp     ,&sols[addr2],solsiz*sizeof(double));
  memcpy(&sols[addr2],&sols[addr1],solsiz*sizeof(double));
  memcpy(&sols[addr1],&soltmp     ,solsiz*sizeof(double));

  /* swap the permutaion table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}


/** Function : renumbering
 *  it modifies the numbering of each node to prevent from cache missing.
 *
 *  - boxVertNbr : number of vertices by box
 *  - mesh : the input mesh which is modified
 *
 *  returning 0 if OK, 1 else
 */
int renumbering(int boxVertNbr, pMesh mesh, pSol sol) {
  pPoint ppt;
  pTetra ptet;
  SCOTCH_Num edgeNbr;
  SCOTCH_Num *vertTab, *edgeTab, *permVrtTab;
  SCOTCH_Graph graf ;
  int    vertNbr, nodeGlbIdx, tetraIdx, ballTetIdx;
  int    i, j, k;
  int    edgeSiz;
  int    *vertOldTab, *permNodTab, nereal, npreal;
  int    *adja,iadr;


  /* Computing the number of vertices and a contiguous tabular of vertices */
  vertNbr = 0;

  ADD_MEM(mesh,(mesh->ne+1)*sizeof(int),"vertOldTab",return(1));
  SAFE_CALLOC(vertOldTab,mesh->ne+1,int);

  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {

    /* Testing if the tetra exists */
    if (!mesh->tetra[tetraIdx].v[0]) continue;
    vertOldTab[tetraIdx] = ++vertNbr;
  }

  if ( vertNbr/2 < BOXSIZE ) {
    /* not enough tetra to renum */
    DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
    return(1);
  }

  /* Allocating memory to compute adjacency lists */
  ADD_MEM(mesh,(vertNbr+2)*sizeof(SCOTCH_Num),"vertTab",
          DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
          return(1));
  SAFE_CALLOC(vertTab,vertNbr+2,SCOTCH_Num);

  if (!memset(vertTab, ~0, sizeof(SCOTCH_Num)*(vertNbr + 2))) {
    perror("  ## Memory problem: memset");
    DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
    DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
    return 1;
  }

  edgeNbr = 1;
  /* Euler-Poincare formulae edgeSiz = 20*mesh->np~4*mesh->ne;
     (2*(12*mesh->np (triangles)-2*mesh->np (boundary triangles))) */
  edgeSiz = vertNbr*4;

  ADD_MEM(mesh,edgeSiz*sizeof(SCOTCH_Num),"edgeTab",
          DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
          DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
          return(1));
  SAFE_CALLOC(edgeTab,edgeSiz,SCOTCH_Num);


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

      /* Testing if edgeTab memory is enough */
      if (edgeNbr >= edgeSiz) {
        ADD_MEM(mesh,0.2*sizeof(SCOTCH_Num),"edgeTab",
                DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
                DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
                return(1));
        edgeSiz *= 1.2;
        SAFE_REALLOC(edgeTab,edgeSiz,SCOTCH_Num);
      }

      edgeTab[edgeNbr++] = vertOldTab[ballTetIdx];
    }
  }
  vertTab[vertNbr+1] = edgeNbr;
  edgeNbr--;

  /* free adjacents to gain memory space */
  DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

  /* Building the graph by calling Scotch functions */
  SCOTCH_graphInit(&graf) ;
  CHECK_SCOTCH(SCOTCH_graphBuild(&graf, (SCOTCH_Num) 1, vertNbr, vertTab+1,
                                 NULL, NULL, NULL, edgeNbr, edgeTab+1, NULL),
               "scotch_graphbuild", 0) ;
#ifndef NDEBUG
  /* don't check in release mode */
  CHECK_SCOTCH(SCOTCH_graphCheck(&graf), "scotch_graphcheck", 0);
#endif

  ADD_MEM(mesh,(vertNbr+1)*sizeof(SCOTCH_Num),"permVrtTab",
          DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
          DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
          DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
          if( !hashTetra(mesh,1) ) return(0);
          return(1));
  SAFE_CALLOC(permVrtTab,vertNbr+1,SCOTCH_Num);

  CHECK_SCOTCH(kPartBoxCompute(graf, vertNbr, boxVertNbr, permVrtTab, mesh),
               "boxCompute", 0);

  SCOTCH_graphExit(&graf) ;

  DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
  DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));

  /* Computing the new point list and modifying the adja strcuture */
  ADD_MEM(mesh,(mesh->np+1)*sizeof(int),"permNodTab",
          DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
          DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
          DEL_MEM(mesh,permVrtTab,(vertNbr+1)*sizeof(SCOTCH_Num));
          DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
          if( !hashTetra(mesh,1) ) return(0);
          return(1));
  SAFE_CALLOC(permNodTab,mesh->np+1,int);

  nereal = 0;
  npreal = 0;

  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {
    ptet = &mesh->tetra[tetraIdx];

    /* Testing if the tetra exists */
    if (!ptet->v[0]) continue;

    nereal++;

    for(j = 0 ; j <= 3 ; j++) {

      nodeGlbIdx = ptet->v[j];

      if ( permNodTab[nodeGlbIdx] ) continue;

      ppt = &mesh->point[nodeGlbIdx];

      if ( !(ppt->tag & MG_NUL) )
        /* Building the new point list */
        permNodTab[nodeGlbIdx] = ++npreal;
    }
  }

  /* Create the final permutation table for tetras (stored in vertOldTab) and *
     modify the numbering of the nodes of each tetra */
  for( tetraIdx = 1; tetraIdx < mesh->ne + 1; tetraIdx++) {
    if ( !mesh->tetra[tetraIdx].v[0] )  continue;
    vertOldTab[tetraIdx] = permVrtTab[vertOldTab[tetraIdx]];
    for(j = 0 ; j <= 3 ; j++) {
      mesh->tetra[tetraIdx].v[j] = permNodTab[mesh->tetra[tetraIdx].v[j]];
    }
  }
  DEL_MEM(mesh,permVrtTab,(vertNbr+1)*sizeof(SCOTCH_Num));

  /* Permute nodes and sol */
  for (j=1; j<= mesh->np; j++) {
    while ( permNodTab[j] != j && permNodTab[j] )
      swapNod(mesh->point,sol->m,permNodTab,j,permNodTab[j],sol->size);
  }
  DEL_MEM(mesh,permNodTab,(mesh->np+1)*sizeof(int));

  /* Permute tetrahedras */
  for (j=1; j<= mesh->ne; j++) {
    while ( vertOldTab[j] != j && vertOldTab[j] )
      swapTet(mesh->tetra/*,mesh->adja*/,vertOldTab,j,vertOldTab[j]);
  }
  DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));

  mesh->ne = nereal;
  mesh->np = npreal;

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

  if( !hashTetra(mesh,0) ) return(0);

  return 1;
}
#endif
