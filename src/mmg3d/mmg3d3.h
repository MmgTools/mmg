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
 * \file mmg3d/mmg3d3.h
 * \brief Common functions and structures between MMG3D5 ans ELASTIC libraries.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#define _MMG5_DEGTOL    0.75
#define _MMG5_DISPREF   0
#define _LS_REFDIR      21
#define _LS_Ver         (1<<0)
#define _LS_Tri         (1<<2)
#define _LS_LAMBDA      10.0e5
#define _LS_MU          8.2e5

enum{ Dirichlet=1, Load};

/* Auxiliary mesh structures, compatible with the elastic library */
typedef struct {
  double   c[3];
  int      ref,old;
} LS_Point;
typedef LS_Point * LS_pPoint;

typedef struct {
  int     v[3],ref;
} LS_Edge;
typedef LS_Edge * LS_pEdge;

typedef struct {
  int     v[6],ref;
} LS_Tria;
typedef LS_Tria * LS_pTria;

typedef struct {
  int     v[10],ref;
} LS_Tetra;
typedef LS_Tetra * LS_pTetra;

typedef struct {
  double   u[3];
  int      ref;
  char     typ,elt,att;
} LS_Cl;
typedef LS_Cl * LS_pCl;

typedef struct {
  double   lambda,mu;
  int      ref;
} LS_Mat;
typedef LS_Mat * LS_pMat;

typedef struct {
  int       dim,ver,np,na,ne,nit,iter,nbcl,nmat;
  double   *u,*bc,*u0,err;
  char     *namein,*nameout,cltyp;
  LS_pCl   cl;
  LS_pMat  mat;
} LS_Sol;
typedef LS_Sol * LS_pSol;

typedef struct {
  int      np,np2,na,nt,ne,npi,nai,nti,nei,hmax,hcur,ver,dim,mark;
  char    *name;
  
  LS_pPoint   point;
  LS_pEdge    edge;
  LS_pTria    tria;
  LS_pTetra   tetra;
} LS_Mesh;
typedef LS_Mesh * LS_pMesh;