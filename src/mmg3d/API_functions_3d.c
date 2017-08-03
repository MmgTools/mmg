/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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
 * \file mmg3d/API_functions_3d.c
 * \brief C API functions definitions for MMG3D library.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning Use the MMG3D_ prefix: MMG5_ prefix will became obsolete soon...
 *
 * \note This file contains some internal functions for the API, see
 * the \ref mmg3d/libmmg3d.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * C API for MMG3D library.
 *
 */

#include "mmg3d.h"

int MMG3D_Init_mesh(const int starter,...) {
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = _MMG3D_Init_mesh_var(argptr);

  va_end(argptr);

  return ier;
}
void MMG3D_Init_fileNames(MMG5_pMesh mesh,MMG5_pSol sol
  ) {

  MMG5_Init_fileNames(mesh,sol);
  return;
}

int MMG3D_Set_inputMeshName(MMG5_pMesh mesh, const char* meshin) {

  return(MMG5_Set_inputMeshName(mesh,meshin));
}

int MMG3D_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solin) {
  return(MMG5_Set_inputSolName(mesh,sol,solin));
}

int MMG3D_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout) {

  return(MMG5_Set_outputMeshName(mesh,meshout));
}

int MMG3D_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solout) {
  return(MMG5_Set_outputSolName(mesh,sol,solout));
}

void MMG3D_Init_parameters(MMG5_pMesh mesh) {

  /* Init common parameters for mmgs and mmg3d. */
  _MMG5_Init_parameters(mesh);

  /* default values for integers */
  /* MMG3D_IPARAM_lag = -1 */
  mesh->info.lag      = -1;
  /* MMG3D_IPARAM_nofem = 0 */
  mesh->info.fem      = 1;
  /* MMG3D_IPARAM_optim = 0 */
  mesh->info.optim    =  0;
  /* MMG3D_IPARAM_optimLES = 0 */
  mesh->info.optimLES  =  0;
  /* MMG3D_IPARAM_nosurf = 0 */
  mesh->info.nosurf   =  0;  /* [0/1]    ,avoid/allow surface modifications */
#ifdef USE_SCOTCH
  mesh->info.renum    = 1;   /* [1/0]    , Turn on/off the renumbering using SCOTCH; */
#else
  mesh->info.renum    = 0;   /* [0]    , Turn on/off the renumbering using SCOTCH; */
#endif

  /* default values for doubles */
  mesh->info.ls       = 0.0;      /* level set value */


#ifndef PATTERN
  /* MMG3D_IPARAM_octree = 64 */
  mesh->info.octree = 32;
#endif
}

int MMG3D_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, int np, int typSol) {

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) && sol->m )
    fprintf(stderr,"\n  ## Warning: %s: old solution deletion.\n",__func__);

  if ( typEntity != MMG5_Vertex ) {
    fprintf(stderr,"\n  ## Error: %s: mmg3d need a solution imposed on vertices.\n",
            __func__);
    return(0);
  }
  if ( typSol == MMG5_Scalar ) {
    sol->size = 1;
  }
  else if ( typSol == MMG5_Vector ) {
    sol->size = 3;
  }
  else if ( typSol == MMG5_Tensor ) {
    sol->size = 6;
  }
  else {
    fprintf(stderr,"\n  ## Error: %s: type of solution not yet implemented.\n",
            __func__);
    return(0);
  }

  sol->dim = 3;
  if ( np ) {
    sol->np  = np;
    sol->npi = np;
    if ( sol->m )
      _MMG5_DEL_MEM(mesh,sol->m,(sol->size*(sol->npmax+1))*sizeof(double));

    sol->npmax = mesh->npmax;
    _MMG5_ADD_MEM(mesh,(sol->size*(sol->npmax+1))*sizeof(double),"initial solution",
                  return 0);
    _MMG5_SAFE_CALLOC(sol->m,(sol->size*(sol->npmax+1)),double,0);
  }
  return(1);
}

int MMG3D_Set_meshSize(MMG5_pMesh mesh, int np, int ne, int nprism,
                       int nt, int nquad, int na ) {
  int k;

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) &&
       ( mesh->point || mesh->tria || mesh->tetra || mesh->edge) )
    fprintf(stderr,"\n  ## Warning: %s: old mesh deletion.\n",__func__);

  if ( !np ) {
    fprintf(stderr,"  ** MISSING DATA:\n");
    fprintf(stderr,"     Your mesh must contains at least points.\n");
    return(0);
  }
  if ( !ne && (mesh->info.imprim > 4 || mesh->info.ddebug) ) {
    fprintf(stderr,"  ** WARNING:\n");
    fprintf(stderr,"     Your mesh don't contains tetrahedra.\n");
  }
  if ( mesh->point )
    _MMG5_DEL_MEM(mesh,mesh->point,(mesh->npmax+1)*sizeof(MMG5_Point));
  if ( mesh->tetra )
    _MMG5_DEL_MEM(mesh,mesh->tetra,(mesh->nemax+1)*sizeof(MMG5_Tetra));
  if ( mesh->prism )
    _MMG5_DEL_MEM(mesh,mesh->prism,(mesh->nprism+1)*sizeof(MMG5_Prism));
  if ( mesh->tria )
    _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(MMG5_Tria));
  if ( mesh->quadra )
    _MMG5_DEL_MEM(mesh,mesh->quadra,(mesh->nquad+1)*sizeof(MMG5_Quad));
  if ( mesh->edge )
    _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(MMG5_Edge));

  mesh->np  = np;
  mesh->ne  = ne;
  mesh->nt  = nt;
  mesh->na  = na;
  mesh->nprism = nprism;
  mesh->nquad  = nquad;

  mesh->npi = mesh->np;
  mesh->nei = mesh->ne;
  mesh->nti = mesh->nt;
  mesh->nai = mesh->na;

  /*tester si -m definie : renvoie 0 si pas ok et met la taille min dans info.mem */
  if( mesh->info.mem > 0) {
    if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne)) {
      if ( !_MMG3D_memOption(mesh) )  return 0;
    } else if(mesh->info.mem < 39) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory  %d\n",__func__,
              mesh->info.mem);
      return(0);
    }
  } else {
    mesh->npmax = MG_MAX(1.5*mesh->np,_MMG3D_NPMAX);
    mesh->nemax = MG_MAX(1.5*mesh->ne,_MMG3D_NEMAX);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,_MMG3D_NTMAX);

  }
  _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                return 0);
  _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point,0);


  _MMG5_ADD_MEM(mesh,(mesh->nemax+1)*sizeof(MMG5_Tetra),"initial tetrahedra",
                return 0);
  _MMG5_SAFE_CALLOC(mesh->tetra,mesh->nemax+1,MMG5_Tetra,0);


  if ( mesh->nprism ) {
    _MMG5_ADD_MEM(mesh,(mesh->nprism+1)*sizeof(MMG5_Prism),"initial prisms",return(0));
    _MMG5_SAFE_CALLOC(mesh->prism,(mesh->nprism+1),MMG5_Prism,0);
  }

  if ( mesh->nt ) {
    _MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(MMG5_Tria),"initial triangles",return(0));
    _MMG5_SAFE_CALLOC(mesh->tria,mesh->nt+1,MMG5_Tria,0);
  }

  if ( mesh->nquad ) {
    _MMG5_ADD_MEM(mesh,(mesh->nquad+1)*sizeof(MMG5_Quad),"initial quadrilaterals",return(0));
    _MMG5_SAFE_CALLOC(mesh->quadra,(mesh->nquad+1),MMG5_Quad,0);
  }

  mesh->namax = mesh->na;
  if ( mesh->na ) {
    _MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"initial edges",return(0));
    _MMG5_SAFE_CALLOC(mesh->edge,(mesh->na+1),MMG5_Edge,0);
  }

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->ne + 1;
  for (k=mesh->npnil; k<mesh->npmax-1; k++) {
    mesh->point[k].tmp  = k+1;
  }
  for (k=mesh->nenil; k<mesh->nemax-1; k++) {
    mesh->tetra[k].v[3] = k+1;
  }

  /* stats */
  if ( abs(mesh->info.imprim) > 6 ) {
    fprintf(stdout,"     NUMBER OF VERTICES       %8d\n",mesh->np);
    if ( mesh->na ) {
      fprintf(stdout,"     NUMBER OF EDGES          %8d\n",mesh->na);
    }
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES      %8d\n",mesh->nt);
    if ( mesh->nquad )
      fprintf(stdout,"     NUMBER OF QUADRILATERALS %8d\n",mesh->nquad);

    fprintf(stdout,"     NUMBER OF TETRAHEDRA     %8d\n",mesh->ne);

    if ( mesh->nprism )
      fprintf(stdout,"     NUMBER OF PRISMS         %8d\n",mesh->nprism);
  }

  return(1);
}

int MMG3D_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, int* np, int* typSol) {

  if ( typEntity != NULL )
    *typEntity = MMG5_Vertex;

  if ( typSol != NULL ) {
    if ( sol->size == 1 )
      *typSol    = MMG5_Scalar;
    else if ( sol->size == 3 )
      *typSol    = MMG5_Vector;
    else if ( sol->size == 6 )
      *typSol    = MMG5_Tensor;
    else
      *typSol    = MMG5_Notype;
  }

  assert( (!sol->np) || (sol->np == mesh->np));

  if ( np != NULL )
    *np = sol->np;

  return(1);
}

int MMG3D_Get_meshSize(MMG5_pMesh mesh, int* np, int* ne, int* nprism,
                       int* nt, int * nquad, int* na) {

  if ( np != NULL )
    *np = mesh->np;
  if ( ne != NULL )
    *ne = mesh->ne;
  if ( nprism != NULL )
    *nprism = mesh->nprism;
  if ( nt != NULL )
    *nt = mesh->nt;
  if ( nquad != NULL )
    *nquad = mesh->nquad;
  if ( na != NULL )
    *na = mesh->na;

  return(1);
}

int MMG3D_Set_vertex(MMG5_pMesh mesh, double c0, double c1, double c2, int ref, int pos) {

  if ( !mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: you must set the number of points with the",
            __func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting vertices in mesh.\n");
    return(0);
  }

  if ( pos > mesh->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate a new point.\n",__func__);
    fprintf(stderr,"    max number of points: %d\n",mesh->npmax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new vertex at position %d.",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of vertices: %d\n",mesh->np);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the vertex.\n");
    return(0);
  }

  mesh->point[pos].c[0] = c0;
  mesh->point[pos].c[1] = c1;
  mesh->point[pos].c[2] = c2;
  mesh->point[pos].ref  = ref;
  mesh->point[pos].tag  = MG_NUL;
  mesh->point[pos].flag = 0;
  mesh->point[pos].tmp = 0;

  return(1);
}

int MMG3D_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, int* ref,
                     int* isCorner, int* isRequired) {

  if ( mesh->npi == mesh->np ) {
    mesh->npi = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
        __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_vertex function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %d\n ",mesh->np);
    }
  }

  mesh->npi++;

  if ( mesh->npi > mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get point.\n",__func__);
    fprintf(stderr,"     The number of call of MMG3D_Get_vertex function");
    fprintf(stderr," can not exceed the number of points: %d\n ",mesh->np);
    return(0);
  }

  *c0  = mesh->point[mesh->npi].c[0];
  *c1  = mesh->point[mesh->npi].c[1];
  *c2  = mesh->point[mesh->npi].c[2];
  if ( ref != NULL )
    *ref = mesh->point[mesh->npi].ref;

  if ( isCorner != NULL ) {
    if ( mesh->point[mesh->npi].tag & MG_CRN )
      *isCorner = 1;
    else
      *isCorner = 0;
  }

  if ( isRequired != NULL ) {
    if ( mesh->point[mesh->npi].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return(1);
}

int  MMG3D_Set_vertices(MMG5_pMesh mesh, double *vertices,int *refs) {

  MMG5_pPoint ppt;
  int i,j;

  /*coordinates vertices*/
  for (i=1;i<=mesh->np;i++)
  {
    ppt = &mesh->point[i];

    j = (i-1)*3;
    ppt->c[0]  = vertices[j];
    ppt->c[1]  = vertices[j+1];
    ppt->c[2]  = vertices[j+2];

    ppt->tag = MG_NUL;
    ppt->flag = 0;
    ppt->tmp = 0;

    if ( refs != NULL )
      ppt->ref   = refs[i-1];
  }

  return 1;
}


int  MMG3D_Get_vertices(MMG5_pMesh mesh, double* vertices, int* refs,
                        int* areCorners, int* areRequired) {
  MMG5_pPoint ppt;
  int i,j;

  for (i=1;i<=mesh->np;i++)
  {
    ppt = &mesh->point[i];

    j = (i-1)*3;
    vertices[j] = ppt->c[0];
    vertices[j+1] = ppt->c[1];
    vertices[j+2] = ppt->c[2];

    j = i-1;
    if ( refs != NULL )
      refs[j] = ppt->ref;

    if ( areCorners !=NULL ) {
      if ( ppt->tag & MG_CRN )
        areCorners[j] = 1;
      else
        areCorners[j] = 0;
    }

    if ( areRequired != NULL ) {
      if ( ppt->tag & MG_REQ )
        areRequired[j] = 1;
      else
        areRequired[j] = 0;
    }
  }

  return 1;
}

int MMG3D_Set_tetrahedron(MMG5_pMesh mesh, int v0, int v1, int v2, int v3, int ref, int pos) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  double vol;
  int    aux,j, ip;

  if ( !mesh->ne ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of elements with the",
      __func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting elements in mesh\n");
    return(0);
  }

  if ( pos > mesh->nemax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate a new element.\n",
      __func__);
    fprintf(stderr,"    max number of element: %d\n",mesh->nemax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->ne ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new tetrahedron at position %d.",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of tetrahedron: %d\n",mesh->ne);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the tetrahedron.\n");
    return(0);
  }

  pt = &mesh->tetra[pos];
  pt->v[0] = v0;
  pt->v[1] = v1;
  pt->v[2] = v2;
  pt->v[3] = v3;
  pt->ref  = abs(ref);

  mesh->point[pt->v[0]].tag &= ~MG_NUL;
  mesh->point[pt->v[1]].tag &= ~MG_NUL;
  mesh->point[pt->v[2]].tag &= ~MG_NUL;
  mesh->point[pt->v[3]].tag &= ~MG_NUL;

  vol = _MMG5_orvol(mesh->point,pt->v);
  if ( vol <= _MMG5_EPSD2 ) {
    fprintf(stderr,"\n  ## Error: %s: tetrahedron %d has volume null.\n",
            __func__,pos);
    for ( ip=0; ip<4; ip++ ) {
      ppt = &mesh->point[pt->v[ip]];
      for ( j=0; j<3; j++ ) {
        if ( fabs(ppt->c[j])>0. ) {
          fprintf(stderr," Check that you don't have a sliver tetrahedron.\n");
          return(0);
        }
      }
    }
    fprintf(stderr,"  All vertices have zero coordinates.");
    fprintf(stderr," Check that you have set the vertices before the tetrahedra.\n");
    return(0);
  }
  else if ( vol < 0.0 ) {
    /* Possibly switch 2 vertices number so that each tet is positively oriented */
    aux = pt->v[2];
    pt->v[2] = pt->v[3];
    pt->v[3] = aux;
    /* mesh->xt temporary used to count reoriented tetra */
    mesh->xt++;
  }

  return(1);
}

int MMG3D_Get_tetrahedron(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* v3,
                          int* ref, int* isRequired) {

  if ( mesh->nei == mesh->ne ) {
    mesh->nei = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of"
              " tetrahedra.\n",__func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_tetrahedron function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of tetrahedron: %d\n ",mesh->ne);
    }
  }

  mesh->nei++;

  if ( mesh->nei > mesh->ne ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get tetra.\n",__func__);
    fprintf(stderr,"    The number of call of MMG3D_Get_tetrahedron function");
    fprintf(stderr," can not exceed the number of tetra: %d\n ",mesh->ne);
    return(0);
  }

  *v0  = mesh->tetra[mesh->nei].v[0];
  *v1  = mesh->tetra[mesh->nei].v[1];
  *v2  = mesh->tetra[mesh->nei].v[2];
  *v3  = mesh->tetra[mesh->nei].v[3];
  if ( ref != NULL ) {
    *ref = mesh->tetra[mesh->nei].ref;
  }

  if ( isRequired != NULL ) {
    if ( mesh->tetra[mesh->nei].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return(1);
}

int  MMG3D_Set_tetrahedra(MMG5_pMesh mesh, int *tetra, int *refs) {
  MMG5_pPoint ppt;
  MMG5_pTetra pt;
  double     vol;
  int        i,ip,aux, j;

  mesh->xp = 0;
  for (i=1;i<=mesh->ne;i++)
  {
    j = (i-1)*4;
    pt = &mesh->tetra[i];
    pt->v[0]  = tetra[j];
    pt->v[1]  = tetra[j+1];
    pt->v[2]  = tetra[j+2];
    pt->v[3]  = tetra[j+3];

    if ( refs != NULL )
      pt->ref   = abs(refs[i-1]);

    mesh->point[pt->v[0]].tag &= ~MG_NUL;
    mesh->point[pt->v[1]].tag &= ~MG_NUL;
    mesh->point[pt->v[2]].tag &= ~MG_NUL;
    mesh->point[pt->v[3]].tag &= ~MG_NUL;

    vol = _MMG5_orvol(mesh->point,pt->v);

    if ( vol <= _MMG5_EPSD2 ) {
      fprintf(stderr,"\n  ## Error: %s: tetrahedron %d has volume null.\n",
              __func__,i);

      for ( ip=0; ip<4; ip++ ) {
        ppt = &mesh->point[pt->v[ip]];
        for ( j=0; j<3; j++ ) {
          if ( fabs(ppt->c[j])>0. ) {
            fprintf(stderr," Check that you don't have a sliver tetrahedron.\n");
            return(0);
          }
        }
      }

      fprintf(stderr,"  All vertices have zero coordinates.");
      fprintf(stderr," Check that you have set the vertices before the tetrahedra.\n");
      return(0);
    }
    else if ( vol < 0.0 ) {
      /* Possibly switch 2 vertices number so that each tet is positively oriented */
      aux = pt->v[2];
      pt->v[2] = pt->v[3];
      pt->v[3] = aux;

      ++mesh->xp;
    }
  }

  return 1;
}

int  MMG3D_Get_tetrahedra(MMG5_pMesh mesh, int *tetra, int *refs, int * areRequired) {
  MMG5_pTetra pt;
  int         i, j;

  for (i=1;i<=mesh->ne;i++)
  {
    j = (i-1)*4;
    pt = &mesh->tetra[i];
    tetra[j]   = pt->v[0];
    tetra[j+2] = pt->v[1];
    tetra[j+1] = pt->v[2];
    tetra[j+3] = pt->v[3];
    if ( refs!=NULL )
      refs[i-1]  = pt->ref ;
    if ( areRequired != NULL ) {
      if ( pt->tag & MG_REQ )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }
  return 1;
}

int MMG3D_Set_prism(MMG5_pMesh mesh, int v0, int v1, int v2,
                    int v3, int v4, int v5, int ref, int pos) {
  MMG5_pPrism pp;

  if ( !mesh->nprism ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of prisms with the",
      __func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting elements in mesh\n");
    return(0);
  }

  if ( pos > mesh->nprism ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new prism at position %d.",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of prism: %d\n",mesh->nprism);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the prism.\n");
    return(0);
  }

  pp = &mesh->prism[pos];
  pp->v[0] = v0;
  pp->v[1] = v1;
  pp->v[2] = v2;
  pp->v[3] = v3;
  pp->v[4] = v4;
  pp->v[5] = v5;
  pp->ref  = ref;

  mesh->point[pp->v[0]].tag &= ~MG_NUL;
  mesh->point[pp->v[1]].tag &= ~MG_NUL;
  mesh->point[pp->v[2]].tag &= ~MG_NUL;
  mesh->point[pp->v[3]].tag &= ~MG_NUL;
  mesh->point[pp->v[4]].tag &= ~MG_NUL;
  mesh->point[pp->v[5]].tag &= ~MG_NUL;


  return(1);
}

int MMG3D_Get_prism(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* v3,
                    int* v4, int* v5, int* ref, int* isRequired) {
  static int npri = 0;

  if ( npri == mesh->nprism ) {
    npri = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of prisms.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_prism function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of prisms: %d\n ",mesh->nprism);
    }
  }

  ++npri;

  if ( npri > mesh->nprism ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get prism.\n",__func__);
    fprintf(stderr,"    The number of call of MMG3D_Get_prism function");
    fprintf(stderr," can not exceed the number of prism: %d\n ",mesh->nprism);
    return(0);
  }

  *v0  = mesh->prism[npri].v[0];
  *v1  = mesh->prism[npri].v[1];
  *v2  = mesh->prism[npri].v[2];
  *v3  = mesh->prism[npri].v[3];
  *v4  = mesh->prism[npri].v[4];
  *v5  = mesh->prism[npri].v[5];

  if ( ref != NULL ) {
    *ref = mesh->prism[npri].ref;
  }

  if ( isRequired != NULL ) {
    if ( mesh->prism[npri].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return(1);
}

int  MMG3D_Set_prisms(MMG5_pMesh mesh, int *prisms, int *refs) {
  MMG5_pPrism pp;
  int         i,j;

  for (i=1;i<=mesh->nprism;i++)
  {
    j = (i-1)*6;
    pp = &mesh->prism[i];
    pp->v[0]  = prisms[j];
    pp->v[1]  = prisms[j+1];
    pp->v[2]  = prisms[j+2];
    pp->v[3]  = prisms[j+3];
    pp->v[4]  = prisms[j+4];
    pp->v[5]  = prisms[j+5];

    if ( refs != NULL )
      pp->ref   = refs[i-1];

    mesh->point[pp->v[0]].tag &= ~MG_NUL;
    mesh->point[pp->v[1]].tag &= ~MG_NUL;
    mesh->point[pp->v[2]].tag &= ~MG_NUL;
    mesh->point[pp->v[3]].tag &= ~MG_NUL;
    mesh->point[pp->v[4]].tag &= ~MG_NUL;
    mesh->point[pp->v[5]].tag &= ~MG_NUL;

  }

  return 1;
}

int  MMG3D_Get_prisms(MMG5_pMesh mesh, int *prisms, int *refs, int * areRequired) {
  MMG5_pPrism pp;
  int         i, j;

  for (i=1;i<=mesh->nprism;i++)
  {
    j = (i-1)*6;
    pp = &mesh->prism[i];
    prisms[j]   = pp->v[0];
    prisms[j+2] = pp->v[1];
    prisms[j+1] = pp->v[2];
    prisms[j+3] = pp->v[3];
    prisms[j+4] = pp->v[4];
    prisms[j+5] = pp->v[5];

    if ( refs!=NULL )
      refs[i-1]  = pp->ref ;
    if ( areRequired != NULL ) {
      if ( pp->tag & MG_REQ )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }
  return 1;
}



int MMG3D_Set_triangle(MMG5_pMesh mesh, int v0, int v1, int v2, int ref,int pos) {

  if ( !mesh->nt ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of triangles"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting triangles in mesh\n");
    return(0);
  }

  if ( pos > mesh->ntmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate a new triangle.\n",
            __func__);
    fprintf(stderr,"    max number of triangle: %d\n",mesh->ntmax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->nt ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new triangle at"
            " position %d.",__func__,pos);
    fprintf(stderr," Overflow of the given number of triangles: %d\n",mesh->nt);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the triangle.\n");
    return(0);
  }

  mesh->tria[pos].v[0] = v0;
  mesh->tria[pos].v[1] = v1;
  mesh->tria[pos].v[2] = v2;
  mesh->tria[pos].ref  = ref;

  return(1);
}

int MMG3D_Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref
                       ,int* isRequired) {
  MMG5_pTria  ptt;

  if ( mesh->nti == mesh->nt ) {
    mesh->nti = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of"
              " triangles.\n",__func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_triangle function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of triangles: %d\n ",mesh->nt);
    }
  }

  mesh->nti++;

  if ( mesh->nti > mesh->nt ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get triangle.\n",__func__);
    fprintf(stderr,"    The number of call of MMG3D_Get_triangle function");
    fprintf(stderr," can not exceed the number of triangles: %d\n ",mesh->nt);
    return(0);
  }

  ptt = &mesh->tria[mesh->nti];
  *v0  = ptt->v[0];
  *v1  = ptt->v[1];
  *v2  = ptt->v[2];
  if ( ref != NULL )
    *ref = ptt->ref;

  if ( isRequired != NULL ) {
    if ( (ptt->tag[0] & MG_REQ) && (ptt->tag[1] & MG_REQ) &&
         (ptt->tag[2] & MG_REQ) )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return(1);
}

int  MMG3D_Set_triangles(MMG5_pMesh mesh, int *tria, int *refs) {
  MMG5_pTria ptt;
  int         i, j;

  for (i=1;i<=mesh->nt;i++)
  {
    j = (i-1)*3;
    ptt = &mesh->tria[i];
    ptt->v[0] = tria[j]  ;
    ptt->v[1] = tria[j+2];
    ptt->v[2] = tria[j+1];
    if ( refs != NULL )
      ptt->ref  = refs[i-1];
  }
  return 1;
}

int  MMG3D_Get_triangles(MMG5_pMesh mesh, int *tria, int *refs, int *areRequired) {
  MMG5_pTria ptt;
  int         i, j;

  for (i=1;i<=mesh->nt;i++)
  {
    j = (i-1)*3;
    ptt = &mesh->tria[i];
    tria[j]   = ptt->v[0];
    tria[j+2] = ptt->v[1];
    tria[j+1] = ptt->v[2];

    if ( refs!=NULL )
      refs[i-1]  = ptt->ref ;
    if ( areRequired != NULL ) {
      if ( (ptt->tag[0] & MG_REQ) && (ptt->tag[1] & MG_REQ) &&
           (ptt->tag[2] & MG_REQ) )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }
  return 1;
}

int MMG3D_Set_quadrilateral(MMG5_pMesh mesh, int v0, int v1, int v2, int v3,
                         int ref,int pos) {

  if ( !mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of quadrilaterals"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting quadrilaterals in mesh\n");
    return(0);
  }

  if ( pos > mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new quadrilateral"
            " at position %d.",__func__,pos);
    fprintf(stderr," Overflow of the given number of quadrilaterals: %d\n",mesh->nquad);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the quadrilateral.\n");
    return(0);
  }

  mesh->quadra[pos].v[0] = v0;
  mesh->quadra[pos].v[1] = v1;
  mesh->quadra[pos].v[2] = v2;
  mesh->quadra[pos].v[3] = v3;
  mesh->quadra[pos].ref  = ref;

  return(1);
}

int MMG3D_Get_quadrilateral(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* v3,
                       int* ref,int* isRequired) {
  MMG5_pQuad  pq;
  static int nqi = 0;

  if ( nqi == mesh->nquad ) {
    nqi = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter"
              " of quadrilaterals.\n",__func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_quadrilateral function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of quadrilaterals: %d\n ",mesh->nquad);
    }
  }

  nqi++;

  if ( nqi > mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get quadrilateral.\n",__func__);
    fprintf(stderr,"    The number of call of MMG3D_Get_quadrilateral function");
    fprintf(stderr," can not exceed the number of quadrilaterals: %d\n ",mesh->nquad);
    return(0);
  }

  pq = &mesh->quadra[nqi];
  *v0  = pq->v[0];
  *v1  = pq->v[1];
  *v2  = pq->v[2];
  *v3  = pq->v[3];
  if ( ref != NULL )
    *ref = pq->ref;

  if ( isRequired != NULL ) {
    if ( (pq->tag[0] & MG_REQ) && (pq->tag[1] & MG_REQ) &&
         (pq->tag[2] & MG_REQ) && (pq->tag[3] & MG_REQ))
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return(1);
}

int  MMG3D_Set_quadrilaterals(MMG5_pMesh mesh, int *quads, int *refs) {
  MMG5_pQuad  pq;
  int         i, j;

  for (i=1;i<=mesh->nquad;i++)
  {
    j = (i-1)*4;
    pq = &mesh->quadra[i];
    pq->v[0] = quads[j]  ;
    pq->v[1] = quads[j+1];
    pq->v[2] = quads[j+2];
    pq->v[3] = quads[j+3];
    if ( refs != NULL )
      pq->ref  = refs[i-1];
  }
  return 1;
}

int  MMG3D_Get_quadrilaterals(MMG5_pMesh mesh, int *quads, int *refs, int *areRequired) {
  MMG5_pQuad  pq;
  int         i, j;

  for (i=1;i<=mesh->nquad;i++)
  {
    j = (i-1)*4;
    pq = &mesh->quadra[i];
    quads[j]   = pq->v[0];
    quads[j+1] = pq->v[1];
    quads[j+2] = pq->v[2];
    quads[j+3] = pq->v[3];

    if ( refs!=NULL )
      refs[i-1]  = pq->ref ;
    if ( areRequired != NULL ) {
      if ( (pq->tag[0] & MG_REQ) && (pq->tag[1] & MG_REQ) &&
           (pq->tag[2] & MG_REQ) && (pq->tag[3] & MG_REQ) )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }
  return 1;
}

int MMG3D_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref, int pos) {

  if ( !mesh->na ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of edges with"
            " the",__func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting edges in mesh\n");
    return(0);
  }
  if ( pos > mesh->namax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate a new edge.\n",
      __func__);
    fprintf(stderr,"    max number of edge: %d\n",mesh->namax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }
  if ( pos > mesh->na ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new edge at position %d.",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of edges: %d\n",mesh->na);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the edge.\n");
    return(0);
  }

  mesh->edge[pos].a = v0;
  mesh->edge[pos].b = v1;
  mesh->edge[pos].ref  = ref;
  mesh->edge[pos].tag |= MG_REF;

  return(1);
}

int MMG3D_Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref
                   ,int* isRidge, int* isRequired) {

  if ( mesh->nai == mesh->na ) {
    mesh->nai = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of edges.\n",
        __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_edge function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of edges: %d\n ",mesh->na);
    }
  }

  mesh->nai++;

  if ( mesh->nai > mesh->na ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get edge.\n",__func__);
    fprintf(stderr,"    The number of call of MMG3D_Get_edge function");
    fprintf(stderr," can not exceed the number of edges: %d\n ",mesh->na);
    return(0);
  }

  *e0  = mesh->edge[mesh->nai].a;
  *e1  = mesh->edge[mesh->nai].b;
  if ( ref!=NULL )
    *ref = mesh->edge[mesh->nai].ref;

  if ( isRidge != NULL ) {
    if ( mesh->edge[mesh->nai].tag & MG_GEO )
      *isRidge = 1;
    else
      *isRidge = 0;
  }

  if ( isRequired != NULL ) {
    if ( mesh->edge[mesh->nai].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return(1);
}

int MMG3D_Set_corner(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag |= MG_CRN;
  return(1);
}

int MMG3D_Set_requiredVertex(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag |= MG_REQ;
  return(1);
}

int MMG3D_Set_requiredTetrahedron(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->ne );
  mesh->tetra[k].tag |= MG_REQ;
  return(1);
}

int MMG3D_Set_requiredTetrahedra(MMG5_pMesh mesh, int *reqIdx, int nreq) {
  int k;

  for ( k=0; k<nreq; ++k ){
    mesh->tetra[reqIdx[k]].tag |= MG_REQ;
  }

  return(1);
}

int MMG3D_Set_requiredTriangle(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->nt );
  mesh->tria[k].tag[0] |= MG_REQ;
  mesh->tria[k].tag[1] |= MG_REQ;
  mesh->tria[k].tag[2] |= MG_REQ;
  return(1);
}

int MMG3D_Set_requiredTriangles(MMG5_pMesh mesh, int* reqIdx, int nreq) {
  int k;

  for ( k=0; k<nreq; ++k ){
    mesh->tria[reqIdx[k]].tag[0] |= MG_REQ;
    mesh->tria[reqIdx[k]].tag[1] |= MG_REQ;
    mesh->tria[reqIdx[k]].tag[2] |= MG_REQ;
  }
  return(1);
}

int MMG3D_Set_parallelTriangle(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->nt );
  mesh->tria[k].tag[0] |= MG_PARBDY;
  mesh->tria[k].tag[1] |= MG_PARBDY;
  mesh->tria[k].tag[2] |= MG_PARBDY;
  return(1);
}

int MMG3D_Set_parallelTriangles(MMG5_pMesh mesh, int* parIdx, int npar) {
  int k;

  for ( k=0; k<npar; ++k ){
    mesh->tria[parIdx[k]].tag[0] |= MG_PARBDY;
    mesh->tria[parIdx[k]].tag[1] |= MG_PARBDY;
    mesh->tria[parIdx[k]].tag[2] |= MG_PARBDY;
  }
  return(1);
}

int MMG3D_Set_ridge(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->na );
  mesh->edge[k].tag |= MG_GEO;
  return(1);
}

int MMG3D_Set_requiredEdge(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->na );
  mesh->edge[k].tag |= MG_REQ;
  return(1);
}

int MMG3D_Set_normalAtVertex(MMG5_pMesh mesh, int k, double n0, double n1, double n2) {

  assert ( k <= mesh->np );
  mesh->point[k].n[0] = n0;
  mesh->point[k].n[1] = n1;
  mesh->point[k].n[2] = n2;

  ++mesh->nc1;

  return(1);
}

int MMG3D_Get_normalAtVertex(MMG5_pMesh mesh, int k, double *n0, double *n1, double *n2) {

  assert ( k <= mesh->np );
  (*n0) = mesh->point[k].n[0];
  (*n1) = mesh->point[k].n[1];
  (*n2) = mesh->point[k].n[2];

  return(1);
}

int MMG3D_Set_scalarSol(MMG5_pSol met, double s, int pos) {

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of"
            " solution with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }
  if ( pos < 1 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",__func__);
    fprintf(stderr,"    Minimal index of the solution position must be 1.\n");
    return(0);
  }
  if ( pos >= met->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",__func__);
    fprintf(stderr,"    max number of solutions: %d\n",met->npmax);
    return(0);
  }

  if ( pos > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new solution at"
            " position %d.",__func__,pos);
    fprintf(stderr," Overflow of the given number of solutions: %d\n",met->np);
    fprintf(stderr,"\n  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return(0);
  }

  met->m[pos] = s;
  return(1);
}


int MMG3D_Get_scalarSol(MMG5_pSol met, double* s) {

  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
        __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_scalarSol function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %d\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get solution.\n",__func__);
    fprintf(stderr,"     The number of call of MMG3D_Get_scalarSol function");
    fprintf(stderr," can not exceed the number of points: %d\n ",met->np);
    return(0);
  }

  *s  = met->m[met->npi];

  return(1);
}

int MMG3D_Set_scalarSols(MMG5_pSol met, double *s ) {
  int k;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of solution"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }

  for ( k=0; k<met->np; ++k )
    met->m[k+1] = s[k];

  return(1);
}

int MMG3D_Get_scalarSols(MMG5_pSol met, double* s) {
  int k;

  for ( k=0; k<met->np; ++k )
    s[k]  = met->m[k+1];

  return(1);
}

int MMG3D_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz, int pos) {

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of solution"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }
  if ( pos < 1 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",__func__);
    fprintf(stderr,"    Minimal index of the solution position must be 1.\n");
    return(0);
  }
  if ( pos >= met->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",__func__);
    fprintf(stderr,"    max number of solutions: %d\n",met->npmax);
    return(0);
  }

  if ( pos > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new solution at"
            " position %d.",__func__,pos);
    fprintf(stderr," Overflow of the given number of solutions: %d\n",met->np);
    fprintf(stderr,"\n  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return(0);
  }

  met->m[3*pos]   = vx;
  met->m[3*pos+1] = vy;
  met->m[3*pos+2] = vz;

  return(1);
}


int MMG3D_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz) {

  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_vectorSol function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %d\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get solution.\n",__func__);
    fprintf(stderr,"     The number of call of MMG3D_Get_vectorSol function");
    fprintf(stderr," can not exceed the number of points: %d\n ",met->np);
    return(0);
  }

  *vx  = met->m[3*met->npi];
  *vy  = met->m[3*met->npi+1];
  *vz  = met->m[3*met->npi+2];

  return(1);
}

int MMG3D_Set_vectorSols(MMG5_pSol met, double *sols) {
  double *m;
  int k,j;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of solution"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }

  for ( k=0; k<met->np; ++k ) {
    j = 3*k;
    m = &met->m[j+3];
    m[0] = sols[j];
    m[1] = sols[j+1];
    m[2] = sols[j+2];
  }

  return(1);
}

int MMG3D_Get_vectorSols(MMG5_pSol met, double* sols) {
  double *m;
  int k, j;

  for ( k=0; k<met->np; ++k ) {
    j = 3*k;
    m = &met->m[j+3];
    sols[j]   = m[0];
    sols[j+1] = m[1];
    sols[j+2] = m[2];
  }

  return(1);
}

int MMG3D_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
                        double m22,double m23, double m33, int pos) {

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of solution"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }
  if ( pos < 1 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",
      __func__);
    fprintf(stderr,"    Minimal index of the solution position must be 1.\n");
    return(0);
  }
  if ( pos >= met->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",__func__);
    fprintf(stderr,"    max number of solutions: %d\n",met->npmax);
    return(0);
  }

  if ( pos > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new solution at "
            "position %d.",__func__,pos);
    fprintf(stderr," Overflow of the given number of solutions: %d\n",met->np);
    fprintf(stderr,"\n  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return(0);
  }

  met->m[6*pos]   = m11;
  met->m[6*pos+1] = m12;
  met->m[6*pos+2] = m13;
  met->m[6*pos+3] = m22;
  met->m[6*pos+4] = m23;
  met->m[6*pos+5] = m33;

  return(1);
}


int MMG3D_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
                        double *m22,double *m23, double *m33) {

  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
        __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_tensorSol function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %d\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get solution.\n",__func__);
    fprintf(stderr,"     The number of call of MMG3D_Get_tensorSol function");
    fprintf(stderr," can not exceed the number of points: %d\n ",met->np);
    return(0);
  }

  *m11 = met->m[6*met->npi];
  *m12 = met->m[6*met->npi+1];
  *m13 = met->m[6*met->npi+2];
  *m22 = met->m[6*met->npi+3];
  *m13 = met->m[6*met->npi+4];
  *m33 = met->m[6*met->npi+5];

  return(1);
}

int MMG3D_Set_tensorSols(MMG5_pSol met, double *sols) {
  double *m;
  int k,j;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of"
            " solution with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }

  for ( k=0; k<met->np; ++k ) {
    j = 6*k;
    m = &met->m[j+6];

    m[0] = sols[j];
    m[1] = sols[j+1];
    m[2] = sols[j+2];
    m[3] = sols[j+3];
    m[4] = sols[j+4];
    m[5] = sols[j+5];
  }
  return(1);
}

int MMG3D_Get_tensorSols(MMG5_pSol met, double *sols) {
  double *m;
  int k,j;

  for ( k=0; k<met->np; ++k ) {
    j = 6*k;
    m = &met->m[j+6];

    sols[j]   = m[0];
    sols[j+1] = m[1];
    sols[j+2] = m[2];
    sols[j+3] = m[3];
    sols[j+4] = m[4];
    sols[j+5] = m[5];
  }

  return(1);
}

void MMG3D_Set_handGivenMesh(MMG5_pMesh mesh) {
  int k, aux;

  /* Possibly switch 2 vertices number so that each tet is positively oriented */
  for (k=1; k<=mesh->ne; k++) {
    if ( _MMG5_orvol(mesh->point,mesh->tetra[k].v) < 0.0 ) {
      /* mesh->xt temporary used to count reoriented tetra */
      mesh->xt++;
      aux = mesh->tetra[k].v[2];
      mesh->tetra[k].v[2] = mesh->tetra[k].v[3];
      mesh->tetra[k].v[3] = aux;
    }
  }
  return;
}

int MMG3D_Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met) {

  if ( (mesh->npi != mesh->np) || (mesh->nei != mesh->ne) ) {
    fprintf(stderr,"\n  ## Error: %s: if you don't use the MMG3D_loadMesh"
            " function,",__func__);
    fprintf(stderr," you must call the MMG3D_Set_meshSize function to have a");
    fprintf(stderr," valid mesh.\n");
    fprintf(stderr," Missing datas.\n");
    return(0);
  }

  if ( met->npi != met->np ) {
    fprintf(stderr,"\n  ## Error: %s: if you don't use the MMG3D_loadSol"
            " function,",__func__);
    fprintf(stderr," you must call the MMG3D_Set_solSize function to have a");
    fprintf(stderr," valid solution.\n");
    fprintf(stderr," Missing datas.\n");
    return(0);
  }

  /*  Check mesh data */
  if ( mesh->info.ddebug ) {
    if ( (!mesh->np) || (!mesh->point) ||
         (!mesh->ne) || (!mesh->tetra) ) {
      fprintf(stderr,"  ** MISSING DATA.\n");
      fprintf(stderr," Check that your mesh contains points and tetrahedra.\n");
      fprintf(stderr," Exit program.\n");
      return(0);
    }
  }

  if ( mesh->dim != 3 ) {
    fprintf(stderr,"  ** 3 DIMENSIONAL MESH NEEDED. Exit program.\n");
    return(0);
  }
  if ( met->dim != 3 ) {
    fprintf(stderr,"  ** WRONG DIMENSION FOR METRIC. Exit program.\n");
    return(0);
  }
  if ( !mesh->ver )  mesh->ver = 2;
  if ( !met ->ver )  met ->ver = 2;

  return(1);
}

static inline
int _MMG3D_skipIso(MMG5_pMesh mesh) {
  MMG5_pTria  ptt,ptt1;
  MMG5_pEdge  pa,pa1;
  int    k;

  if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
    fprintf(stderr,"\n  ## Warning: %s: skip of all entites with %d reference.\n",
            __func__,MG_ISO);

  /* Skip triangles with MG_ISO refs */
  k = 1;
  do {
    ptt = &mesh->tria[k];
    if ( abs(ptt->ref) != MG_ISO ) continue;
    /* here ptt is the first tri of mesh->tria that we want to delete */
    do {
      ptt1 = &mesh->tria[mesh->nti];
    }
    while( (abs(ptt1->ref) == MG_ISO) && (k <= --mesh->nti) );

    if ( abs(ptt1->ref) != MG_ISO )
      /* ptt1 is the last tri of mesh->tria that we want to keep */
      memcpy(ptt,ptt1,sizeof(MMG5_Tria));
  } while( ++k <= mesh->nti );

  if ( mesh->nti < mesh->nt ) {
    if( !mesh->nti )
      _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(MMG5_Tria));
    else {
      _MMG5_ADD_MEM(mesh,mesh->nti-mesh->nt,"triangles",return(0));
      _MMG5_SAFE_RECALLOC(mesh->tria,mesh->nt+1,(mesh->nti+1),MMG5_Tria,"triangles",0);
    }
    mesh->nt = mesh->nti;
  }

  /* Skip edges with MG_ISO refs */
  k = 1;
  do {
    pa = &mesh->edge[k];
    if ( abs(pa->ref) != MG_ISO ) {
      pa->ref = abs(pa->ref);
      continue;
    }
    /* here pa is the first edge of mesh->edge that we want to delete */
    do {
      pa1 = &mesh->edge[mesh->nai];
    }
    while( (abs(pa1->ref) == MG_ISO) && (k <= --mesh->nai) );

    if ( abs(pa1->ref) != MG_ISO ) {
      /* pa1 is the last edge of mesh->edge that we want to keep */
      memcpy(pa,pa1,sizeof(MMG5_Edge));
      pa1->ref = abs(pa1->ref);
    }
  } while( ++k <= mesh->nai );

  if ( mesh->nai < mesh->na ) {
    if( !mesh->nai )
      _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->nai+1)*sizeof(MMG5_Edge));
    else {
      _MMG5_ADD_MEM(mesh,mesh->nai-mesh->na,"Edges",return(0));
      _MMG5_SAFE_RECALLOC(mesh->edge,mesh->na+1,(mesh->nai+1),MMG5_Edge,"edges",0);
    }
    mesh->na = mesh->nai;
  }

  /* delete tetrahedra references */
  for (k=1; k<=mesh->ne; k++) {
    mesh->tetra[k].ref = 0;
  }
  return(1);
}

int MMG3D_Set_iparameter(MMG5_pMesh mesh, MMG5_pSol sol, int iparam,int val){
  int k;

  switch ( iparam ) {
    /* Integer parameters */
  case MMG3D_IPARAM_verbose :
    mesh->info.imprim   = val;
    break;
  case MMG3D_IPARAM_mem :
    if ( val <= 0 ) {
      fprintf(stderr,"\n  ## Warning: %s: maximal memory authorized must be"
              " strictly positive.\n",__func__);
      fprintf(stderr,"  Reset to default value.\n");
    }
    else
      mesh->info.mem      = val;
    if ( !_MMG3D_memOption(mesh) )  return 0;
    break;
#ifndef PATTERN
  case MMG3D_IPARAM_octree :
    mesh->info.octree   = val;
    break;
#endif
  case MMG3D_IPARAM_debug :
    mesh->info.ddebug   = val;
    break;
  case MMG3D_IPARAM_angle :
    /* free table that may contains old ridges */
    if ( mesh->htab.geom )
      _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));
    if ( mesh->xpoint )
      _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
    if ( mesh->xtetra )
      _MMG5_DEL_MEM(mesh,mesh->xtetra,(mesh->xtmax+1)*sizeof(MMG5_xTetra));
    if ( !val )
      mesh->info.dhd    = -1.;
    else {
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stderr,"\n  ## Warning: %s: angle detection parameter set"
                " to default value\n",__func__);
      mesh->info.dhd    = _MMG5_ANGEDG;
    }
    break;
  case MMG3D_IPARAM_nofem :
    mesh->info.fem      = (val==1)? 0 : 1;
    break;
  case MMG3D_IPARAM_opnbdy :
    mesh->info.opnbdy = val;
    break;
  case MMG3D_IPARAM_iso :
    mesh->info.iso      = val;
    if ( mesh->info.iso )
      if ( mesh->nt && !_MMG3D_skipIso(mesh) )
        return 0;
    break;
  case MMG3D_IPARAM_lag :
#ifdef USE_ELAS
    if ( val < 0 || val > 2 )
      return 0;
    mesh->info.lag = val;
#else
    fprintf(stderr,"\n  ## Error: %s"
            " \"lagrangian motion\" option unavailable (-lag):\n"
            " set the USE_ELAS CMake's flag to ON when compiling the mmg3d"
            " library to enable this feature.\n",__func__);
    return(0);
#endif
    break;
  case MMG3D_IPARAM_optim :
    mesh->info.optim = val;
    break;
  case MMG3D_IPARAM_optimLES :
    mesh->info.optimLES = val;
    break;
  case MMG3D_IPARAM_noinsert :
    mesh->info.noinsert = val;
    break;
  case MMG3D_IPARAM_noswap :
    mesh->info.noswap   = val;
    break;
  case MMG3D_IPARAM_nomove :
    mesh->info.nomove   = val;
    break;
  case MMG3D_IPARAM_nosurf :
    mesh->info.nosurf   = val;
    break;
  case MMG3D_IPARAM_numberOfLocalParam :
    if ( mesh->info.par ) {
      _MMG5_DEL_MEM(mesh,mesh->info.par,mesh->info.npar*sizeof(MMG5_Par));
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stderr,"\n  ## Warning: %s: new local parameter values\n",__func__);
    }
    mesh->info.npar   = val;
    mesh->info.npari  = 0;
    mesh->info.parTyp = 0;

    _MMG5_ADD_MEM(mesh,mesh->info.npar*sizeof(MMG5_Par),"parameters",
                  printf("  Exit program.\n");
                  return 0);
    _MMG5_SAFE_CALLOC(mesh->info.par,mesh->info.npar,MMG5_Par,0);

    for (k=0; k<mesh->info.npar; k++) {
      mesh->info.par[k].elt   = MMG5_Noentity;
      mesh->info.par[k].ref   = INT_MAX;
      mesh->info.par[k].hausd = mesh->info.hausd;
      mesh->info.par[k].hmin  = mesh->info.hmin;
      mesh->info.par[k].hmax  = mesh->info.hmax;
    }

    break;
#ifdef USE_SCOTCH
  case MMG3D_IPARAM_renum :
    mesh->info.renum    = val;
    break;
#endif
  case MMG3D_IPARAM_anisosize :
    if ( !MMG3D_Set_solSize(mesh,sol,MMG5_Vertex,0,MMG5_Tensor) )
      return 0;

  default :
    fprintf(stderr,"\n  ## Error: %s: unknown type of parameter\n",__func__);
    return(0);
  }
  /* other options */

  return(1);
}

int MMG3D_Get_iparameter(MMG5_pMesh mesh, int iparam) {

  switch ( iparam ) {
    /* Integer parameters */
  case MMG3D_IPARAM_verbose :
    return ( mesh->info.imprim );
    break;
  case MMG3D_IPARAM_mem :
    return ( mesh->info.mem );
    break;
#ifndef PATTERN
  case MMG3D_IPARAM_octree :
    return ( mesh->info.octree );
    break;
#endif
  case MMG3D_IPARAM_debug :
    return ( mesh->info.ddebug );
    break;
  case MMG3D_IPARAM_angle :
    if ( mesh->info.dhd <= 0. ) {
      return ( 0 );
    }
    else {
      return ( 1 );
    }
    break;
  case MMG3D_IPARAM_iso :
    return ( mesh->info.iso );
    break;
  case MMG3D_IPARAM_lag :
    return ( mesh->info.lag );
    break;
  case MMG3D_IPARAM_noinsert :
    return ( mesh->info.noinsert );
    break;
  case MMG3D_IPARAM_noswap :
    return ( mesh->info.noswap );
    break;
  case MMG3D_IPARAM_nomove :
    return ( mesh->info.nomove );
    break;
  case MMG3D_IPARAM_nosurf :
    return ( mesh->info.nosurf );
    break;
  case MMG3D_IPARAM_numberOfLocalParam :
    return ( mesh->info.npar );
    break;
#ifdef USE_SCOTCH
  case MMG3D_IPARAM_renum :
    return ( mesh->info.renum );
    break;
#endif
  default :
    fprintf(stderr,"\n  ## Error: %s: unknown type of parameter\n",__func__);
    return 0;
  }
}

int MMG3D_Set_dparameter(MMG5_pMesh mesh, MMG5_pSol sol, int dparam, double val){

  switch ( dparam ) {
    /* double parameters */
  case MMG3D_DPARAM_angleDetection :
    mesh->info.dhd = val;
    mesh->info.dhd = MG_MAX(0.0, MG_MIN(180.0,mesh->info.dhd));
    mesh->info.dhd = cos(mesh->info.dhd*M_PI/180.0);
    break;
  case MMG3D_DPARAM_hmin :
    mesh->info.hmin     = val;
    break;
  case MMG3D_DPARAM_hmax :
    mesh->info.hmax     = val;
    break;
  case MMG3D_DPARAM_hsiz :
    mesh->info.hsiz     = val;
    break;
  case MMG3D_DPARAM_hgrad :
    mesh->info.hgrad    = val;
    if ( mesh->info.hgrad < 0.0 )
      mesh->info.hgrad = -1.0;
    else
      mesh->info.hgrad = log(mesh->info.hgrad);
    break;
  case MMG3D_DPARAM_hausd :
    if ( val <=0 ) {
      fprintf(stderr,"\n  ## Error: %s: hausdorff number must be strictly"
              " positive.\n",__func__);
      return(0);
    }
    else
      mesh->info.hausd    = val;
    break;
  case MMG3D_DPARAM_ls :
    mesh->info.ls       = val;
    break;
  default :
    fprintf(stderr,"\n  ## Error: %s: unknown type of parameter\n", __func__);
    return(0);
  }
  return(1);
}

int MMG3D_Set_localParameter(MMG5_pMesh mesh,MMG5_pSol sol, int typ, int ref,
                             double hmin,double hmax,double hausd){
  MMG5_pPar par;
  int k;

  if ( !mesh->info.npar ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of local"
            " parameters",__func__);
    fprintf(stderr," with the MMG3D_Set_iparameters function before setting");
    fprintf(stderr," values in local parameters structure. \n");
    return(0);
  }
  if ( mesh->info.npari >= mesh->info.npar ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new local parameter.\n",
            __func__);
    fprintf(stderr,"    max number of local parameters: %d\n",mesh->info.npar);
    return(0);
  }
  if ( typ != MMG5_Triangle && typ != MMG5_Tetrahedron ) {
    fprintf(stderr,"\n  ## Warning: %s: you must apply your local parameters",
      __func__);
    fprintf(stderr," on triangles (MMG5_Triangle or %d) or tetrahedron"
            " (MMG5_Tetrahedron or %d).\n",MMG5_Triangle,MMG5_Tetrahedron);
    fprintf(stderr,"\n  ## Unknown type of entity: ignored.\n");
    return(0);
  }
  if ( ref < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative references are not allowed.\n",
      __func__);
    return(0);
  }

  for (k=0; k<mesh->info.npari; k++) {
    par = &mesh->info.par[k];

    if ( par->elt == typ && par->ref == ref ) {
      par->hausd = hausd;
      par->hmin  = hmin;
      par->hmax  = hmax;
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
          fprintf(stderr,"\n  ## Warning: %s: new parameters (hausd, hmin and hmax)",
            __func__);
          fprintf(stderr," for entities of type %d and of ref %d\n",typ,ref);
      }
      return 1;
    }
  }

  mesh->info.par[mesh->info.npari].elt   = typ;
  mesh->info.par[mesh->info.npari].ref   = ref;
  mesh->info.par[mesh->info.npari].hmin  = hmin;
  mesh->info.par[mesh->info.npari].hmax  = hmax;
  mesh->info.par[mesh->info.npari].hausd = hausd;

  switch ( typ )
  {
  case ( MMG5_Vertex ):
    mesh->info.parTyp |= MG_Vert;
    break;
  case ( MMG5_Triangle ):
    mesh->info.parTyp |= MG_Tria;
    break;
  case ( MMG5_Tetrahedron ):
    mesh->info.parTyp |= MG_Tetra;
    break;
  }

  mesh->info.npari++;

  return(1);
}

int MMG3D_Free_all(const int starter,...)
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = _MMG3D_Free_all_var(argptr);

  va_end(argptr);

  return ier;
}

int MMG3D_Free_structures(const int starter,...)
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = _MMG3D_Free_structures_var(argptr);

  va_end(argptr);

  return ier;
}

int MMG3D_Free_names(const int starter,...)
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = _MMG3D_Free_names_var(argptr);

  va_end(argptr);

  return ier;
}
