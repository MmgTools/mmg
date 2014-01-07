 /**
 *
 * Written by Cecile Dobrzynski (IMB), Charles Dapogny and Pascal Frey (LJLL)
 * Copyright (c) 2004- IMB/LJLL.
 * All rights reserved.
 *
 * mmg3dlib(int *options_i,double *options_d ): to use mmg3d via a library
 *
 * Integers parameters:
 *    MMG5_IPARAM_verbose           = [-10..10] , Tune level of verbosity;
 *    MMG5_IPARAM_mem               = [n/-1]    , Set memory size to n Mbytes/keep the default value;
 *    MMG5_IPARAM_debug             = [1/0]     , Turn on/off debug mode;
 *    MMG5_IPARAM_angle             = [1/0]     , Turn on/off angle detection;
 *    MMG5_IPARAM_iso               = [1/0]     , Turn on/off levelset meshing;
 *    MMG5_IPARAM_noinsert          = [1/0]     , avoid/allow point insertion/deletion;
 *    MMG5_IPARAM_noswap            = [1/0]     , avoid/allow edge or face flipping;
 *    MMG5_IPARAM_nomove            = [1/0]     , avoid/allow point relocation;
 *    MMG5_IPARAM_numberOflocalParam= [n]       , number of local parameters;
 *    MMG5_IPARAM_renum             = [1/0]     , Turn on/off the renumbering using SCOTCH;
 *    MMG5_IPARAM_sing              = [1/0]     , Turn on/off the insertion of singularities
 *                                        (need to compile with -DSINGUL flag);
 * Doble parameters:
 *    MMG5_DPARAM_dhd   = [val]     , angle detection;
 *    MMG5_DPARAM_hmin  = [val]     , minimal mesh size;
 *    MMG5_DPARAM_hmax  = [val]     , maximal mesh size;
 *    MMG5_DPARAM_hausd = [val]     , control Hausdorff distance;
 *    MMG5_DPARAM_hgrad = [val]     , control gradation;
 *    MMG5_DPARAM_ls    = [val]     , level set value;
 **/

#include "mmg3d.h"

/** Allocate the mesh and sol structures and initialize it to default values */
void Init_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol
#ifdef SINGUL
              , MMG5_pSingul *sing
#endif
              ) {

  /* mesh allocation */
  if ( *mesh ) {
    free(*mesh);
    *mesh = NULL;
  }
  *mesh = (MMG5_pMesh)calloc(1,sizeof(MMG5_Mesh));
  if ( !*mesh ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }

  /* sol allocation */
  if ( *sol ) {
    free(*sol);
    *sol = NULL;
  }
  *sol           = (MMG5_pSol)calloc(1,sizeof(MMG5_Sol));
  if ( !*sol ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }

#ifdef SINGUL
  /* singul allocation */
  if ( *sing ) {
    free(*sing);
    *sing = NULL;
  }
  *sing          = (MMG5_pSingul)calloc(1,sizeof(MMG5_Singul));
  if ( !*sing ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
#endif

  (*mesh)->dim = 3;
  (*mesh)->ver = 2;
  (*sol)->dim  = 3;
  (*sol)->ver  = 2;
  (*sol)->size = 1;

  MMG5_Init_parameters(*mesh);

  /* file name initialisations */
  MMG5_Set_inputMeshName(*mesh,"");
  MMG5_Set_outputMeshName(*mesh,"");

  MMG5_Set_inputSolName(*mesh,*sol,"");
  MMG5_Set_outputSolName(*mesh,*sol,"");

#ifdef SINGUL
  MMG5_Set_inputSingulName(*mesh,*sing,"");
#endif

  return;
}

/** Initialization of parameters stored in the Info structure */
void Init_parameters(pMesh mesh) {

  /* default values for integers */
  mesh->info.imprim   = -99; /**< [-10..10],Tune level of imprim */
  mesh->info.mem      = -1;  /**< [n/-1]   ,Set memory size to n Mbytes/keep the default value */
  mesh->info.ddebug   =  0;  /**< [0/1]    ,Turn on/off debug mode */
  mesh->info.iso      =  0;  /**< [0/1]    ,Turn on/off levelset meshing */
  mesh->info.noinsert =  0;  /**< [0/1]    ,avoid/allow point insertion/deletion */
  mesh->info.noswap   =  0;  /**< [0/1]    ,avoid/allow edge or face flipping */
  mesh->info.nomove   =  0;  /**< [0/1]    ,avoid/allow point relocation */
  mesh->info.npar     =  0;  /**< [n]      ,number of local parameters */
#ifdef USE_SCOTCH
  mesh->info.renum    = 1;   /**< [1/0]    , Turn on/off the renumbering using SCOTCH; */
#else
  mesh->info.renum    = 0;   /**< [1/0]    , Turn on/off the renumbering using SCOTCH; */
#endif
  mesh->info.sing     =  0;  /**< [0/1]    ,preserve internal singularities */

  /* default values for doubles */
  mesh->info.dhd      = ANGEDG;   /**< angle detection; */
  mesh->info.hmin     = 0.0;      /**< minimal mesh size; */
  mesh->info.hmax     = FLT_MAX;  /**< maximal mesh size; */
  mesh->info.hausd    = 0.01;     /**< control Hausdorff */
  mesh->info.hgrad    = 0.1;      /**< control gradation; */
  mesh->info.ls       = 0.0;      /**< level set value */
}

/** Set the name of input mesh */
int Set_inputMeshName(MMG5_pMesh mesh, char* meshin) {

  if ( mesh->namein ) {
    free ( mesh->namein );
    mesh->namein = NULL;
  }
  if ( strlen(meshin) ) {
    mesh->namein = (char *)calloc(strlen(meshin)+1,sizeof(char));
    if ( !mesh->namein ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    strcpy(mesh->namein,meshin);
  }
  else {
    mesh->namein = (char *)calloc(128,sizeof(char));
    if ( !mesh->namein ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    strcpy(mesh->namein,"mesh.mesh");
    if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: no name given for input mesh.\n");
      fprintf(stdout,"     Use of default value \"mesh.mesh\".\n");
    }
  }

  return(1);
}

/** Set the name of input sol */
int Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solin) {
  char *ptr, *ptr2;

  if ( sol->namein ) {
    free ( sol->namein );
    sol->namein = NULL;
  }
  if ( strlen(solin) ) {
    sol->namein = (char *)calloc(strlen(solin)+1,sizeof(char));
    if ( !sol->namein ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    strcpy(sol->namein,solin);
  }
  else {
    if ( strlen(mesh->namein) ) {
      sol->namein = (char *)calloc( strlen(mesh->namein)+1,sizeof(char));
      if ( !sol->namein ) {
        perror("  ## Memory problem: calloc");
        exit(EXIT_FAILURE);
      }
      strcpy(sol->namein,mesh->namein);

      ptr = strstr(sol->namein,".mesh");
      if ( ptr )
        /* the sol file is renamed with the meshfile without extension */
        *ptr = '\0';
    }
    else {
      sol->namein = (char *)calloc(128,sizeof(char));
      if ( !sol->namein ) {
        perror("  ## Memory problem: calloc");
        exit(EXIT_FAILURE);
      }
      strcpy(sol->namein,"mesh.sol");
    }
  }
  return(1);
}

/** Set the name of output mesh */
int Set_outputMeshName(MMG5_pMesh mesh, char* meshout) {
  char *ptr,*ptr2;

  if ( mesh->nameout ) {
    free ( mesh->nameout );
    mesh->nameout = NULL;
  }
  if ( strlen(meshout) ) {
    mesh->nameout = (char *)calloc(strlen(meshout)+1,sizeof(char));
    if ( !mesh->nameout ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    strcpy(mesh->nameout,meshout);
  }
  else {
    if ( strlen(mesh->namein) ) {
      mesh->nameout = (char *)calloc( strlen(mesh->namein)+8,sizeof(char));
      if ( !mesh->nameout ) {
        perror("  ## Memory problem: calloc");
        exit(EXIT_FAILURE);
      }
      strcpy(mesh->nameout,mesh->namein);
      ptr = strstr(mesh->nameout,".mesh");
      if ( !ptr ) {
        /* filename without extension */
        strcat(mesh->nameout,".o.mesh");
      }
      else {
        ptr2 = strstr(mesh->nameout,".meshb");
        if ( ptr2 ) {
          /* filename with .meshb extention */
          *ptr2 = '\0';
          strcat(mesh->nameout,".o.meshb");
        }
        else {
          /* filename with .mesh extention */
          *ptr = '\0';
          strcat(mesh->nameout,".o.mesh");
        }
      }
    }
    else {
      mesh->nameout = (char *)calloc(128,sizeof(char));
      if ( !mesh->nameout ) {
        perror("  ## Memory problem: calloc");
        exit(EXIT_FAILURE);
      }
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
        fprintf(stdout,"  ## Warning: no name given for output mesh.\n");
        fprintf(stdout,"     Use of default value \"mesh.o.mesh\".\n");
      }
      strcpy(mesh->nameout,"mesh.o.mesh");
    }
  }
  return(1);
}

/** Set the name of output sol */
int Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, char* solout) {
  char *ptr;

  if ( sol->nameout ) {
    free ( sol->nameout );
    sol->nameout = NULL;
  }
  if ( strlen(solout) ) {
    sol->nameout = (char *)calloc(strlen(solout)+1,sizeof(char));
    if ( !sol->nameout ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    strcpy(sol->nameout,solout);
  }
  else {
    if ( strlen(mesh->nameout) ) {
      sol->nameout = (char *)calloc( strlen(mesh->nameout)+1,sizeof(char));
      if ( !sol->nameout ) {
        perror("  ## Memory problem: calloc");
        exit(EXIT_FAILURE);
      }
      strcpy(sol->nameout,mesh->nameout);

      ptr = strstr(sol->nameout,".mesh");
      if ( ptr )
        /* the sol file is renamed with the meshfile without extension */
        *ptr = '\0';
    }
    else {
      fprintf(stdout,"  ## Error: no name for output mesh. please, use");
      fprintf(stdout," the MMG5_Set_outputMeshName to set the mesh name.\n");
      return(0);
    }
  }
  return(1);
}

#ifdef SINGUL
/** Set the name of input singularities file */
int Set_inputSingulName(MMG5_pMesh mesh,MMG5_pSingul sing, char* singin) {

  if ( sing->namein ) {
    free ( sing->namein );
    sing->namein = NULL;
  }
  if ( strlen(singin) ) {
    sing->namein = (char *)calloc(strlen(singin)+1,sizeof(char));
    if ( !sing->namein ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    strcpy(sing->namein,singin);
  }
  else {
    sing->namein = (char *)calloc(128,sizeof(char));
    if ( !sing->namein ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    strcpy(sing->namein,"singularities.mesh");
    if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: no name given for input singularities.\n");
      fprintf(stdout,"     Use of default value \"singularities.mesh\".\n");
    }
  }
  return(1);
}
#endif


/** Set the solution number, dimension and type */
int Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, int np, int typSol) {

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) && sol->m )
    fprintf(stdout,"  ## Warning: reset of solution\n");

  if ( typEntity != MMG5_Vertex ) {
    fprintf(stdout,"  ## Error: MMG3D5 need a solution imposed on vertices\n");
    return(0);
  }
  if ( typSol != MMG5_Scalar ) {
    fprintf(stdout,"  ## Error: anisotropic adaptation not yet implemented\n");
    return(0);
  }
  else sol->size = 1;

  sol->dim = 3;
  sol->npmax = mesh->npmax;
  if ( np ) {
    sol->np = np;
    if ( sol->m ) {
      free(sol->m);
      sol->m = NULL;
    }
    sol->m    = (double*)calloc(sol->npmax+1,sol->size*sizeof(double));
    if ( !sol->m ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
  }
  return(1);
}

/** Set the number of vertices, tetrahedra, triangles and edges of the mesh
    and allocate the associated tables. If call twice, reset the whole mesh to
    realloc it at the new size */
int Set_meshSize(MMG5_pMesh mesh, int np, int ne, int nt, int na) {
  int k;

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) &&
       ( mesh->point || mesh->tria || mesh->tetra || mesh->edge) )
    fprintf(stdout,"  ## Warning: reset of mesh\n");

  mesh->np  = np;
  mesh->ne  = ne;
  mesh->nt  = nt;
  mesh->na  = na;
  mesh->npi = 0;
  mesh->nei = 0;
  mesh->nti = 0;
  mesh->nai = 0;

  mesh->npmax = MG_MAX(1.5*mesh->np,NPMAX);
  mesh->ntmax = MG_MAX(1.5*mesh->nt,NTMAX);
  mesh->nemax = MG_MAX(1.5*mesh->ne,NEMAX);
  mesh->namax = mesh->na;

  if ( mesh->point ) {
    free(mesh->point);
    mesh->point = NULL;
  }
  mesh->point = (MMG5_pPoint)calloc(mesh->npmax+1,sizeof(Point));
  if ( !mesh->point ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }

  if ( mesh->tetra ) {
    free(mesh->tetra);
    mesh->tetra = NULL;
  }
  mesh->tetra = (MMG5_pTetra)calloc(mesh->nemax+1,sizeof(Tetra));
  if ( !mesh->tetra ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }

  if ( mesh->tria ) {
    free(mesh->tria);
    mesh->tria = NULL;
  }
  if ( mesh->nt ) {
    mesh->tria  = (MMG5_pTria)calloc(mesh->ntmax+1,sizeof(Tria));
    if ( !mesh->tria ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
  }

  if ( mesh->edge ) {
    free(mesh->edge);
    mesh->edge = NULL;
  }
  if ( mesh->na ) {
    mesh->edge  = (MMG5_pEdge)calloc(mesh->namax+1,sizeof(Edge));
    if ( !mesh->edge ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
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

  return(1);
}

#ifdef SINGUL
/** Set the number of singular vertices and edges and allocate the associated tables. */
int Set_singulSize(MMG5_pSingul sing, int np, int na) {

  if ( sing->point || sing->edge )
    fprintf(stdout,"  ## Warning: reset of singularites\n");

  sing->ns = np;
  sing->na = na;

  if ( sing->ns ) {
    if ( sing->point ) {
      free(sing->point);
      sing->point = NULL;
    }
    sing->point  = (MMG5_psPoint)calloc(sing->ns+1,sizeof(sPoint));
    if ( !sing->point ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
  }

  if ( sing->na ) {
    if ( sing->edge ) {
      free(sing->edge);
      sing->edge = NULL;
    }
    sing->edge  = (MMG5_pEdge)calloc(sing->na+1,sizeof(Edge));
    if ( !sing->edge ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
  }
  return(1);
}
#endif

/** Get the solution number, dimension and type */
int Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, int* np, int* typSol) {
  int k;

  *typEntity = MMG5_Vertex;
  *typSol    = sol->size;

  if ( (sol->np != mesh->np) || (sol->npi==0) ) {
   fprintf(stdout,"  ## Error: You must first call the");
   fprintf(stdout," MMG5_Get_meshSize function before getting the size of solution\n");
   return(0);
  }

  *np = 0;
  for (k=sol->np; k>=1; k--) {
    if ( MG_VOK(&mesh->point[k]) ) {
      if ( mesh->point[mesh->npi].tmp != ++(*np) ) {
        fprintf(stdout,"  ## Error: problem when counting the solutions");
        fprintf(stdout,"  ## Check that you call the MMG5_Get_meshSize");
        fprintf(stdout," function before calling MMG5_Get_solSize.\n");
        return(0);
      }
    }
  }

  assert( ((*np) == sol->np) && (*np) == mesh->np );

  return(1);
}

/** Get the number of vertices, tetrahedra, triangles and edges of the mesh. */
int Get_meshSize(MMG5_pMesh mesh, MMG5_pSol sol, int* np, int* ne, int* nt, int* na) {
  int    k, i;
  pPoint ppt;
  pTetra pt;
  hgeom  *ph;

  for (k=mesh->np; k>=1; k--) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) mesh->point[k].tmp = ++(*np);
  }
  mesh->npi = mesh->np;
  sol->npi  = mesh->np;
  mesh->np  = *np;
  sol->np   = mesh->np;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( MG_EOK(pt) ) (*ne)++;
  }
  mesh->nei = mesh->ne;
  mesh->ne  = *ne;

  mesh->nt = 0;
  if ( mesh->tria ){
    free(mesh->tria);
    mesh->tria=NULL;
  }
  chkNumberOfTri(mesh);
  mesh->nti = mesh->nt;
  *nt       = mesh->nt;

  if ( bdryTria(mesh) ) {
    free(mesh->adjt);
    mesh->adjt=NULL;
    free(mesh->adja);
    mesh->adja = NULL;

    /* build hash table for edges */
    if ( mesh->htab.geom ) {
      free(mesh->htab.geom);
      mesh->htab.geom=NULL;
    }
    /* in the wost case (all edges are marked), we will have around 1 edge per *
     * triangle (we count edges only one time) */
    *na = 0;
    if ( hNew(&mesh->htab,mesh->nt,3*(mesh->nt),0) ) {
      for (k=1; k<=mesh->ne; k++) {
        pt   = &mesh->tetra[k];
        if ( MG_EOK(pt) &&  pt->xt ) {
          for (i=0; i<6; i++) {
            if ( mesh->xtetra[pt->xt].edg[i] ||
                 ( MG_EDG(mesh->xtetra[pt->xt].tag[i] ) ||
                   (mesh->xtetra[pt->xt].tag[i] & MG_REQ) ) )
              hEdge(mesh,pt->v[iare[i][0]],pt->v[iare[i][1]],
                    mesh->xtetra[pt->xt].edg[i],mesh->xtetra[pt->xt].tag[i]);
          }
        }
      }
      /* edges + ridges + required edges */
      for (k=0; k<=mesh->htab.max; k++) {
        ph = &mesh->htab.geom[k];
        if ( !ph->a )  continue;
        (*na)++;
      }
      mesh->na = (*na);
      (*na) = 0;
      if ( mesh->na ) {
        if ( mesh->edge ) {
          free(mesh->edge);
          mesh->edge = NULL;
        }
        mesh->edge = (pEdge)calloc(mesh->na+1,sizeof(Edge));

        for (k=0; k<=mesh->htab.max; k++) {
          ph = &mesh->htab.geom[k];
          if ( !ph->a )  continue;
          (*na)++;
          mesh->edge[*na].a   = mesh->point[ph->a].tmp;
          mesh->edge[*na].b   = mesh->point[ph->b].tmp;
          mesh->edge[*na].ref = ph->ref;
          mesh->edge[*na].tag = ph->tag;
        }
      }
      free(mesh->htab.geom);
      mesh->htab.geom = NULL;
    }
  }
  assert ( (*na) == mesh->na);
  mesh->nai = mesh->na;

  return(1);
}

/** Set vertex of coordinates c0,c1,c2 and reference ref in mesh structure */
int Set_vertex(MMG5_pMesh mesh, double c0, double c1, double c2, int ref) {

  if ( !mesh->np ) {
   fprintf(stdout,"  ## Error: You must set the number of points with the");
   fprintf(stdout," MMG5_Set_meshSize function before setting vertices in mesh\n");
   return(0);
  }

  if ( ++mesh->npi > mesh->npmax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new point.\n");
    fprintf(stdout,"    max number of points: %d\n",mesh->npmax);
    fprintf(stdout,"  ## Check the mesh size or ");
    fprintf(stdout,"increase the allocated memory with the -m option.\n");
    return(0);
  }
  if ( mesh->npi > mesh->np ) {
    if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: to much vertices given. Number of vertices");
      fprintf(stdout," modified: old number %d, new number %d.\n",
              mesh->np,mesh->npi);
    }
    mesh->np = mesh->npi;
  }

  mesh->point[mesh->npi].c[0] = c0;
  mesh->point[mesh->npi].c[1] = c1;
  mesh->point[mesh->npi].c[2] = c2;
  mesh->point[mesh->npi].ref  = ref;
  mesh->point[mesh->npi].tag  = MG_NUL;
  mesh->point[mesh->npi].flag = 0;

  return(1);
}

/** Get coordinates c0,c1,c2 and reference ref of next vertex of mesh  */
int Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, int* ref,
               int* isCorner, int* isRequired) {
  int tmp;

  tmp = 0;

  while ( mesh->npi && !MG_VOK((&mesh->point)[mesh->npi]) )
    {
      mesh->npi--;
    }

  if ( mesh->npi <= 0 ) {
    fprintf(stdout,"  ## Error: unable to get point.\n");
    fprintf(stdout,"     The number of call of MMG5_Get_vertex function");
    fprintf(stdout," can not exceed the number of points: %d\n ",mesh->np);
    return(0);
  }

  if ( mesh->point[mesh->npi].tmp != ++tmp ) {
    if ( tmp==1 ) {
      fprintf(stdout,"  ## Error: problem when getting the");
      fprintf(stdout," %d^th point.\n",tmp);
      fprintf(stdout,"  ## Check that you call the MMG5_Get_meshSize");
      fprintf(stdout," function before calling MMG5_Get_Vertex.\n");
      return(0);
    }
    else
      {
        fprintf(stdout,"  ## Error: problem when getting the");
        fprintf(stdout," %d^th point.\n",tmp);
        return(0);
      }
  }

  *c0  = mesh->point[mesh->npi].c[0];
  *c1  = mesh->point[mesh->npi].c[1];
  *c2  = mesh->point[mesh->npi].c[2];
  *ref = mesh->point[mesh->npi].ref;
  if ( mesh->point[mesh->npi].tag & MG_CRN )
    *isCorner = 1;
  else
    *isCorner = 0;

  if ( mesh->point[mesh->npi].tag & MG_REQ )
    *isRequired = 1;
  else
    *isRequired = 0;

  mesh->npi--;

  return(1);
}

/** Set tetrahedra of vertices v0,v1,v2,v3 and reference ref in mesh structure */
int Set_tetrahedra(MMG5_pMesh mesh, int v0, int v1, int v2, int v3, int ref) {
  pTetra pt;
  double aux;

  if ( !mesh->ne ) {
   fprintf(stdout,"  ## Error: You must set the number of elements with the");
   fprintf(stdout," MMG5_Set_meshSize function before setting elements in mesh\n");
   return(0);
  }
  if ( mesh->npi < mesh->np ) {
   fprintf(stdout,"  ## Error: You must set the vertices before the tetrahedras.\n");
   return(0);
  }
  if ( ++mesh->nei > mesh->nemax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new element.\n");
    fprintf(stdout,"    max number of element: %d\n",mesh->nemax);
    fprintf(stdout,"  ## Check the mesh size or ");
    fprintf(stdout,"increase the allocated memory with the -m option.\n");
    return(0);
  }
  if ( mesh->nei > mesh->ne ) {
    if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: to much tetra given. Number of tetras");
      fprintf(stdout," modified: old number %d, new number %d.\n",
              mesh->ne,mesh->nei);
    }
    mesh->ne = mesh->nei;
  }

  pt = &mesh->tetra[mesh->nei];
  pt->v[0] = v0;
  pt->v[1] = v1;
  pt->v[2] = v2;
  pt->v[3] = v3;
  pt->ref  = ref;

  mesh->point[pt->v[0]].tag &= ~MG_NUL;
  mesh->point[pt->v[1]].tag &= ~MG_NUL;
  mesh->point[pt->v[2]].tag &= ~MG_NUL;
  mesh->point[pt->v[3]].tag &= ~MG_NUL;

  /* Possibly switch 2 vertices number so that each tet is positively oriented */
  if ( orvol(mesh->point,pt->v) < 0.0 ) {
    aux = pt->v[2];
    pt->v[2] = pt->v[3];
    pt->v[3] = aux;
  }

  return(1);
}

/** Get vertices v0,v1,v2,v3 and reference ref of next tetra of mesh  */
int Get_tetrahedra(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* v3,
                   int* ref, int* isRequired) {

  if ( mesh->npi ) {
    fprintf(stdout,"  ## Error: You must get all the vertices");
    fprintf(stdout," using the MMG5_Get_Vertex function");
    fprintf(stdout," before getting tetrahedra.\n");
    return(0);
  }

  while ( mesh->nei && !MG_EOK((&mesh->tetra)[mesh->nei]) )
    {
      mesh->npi--;
    }

  if ( mesh->nei <= 0 ) {
    fprintf(stdout,"  ## Error: unable to get tetra.\n");
    fprintf(stdout,"    The number of call of MMG5_Get_tetrahedra function");
    fprintf(stdout," can not exceed the number of tetra: %d\n ",mesh->ne);
    return(0);
  }

  *v0  = mesh->point[mesh->tetra[mesh->nei].v[0]].tmp;
  *v1  = mesh->point[mesh->tetra[mesh->nei].v[1]].tmp;
  *v2  = mesh->point[mesh->tetra[mesh->nei].v[2]].tmp;
  *v3  = mesh->point[mesh->tetra[mesh->nei].v[3]].tmp;
  *ref = mesh->tetra[mesh->nei].ref;

  if ( mesh->tetra[mesh->nei].tag & MG_REQ )
    *isRequired = 1;
  else
    *isRequired = 0;

  mesh->nei--;

  return(1);
}

/** Set triangle of vertices v0,v1,v2 and reference ref in mesh structure */
int Set_triangle(MMG5_pMesh mesh, int v0, int v1, int v2, int ref) {

  if ( !mesh->nt ) {
   fprintf(stdout,"  ## Error: You must set the number of triangles with the");
   fprintf(stdout," MMG5_Set_meshSize function before setting triangles in mesh\n");
   return(0);
  }
  if ( mesh->npi < mesh->np ) {
   fprintf(stdout,"  ## Error: You must set the vertices before the triangles.\n");
   return(0);
  }
  if ( ++mesh->nti > mesh->ntmax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new triangle.\n");
    fprintf(stdout,"    max number of triangle: %d\n",mesh->ntmax);
    fprintf(stdout,"  ## Check the mesh size or ");
    fprintf(stdout,"increase the allocated memory with the -m option.\n");
    return(0);
  }
  if ( mesh->nti > mesh->nt ) {
    if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: to much triangles given. Number of triangles");
      fprintf(stdout," modified: old number %d, new number %d.\n",
              mesh->nt,mesh->nti);
    }
    mesh->nt = mesh->nti;
  }

  mesh->tria[mesh->nti].v[0] = v0;
  mesh->tria[mesh->nti].v[1] = v1;
  mesh->tria[mesh->nti].v[2] = v2;
  mesh->tria[mesh->nti].ref  = ref;

  return(1);
}

/** Get vertices v0,v1,v2 and reference ref of next triangle of mesh  */
int Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref
                 ,int* isRequired) {
  pTria  ptt;

  if ( mesh->npi ) {
    fprintf(stdout,"  ## Error: You must get all the vertices");
    fprintf(stdout," using the MMG5_Get_Vertex function");
    fprintf(stdout," before getting triangles.\n");
    return(0);
  }

  if ( mesh->nti <= 0 ) {
    fprintf(stdout,"  ## Error: unable to get triangle.\n");
    fprintf(stdout,"    The number of call of MMG5_Get_triangle function");
    fprintf(stdout," can not exceed the number of triangles: %d\n ",mesh->nt);
    return(0);
  }

  ptt = &mesh->tria[mesh->nti];
  *v0  = mesh->point[ptt->v[0]].tmp;
  *v1  = mesh->point[ptt->v[1]].tmp;
  *v2  = mesh->point[ptt->v[2]].tmp;
  *ref = ptt->ref;

  if ( (ptt->tag[0] & MG_REQ) && (ptt->tag[1] & MG_REQ) &&
       (ptt->tag[2] & MG_REQ) )
    *isRequired = 1;
  else
    *isRequired = 0;

  mesh->nti--;

  return(1);
}

/** Set edges of extremities v0,v1 and reference ref in mesh structure */
int Set_edges(MMG5_pMesh mesh, int v0, int v1, int ref) {

  if ( !mesh->na ) {
   fprintf(stdout,"  ## Error: You must set the number of edges with the");
   fprintf(stdout," MMG5_Set_meshSize function before setting edges in mesh\n");
   return(0);
  }
  if ( ++mesh->nai > mesh->namax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new edge.\n");
    fprintf(stdout,"    max number of edge: %d\n",mesh->namax);
    fprintf(stdout,"  ## Check the mesh size or ");
    fprintf(stdout,"increase the allocated memory with the -m option.\n");
    return(0);
  }
  if ( mesh->nai > mesh->na ) {
    if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: to much edges given. Number of edges");
      fprintf(stdout," modified: old number %d, new number %d.\n",
              mesh->na,mesh->nai);
    }
    mesh->na = mesh->nai;
  }

  mesh->edge[mesh->nai].a = v0;
  mesh->edge[mesh->nai].b = v1;
  mesh->edge[mesh->nai].ref  = ref;

  return(1);
}

/** Get extremities e0,e1 and reference ref of next edge of mesh  */
int Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref
             ,int* isRidge, int* isRequired) {

 if ( mesh->npi ) {
    fprintf(stdout,"  ## Error: You must get all the vertices");
    fprintf(stdout," using the MMG5_Get_Vertex function");
    fprintf(stdout," before getting edges.\n");
    return(0);
  }

  if ( mesh->nai <= 0 ) {
    fprintf(stdout,"  ## Error: unable to get edge.\n");
    fprintf(stdout,"    The number of call of MMG5_Get_edge function");
    fprintf(stdout," can not exceed the number of edges: %d\n ",mesh->na);
    return(0);
  }

  *e0  = mesh->point[mesh->edge[mesh->nai].a].tmp;
  *e1  = mesh->point[mesh->edge[mesh->nai].b].tmp;
  *ref = mesh->edge[mesh->nai].ref;

  if ( mesh->edge[mesh->nai].tag & MG_GEO )
    *isRidge = 1;
  else
    *isRidge = 0;

  if ( mesh->edge[mesh->nai].tag & MG_REQ )
    *isRequired = 1;
  else
    *isRequired = 0;

  mesh->nai--;

  return(1);
}

/** Set corner at point k  */
int Set_corner(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag |= MG_CRN;
  return(1);
}

/** Set point k as required  */
int Set_requiredVertex(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag |= MG_REQ;
  return(1);
}

/** Set element k as required  */
int Set_requiredTetrahedra(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->ne );
  mesh->tetra[k].tag |= MG_REQ;
  return(1);
}

/** Set triangle k as required  */
int Set_requiredTriangle(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->nt );
  mesh->tria[k].tag[0] |= MG_REQ;
  mesh->tria[k].tag[1] |= MG_REQ;
  mesh->tria[k].tag[2] |= MG_REQ;
  return(1);
}

/** Set ridge at edge k  */
int Set_ridge(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->na );
  mesh->edge[k].tag |= MG_GEO;
  return(1);
}

/** Set edge k as required  */
int Set_requiredEdge(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->na );
  mesh->edge[k].tag |= MG_REQ;
  return(1);
}

/** Set scalar value s in solution structure */
int Set_scalarSol(MMG5_pSol met, double s) {

  if ( !met->np ) {
   fprintf(stdout,"  ## Error: You must set the number of solution with the");
   fprintf(stdout," MMG5_Set_solSize function before setting values");
   fprintf(stdout," in solution structure \n");
   return(0);
  }
  if ( met->npi >= met->npmax ) {
    fprintf(stdout,"  ## Error: unable to set a new solution.\n");
    fprintf(stdout,"    max number of solutions: %d\n",met->npmax);
    return(0);
  }

  met->npi++;
  met->m[met->npi] = s;
  return(1);
}

/** Get solution s of next vertex of mesh  */
int Get_scalarSol(MMG5_pMesh mesh, MMG5_pSol met, double* s) {
  int tmp;

  tmp = 0;

  while ( met->npi && !MG_VOK((&mesh->point)[met->npi]) )
    {
      met->npi--;
    }

  if ( met->npi <= 0 ) {
    fprintf(stdout,"  ## Error: unable to get solution.\n");
    fprintf(stdout,"     The number of call of MMG5_Get_scalarSol function");
    fprintf(stdout," can not exceed the number of points: %d\n ",mesh->np);
    return(0);
  }

  if ( mesh->point[met->npi].tmp != ++tmp ) {
    fprintf(stdout,"  ## Error: problem when getting the");
    fprintf(stdout," %d^th solution.\n",tmp);
    fprintf(stdout,"  ## Check that you call the MMG5_Get_meshSize");
    fprintf(stdout," function before calling MMG5_Get_scalarSol.\n");
    return(0);
  }

  *s  = met->m[met->npi];

  met->npi--;

  return(1);
}

#ifdef SINGUL
/** Set singular point of coordinates c0,c1,c2 in the singularities structure */
int Set_singulVertex(MMG5_pSingul sing, double c0, double c1, double c2, int typ) {

  if ( !sing->ns ) {
    fprintf(stdout,"  ## Error: You must set the number of singular vertex with the");
    fprintf(stdout," MMG5_Set_singulSize function before setting values");
    fprintf(stdout," in singularities structure. \n");
    return(0);
  }
  if ( sing->nsi >= sing->ns ) {
    fprintf(stdout,"  ## Error: unable to set a new singularity.\n");
    fprintf(stdout,"    max number of singular vertices: %d\n",sing->ns);
    return(0);
  }

  sing->nsi++;
  sing->point[sing->nsi].c[0] = c0;
  sing->point[sing->nsi].c[1] = c1;
  sing->point[sing->nsi].c[2] = c2;
  return(1);
}


/** Set singular edge of extremities v0,v1 and reference ref
    in the singularities structure */
int Set_singulEdge(MMG5_pSingul sing, int v0, int v1, int ref) {

  if ( !sing->na ) {
    fprintf(stdout,"  ## Error: You must set the number of singular edges with the");
    fprintf(stdout," MMG5_Set_singulSize function before setting values");
    fprintf(stdout," in singularities structure \n");
    return(0);
  }
  if ( sing->nai >= sing->na ) {
    fprintf(stdout,"  ## Error: unable to set a new singularity.\n");
    fprintf(stdout,"    max number of singular edges: %d\n",sing->na);
    return(0);
  }
  if ( (v0 > sing->ns)||(v1 > sing->ns) ) {
    fprintf(stdout,"  ## Error: edge extremity overflow number ov vertices.\n");
    fprintf(stdout,"    max number of singular vertices: %d\n",sing->ns);
    fprintf(stdout,"  ## Tou must insert all singular vertices before edges\n");
    return(0);
  }

  sing->nai++;
  sing->edge[sing->nai].a   = v0;
  sing->edge[sing->nai].b   = v1;
  sing->edge[sing->nai].ref = ref;
  return(1);
}

/** Set corner at singular vertex k  */
int Set_singulCorner(MMG5_pSingul sing, int k) {
  assert ( k <= sing->ns );
  sing->point[k].tag |= MG_CRN;
  return(1);
}

/** Set required vertex at singular vertex k  */
int Set_singulRequiredVertex(MMG5_pSingul sing, int k) {
  assert ( k <= sing->ns );
  sing->point[k].tag |= MG_REQ;
  return(1);
}

/** Set required edge at singular edge k  */
int Set_singulRidge(MMG5_pSingul sing, int k) {
  assert ( k <= sing->na );
  sing->edge[k].tag |= MG_GEO;
  return(1);
}

/** Set required edge at singular edge k  */
int Set_singulRequiredEdge(MMG5_pSingul sing, int k) {
  assert ( k <= sing->na );
  sing->edge[k].tag |= MG_REQ;
  return(1);
}
#endif


/** Check if the number of given entities match with mesh and sol size (not mandatory)
    and check mesh datas  */
int Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met
#ifdef SINGUL
                 ,MMG5_pSingul sing
#endif
                 ) {

  if ( (mesh->npi != mesh->np) || (mesh->nei != mesh->ne) ||
       (mesh->nti != mesh->nt) || (mesh->nai != mesh->na) ) {
    fprintf(stdout,"  ## Error: problem in the mesh structure\n");
    return(0);
  }
  if ( met->npi != met->np ) {
    fprintf(stdout,"  ## Error: problem in the solution structure\n");
    return(0);
  }

#ifdef SINGUL
  if ( (sing->nsi != sing->ns) || (sing->nai != sing->na) ) {
    fprintf(stdout,"  ## Error: problem in the singularities structure\n");
    return(0);
  }
#endif

  /*  Check mesh data */
  if ( mesh->info.ddebug ) {
    if ( (!mesh->np) || (!mesh->point) ||
         (!mesh->ne) || (!mesh->tetra) ) {
      fprintf(stdout,"  ** MISSING DATA. Exit program.\n");
      return(0);
    }
  }

  if ( mesh->dim != 3 ) {
    fprintf(stdout,"  ** 3 DIMENSIONAL MESH NEEDED. Exit program.\n");
    return(0);
  }
  if ( met->dim != 3 ) {
    fprintf(stdout,"  ** WRONG DIMENSION FOR METRIC. Exit program.\n");
    return(0);
  }
  if ( !mesh->ver )  mesh->ver = 2;
  if ( !met ->ver )  met ->ver = 2;

  return(1);
}

/** realloc memory to satisfy the info.mem option */
static inline
int reallocMemory(MMG5_pMesh mesh) {
  int    k,old_ntmax, old_npmax, old_nemax;
  void   *ptmp;
  pPoint point;
  pTetra tetra;
  pTria  tria;

  old_npmax = mesh->npmax;
  old_ntmax = mesh->ntmax;
  old_nemax = mesh->nemax;

  memRepartition(mesh);

  if ( mesh->npmax != old_npmax ) {
    point = (pPoint)realloc(mesh->point,(mesh->npmax+1)*sizeof(Point));
    if ( !point ){
      perror("  ## Memory problem: realloc");
      free(mesh->point);
      exit(EXIT_FAILURE);
    }
    mesh->point = point;

    /* keep track of empty links */
    if ( old_npmax < mesh->npmax ) {
      for (k=old_npmax; k<mesh->npmax-1; k++)
        mesh->point[k].tmp  = k+1;
    }
    else
      mesh->point[mesh->npmax].tmp  = 0;
  }

  if ( mesh->nemax != old_nemax ) {
    tetra = (pTetra)realloc(mesh->tetra,(mesh->nemax+1)*sizeof(Tetra));
    if ( !tetra ){
      free(mesh->tetra);
      perror("  ## Memory problem: realloc");
      exit(EXIT_FAILURE);
    }
    mesh->tetra = tetra;

    /* keep track of empty links */
    if ( old_nemax < mesh->nemax ) {
      for (k=old_nemax; k<mesh->nemax-1; k++)
        mesh->tetra[k].v[3] = k+1;
    }
    else
      mesh->tetra[mesh->nemax].v[3] = 0;
  }

  if ( mesh->ntmax != old_ntmax ) {
    if ( mesh->nt ) {
      tria = (pTria)realloc(mesh->tria,(mesh->ntmax+1)*sizeof(Tria));
      if ( !mesh->tria ){
        free(mesh->tria);
        perror("  ## Memory problem: realloc");
        exit(EXIT_FAILURE);
      }
      mesh->tria = tria;
    }
  }
  return(1);
}

/** skip the MG_ISO references in an input mesh */
static inline
int skipIso(MMG5_pMesh mesh) {
  pTria  ptt,ptt1;
  pEdge  pa,pa1;
  int    k;

  if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
    fprintf(stdout,"  ## Warning: skip of all entites with %d reference.\n",MG_ISO);

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
      memcpy(ptt,ptt1,sizeof(Tria));
  } while( ++k <= mesh->nti );

  if ( mesh->nti < mesh->nt ) {
    mesh->nt = mesh->nti;
    if( !mesh->nt ){
      free(mesh->tria);
      mesh->tria=NULL;
    }
    else {
      ptt = (pTria)realloc(mesh->tria,(mesh->nt+1)*sizeof(Tria));
      if ( !ptt ) {
        free(mesh->tria);
        perror("  ## Memory problem: realloc");
        exit(EXIT_FAILURE);
      }
      mesh->tria = ptt;
    }
  }

  /* Skip edges with MG_ISO refs */
  k = 1;
  do {
    pa = &mesh->edge[k];
    if ( abs(pa->ref) != MG_ISO ) {
      pa->ref = abs(pa->ref);
      continue;
    }
    /* here pa is the first tri of mesh->edge that we want to delete */
    do {
      pa1 = &mesh->edge[mesh->nai];
    }
    while( (abs(pa1->ref) == MG_ISO) && (k <= --mesh->nai) );

    if ( abs(pa1->ref) != MG_ISO ) {
      /* pa1 is the last tri of mesh->edge that we want to keep */
      memcpy(pa,pa1,sizeof(Edge));
      pa1->ref = abs(pa1->ref);
    }
  } while( ++k <= mesh->nai );

  if ( mesh->nai < mesh->na ) {
    mesh->na = mesh->nai;
    if( !mesh->na ){
      free(mesh->edge);
      mesh->edge=NULL;
    }
    else {
      pa = (pEdge)realloc(mesh->edge,(mesh->na+1)*sizeof(Edge));
      if ( !pa ) {
        free(mesh->edge);
        perror("  ## Memory problem: realloc");
        exit(EXIT_FAILURE);
      }
      mesh->edge = pa;
    }
  }

  /* delete tetrahedra references */
  for (k=1; k<=mesh->ne; k++)  mesh->tetra[k].ref = 0;
  return(1);
}

/** Set integer parameter iparam at value val */
int Set_iparameters(MMG5_pMesh mesh, int iparam, int val){
  int k;

  switch ( iparam ) {
    /* Integer parameters */
  case MMG5_IPARAM_verbose :
    mesh->info.imprim   = val;
    break;
  case MMG5_IPARAM_mem :
    mesh->info.mem      = val;
    if ( !reallocMemory(mesh) )  exit(EXIT_FAILURE);
    break;
  case MMG5_IPARAM_debug :
    mesh->info.ddebug   = val;
    break;
  case MMG5_IPARAM_angle :
    if ( !val )
      mesh->info.dhd    = -1.;
    else {
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stdout,"  ## Warning: angle detection parameter set to default value\n");
      mesh->info.dhd    = ANGEDG;
    }
    break;
  case MMG5_IPARAM_iso :
    mesh->info.iso      = val;
    if ( mesh->info.iso )
      if ( !skipIso(mesh) )
        exit(EXIT_FAILURE);
    break;
  case MMG5_IPARAM_noinsert :
    mesh->info.noinsert = val;
    break;
  case MMG5_IPARAM_noswap :
    mesh->info.noswap   = val;
    break;
  case MMG5_IPARAM_nomove :
    mesh->info.nomove   = val;
    break;
  case MMG5_IPARAM_numberOfLocalParam :
    mesh->info.npar     = val;
    mesh->info.npari = 0;
    if ( mesh->info.par ) {
      free(mesh->info.par);
      mesh->info.par   = NULL;
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stdout,"  ## Warning: reset of local parameter values\n");
    }
    mesh->info.par = (pPar)calloc(mesh->info.npar,sizeof(Par));
    if ( !mesh->info.par ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    for (k=0; k<mesh->info.npar; k++) {
      mesh->info.par[k].elt   = MMG5_Noentity;
      mesh->info.par[k].ref   = INT_MAX;
      mesh->info.par[k].hausd = mesh->info.hausd;
    }

    break;
#ifdef USE_SCOTCH
  case MMG5_IPARAM_renum :
    mesh->info.renum    = val;
    break;
#endif
#ifdef SINGUL
  case MMG5_IPARAM_sing :
    mesh->info.sing     = val;
    break;
#endif
  default :
    fprintf(stdout,"  ## Error: unknown type of parameter\n");
    return(0);
  }
  /* other options */
  mesh->info.fem      = 0;
  return(1);
}

/** Set double parameter dparam at value val */
int Set_dparameters(MMG5_pMesh mesh, int dparam, double val){

  switch ( dparam ) {
    /* double parameters */
  case MMG5_DPARAM_angleDetection :
    mesh->info.dhd = val;
    mesh->info.dhd = MG_MAX(0.0, MG_MIN(180.0,mesh->info.dhd));
    mesh->info.dhd = cos(mesh->info.dhd*M_PI/180.0);
    break;
  case MMG5_DPARAM_hmin :
    mesh->info.hmin     = val;
    break;
  case MMG5_DPARAM_hmax :
    mesh->info.hmax     = val;
    break;
  case MMG5_DPARAM_hgrad :
    mesh->info.hgrad    = val;
    if ( mesh->info.hgrad < 0.0 )
      mesh->info.hgrad = -1.0;
    else
      mesh->info.hgrad = log(mesh->info.hgrad);
    break;
  case MMG5_DPARAM_hausd :
    mesh->info.hausd    = val;
    break;
  case MMG5_DPARAM_ls :
    mesh->info.ls       = val;
    break;
  default :
    fprintf(stdout,"  ## Error: unknown type of parameter\n");
    return(0);
  }
  return(1);
}

/** Set local parameters: set the hausdorff value at val for all elements
    of type typ and reference ref  */
int Set_localParameters(MMG5_pMesh mesh, int typ, int ref, double val){
  int k;

  if ( !mesh->info.npar ) {
    fprintf(stdout,"  ## Error: You must set the number of local parameters");
    fprintf(stdout," with the MMG5_Set_iparameters function before setting");
    fprintf(stdout," values in local parameters structure. \n");
    return(0);
  }
  if ( mesh->info.npari > mesh->info.npar ) {
    fprintf(stdout,"  ## Error: unable to set a new local parameter.\n");
    fprintf(stdout,"    max number of local parameters: %d\n",mesh->info.npar);
    return(0);
  }

  switch ( typ ) {
    /* double parameters */
  case MMG5_Triangle :
    for (k=0; k<mesh->info.npari; k++) {
      if ( mesh->info.par[k].ref == ref ) {
        mesh->info.par[k].hausd = val;
        if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
          fprintf(stdout,"  ## Warning: new hausdorff value for triangles");
          fprintf(stdout," of ref %d\n",ref);
        }
        return(1);
      }
    }
    if ( mesh->info.npari == mesh->info.npar ) {
      fprintf(stdout,"  ## Error: unable to set a new local parameter.\n");
      fprintf(stdout,"    max number of local parameters: %d\n",mesh->info.npar);
      return(0);
    }
    mesh->info.par[mesh->info.npari].elt   = typ;
    mesh->info.par[mesh->info.npari].ref   = ref;
    mesh->info.par[mesh->info.npari].hausd = val;
    mesh->info.npari++;
    break;
  default :
    printf("plop %d\n",typ);
    fprintf(stdout,"  ## Warning: you must apply local hausdorff number");
    fprintf(stdout," on triangles (MMG5_Triangle or %d).\n",MMG5_Triangle);
    fprintf(stdout,"  ## Ignored.\n");
    return(1);
  }

  return(1);
}


/** File names deallocations before return */
void Free_names(pMesh mesh,pSol met
#ifdef SINGUL
               ,pSingul singul
#endif
               ){
  /* mesh */
  free(mesh->nameout);
  mesh->nameout = NULL;
  free(mesh->namein);
  mesh->namein = NULL;

  /* met */
  if ( met->namein ) {
    free(met->namein);
    met->namein = NULL;
  }
  if ( met->nameout ) {
    free(met->nameout);
    met->nameout = NULL;
  }

#ifdef SINGUL
  /* singul */
  if ( mesh->info.sing ) {
    if ( singul->namein ) {
      free(singul->namein);
      singul->namein=NULL;
    }
   }
#endif
}

/** Structure deallocations before return */
void Free_structures(pMesh mesh,pSol met
#ifdef SINGUL
             ,pSingul singul
#endif
             ){

#ifdef SINGUL
  Free_names(mesh,met,singul);
#else
  Free_names(mesh,met);
#endif

  /* mesh */
  free(mesh->point);
  mesh->point = NULL;
  free(mesh->tetra);
  mesh->tetra = NULL;
  free(mesh->adja);
  mesh->adja = NULL;
  if ( mesh->xpoint ) {
    free(mesh->xpoint);
    mesh->xpoint = NULL;
  }
  if ( mesh->htab.geom ) {
    free(mesh->htab.geom);
    mesh->htab.geom = NULL;
  }

  if ( mesh->tria ) {
    free(mesh->tria);
    mesh->tria = NULL;
  }
  if ( mesh->xtetra ) {
    free(mesh->xtetra);
    mesh->xtetra = NULL;
  }

  /* met */
  if ( /*!mesh->info.iso &&*/ met->m ) {
    free(met->m);
    met->m = NULL;
  }

  /* mesh->info */
  if ( mesh->info.npar && mesh->info.par ) {
    free(mesh->info.par);
    mesh->info.par = NULL;
  }

#ifdef SINGUL
  /* singul */
  if ( singul->ns ) {
    free(singul->point);
    singul->point=NULL;
  }
  if ( singul->na ) {
    free(singul->edge);
    singul->edge=NULL;
  }
#endif
}
