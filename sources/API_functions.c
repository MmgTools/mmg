 /**
 *
 * Written by Cecile Dobrzynski (IMB), Charles Dapogny,
 * Pascal Frey (LJLL) and Algiane Froehly
 * Copyright (c) 2004- IMB/LJLL.
 * All rights reserved.
 *
 * Integers parameters:
 *    MMG5_IPARAM_verbose           = [-10..10] , Tune level of verbosity;
 *    MMG5_IPARAM_mem               = [n/-1]    , Set maximal memory size to n Mbytes/keep the default value;
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
 *    MMG5_IPARAM_bucket            = [val]     , Specify the size of the bucket per dimension (Delaunay)
 *                                        (need to compile with PATTERN=NO);
 * Double parameters:
 *    MMG5_DPARAM_angleDetection   = [val]     , angle detection;
 *    MMG5_DPARAM_hmin             = [val]     , minimal mesh size;
 *    MMG5_DPARAM_hmax             = [val]     , maximal mesh size;
 *    MMG5_DPARAM_hausd = [val]     , control global Hausdorff distance
 *                                    (on all the boundary surfaces of the mesh);
 *    MMG5_DPARAM_hgrad            = [val]     , control gradation;
 *    MMG5_DPARAM_ls               = [val]     , level set value;
 **/

#include "mmg3d.h"


/** Allocate the mesh and sol structures */
static inline
void Alloc_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol
#ifdef SINGUL
    , MMG5_pSingul *sing
#endif
    ) {

  /* mesh allocation */
  if ( *mesh )  SAFE_FREE(*mesh);
  SAFE_CALLOC(*mesh,1,MMG5_Mesh);

  /* sol allocation */
  if ( *sol )  DEL_MEM(*mesh,*sol,sizeof(MMG5_Sol));
  SAFE_CALLOC(*sol,1,MMG5_Sol);

#ifdef SINGUL
  /* singul allocation */
  if ( *sing )  DEL_MEM(*mesh,*sing,sizeof(MMG5_Singul));
  SAFE_CALLOC(*sing,1,MMG5_Singul);

#endif
  return;
}

/** Initialization of mesh and sol structures */
static inline
void Init_woalloc_mesh(MMG5_pMesh mesh, MMG5_pSol sol
#ifdef SINGUL
              , MMG5_pSingul sing
#endif
              ) {

  (mesh)->dim = 3;
  (mesh)->ver = 2;
  (sol)->dim  = 3;
  (sol)->ver  = 2;
  (sol)->size = 1;

  /* Default parameters values */
  MMG5_Init_parameters(mesh);

  /* Default vaules for file names */
#ifndef SINGUL
  MMG5_Init_fileNames(mesh,sol);
#else
  MMG5_Init_fileNames(mesh,sol,sing);
#endif

  return;
}

/** Allocate the mesh and sol structures and initialize it to default values */
void Init_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol
#ifdef SINGUL
              , MMG5_pSingul *sing
#endif
              ) {

#ifndef SINGUL
  /* allocations */
  Alloc_mesh(mesh,sol);
  /* initialisations */
  Init_woalloc_mesh(*mesh,*sol);
#else
  Alloc_mesh(mesh,sol,sing);
  /* initialisations */
  Init_woalloc_mesh(*mesh,*sol,*sing);
#endif

  return;
}

/** Initialization of parameters stored in the Info structure */
void Init_parameters(pMesh mesh) {

  /* default values for integers */
  mesh->info.imprim   =  4; /**< [-10..10],Tune level of imprim */
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

  /* initial value for memMax and gap */
  mesh->gap = 0.2;
  mesh->memMax = memSize();
  if ( mesh->memMax )
  /* maximal memory = 50% of total physical memory */
    mesh->memMax = mesh->memMax*50/100;
  else {
    /* default value = 800 Mo */
    printf("  Maximum memory set to default value: %d Mo.\n",MEMMAX);
    mesh->memMax = MEMMAX << 20;
  }

#ifndef PATTERN
  mesh->info.bucket = 64;
#endif
}

/** default values for file names */
void Init_fileNames(pMesh mesh,pSol sol
#ifdef SINGUL
                    ,pSingul sing
#endif
                    ) {
  MMG5_Set_inputMeshName(mesh,"");
  MMG5_Set_outputMeshName(mesh,"");

  MMG5_Set_inputSolName(mesh,sol,"");
  MMG5_Set_outputSolName(mesh,sol,"");

#ifdef SINGUL
  MMG5_Set_inputSingulName(mesh,sing,"");
#endif
  return;
}


/** Set the name of input mesh */
int Set_inputMeshName(MMG5_pMesh mesh, char* meshin) {

  if ( mesh->namein )
    DEL_MEM(mesh,mesh->namein,(strlen(mesh->namein)+1)*sizeof(char));

  if ( strlen(meshin) ) {
    ADD_MEM(mesh,(strlen(meshin)+1)*sizeof(char),"input mesh name",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(mesh->namein,strlen(meshin)+1,char);
    strcpy(mesh->namein,meshin);
  }
  else {
    ADD_MEM(mesh,10*sizeof(char),"input mesh name",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(mesh->namein,10,char);
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
  char *ptr;

  if ( sol->namein )
    DEL_MEM(mesh,sol->namein,(strlen(sol->namein)+1)*sizeof(char));
  if ( strlen(solin) ) {
    ADD_MEM(mesh,(strlen(solin)+1)*sizeof(char),"input sol name",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(sol->namein,strlen(solin)+1,char);
    strcpy(sol->namein,solin);
  }
  else {
    if ( strlen(mesh->namein) ) {
      SAFE_CALLOC(sol->namein,strlen(mesh->namein)+1,char);
      strcpy(sol->namein,mesh->namein);
      ptr = strstr(sol->namein,".mesh");
      if ( ptr ) {
        /* the sol file is renamed with the meshfile without extension */
        *ptr = '\0';
        SAFE_REALLOC(sol->namein,(strlen(sol->namein)+1),char,"input sol name");
      }
      ADD_MEM(mesh,(strlen(sol->namein)+1)*sizeof(char),"input sol name",
              printf("  Exit program.\n");
              exit(EXIT_FAILURE));
    }
    else {
      ADD_MEM(mesh,9*sizeof(char),"input sol name",
              printf("  Exit program.\n");
              exit(EXIT_FAILURE));
      SAFE_CALLOC(sol->namein,9,char);
      strcpy(sol->namein,"mesh.sol");
    }
  }
  return(1);
}

/** Set the name of output mesh */
int Set_outputMeshName(MMG5_pMesh mesh, char* meshout) {
  char *ptr;

  if ( mesh->nameout )
    DEL_MEM(mesh,mesh->nameout,(strlen(mesh->nameout)+1)*sizeof(char));

  if ( strlen(meshout) ) {
    ADD_MEM(mesh,(strlen(meshout)+1)*sizeof(char),"output mesh name",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(mesh->nameout,strlen(meshout)+1,char);
    strcpy(mesh->nameout,meshout);
  }
  else {
    if ( strlen(mesh->namein) ) {
      ADD_MEM(mesh,(strlen(mesh->namein)+3)*sizeof(char),"output mesh name",
              printf("  Exit program.\n");
              exit(EXIT_FAILURE));
      SAFE_CALLOC(mesh->nameout,strlen(mesh->namein)+3,char);
      strcpy(mesh->nameout,mesh->namein);
      ptr = strstr(mesh->nameout,".mesh");
      if ( !ptr ) {
        /* filename without extension */
        strcat(mesh->nameout,".o");
      }
      else {
        *ptr = '\0';
        strcat(mesh->nameout,".o.mesh");
      }
      ptr = strstr(mesh->namein,".meshb");
      if ( ptr ) {
        /* filename with .meshb extention */
        strcat(mesh->nameout,"b");
      }
    }
    else {
      ADD_MEM(mesh,7*sizeof(char),"output mesh name",
              printf("  Exit program.\n");
              exit(EXIT_FAILURE));
      SAFE_CALLOC(mesh->nameout,7,char);
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

  if ( sol->nameout )
    DEL_MEM(mesh,sol->nameout,(strlen(sol->nameout)+1)*sizeof(char));

  if ( strlen(solout) ) {
    ADD_MEM(mesh,(strlen(solout)+1)*sizeof(char),"output sol name",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(sol->nameout,strlen(solout)+1,char);
    strcpy(sol->nameout,solout);
  }
  else {
    if ( strlen(mesh->nameout) ) {
      ptr = strstr(mesh->nameout,".mesh");
      if ( ptr )
        SAFE_CALLOC(sol->nameout,strlen(mesh->nameout)+1,char);
      else
        SAFE_CALLOC(sol->nameout,strlen(mesh->nameout)+5,char);
      strcpy(sol->nameout,mesh->nameout);
      ptr = strstr(sol->nameout,".mesh");
      if ( ptr )
        /* the sol file is renamed with the meshfile without extension */
        *ptr = '\0';
      strcat(sol->nameout,".sol");
      ADD_MEM(mesh,(strlen(sol->nameout)+1)*sizeof(char),"output sol name",
              printf("  Exit program.\n");
              exit(EXIT_FAILURE));
      SAFE_REALLOC(sol->nameout,(strlen(sol->nameout)+1),char,"output sol name");
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

  if ( sing->namein )
    DEL_MEM(mesh,sing->namein,(strlen(sing->namein)+1)*sizeof(char));

  if ( strlen(singin) ) {
    ADD_MEM(mesh,(strlen(singin)+1)*sizeof(char),
            "input singularities file name",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(sing->namein,strlen(singin)+1,char);
    strcpy(sing->namein,singin);
  }
  else {
    ADD_MEM(mesh,19*sizeof(char),"input singularities file name",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(sing->namein,19,char);
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
    fprintf(stdout,"  ## Warning: new solution\n");

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
  if ( np ) {
    sol->np  = np;
    sol->npi = np;
    if ( sol->m )
      DEL_MEM(mesh,sol->m,(sol->size*sol->npmax+1)*sizeof(double));

    sol->npmax = mesh->npmax;
    ADD_MEM(mesh,(sol->size*sol->npmax+1)*sizeof(double),"initial solution",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(sol->m,(sol->npmax*sol->size+1),double);
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
    fprintf(stdout,"  ## Warning: new mesh\n");

  mesh->np  = np;
  mesh->ne  = ne;
  mesh->nt  = nt;
  mesh->na  = na;
  mesh->npi = mesh->np;
  mesh->nei = mesh->ne;
  mesh->nti = mesh->nt;
  mesh->nai = mesh->na;

  if ( mesh->point )
    DEL_MEM(mesh,mesh->point,(mesh->npmax+1)*sizeof(Point));

  mesh->npmax = MG_MAX(1.5*mesh->np,NPMAX);
  ADD_MEM(mesh,(mesh->npmax+1)*sizeof(Point),"initial vertices",
          printf("  Exit program.\n");
          exit(EXIT_FAILURE));
  SAFE_CALLOC(mesh->point,mesh->npmax+1,Point);

  if ( mesh->tetra )
    DEL_MEM(mesh,mesh->tetra,(mesh->nemax+1)*sizeof(Tetra));

  mesh->nemax = MG_MAX(1.5*mesh->ne,NEMAX);
  ADD_MEM(mesh,(mesh->nemax+1)*sizeof(Tetra),"initial tetrahedra",
          printf("  Exit program.\n");
          exit(EXIT_FAILURE));
  SAFE_CALLOC(mesh->tetra,mesh->nemax+1,Tetra);

  if ( mesh->tria )
    DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(Tria));

  mesh->ntmax = MG_MAX(1.5*mesh->nt,NTMAX);
  if ( mesh->nt ) {
    ADD_MEM(mesh,(mesh->nt+1)*sizeof(Tria),"initial triangles",return(0));
    SAFE_CALLOC(mesh->tria,mesh->nt+1,Tria);
  }

  if ( mesh->edge )
    DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(Edge));

  mesh->namax = mesh->na;
  if ( mesh->na ) {
    ADD_MEM(mesh,(mesh->na+1)*sizeof(Edge),"initial edges",return(0));
    SAFE_CALLOC(mesh->edge,(mesh->na+1),Edge);
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
    fprintf(stdout,"     NUMBER OF VERTICES     %8d\n",mesh->np);
    if ( mesh->na ) {
      fprintf(stdout,"     NUMBER OF EDGES        %8d\n",mesh->na);
    }
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES    %8d\n",mesh->nt);
    fprintf(stdout,"     NUMBER OF ELEMENTS     %8d\n",mesh->ne);
  }
  return(1);
}

#ifdef SINGUL
/** Set the number of singular vertices and edges and allocate
    the associated tables. */
int Set_singulSize(MMG5_pMesh mesh,MMG5_pSingul sing, int np, int na) {

  if ( sing->point || sing->edge )
    fprintf(stdout,"  ## Warning: new singularites\n");

  sing->ns = np;
  sing->na = na;

  if ( sing->ns ) {
    if ( sing->point )
      DEL_MEM(mesh,sing->point,(sing->ns+1)*sizeof(sPoint));

    ADD_MEM(mesh,(sing->ns+1)*sizeof(sPoint),"vertex singularities",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(sing->point,sing->ns+1,sPoint);
  }

  if ( sing->na ) {
    if ( sing->edge )
      DEL_MEM(mesh,sing->edge,(sing->na+1)*sizeof(Edge));

    ADD_MEM(mesh,(sing->na+1)*sizeof(Edge),"edge singularities",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(sing->edge,sing->na+1,Edge);
  }
  return(1);
}
#endif

/** Get the solution number, dimension and type */
int Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, int* np, int* typSol) {

  *typEntity = MMG5_Vertex;
  *typSol    = sol->size;

  assert(sol->np = mesh->np);

  *np = mesh->np;
  sol->npi = 0;

  return(1);
}

/** Get the number of vertices, tetrahedra, triangles and edges of the mesh. */
int Get_meshSize(MMG5_pMesh mesh, int* np, int* ne, int* nt, int* na) {

  if ( np != NULL )
    *np = mesh->np;
  if ( ne != NULL )
    *ne = mesh->ne;
  if ( nt != NULL )
    *nt = mesh->nt;
  if ( na != NULL )
    *na = mesh->na;

  mesh->npi = 0;
  mesh->nei = 0;
  mesh->nti = 0;
  mesh->nai = 0;

  return(1);
}

/** Set vertex of coordinates c0,c1,c2 and reference ref at position pos
    in mesh structure */
int Set_vertex(MMG5_pMesh mesh, double c0, double c1, double c2, int ref, int pos) {

  if ( !mesh->np ) {
   fprintf(stdout,"  ## Error: You must set the number of points with the");
   fprintf(stdout," MMG5_Set_meshSize function before setting vertices in mesh\n");
   return(0);
  }

  if ( pos > mesh->npmax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new point.\n");
    fprintf(stdout,"    max number of points: %d\n",mesh->npmax);
    INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->np ) {
    fprintf(stdout,"  ## Error: attempt to set new vertex at position %d.",pos);
    fprintf(stdout," Overflow of the given number of vertices: %d\n",mesh->np);
    fprintf(stdout,"  ## Check the mesh size, its compactness or the position");
    fprintf(stdout," of the vertex.\n");
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

/** Get coordinates c0,c1,c2 and reference ref of next vertex of mesh  */
int Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, int* ref,
               int* isCorner, int* isRequired) {

  mesh->npi++;

  if ( mesh->npi > mesh->np ) {
    fprintf(stdout,"  ## Error: unable to get point.\n");
    fprintf(stdout,"     The number of call of MMG5_Get_vertex function");
    fprintf(stdout," can not exceed the number of points: %d\n ",mesh->np);
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

/** Set tetrahedra of vertices v0,v1,v2,v3 and reference ref at position pos
    in mesh structure */
int Set_tetrahedra(MMG5_pMesh mesh, int v0, int v1, int v2, int v3, int ref, int pos) {
  pTetra pt;
  pPoint ppt;
  double aux, vol;
  int    j, ip;

  if ( !mesh->ne ) {
   fprintf(stdout,"  ## Error: You must set the number of elements with the");
   fprintf(stdout," MMG5_Set_meshSize function before setting elements in mesh\n");
   return(0);
  }

  if ( pos > mesh->nemax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new element.\n");
    fprintf(stdout,"    max number of element: %d\n",mesh->nemax);
    INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->ne ) {
    fprintf(stdout,"  ## Error: attempt to set new tetrahedra at position %d.",pos);
    fprintf(stdout," Overflow of the given number of tetrahedra: %d\n",mesh->ne);
    fprintf(stdout,"  ## Check the mesh size, its compactness or the position");
    fprintf(stdout," of the tetrahedra.\n");
    return(0);
  }

  pt = &mesh->tetra[pos];
  pt->v[0] = v0;
  pt->v[1] = v1;
  pt->v[2] = v2;
  pt->v[3] = v3;
  pt->ref  = ref;

  mesh->point[pt->v[0]].tag &= ~MG_NUL;
  mesh->point[pt->v[1]].tag &= ~MG_NUL;
  mesh->point[pt->v[2]].tag &= ~MG_NUL;
  mesh->point[pt->v[3]].tag &= ~MG_NUL;

  vol = orvol(mesh->point,pt->v);
  if ( vol == 0.0 ) {
    fprintf(stdout,"  ## Error: tetrahedra volume is null.\n");
    for ( ip=0; ip<4; ip++ ) {
      ppt = &mesh->point[pt->v[ip]];
      for ( j=0; j<3; j++ ) {
        if ( fabs(ppt->c[j])>0. ) {
          fprintf(stdout," Check that you don't have a sliver tetrahedra.\n");
          return(0);
        }
      }
    }
    fprintf(stdout,"  All vertices have zero coordinates.");
    fprintf(stdout," Check that you have set the vertices before the tetrahedras.\n");
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

  pt->qual = orcal(mesh,pos);

  return(1);
}

/** Get vertices v0,v1,v2,v3 and reference ref of next tetra of mesh  */
int Get_tetrahedra(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* v3,
                   int* ref, int* isRequired) {

  mesh->nei++;

  if ( mesh->nei > mesh->ne ) {
    fprintf(stdout,"  ## Error: unable to get tetra.\n");
    fprintf(stdout,"    The number of call of MMG5_Get_tetrahedra function");
    fprintf(stdout," can not exceed the number of tetra: %d\n ",mesh->ne);
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

/** Set triangle of vertices v0,v1,v2 and reference ref at position pos
    in mesh structure */
int Set_triangle(MMG5_pMesh mesh, int v0, int v1, int v2, int ref,int pos) {

  if ( !mesh->nt ) {
   fprintf(stdout,"  ## Error: You must set the number of triangles with the");
   fprintf(stdout," MMG5_Set_meshSize function before setting triangles in mesh\n");
   return(0);
  }

  if ( pos > mesh->ntmax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new triangle.\n");
    fprintf(stdout,"    max number of triangle: %d\n",mesh->ntmax);
    INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->nt ) {
    fprintf(stdout,"  ## Error: attempt to set new triangle at position %d.",pos);
    fprintf(stdout," Overflow of the given number of triangles: %d\n",mesh->nt);
    fprintf(stdout,"  ## Check the mesh size, its compactness or the position");
    fprintf(stdout," of the triangle.\n");
    return(0);
  }

  mesh->tria[pos].v[0] = v0;
  mesh->tria[pos].v[1] = v1;
  mesh->tria[pos].v[2] = v2;
  mesh->tria[pos].ref  = ref;

  return(1);
}

/** Get vertices v0,v1,v2 and reference ref of next triangle of mesh  */
int Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref
                 ,int* isRequired) {
  pTria  ptt;

  mesh->nti++;

  if ( mesh->nti > mesh->nt ) {
    fprintf(stdout,"  ## Error: unable to get triangle.\n");
    fprintf(stdout,"    The number of call of MMG5_Get_triangle function");
    fprintf(stdout," can not exceed the number of triangles: %d\n ",mesh->nt);
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

/** Set edges of extremities v0,v1 and reference ref at position pos
    in mesh structure */
int Set_edges(MMG5_pMesh mesh, int v0, int v1, int ref, int pos) {

  if ( !mesh->na ) {
   fprintf(stdout,"  ## Error: You must set the number of edges with the");
   fprintf(stdout," MMG5_Set_meshSize function before setting edges in mesh\n");
   return(0);
  }
  if ( pos > mesh->namax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new edge.\n");
    fprintf(stdout,"    max number of edge: %d\n",mesh->namax);
    INCREASE_MEM_MESSAGE();
    return(0);
  }
  if ( pos > mesh->na ) {
    fprintf(stdout,"  ## Error: attempt to set new edge at position %d.",pos);
    fprintf(stdout," Overflow of the given number of edges: %d\n",mesh->na);
    fprintf(stdout,"  ## Check the mesh size, its compactness or the position");
    fprintf(stdout," of the edge.\n");
    return(0);
  }

  mesh->edge[pos].a = v0;
  mesh->edge[pos].b = v1;
  mesh->edge[pos].ref  = ref;
  mesh->edge[pos].tag |= MG_REF;

  return(1);
}

/** Get extremities e0,e1 and reference ref of next edge of mesh  */
int Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref
             ,int* isRidge, int* isRequired) {

  mesh->nai++;

  if ( mesh->nai > mesh->na ) {
    fprintf(stdout,"  ## Error: unable to get edge.\n");
    fprintf(stdout,"    The number of call of MMG5_Get_edge function");
    fprintf(stdout," can not exceed the number of edges: %d\n ",mesh->na);
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

/** Set scalar value s at position pos in solution structure */
int Set_scalarSol(MMG5_pSol met, double s, int pos) {

  if ( !met->np ) {
   fprintf(stdout,"  ## Error: You must set the number of solution with the");
   fprintf(stdout," MMG5_Set_solSize function before setting values");
   fprintf(stdout," in solution structure \n");
   return(0);
  }

  if ( pos >= met->npmax ) {
    fprintf(stdout,"  ## Error: unable to set a new solution.\n");
    fprintf(stdout,"    max number of solutions: %d\n",met->npmax);
    return(0);
  }

  if ( pos > met->np ) {
    fprintf(stdout,"  ## Error: attempt to set new solution at position %d.",pos);
    fprintf(stdout," Overflow of the given number of solutions: %d\n",met->np);
    fprintf(stdout,"  ## Check the solution size, its compactness or the position");
    fprintf(stdout," of the solution.\n");
    return(0);
  }

  met->m[pos] = s;
  return(1);
}

/** Get solution s of next vertex of mesh  */
int Get_scalarSol(MMG5_pSol met, double* s) {

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stdout,"  ## Error: unable to get solution.\n");
    fprintf(stdout,"     The number of call of MMG5_Get_scalarSol function");
    fprintf(stdout," can not exceed the number of points: %d\n ",met->np);
    return(0);
  }

  *s  = met->m[met->npi];

  return(1);
}

#ifdef SINGUL
/** Set singular point of coordinates c0,c1,c2 at position pos
    in the singularities structure */
int Set_singulVertex(MMG5_pSingul sing, double c0, double c1,
                     double c2, int typ, int pos) {

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
  if ( pos > sing->ns ) {
    fprintf(stdout,"  ## Error: attempt to set new singular");
    fprintf(stdout," vertex at position %d.",pos);
    fprintf(stdout," Overflow of the given number of sigular vertices: %d\n",sing->ns);
    fprintf(stdout,"  ## Check the singular mesh size, its compactness");
    fprintf(stdout," or the position of the singular vertex.\n");
    return(0);
  }

  sing->nsi++;
  sing->point[pos].c[0] = c0;
  sing->point[pos].c[1] = c1;
  sing->point[pos].c[2] = c2;
  return(1);
}


/** Set singular edge of extremities v0,v1 and reference ref
    at position pos in the singularities structure */
int Set_singulEdge(MMG5_pSingul sing, int v0, int v1, int ref, int pos) {

  if ( !sing->na ) {
    fprintf(stdout,"  ## Error: You must set the number of singular edges with the");
    fprintf(stdout," MMG5_Set_singulSize function before setting values");
    fprintf(stdout," in singularities structure \n");
    return(0);
  }
  if ( pos >= sing->na ) {
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
  if ( pos > sing->na ) {
    fprintf(stdout,"  ## Error: attempt to set new singular");
    fprintf(stdout," edge at position %d.",pos);
    fprintf(stdout," Overflow of the given number of singular edges: %d\n",sing->na);
    fprintf(stdout,"  ## Check the singular mesh size, its compactness");
    fprintf(stdout," or the position of the edge.\n");
    return(0);
  }

  sing->edge[pos].a   = v0;
  sing->edge[pos].b   = v1;
  sing->edge[pos].ref = ref;
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

/** To mark as ended a mesh given without using the API functions
    (for example, mesh given by mesh->point[i] = 0 ...)*/
void Set_handGivenMesh(MMG5_pMesh mesh) {
  int k, aux;

  /* Possibly switch 2 vertices number so that each tet is positively oriented */
  for (k=1; k<=mesh->ne; k++) {
    if ( orvol(mesh->point,mesh->tetra[k].v) < 0.0 ) {
      /* mesh->xt temporary used to count reoriented tetra */
      mesh->xt++;
      aux = mesh->tetra[k].v[2];
      mesh->tetra[k].v[2] = mesh->tetra[k].v[3];
      mesh->tetra[k].v[3] = aux;
    }
  }
  return;
}

/** Check if the number of given entities match with mesh and sol size (not mandatory)
    and check mesh datas  */
int Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met) {

  if ( (mesh->npi != mesh->np) || (mesh->nei != mesh->ne) ) {
    fprintf(stdout,"  ## Error: if you don't use the MMG5_loadMesh function,");
    fprintf(stdout," you must call the MMG5_Set_meshSize function to have a");
    fprintf(stdout," valid mesh.\n");
    fprintf(stdout," Missing datas.\n");
    return(0);
  }

  if ( met->npi != met->np ) {
    fprintf(stdout,"  ## Error: if you don't use the MMG5_loadMet function,");
    fprintf(stdout," you must call the MMG5_Set_solSize function to have a");
    fprintf(stdout," valid solution.\n");
    fprintf(stdout," Missing datas.\n");
    return(0);
  }

  /*  Check mesh data */
  if ( mesh->info.ddebug ) {
    if ( (!mesh->np) || (!mesh->point) ||
         (!mesh->ne) || (!mesh->tetra) ) {
      fprintf(stdout,"  ** MISSING DATA.\n");
      fprintf(stdout," Check that your mesh contains points and tetrahedra.\n");
      fprintf(stdout," Exit program.\n");
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
    if( !mesh->nti )
      DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(Tria));
    else {
      ADD_MEM(mesh,mesh->nti-mesh->nt,"triangles",return(0));
      SAFE_RECALLOC(mesh->tria,mesh->nt+1,(mesh->nti+1),Tria,"triangles");
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
      memcpy(pa,pa1,sizeof(Edge));
      pa1->ref = abs(pa1->ref);
    }
  } while( ++k <= mesh->nai );

  if ( mesh->nai < mesh->na ) {
    if( !mesh->nai )
      DEL_MEM(mesh,mesh->edge,(mesh->nai+1)*sizeof(Edge));
    else {
      ADD_MEM(mesh,mesh->nai-mesh->na,"Edges",return(0));
      SAFE_RECALLOC(mesh->edge,mesh->na+1,(mesh->nai+1),Edge,"edges");
    }
    mesh->na = mesh->nai;
  }

  /* delete tetrahedra references */
  for (k=1; k<=mesh->ne; k++) {
    mesh->tetra[k].ref = 0;
  }
  return(1);
}

/** Set integer parameter iparam at value val */
int Set_iparameters(MMG5_pMesh mesh, MMG5_pSol sol, int iparam, int val){
  int k;

  switch ( iparam ) {
    /* Integer parameters */
  case MMG5_IPARAM_verbose :
    mesh->info.imprim   = val;
    break;
  case MMG5_IPARAM_mem :
    if ( val <= 0 ) {
      fprintf(stdout,"  ## Warning: maximal memory authorized must be strictly positive.\n");
      fprintf(stdout,"  Reset to default value.\n");
    }
    else
      mesh->info.mem      = val;
    memOption(mesh);
    break;
#ifndef PATTERN
  case MMG5_IPARAM_bucket :
    mesh->info.bucket   = val;
    break;
#endif
  case MMG5_IPARAM_debug :
    mesh->info.ddebug   = val;
    break;
  case MMG5_IPARAM_angle :
    /* free table that may contains old ridges */
    if ( mesh->htab.geom )
      DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(hgeom));
    if ( mesh->xpoint )
      DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(xPoint));
    if ( mesh->xtetra )
      DEL_MEM(mesh,mesh->xtetra,(mesh->xtmax+1)*sizeof(xTetra));
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
      if ( mesh->nt && !skipIso(mesh) )
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
    if ( mesh->info.par ) {
      DEL_MEM(mesh,mesh->info.par,mesh->info.npar*sizeof(Par));
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stdout,"  ## Warning: new local parameter values\n");
    }
    mesh->info.npar  = val;
    mesh->info.npari = 0;
    ADD_MEM(mesh,mesh->info.npar*sizeof(Par),"parameters",
            printf("  Exit program.\n");
            exit(EXIT_FAILURE));
    SAFE_CALLOC(mesh->info.par,mesh->info.npar,Par);

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
int Set_dparameters(MMG5_pMesh mesh, MMG5_pSol sol, int dparam, double val){

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
    if ( val <=0 ) {
      fprintf(stdout,"  ## Warning: hausdorff number must be strictly positive.\n");
      fprintf(stdout,"  Reset to default value.\n");
    }
    else
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
int Set_localParameters(MMG5_pMesh mesh,MMG5_pSol sol, int typ, int ref, double val){
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
  if ( mesh->nameout ) {
    DEL_MEM(mesh,mesh->nameout,(strlen(mesh->nameout)+1)*sizeof(char));
  }

  if ( mesh->namein ) {
    DEL_MEM(mesh,mesh->namein,(strlen(mesh->namein)+1)*sizeof(char));
  }

  /* met */
  if ( met ) {
    if ( met->namein ) {
      DEL_MEM(mesh,met->namein,(strlen(met->namein)+1)*sizeof(char));
    }

    if ( met->nameout ) {
      DEL_MEM(mesh,met->nameout,(strlen(met->nameout)+1)*sizeof(char));
    }
  }
#ifdef SINGUL
  /* singul */
  if ( singul->namein ) {
    DEL_MEM(mesh,singul->namein,(strlen(singul->namein)+1)*sizeof(char));
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
  if ( mesh->point )
    DEL_MEM(mesh,mesh->point,(mesh->npmax+1)*sizeof(Point));

  if ( mesh->tetra )
    DEL_MEM(mesh,mesh->tetra,(mesh->nemax+1)*sizeof(Tetra));

  if ( mesh->edge )
    DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(Edge));

  if ( mesh->adja )
    DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

  if ( mesh->xpoint )
    DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(xPoint));

  if ( mesh->htab.geom )
    DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(hgeom));

  if ( mesh->tria )
    DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(Tria));

  if ( mesh->xtetra )
    DEL_MEM(mesh,mesh->xtetra,(mesh->xtmax+1)*sizeof(xTetra));

  /* met */
  if ( /*!mesh->info.iso &&*/ met && met->m )
    DEL_MEM(mesh,met->m,(met->size*met->npmax+1)*sizeof(double));

  /* mesh->info */
  if ( mesh->info.npar && mesh->info.par )
    DEL_MEM(mesh,mesh->info.par,mesh->info.npar*sizeof(Par));

#ifdef SINGUL
  /* singul */
  if ( singul->point )
    DEL_MEM(mesh,singul->point,(singul->ns+1)*sizeof(sPoint));
  if ( singul->edge )
    DEL_MEM(mesh,singul->edge,(singul->na+1)*sizeof(Edge));
#endif

  if ( mesh->info.imprim>6 || mesh->info.ddebug )
    printf("  MEMORY USED AT END (bytes) %lld\n",mesh->memCur);
}
