#include "mmg2d.h"
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Allocate the mesh and solutions structures at \a MMG2D format.
 *
 */
static inline
void MMG2_Alloc_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol
  ) {

  /* mesh allocation */
  if ( *mesh )  _MMG5_SAFE_FREE(*mesh);
  _MMG5_SAFE_CALLOC(*mesh,1,MMG5_Mesh);

  /* sol allocation */
  if ( *sol )  _MMG5_DEL_MEM(*mesh,*sol,sizeof(MMG5_Sol));
  _MMG5_SAFE_CALLOC(*sol,1,MMG5_Sol);

  return;
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Initialization of mesh and solution structures to their default
 * values (default names, versions, dimensions...).
 *
 */
static inline
void MMG2_Init_woalloc_mesh(MMG5_pMesh mesh, MMG5_pSol sol
  ) {

  _MMG2_Set_commonFunc();

  (mesh)->dim = 2;
  (mesh)->ver = 2;
  (sol)->dim  = 2;
  (sol)->ver  = 2;
  (sol)->size = 1;

  /* Default parameters values */
  _MMG2_Init_parameters(mesh);

  /* Default vaules for file names */
  MMG5_Init_fileNames(mesh,sol);

  return;
}
/**
 * \param mesh pointer toward a pointer toward the mesh structure.
 * \param sol pointer toward a pointer toward the sol structure.
 *
 * Allocate the mesh and solution structures and initialize it to
 * their default values.
 *
 */
void MMG2_Init_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol
  ) {

  /* allocations */
  MMG2_Alloc_mesh(mesh,sol);
  /* initialisations */
  MMG2_Init_woalloc_mesh(*mesh,*sol);
  return;
}
/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 */
void _MMG2_Init_parameters(MMG5_pMesh mesh) {

  /* Init common parameters for mmg2d, mmgs and mmg3d. */
  _MMG5_mmgInit_parameters(mesh);

 /* default values for integers */
  /** MMG5_IPARAM_iso = 0 */
  mesh->info.iso      =  0;  /* [0/1]    ,Turn on/off levelset meshing */
  /** MMG5_IPARAM_lag = -1 */
  mesh->info.lag      = -1;
  /** MMG5_IPARAM_optim = 0 */
  mesh->info.optim    =  0;
  /** MMG5_IPARAM_nosurf = 0 */
  mesh->info.nosurf   =  0;  /* [0/1]    ,avoid/allow surface modifications */

  mesh->info.renum    = 0;   /* [0]    , Turn on/off the renumbering using SCOTCH; */
  mesh->info.nreg    = 0;
  /* default values for doubles */
  mesh->info.ls       = 0.0;      /* level set value */
  mesh->info.hgrad    = 1.3;      /* control gradation; */

  mesh->info.dhd  = 135.;

  //mesh->info.imprim = -7;

  /** MMG5_IPARAM_bucket = 64 */
  mesh->info.bucket = 64;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param iparam integer parameter to set (see \a MMG5_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 */
int MMG5_Set_iparameter(MMG5_pMesh mesh, MMG5_pSol sol, int iparam, int val){
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
    _MMG5_memOption(mesh);
    if(mesh->np && (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt )) {
      return(0);
    } else if(mesh->info.mem < 39)
      return(0);
    break;
  case MMG5_IPARAM_bucket :
    mesh->info.bucket   = val;
    break;
  case MMG5_IPARAM_debug :
    mesh->info.ddebug   = val;
    break;
  case MMG5_IPARAM_angle :
#warning dhd
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
        fprintf(stdout,"  ## Warning: angle detection parameter set to default value\n");
      mesh->info.dhd    = _MMG5_ANGEDG;
    }
    break;
  case MMG5_IPARAM_iso :
    mesh->info.iso      = val;
#warning : iso case
    // if ( mesh->info.iso )
    //if ( mesh->nt && !MMG5_skipIso(mesh) )
    //  exit(EXIT_FAILURE);
    break;
  case MMG5_IPARAM_lag :
    if ( val < 0 || val > 2 )
      exit(EXIT_FAILURE);
    mesh->info.lag = val;
    break;
  case MMG5_IPARAM_msh :
    mesh->info.nreg = val;
    break;
  case MMG5_IPARAM_numsubdomain :
    mesh->info.renum = val;
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
  default :
    fprintf(stdout,"  ## Error: unknown type of parameter\n");
    return(0);
  }
  /* other options */
  mesh->info.fem      = 0;
  return(1);
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param dparam double parameter to set (see \a MMG5_Param structure).
 * \val value of the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
 *
 */
int MMG5_Set_dparameter(MMG5_pMesh mesh, MMG5_pSol sol, int dparam, double val){

  switch ( dparam ) {
    /* double parameters */
  case MMG5_DPARAM_angleDetection :
    mesh->info.dhd = val;
    mesh->info.dhd = MG_MAX(0.0, MG_MIN(180.0,mesh->info.dhd));
    mesh->info.dhd = 180. - mesh->info.dhd;
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
    break;
  case MMG5_DPARAM_hausd :
    if ( val <=0 ) {
      fprintf(stdout,"  ## Error: hausdorff number must be strictly positive.\n");
      return(0);
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



/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param nt number of triangles.
 * \param na number of edges.
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of vertices, tetrahedra, triangles and edges of the
 * mesh and allocate the associated tables. If call twice, reset the
 * whole mesh to realloc it at the new size
 *
 */
int MMG5_Set_meshSize(MMG5_pMesh mesh, int np, int nt, int na) {
  int k;

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) &&
       ( mesh->point || mesh->tria || mesh->edge) )
    fprintf(stdout,"  ## Warning: new mesh\n");

  mesh->np  = np;
  mesh->nt  = nt;
  mesh->na  = na;
  mesh->npi = mesh->np;
  mesh->nti = mesh->nt;
  mesh->nai = mesh->na;

  if ( mesh->point )
    _MMG5_DEL_MEM(mesh,mesh->point,(mesh->npmax+1)*sizeof(MMG5_Point));
  if ( mesh->tria )
    _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->ntmax+1)*sizeof(MMG5_Tria));
  if ( mesh->edge )
    _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->namax+1)*sizeof(MMG5_Edge));

  /*tester si -m definie : renvoie 0 si pas ok et met la taille min dans info.mem */
  if( mesh->info.mem > 0) {
    if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->namax < mesh->na) ) {
      _MMG5_memOption(mesh);
      //     printf("pas de pbs ? %d %d %d %d %d %d -- %d\n",mesh->npmax,mesh->np,
      //     mesh->ntmax,mesh->nt,mesh->nemax,mesh->ne,mesh->info.mem);
      if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt)) {
        fprintf(stdout,"mem insuffisante np : %d %d nt : %d %d \n"
                ,mesh->npmax,mesh->np,
                mesh->ntmax,mesh->nt);
        return(0);
      }
      else
        return(1);
    } else if(mesh->info.mem < 39) {
      printf("mem insuffisante %d\n",mesh->info.mem);
      return(0);
    }
  } else {
    mesh->npmax = MG_MAX(1.5*mesh->np,_MMG5_NPMAX);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,_MMG5_NEMAX);

  }
  _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point);

  _MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),"initial triangles",return(0));
  _MMG5_SAFE_CALLOC(mesh->tria,mesh->ntmax+1,MMG5_Tria);

  mesh->namax =  MG_MAX(mesh->na,_MMG5_NEDMAX);
  _MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"initial edges",return(0));
  _MMG5_SAFE_CALLOC(mesh->edge,(mesh->namax+1),MMG5_Edge);

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->nt + 1;
  mesh->nanil = mesh->na + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++) {
    mesh->point[k].tmp  = k+1;
  }
  for (k=mesh->nenil; k<mesh->ntmax-1; k++) {
    mesh->tria[k].v[2] = k+1;
  }
  for (k=mesh->nanil; k<mesh->namax-1; k++) {
    mesh->edge[k].b = k+1;
  }
   
  /* memory alloc */
  _MMG5_ADD_MEM(mesh,(3*mesh->ntmax+5)*sizeof(int),"adjacency table",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,int);

  /* stats */
  if ( abs(mesh->info.imprim) > 6 ) {
    fprintf(stdout,"     NUMBER OF VERTICES     %8d\n",mesh->np);
    if ( mesh->na ) {
      fprintf(stdout,"     NUMBER OF EDGES        %8d\n",mesh->na);
    }
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES    %8d\n",mesh->nt);
  }
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typEntity type of solutions entities (vertices, triangles...).
 * \param np number of solutions.
 * \param typSol type of solution (scalar, vectorial...).
 * \return 0 if failed, 1 otherwise.
 *
 * Set the solution number, dimension and type.
 *
 */
int MMG5_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, int np, int typSol) {

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) && sol->m )
    fprintf(stdout,"  ## Warning: new solution\n");

  if ( typEntity != MMG5_Vertex ) {
    fprintf(stdout,"  ## Error: MMG2D5 need a solution imposed on vertices\n");
    return(0);
  }
  if ( typSol == MMG5_Scalar ) {
    sol->size = 1;
  }
  else if ( typSol == MMG5_Tensor ) {
    sol->size = 6;
  }
  else {
    fprintf(stdout,"  ## Error: type of solution not yet implemented\n");
    return(0);
  }

  sol->dim = 2;
  if ( np ) {
    sol->np  = np;
    sol->npi = np;
    if ( sol->m )
      _MMG5_DEL_MEM(mesh,sol->m,(sol->size*(sol->npmax+1))*sizeof(double));

    sol->npmax = mesh->npmax;
    _MMG5_ADD_MEM(mesh,(sol->size*(sol->npmax+1))*sizeof(double),"initial solution",
                  printf("  Exit program.\n");
                  exit(EXIT_FAILURE));
    _MMG5_SAFE_CALLOC(sol->m,(sol->size*(sol->npmax+1)),double);
  }
  return(1);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param np pointer toward the number of vertices.
 * \param nt pointer toward the number of triangles.
 * \param na pointer toward the number of edges.
 * \return 1.
 *
 * Get the number of vertices, triangles and edges of the mesh.
 *
 */
int MMG5_Get_meshSize(MMG5_pMesh mesh, int* np, int* nt, int* na) {

  if ( np != NULL )
    *np = mesh->np;
  if ( nt != NULL )
    *nt = mesh->nt;
  if ( na != NULL )
    *na = mesh->na;

  return(1);
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param c0 coordinate of the point along the first dimension.
 * \param c1 coordinate of the point along the second dimension.
 * \param ref point reference.
 * \param pos position of the point in the mesh.
 * \return 1.
 *
 * Set vertex of coordinates \a c0, \a c1 and reference \a ref
 * at position \a pos in mesh structure
 *
 */
int MMG5_Set_vertex(MMG5_pMesh mesh, double c0, double c1, int ref, int pos) {

  if ( !mesh->np ) {
    fprintf(stdout,"  ## Error: You must set the number of points with the");
    fprintf(stdout," MMG5_Set_meshSize function before setting vertices in mesh\n");
    return(0);
  }

  if ( pos > mesh->npmax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new point.\n");
    fprintf(stdout,"    max number of points: %d\n",mesh->npmax);
    _MMG5_INCREASE_MEM_MESSAGE();
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
  mesh->point[pos].ref  = ref;
  mesh->point[pos].tag  = MG_NUL;
  mesh->point[pos].flag = 0;
  mesh->point[pos].tmp = 0;

  return(1);
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param num integer
 * \param c0 pointer toward the coordinate of the point along the first dimension.
 * \param c1 pointer toward the coordinate of the point along the second dimension.
 * \param ref poiter to the point reference.
 * \param isCorner pointer toward the flag saying if point is corner.
 * \param isCorner pointer toward the flag saying if point is required.
 * \return 1.
 *
 * Get coordinates \a c0, \a c1 and reference \a ref of 
 * vertex num of mesh.
 *
 */
int MMG5_Get_vertex(MMG5_pMesh mesh,int num, double* c0, double* c1, int* ref,
                    int* isCorner, int* isRequired) {

  if ( num > mesh->np ) {
    fprintf(stdout,"  ## Error: unable to get point.\n");
    fprintf(stdout,"     The number %d in MMG5_Get_vertex function",num);
    fprintf(stdout,"  exceed the max number of points: %d\n ",mesh->np);
    return(0);
  }

  *c0  = mesh->point[num].c[0];
  *c1  = mesh->point[num].c[1];
  if ( ref != NULL )
    *ref = mesh->point[num].ref;

  if ( isCorner != NULL ) {
    if ( mesh->point[num].tag & M_CORNER )
      *isCorner = 1;
    else
      *isCorner = 0;
  }

  if ( isRequired != NULL ) {
    if ( mesh->point[num].tag & M_REQUIRED )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of triangle.
 * \param v1 second vertex of triangle.
 * \param v2 third vertex of triangle.
 * \param ref triangle reference.
 * \param pos triangle position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set triangle of vertices \a v0, \a v1, \a v2 and reference \a ref
 * at position \a pos in mesh structure.
 *
 */
int MMG5_Set_triangle(MMG5_pMesh mesh, int v0, int v1, int v2, int ref, int pos) {
  MMG5_pTria pt;
  MMG5_pPoint ppt;
  double  vol;
  int    i,j, ip,tmp;

  if ( !mesh->nt ) {
    fprintf(stdout,"  ## Error: You must set the number of elements with the");
    fprintf(stdout," MMG5_Set_meshSize function before setting elements in mesh\n");
    return(0);
  }

  if ( pos > mesh->ntmax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new element.\n");
    fprintf(stdout,"    max number of element: %d\n",mesh->ntmax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->nt ) {
    fprintf(stdout,"  ## Error: attempt to set new triangle at position %d.",pos);
    fprintf(stdout," Overflow of the given number of triangle: %d\n",mesh->nt);
    fprintf(stdout,"  ## Check the mesh size, its compactness or the position");
    fprintf(stdout," of the triangle.\n");
    return(0);
  }

  pt = &mesh->tria[pos];
  pt->v[0] = v0;
  pt->v[1] = v1;
  pt->v[2] = v2;
  pt->ref  = ref;

  mesh->point[pt->v[0]].tag &= ~MG_NUL;
  mesh->point[pt->v[1]].tag &= ~MG_NUL;
  mesh->point[pt->v[2]].tag &= ~MG_NUL;

  for(i=0 ; i<3 ; i++)
    pt->edg[i] = 0;
      
  vol = MMG2_quickarea(mesh->point[pt->v[0]].c,mesh->point[pt->v[1]].c,
                           mesh->point[pt->v[2]].c);
  if(vol < 0) {
    printf("Tr %d bad oriented\n",pos);
    tmp = pt->v[2];
    pt->v[2] = pt->v[1];
    pt->v[1] = tmp;
    /* mesh->xt temporary used to count reoriented tetra */
    mesh->xt++;
  }

  return(1);
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of edge.
 * \param v1 second vertex of edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set edge of vertices \a v0, \a v1 and reference \a ref
 * at position \a pos in mesh structure.
 *
 */
int MMG5_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref, int pos) {
  MMG5_pEdge pt;
  MMG5_pPoint ppt;
  double  vol;
  int    i,j, ip,tmp;

  if ( !mesh->na ) {
    fprintf(stdout,"  ## Error: You must set the number of elements with the");
    fprintf(stdout," MMG5_Set_meshSize function before setting elements in mesh\n");
    return(0);
  }

  if ( pos > mesh->namax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new element.\n");
    fprintf(stdout,"    max number of element: %d\n",mesh->namax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->na ) {
    fprintf(stdout,"  ## Error: attempt to set new edge at position %d.",pos);
    fprintf(stdout," Overflow of the given number of edge: %d\n",mesh->na);
    fprintf(stdout,"  ## Check the mesh size, its compactness or the position");
    fprintf(stdout," of the edge.\n");
    return(0);
  }

  pt = &mesh->edge[pos];
  pt->a = v0;
  pt->b = v1;
  pt->ref  = ref;

  mesh->point[pt->a].tag &= ~MG_NUL;
  mesh->point[pt->b].tag &= ~MG_NUL;

  return(1);
}
/**
 * \param met pointer toward the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure
 *
 */
int MMG5_Set_scalarSol(MMG5_pSol met, double s, int pos) {

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
/**
 * \param met pointer toward the sol structure.
 * \param s solution symetric tensor value (s11 s12 s22)
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensor value \a s at position \a pos in solution structure
 *
 */
int MMG5_Set_tensorSol(MMG5_pSol met, double* s, int pos) {
  int isol;

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
  isol = (pos-1) * met->size + 1;
  met->m[isol + 0] = s[0];
  met->m[isol + 1] = s[1];
  met->m[isol + 2] = s[2];
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Check if the number of given entities match with mesh and sol size
 * (not mandatory) and check mesh datas.
 *
 */
int MMG5_Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met) {

  if ( (mesh->npi != mesh->np) || (mesh->nti != mesh->nt) ) {
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
         (!mesh->nt)  ) {
      fprintf(stdout,"  ** MISSING DATA.\n");
      fprintf(stdout," Check that your mesh contains points.\n");
      fprintf(stdout," Exit program.\n");
      return(0);
    }
  }

  if ( mesh->dim != 2 ) {
    fprintf(stdout,"  ** 2 DIMENSIONAL MESH NEEDED. Exit program.\n");
    return(0);
  }
  if ( met->dim != 2 ) {
    fprintf(stdout,"  ** WRONG DIMENSION FOR METRIC. Exit program.\n");
    return(0);
  }
  if ( !mesh->ver )  mesh->ver = 2;
  if ( !met ->ver )  met ->ver = 2;

  return(1);
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * Structure deallocations before return.
 *
 */
void MMG5_Free_structures(MMG5_pMesh mesh,MMG5_pSol met
  ){

  MMG5_Free_names(mesh,met);

  /* mesh */
  if ( mesh->point )
    _MMG5_DEL_MEM(mesh,mesh->point,(mesh->npmax+1)*sizeof(MMG5_Point));

  if ( mesh->edge )
    _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(MMG5_Edge));

  if ( mesh->adja )
    _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

  if ( mesh->xpoint )
    _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));

  if ( mesh->tria )
    _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(MMG5_Tria));

  /* met */
  if ( /*!mesh->info.iso &&*/ met && met->m )
    _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));

  if ( mesh->info.imprim>6 || mesh->info.ddebug )
    printf("  MEMORY USED AT END (bytes) %lld\n",mesh->memCur);
}
