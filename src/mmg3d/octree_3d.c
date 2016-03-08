#include "mmg3d.h"

#define LFILT 0.2

#warning fonction to be deleted after debug
void _MMG3D_verifOctreeRec(MMG5_pMesh mesh, _MMG3D_octree_s* q, double* ver, const int no, const int nv)
{
  int i;
  int quadrant;
  int dim = mesh->dim;
  int test = 0;


  if (q->v != NULL)
  {
    for ( i = 0; i<q->nbVer; ++i)
    {
      if (q->v[i] == no)
      {
        test = 1;
        break;
      } 
    }
    assert(test);
   
  }else if (q->branches != NULL)
  {
    quadrant = 0;
    for ( i = 0; i<dim; ++i)
    {
      quadrant += ((double) (ver[i]>0.5))*(1<<i);
      ver[i] -= ((double) (ver[i]>0.5))*0.5;
      ver[i] *= 2;
    }

    _MMG3D_verifOctreeRec(mesh, &(q->branches[quadrant]), ver, no, nv);
  }
}

/**
 * \param q pointer toward the octree cell
 *
 * Initialisation of the octree cell.
 *
 */
void _MMG3D_initOctree_s( _MMG3D_octree_s* q)
{
  q->nbVer = 0;
  q->depth = 0;
  q->v = NULL;
  q->branches = NULL;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward the global octree
 * \param nv number of vertices in the cell subtree
 *
 * Initialisation of the octree cell.
 *
 */
void _MMG3D_initOctree(MMG5_pMesh mesh,_MMG3D_pOctree q, int nv)
{
#warning ajeter
  if ( mesh->info.imprim > 4 ) printf("OCTREE INIT\n");

  q->nv = nv;

  _MMG5_ADD_MEM(mesh,sizeof(_MMG3D_octree_s),"octree structure",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_MALLOC(q->q0,1, _MMG3D_octree_s);

  _MMG3D_initOctree_s(q->q0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward the octree cell
 * \param nv number of vertices in the cell subtree
 *
 * Free the octree cell.
 *
 */
void _MMG3D_freeOctree_s(MMG5_pMesh mesh,_MMG3D_octree_s* q, int nv)
{
  int nbBitsInt,depthMax,dim,i,sizTab,sizBr;

  dim       = mesh->dim;
  sizBr     = 1<<dim;
  nbBitsInt = sizeof(int)*8;
  depthMax  = nbBitsInt/dim - 1;

  if (q->nbVer>nv && q->depth < depthMax )
  {
    for ( i = 0; i<sizBr; i++) {
      _MMG3D_freeOctree_s(mesh,&(q->branches[i]), nv);
    }
    _MMG5_DEL_MEM(mesh,q->branches,sizBr*sizeof(_MMG3D_octree_s));
  }
  if (q->nbVer>0 && q->depth>0 ) {
    if ( q->nbVer<= nv )
      _MMG5_DEL_MEM(mesh,q->v,nv*sizeof(int));
    else {
      if ( q->depth != depthMax ) {
        sizTab = nv;
      }
      else {
        sizTab = (q->nbVer%nv != 0)? 1 : 0;
        sizTab = nv * ((int)(q->nbVer/nv) + sizTab);
      }
      _MMG5_DEL_MEM(mesh,q->v,sizTab*sizeof(int));
    }
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward the global octree.
 *
 * Free the global octree structure.
 *
 */
void _MMG3D_freeOctree(MMG5_pMesh mesh,_MMG3D_pOctree q)
{
  _MMG3D_freeOctree_s(mesh,q->q0, q->nv);
  _MMG5_DEL_MEM(mesh,q->q0,sizeof(_MMG3D_octree_s));
}

/**
 * \param q pointer toward an octree cell.
 * \param ver coordinates of the point.
 * \param dim space dimension (should be 3).
 * \param coordinate coordinate in the octree, should be 0 at first call.
 *
 * Get the integer containing the coordinates.
 *
 */
void _MMG3D_getOctreeCoordinateRec(_MMG3D_octree_s* q, double* ver, int dim, int* coordinate)
{  
  int i;
  int quadrant;
  if (q->branches != NULL)
  {
    quadrant = 0;
    for ( i = 0; i<dim; ++i)
    {
      quadrant += ((double) (ver[i]>0.5))*(1<<i);
      ver[i] -= ((double) (ver[i]>0.5))*0.5;
      ver[i] *= 2;
    }
    (*coordinate) += quadrant<<(dim*q->depth);
    _MMG3D_getOctreeCoordinateRec(&(q->branches[quadrant]), dim, ver, coordinate);
  }
}

/**
 * \param q pointer toward the global octree.
 * \param ver coordinates of the point.
 * \param dim space dimension (should be 3).
 *
 * Get the integer containing the coordinates
 *
 */
int _MMG3D_getOctreeCoordinate(_MMG3D_pOctree q, double* ver, int dim)
{
  double *pt;
  int coordinate;
  coordinate = 0;
  
  _MMG5_SAFE_MALLOC(pt,dim,double);
  memcpy(pt, ver ,dim*sizeof(double));
  
  _MMG3D_getOctreeCoordinateRec(q->q0, dim, pt, &coordinate);

  _MMG5_SAFE_FREE(pt);
  return coordinate;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward the global octree.
 * \param no index of the moved point.
 * \param newVer new coordinates for the moved point.
 * \param oldVer old coordinates for the moved point.
 *
 * Move one point in the octree structure. /!\ the vertex of index \a no
 * can have either the new or the old coordinates in the mesh but all
 * other vertices should have the same coordinates as when they were inserted
 * into the octree. (ie: one move at a time in the mesh and the octree)
 *
 */
void _MMG3D_moveOctree(MMG5_pMesh mesh, _MMG3D_pOctree q, int no, double* newVer, double* oldVer)
{
  int oldCoor, newCoor;
  double *pt;
  int dim;  
  dim = mesh->dim;
    
  _MMG5_SAFE_MALLOC(pt,dim,double);
  
  memcpy(pt, oldVer ,dim*sizeof(double));
  oldCoor = _MMG3D_getOctreeCoordinate(q, oldVer, dim);
  memcpy(pt, newVer ,dim*sizeof(double));
  newCoor = _MMG3D_getOctreeCoordinate(q, mesh->point[no].c, dim);
  
  if (newCoor == oldCoor)
    return;
  else // it could be possible to combine delOctree and addOctree to keep it local...
  {
    _MMG5_SAFE_MALLOC(pt,dim,double);
          
    /* delOctree */
    memcpy(pt, oldVer ,dim*sizeof(double));
    _MMG3D_delOctreeRec(mesh, q->q0, pt , no, q->nv);

    /* addOctree */
    memcpy(pt, newVer ,dim*sizeof(double));
    _MMG3D_addOctreeRec(mesh, q->q0, pt , no, q->nv);
  }
  _MMG5_SAFE_FREE(pt);
}

#warning remove one pointer on qlist ?
/**
 * \param q pointer toward the octree cell.
 * \param centre coordinates of the centre of the current subtree.
 * \param rect rectangle that we want to intersect with the subtree. We define
 * it given: the coordinates of one corner of the rectange and the length of
 * the rectangle in each dimension.
 * \param qlist pointer toward the list of pointer over the sub octrees that
 *  intersect \a rect.
 * \param nv number of vertices in the subtree.
 * \param dim dimension work.
 * \param index number of octree cells that intersect \a rect
 *
 * Count or list the number of octree cells that intersect the rectangle
 * \a rect depending on initialization of value on which \a qlist points on.
 *
 */
void _MMG3D_getListSquareRec(_MMG3D_octree_s* q, double* center, double* rect,
                             _MMG3D_octree_s*** qlist, int nv, int dim, int* index)
{
  double *recttemp;
  double *centertemp;
  double l = 1./(1<<(q->depth+1));

  _MMG5_SAFE_MALLOC(recttemp,2*dim,double);
  _MMG5_SAFE_MALLOC(centertemp,dim,double);

  if (q->branches==NULL)
  {
    if ((*qlist)!= NULL)
      (*qlist)[*index] = q;
    (*index)++;
    //~ fprintf(stdout,"Indexe incrémenté %i\n", *index);
  }else
  {

    if (dim == 2)
    {
      if (rect[0]<center[0]&&rect[1]<center[1]) // branch 0
      {
        recttemp[0] = rect[0];
        recttemp[1] = rect[1];
        recttemp[2] = rect[0]+rect[2] < center[0] ? rect[2]:center[0]-rect[0];
        recttemp[3] = rect[1]+rect[3] < center[1] ? rect[3]:center[1]-rect[1];

        centertemp[0] = center[0]-l/2;
        centertemp[1] = center[1]-l/2;

        _MMG3D_getListSquareRec(&(q->branches[0]), centertemp, recttemp, qlist, nv, dim, index);
      }
      if (rect[0]+rect[2] >center[0] && rect[1]<center[1]) // branch 1
      {
        recttemp[0] = rect[0]<center[0] ? center[0]:rect[0];
        recttemp[1] = rect[1];
        recttemp[2] = rect[0]+rect[2]-recttemp[0];
        recttemp[3] = rect[1]+rect[3] < center[1] ? rect[3]:center[1]-rect[1];

        centertemp[0] = center[0]+l/2;
        centertemp[1] = center[1]-l/2;

        _MMG3D_getListSquareRec(&(q->branches[1]), centertemp, recttemp, qlist, nv, dim, index);
      }
      if (rect[0]<center[0] && rect[1]+rect[3]>center[1]) // branch 2
      {
        recttemp[0] = rect[0];
        recttemp[1] = rect[1]>center[1] ? rect[1]:center[1];
        recttemp[2] = rect[0]+rect[2] < center[0] ? rect[2]:center[0]-rect[0];
        recttemp[3] = rect[1] + rect[3]- recttemp[1];

        centertemp[0] = center[0]-l/2;
        centertemp[1] = center[1]+l/2;

        _MMG3D_getListSquareRec(&(q->branches[2]), centertemp, recttemp, qlist, nv, dim, index);
      }
      if (rect[0]+rect[2] >center[0] && rect[1]+rect[3]>center[1]) // branch 3
      {
        recttemp[0] = rect[0]>center[0] ? rect[0]:center[0];
        recttemp[1] = rect[1]>center[1] ? rect[1]:center[1];
        recttemp[2] = rect[0]+rect[2]-recttemp[0];
        recttemp[3] = rect[1]+rect[3]-recttemp[1];

        centertemp[0] = center[0]+l/2;
        centertemp[1] = center[1]+l/2;

        _MMG3D_getListSquareRec(&(q->branches[3]), centertemp, recttemp, qlist, nv, dim, index);
      }
    }else if (dim == 3)
    {
      if (rect[0]<center[0]&&rect[1]<center[1]) // branch 0
      {
        recttemp[0] = rect[0];
        recttemp[1] = rect[1];

        recttemp[3] = rect[0]+rect[3] < center[0] ? rect[3]:center[0]-rect[0];
        recttemp[4] = rect[1]+rect[4] < center[1] ? rect[4]:center[1]-rect[1];

        centertemp[0] = center[0]-l/2;
        centertemp[1] = center[1]-l/2;
        if (rect[2] < center[2]) // branch 0
        {

          recttemp[2] = rect[2];
          recttemp[5] = rect[2]+rect[5] < center[2] ? rect[5]:center[2]-rect[2];


          centertemp[2] = center[2]-l/2;

          _MMG3D_getListSquareRec(&(q->branches[0]), centertemp, recttemp, qlist, nv, dim, index);
        }
        if (rect[2]+rect[5] > center[2]) // branch 4
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          _MMG3D_getListSquareRec(&(q->branches[4]), centertemp, recttemp, qlist, nv, dim, index);
        }
      }
      if (rect[0]+rect[2] >center[0] && rect[1]<center[1]) // branch 1
      {
        recttemp[0] = rect[0]<center[0] ? center[0]:rect[0];
        recttemp[1] = rect[1];
        recttemp[3] = rect[0]+rect[3]-recttemp[0];
        recttemp[4] = rect[1]+rect[4] < center[1] ? rect[4]:center[1]-rect[1];
        centertemp[0] = center[0]+l/2;
        centertemp[1] = center[1]-l/2;

        if(rect[2] < center[2]) // branch 1
        {
          recttemp[2] = rect[2];
          recttemp[5] = rect[2]+rect[5] < center[2] ? rect[5]:center[2]-rect[2];
          centertemp[2] = center[2]-l/2;

          _MMG3D_getListSquareRec(&(q->branches[1]), centertemp, recttemp, qlist, nv, dim, index);
        }
        if (rect[2]+rect[5] > center[2]) // branch 5
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          _MMG3D_getListSquareRec(&(q->branches[5]), centertemp, recttemp, qlist, nv, dim, index);
        }
      }
      if (rect[0]<center[0] && rect[1]+rect[3]>center[1]) // branch 2
      {
        recttemp[0] = rect[0];
        recttemp[1] = rect[1]>center[1] ? rect[1]:center[1];
        recttemp[3] = rect[0]+rect[3] < center[0] ? rect[3]:center[0]-rect[0];
        recttemp[4] = rect[1] + rect[4]- recttemp[1];

        centertemp[0] = center[0]-l/2;
        centertemp[1] = center[1]+l/2;


        if(rect[2] < center[2]) // branch 2
        {

          recttemp[2] = rect[2];
          recttemp[5] = rect[2] + rect[5]< center[2] ? rect[5]:center[2]-rect[2];

          centertemp[2] = center[2]-l/2;

          _MMG3D_getListSquareRec(&(q->branches[2]), centertemp, recttemp, qlist, nv, dim, index);
        }
        if (rect[2]+rect[5] >= center[2]) // branch 6
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          _MMG3D_getListSquareRec(&(q->branches[6]), centertemp, recttemp, qlist, nv, dim, index);
        }
      }
      if (rect[0]+rect[2] >center[0] && rect[1]+rect[3]>center[1]) // branch 3
      {
        recttemp[0] = rect[0]>center[0] ? rect[0]:center[0];
        recttemp[1] = rect[1]>center[1] ? rect[1]:center[1];
        recttemp[3] = rect[0]+rect[3]-recttemp[0];
        recttemp[4] = rect[1]+rect[4]-recttemp[1];
        centertemp[0] = center[0]+l/2;
        centertemp[1] = center[1]+l/2;

        if(rect[2] < center[2]) // branch 3
        {

          recttemp[2] = rect[2];

          recttemp[5] = rect[2] + rect[5]< center[2] ? rect[5]:center[2]-rect[2];

          centertemp[2] = center[2]-l/2;

          _MMG3D_getListSquareRec(&(q->branches[3]), centertemp, recttemp, qlist, nv, dim, index);
        }
        if (rect[2]+rect[5] > center[2]) // branch 7
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          _MMG3D_getListSquareRec(&(q->branches[7]), centertemp, recttemp, qlist, nv, dim, index);
        }
      }
    }
  }
  _MMG5_SAFE_FREE(recttemp);
  _MMG5_SAFE_FREE(centertemp);

}

/**
 * \param q pointer toward the global octree structure.
 * \param rect rectangle that we want to intersect with the subtree. We define
 * it given: the coordinates of one corner of the rectange and the length of
 * the rectangle in each dimension.
 * \param qlist pointer toward the list of pointer over the sub octrees that
 *  intersect \a rect.
 *
 * \return index, the number of subtrees in the list, 0 if fail.
 *
 * List the number of octree cells that intersect the rectangle \a rect.
 *
 */
int _MMG3D_getListSquare(MMG5_pMesh mesh,_MMG3D_pOctree q, double* rect,
                         _MMG3D_octree_s*** qlist)
{
  double *rect2, *center;
  int    i,index,dim;

  dim = mesh->dim;
  index = 0;

  _MMG5_SAFE_MALLOC(rect2,2*dim,double);
  _MMG5_SAFE_MALLOC(center,dim,double);

  memcpy(rect2, rect, sizeof(double)*dim*2);

  for (i = 0; i < dim; ++i)
    center[i] = 0.5;

  //~ fprintf(stdout,"q nulle ? %i\n", q->branches == NULL);
  assert(!(*qlist));
  _MMG3D_getListSquareRec(q->q0,center, rect2, qlist, q->nv, dim, &index);
  //~ fprintf(stdout, "index : %i\n", index);

  _MMG5_ADD_MEM(mesh,index*sizeof(_MMG3D_octree_s*),"octree cell",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));

  _MMG5_SAFE_MALLOC(*qlist,index,_MMG3D_octree_s*);

  for (i = 0; i<index; i++)
    (*qlist)[i] = NULL;


  index = 0;
  memcpy(rect2, rect, sizeof(double)*dim*2);
  for (i = 0; i < dim; i++)
    center[i] = 0.5;

  _MMG3D_getListSquareRec(q->q0, center, rect2, qlist, q->nv, dim, &index);

  _MMG5_SAFE_FREE(rect2);
  _MMG5_SAFE_FREE(center);

  return index;
}

#warning change the return values to not exit when malloc fail (we can save the mesh or try to continue to mesh)
/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward an octree cell.
 * \param ver vertex coordinates scaled such that the quadrant is [0;1]x[0;1]x[0;1]
 * \param no vertex index in the mesh.
 * \param nv maximum number of points in an octree cell. 
 *
 * Add vertex in the suitable quadrant of the octree. This function is
 * recursively called until we reach the last one. At each step, the vertex
 * coordinates are scaled such as the quadrant is the [0;1]x[0;1]x[0;1] box.
 *
 */
void _MMG3D_addOctreeRec(MMG5_pMesh mesh, _MMG3D_octree_s* q, double* ver,
                           const int no, int nv)
{
  double   *pt;
  int      dim, nbBitsInt,depthMax,i,j,k;
  int      quadrant,sizBr;

  nbBitsInt = sizeof(int)*8;
  dim       = mesh->dim;
  depthMax  = nbBitsInt/dim - 1;
  sizBr     = 1<<dim;

  _MMG5_SAFE_MALLOC(pt,dim,double);

  if ( q->depth < depthMax )
  {
    if (q->nbVer < nv)
    {
      //~ fprintf(stdout, "ok1 %f %f \n", ver[0], ver[1]);
      if (q->v == NULL)
      {
        //~ fprintf(stdout, "ok1bis\n");
        _MMG5_ADD_MEM(mesh,nv*sizeof(int),"octree vertices table",
                      printf("  Exit program.\n");
                      exit(EXIT_FAILURE));
        _MMG5_SAFE_MALLOC(q->v,nv,int);
      }

// on peut faire mieux en triant les vertices ajoutés, ça permet ensuite de
// faire directement descendre les points à leur place lorsque l'on dépasse le
// nombre de points dans la feuille et de rechercher plus rapidement dans la
// feuille
        q->v[q->nbVer] = no;
      // intéressant seulement si nv est grand
        q->nbVer++;
      return;
    }
    else if (q->nbVer == nv && q->branches==NULL)
    {
      /* creation of sub-branch and relocation of vertices in the sub-branches */
      //~ fprintf(stdout, "ok2\n");
      _MMG5_ADD_MEM(mesh,sizBr*sizeof(_MMG3D_octree_s),"octree branches",
                    printf("  Exit program.\n");
                    exit(EXIT_FAILURE));
      _MMG5_SAFE_MALLOC(q->branches,sizBr,_MMG3D_octree_s);

      for ( i = 0; i<sizBr; i++)
      {
        _MMG3D_initOctree_s(&(q->branches[i]));
        q->branches[i].depth = q->depth+1;
      }
      q->nbVer++;
      for (i = 0; i<nv; i++)
      {

        memcpy(pt, mesh->point[q->v[i]].c ,dim*sizeof(double));
        for ( j =0; j < q->depth; j++)
        {
          for (k = 0; k<dim; k++)
          {
            pt[k] -= ((double) (pt[k]>0.5))*0.5;
            pt[k] *= 2;
          }
        }
        //~ addVertexRec(mesh, q, mesh->point[q->v[i]].c, q->v[i]);
        _MMG3D_addOctreeRec(mesh, q, pt, q->v[i],nv);
        q->nbVer--;
      }
      _MMG3D_addOctreeRec(mesh, q, ver, no, nv);
      q->nbVer--;
      _MMG5_DEL_MEM(mesh,q->v,nv*sizeof(int));

    }else
    {
      quadrant = 0;
      for ( i = 0; i<dim; i++)
      {
        quadrant += ((double) (ver[i]>0.5))*(1<<i);
        ver[i] -= ((double) (ver[i]>0.5))*0.5;
        ver[i] *= 2;
      }

      q->nbVer++;
      _MMG3D_addOctreeRec(mesh, &(q->branches[quadrant]), ver, no, nv);
    }
  }else
  {
    if (q->nbVer%nv == 0)
      _MMG5_SAFE_REALLOC(q->v,q->nbVer+nv,int,"octree");

    q->v[q->nbVer] = no;
    q->nbVer++;
  }
  _MMG5_SAFE_FREE(pt);
}

/**
* \param pointer toward the mesh structure
* \param q pointer toward the global octree structure
* \param no index of the point to add to the octree
*
* Add the vertex of index \a no to the octree.
*
*/
void _MMG3D_addOctree(MMG5_pMesh mesh, _MMG3D_pOctree q, const int no)
{
  double *pt;
  int    dim;

  dim = mesh->dim;
  _MMG5_SAFE_MALLOC(pt,dim,double);

  memcpy(pt, mesh->point[no].c ,dim*sizeof(double));
  _MMG3D_addOctreeRec(mesh, q->q0, pt , no, q->nv);
  memcpy(pt, mesh->point[no].c ,dim*sizeof(double));
  #warning to be removed
  _MMG3D_verifOctreeRec(mesh, q->q0,pt,no,q->nv);
  _MMG5_SAFE_FREE(pt);
}

/**
* \param q pointer toward a terminal octree cell (containing vertex)
* \param no index of the point to delete from the octree
*
* Delete the vertex of index \a no from the terminal octree cell.
*
*/
void _MMG3D_delOctreeVertex(_MMG3D_octree_s* q, int indNo)
{
  assert(q->v);
  assert(q->nbVer>indNo);
  memmove(&q->v[indNo],&q->v[indNo+1], (q->nbVer-indNo-1)*sizeof(int));
  --q->nbVer;
}

/**
 * \param q0 pointer toward an octree cell.
 * \param q pointer toward an octree cell.
 * \param dim dimension of the space (=3).
 * \param nv maximum number of points in an octree cell.
 * \param index next index in the array to be filled.
 *
 * Merge sub-branches \a q of \a q0, in their parent \a q0. \a q0 should
 * contain no more than nv vertices.
 *
 */
void _MMG3D_mergeBranchesRec(_MMG3D_octree_s* q0, _MMG3D_octree_s* q, int dim, int nv, int* index)
{
  int i;
  
  if (q->v != NULL)
  {
    if (*index+q->nbVer == nv+1)
      fprintf(stdout, "Index too high %i\n", *index+q->nbVer);
    fprintf(stdout, "Index %i nbVer %i\n", *index,q->nbVer);
    assert(*index+q->nbVer<=nv);
    memcpy(&(q0->v[*index]), q->v, q->nbVer*sizeof(int));
    (*index)+= q->nbVer;
  }else if (q->branches != NULL)
  {
    for (i = 0; i<(1<<dim); ++i)
      _MMG3D_mergeBranchesRec(q0, &(q->branches[i]), dim, nv, index);
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward an octree cell.
 * \param dim dimension of the space (=3)
 * \param nv maximum number of points in an octree cell. 
 *
 * Merge branches that have a parent counting less than nv vertices.
 *
 */
void _MMG3D_mergeBranches(MMG5_pMesh mesh,_MMG3D_octree_s* q, int dim, int nv)
{
  int index;
  int i;
  index = 0;
  assert(q->v);
  assert(q->branches);
  assert(q->nbVer==nv);

  for (i = 0; i<(1<<dim); ++i)
  {
    //~ _MMG3D_printSubArbre(&(q->branches[i]),nv,dim);
    _MMG3D_mergeBranchesRec(q, &(q->branches[i]), dim, nv, &index); 
    _MMG3D_freeOctree_s(mesh,&(q->branches[i]), nv);
  }
  q->nbVer = index;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward an octree cell.
 * \param ver vertex coordinates scaled such that the quadrant is [0;1]x[0;1]x[0;1]
 * \param no vertex index in the mesh.
 * \param nv maximum number of points in an octree cell. 
 *
 * Delete vertex \a no from the octree. This function is recursively 
 * called until we reach the terminal octree cell containing the vertex 
 * \a no. At each step, the vertex coordinates are scaled such as the
 * quadrant is the [0;1]x[0;1]x[0;1] box.
 *
 */
void _MMG3D_delOctreeRec(MMG5_pMesh mesh, _MMG3D_octree_s* q, double* ver, const int no, const int nv)
{
  int i;
  int quadrant;
  int dim = mesh->dim;
  int depthMax;
  int nbBitsInt;
  #warning test to delete
  int test = 0;
  nbBitsInt = sizeof(int)*8;

  depthMax  = nbBitsInt/dim - 1;

  if (q->v != NULL)
  {
    assert(q->v);
    for ( i = 0; i<q->nbVer; ++i)
    {
      if (q->v[i] == no)
      {
        _MMG3D_delOctreeVertex(q, i);
        if ( q->nbVer == 0)
          _MMG5_DEL_MEM(mesh,q->v,nv*sizeof(int));
        test = 1;
        break;
      } 
    }
    if (!test)
    {
      fprintf(stdout,"Le petit %i est attendu à l'accueil\n Accueil : ", no);
      for ( i = 0; i<q->nbVer; ++i)
        fprintf(stdout,"%i ",q->v[i]);
      fprintf(stdout,"\n");
    }else
    {
      fprintf(stdout,"OK :)\n");
    }
      //assert(test);
   
  }else if ( q->nbVer == nv+1)
  {
    quadrant = 0;
    for ( i = 0; i<dim; ++i)
    {
      quadrant += ((double) (ver[i]>0.5))*(1<<i);
      ver[i] -= ((double) (ver[i]>0.5))*0.5;
      ver[i] *= 2;
    }
    --q->nbVer;
    fprintf(stdout,"quadrant %i nbVer %i\n", quadrant,q->branches[quadrant].nbVer);
    #warning calling recursively here is not optimal
    _MMG3D_delOctreeRec(mesh, &(q->branches[quadrant]), ver, no, nv);
    fprintf(stdout,"quadrant %i nbVer %i\n", quadrant,q->branches[quadrant].nbVer);
    
    _MMG5_ADD_MEM(mesh,nv*sizeof(int),"octree vertices table",
                  printf("  Exit program.\n");
                  exit(EXIT_FAILURE));
    _MMG5_SAFE_MALLOC(q->v,nv,int);
    _MMG3D_mergeBranches(mesh,q,dim,nv);
      
  }else if (q->branches != NULL)
  {
    quadrant = 0;
    for ( i = 0; i<dim; ++i)
    {
      quadrant += ((double) (ver[i]>0.5))*(1<<i);
      ver[i] -= ((double) (ver[i]>0.5))*0.5;
      ver[i] *= 2;
    }

    --q->nbVer;
    _MMG3D_delOctreeRec(mesh, &(q->branches[quadrant]), ver, no, nv);
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward the global octree.
 * \param no reference of the vertex to be deleted.
 *
 * Delete the vertex \a no from the octree structure.
 *
 */
void _MMG3D_delOctree(MMG5_pMesh mesh, _MMG3D_pOctree q, const int no)
{
  double *pt;
  int    dim;

  dim = mesh->dim;
  _MMG5_SAFE_MALLOC(pt,dim,double);

  memcpy(pt, mesh->point[no].c ,dim*sizeof(double));
  _MMG3D_delOctreeRec(mesh, q->q0, pt , no, q->nv);

  _MMG5_SAFE_FREE(pt);
}


/**
* \param q pointer toward an octree cell
* \param depth depth of the subtree
* \param nv number of vertices in the subtree
* \param dim dimension in which we work
*
* Print the depth \a depth of the subtree of \a q.
*
* \warning debug function, not safe
*/
void _MMG3D_printArbreDepth(_MMG3D_octree_s* q, int depth, int nv, int dim)
{
  int i;
  if ( q->depth < depth && q->nbVer > nv)
  {

    for (i = 0; i < (1<<dim); i++)
    {
      _MMG3D_printArbreDepth(&(q->branches[i]),depth, nv, dim);
    }
    //~ fprintf(stdout," | ");
  }else if (q->depth == depth)
  {
    fprintf(stdout,"%i ",q->nbVer);
  }
}

/**
* \param q pointer toward the global octree structure
*
* Print the octree.
*
* \warning debug function, not safe
*
*/
void _MMG3D_printArbre(_MMG3D_pOctree q)
{
  int dim;

  dim = 3;
  int i;
  for (i = 0; i<sizeof(int)*8/dim; i++)
  {
    fprintf(stdout,"\n profondeur %i \n", i);
    _MMG3D_printArbreDepth(q->q0, i, q->nv, dim);

  }
  fprintf(stdout,"\n fin \n");
}

/**
* \param q pointer toward an octree cell
* \param nv maximum number of vertices in an octree leaf
* \param dim spacial dimension 
* 
* Print the octree.
*
* \warning debug function, not safe
*
*/
void _MMG3D_printSubArbre(_MMG3D_octree_s* q, int nv, int dim)
{
  int i;
  for (i = 0; i<sizeof(int)*8/dim; i++)
  {
    fprintf(stdout,"\n profondeur %i \n", i);
    _MMG3D_printArbreDepth(q, i, nv, dim);

  }
  fprintf(stdout,"\n fin \n");
}


/**
* \param q pointer toward an octree cell
* \param nv maximum number of vertices in an octree leaf
* \param dim dimension in which we work
*
* Print the memory size of the octree.
*
* \warning debug function, not safe
*/
int _MMG3D_sizeArbreRec(_MMG3D_octree_s* q, int nv, int dim)
{
  int i, sizeBranches;

  sizeBranches = 0;

  if (q->branches != NULL)
  {

    for (i= 0; i <4; i++)
    {
      sizeBranches += _MMG3D_sizeArbreRec(&(q->branches[i]),nv,dim)
        +sizeof(_MMG3D_octree_s)+(1<<dim)*sizeof(_MMG3D_octree_s*);
    }
    return sizeBranches;
  }else if(q->v != NULL)
  {
    return nv*sizeof(int)+sizeof(_MMG3D_octree_s);
  }else
  {
    return sizeof(_MMG3D_octree_s);
  }
}

/**
* \param q pointer toward the global octree structure
*
* Print the octree memory size.
*
* \warning debug function, not safe
*
*/
int _MMG3D_sizeArbre(_MMG3D_pOctree q,int dim)
{
  return _MMG3D_sizeArbreRec(q->q0, q->nv, dim);
}


/**
* \param q pointer toward an octree cell
* \param nv number of vertices in the subtree
* \param dim dimension in which we work
*
* Print the memory size of the octree for point stored with a linked list.
*
* \warning debug function, not safe
*/
static inline
int _MMG3D_sizeArbreLinkRec(_MMG3D_octree_s* q, int nv, int dim)
{
  int sizeBranches,i;

  sizeBranches = 0;

  if (q->branches != NULL)
  {

    for (i= 0; i <4; i++)
    {
      sizeBranches += _MMG3D_sizeArbreLinkRec(&(q->branches[i]), nv, dim)
        +sizeof(_MMG3D_octree_s)+(1<<dim)*sizeof(_MMG3D_octree_s*);
    }
    return sizeBranches;
  }else if(q->v != NULL)
  {
    return sizeof(int)+sizeof(_MMG3D_octree_s);
  }else
  {
    return sizeof(_MMG3D_octree_s);
  }
}

/**
* \param q pointer toward the global octree
*
* Print the memory size of the octree for point stored with a linked list.
*
* \warning debug function, not safe
*/
static inline
int _MMG3D_sizeArbreLink(_MMG3D_pOctree q)
{
  int dim;

  dim = 3;
  return _MMG3D_sizeArbreLinkRec(q->q0, q->nv, dim)+q->q0->nbVer*sizeof(int);
}


static inline
int NearNeighborBrutForce(MMG5_pMesh mesh, int no, double l, int dim, int nb_pt)
{
  double lmin =10;
  int nmin = -1;
  double x[dim], norm;
  int i, j;
  for(i = 1; i <= nb_pt; i++)
  {
    norm = 0;
    if(i != no)
    {
      for (j = 0; j<dim; j++)
      {

        x[j] =  mesh->point[i].c[j] - mesh->point[no].c[j];
      }
      for (j = 0; j<dim; j++)
      {
        norm += x[j]*x[j];
      }
      //~ fprintf(stdout, "c'est la norme %g entre %i et %i\n", norm, i, no);
      if(lmin>norm)
      {
        //~ fprintf(stdout,"Valide\n");
        lmin = norm;
        nmin = i;
      }
    }
  }
  if (sqrt(lmin)<l)
    return nmin;
  else
    return -1;
}

static inline
int NearNeighborSquare(MMG5_pMesh mesh, _MMG3D_pOctree q, int no, double l, int dim)
{
  MMG5_pPoint   ppt,ppt1;
  _MMG3D_octree_s** qlist;
  int ns,nver;
  double rect[2*dim];
  double lmin =10;
  double x,y,z;
  int nmin;
  int i, j;

  ppt = &mesh->point[no];
  rect[0] = ppt->c[0]-l;
  rect[1] = ppt->c[1]-l;
  rect[2] = ppt->c[2]-l;
  rect[3] = 2*l;
  rect[4] = 2*l;
  rect[5] = 2*l;
  qlist = NULL;
  ns = _MMG3D_getListSquare(mesh,q, rect, &qlist);

  
  for (i = 0; i < ns; i++)
  {
    for (j = 0; j<qlist[i]->nbVer; j++)
    {
      nver = qlist[i]->v[j];
      if(nver != no)
      {
        ppt  = &mesh->point[nver];
        ppt1 = &mesh->point[no];

        x = ppt->c[0] - ppt1->c[0];
        y = ppt->c[1] - ppt1->c[1];
        z = ppt->c[2] - ppt1->c[2];
        x = x*x+y*y+z*z;
        if(lmin>x)
        {
          lmin = x;
          nmin = nver;
        }
      }
    }
  }

  _MMG5_SAFE_FREE(qlist);

  if (sqrt(lmin)<l)
    return nmin;
  else
    return -1;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param octree pointer toward the octree structure.
 * \param ip index of point to check.
 *
 * \return 1 if we can insert \a ip, 0 otherwise
 *
 * Check if the vertex \a ip is not too close from another one.
 *
 */
int _MMG3D_octreein_iso(MMG5_pMesh mesh,MMG5_pSol sol,_MMG3D_pOctree octree,int ip) {
  MMG5_pPoint     ppt,pp1;
  _MMG3D_octree_s **lococ;
  double          d2,ux,uy,uz,hpi,hp1,hpi2,methalo[6];
  int             ip1,i,j;
  int             ncells;

  lococ = NULL;
  ppt = &mesh->point[ip];

  hpi = LFILT * sol->m[ip];
  hp1 = hpi*hpi;

  /* methalo is the box that we want to intersect with the octree, thus, the limit
   * of the filter. We give: the coordinates of one of the corner of the box and
   * the length of the box in each direction. */
  methalo[0] = ppt->c[0] - hpi;
  methalo[1] = ppt->c[1] - hpi;
  methalo[2] = ppt->c[2] - hpi;
  methalo[3] = methalo[4] = methalo[5] = 2.*hpi;

  ncells = _MMG3D_getListSquare(mesh,octree, methalo, &lococ);

  /* Check the octree cells */
  for ( i=0; i<ncells; ++i ) {
    for (j=0; j<lococ[i]->nbVer; ++j)
    {
      ip1  = lococ[i]->v[j];
      pp1  = &mesh->point[ip1];


#warning to remove when delOctree will be implemented
      //~ if ( !MG_VOK(pp1) ) continue;

      hpi2 = LFILT * sol->m[ip1];

      ux = pp1->c[0] - ppt->c[0];
      uy = pp1->c[1] - ppt->c[1];
      uz = pp1->c[2] - ppt->c[2];

      d2 = ux*ux + uy*uy + uz*uz;

      if ( d2 < hp1 || d2 < hpi2*hpi2 )  {
        //printf("filtre current %d : %e %e %e %e\n",ip1,d2,hp1,d2,hpi2*hpi2);
        return(0);
      }
    }

  }
  _MMG5_DEL_MEM(mesh,lococ,ncells*sizeof(_MMG3D_octree_s**));
  return(1);
}
