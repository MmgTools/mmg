#include "mmg3d.h"

#define LFILT 0.2

#warning nbVer utilise? to remove?

/* double _MMG3D_distMet(double* x, double* metric, int dim) */
/* { */
/*   double dist = 0; */
/*   int ind = 0; */
/*   for(int i = 0; i< dim; i++) */
/*   { */
/*     for (int j = i; j<dim; j++) */
/*     { */
/*       dist += metric[ind]*x[i]*x[j]; */
/*       ind++; */
/*     } */
/*   } */
/*   return dist; */
/* } */

/* void sorties_vtk(MMG5_pMesh mesh, int* connectionList, int nb_pt, char* fich) */
/* { */
/*   MMG5_pPoint ppt; */
/*   int imax = nb_pt, jmax = nb_pt; */

/*   FILE* fileVtk; */

/*   if (!(fileVtk = fopen(fich,"wb"))) */
/*     fprintf(stdout,"sorties_vtk: unable to open fije %s\n", fich); */

/*   // ----------------------Header------------------------------------- */
/*   fprintf(fileVtk,"# vtk DataFile Version 3.1\n"); */
/*   fprintf(fileVtk,"Mesh nodes\n"); */
/*   fprintf(fileVtk,"ASCII\n"); */
/*   fprintf(fileVtk,"DATASET UNSTRUCTURED_GRID\n\n"); */

/*   //----------------------Nodes position------------------------------ */
/*   fprintf(fileVtk,"POINTS    %i    float\n", nb_pt); */
/*   for (int i = 1; i<=nb_pt; i++) */
/*   { */
/*     ppt = &mesh->point[i]; */
/*     fprintf(fileVtk,"%g    %g    %g\n",ppt->c[0],ppt->c[1],ppt->c[2]); */
/*   } */

/*   fprintf(fileVtk,"\n"); */

/*   //---------------------Elements------------------------------------- */
/*   fprintf(fileVtk,"CELLS    %i    %i\n", nb_pt, 3*nb_pt); */
/*   for (int j = 1; j<=nb_pt; j++) */
/*   { */
/*     fprintf(fileVtk, "2    %i    %i\n", j-1, connectionList[j]-1); */
/*   } */

/*   fprintf(fileVtk,"\n"); */

/*   fprintf(fileVtk,"CELL_TYPES    %i\n",nb_pt); */
/*   for (int j = 1; j<=nb_pt; j++) */
/*   { */
/*     fprintf(fileVtk,"3 "); */
/*   } */
/*   fprintf(fileVtk, "\n\n"); */

/*   //---------------------Data----------------------------------------- */
/*   fprintf(fileVtk,"POINT_DATA    %i\n", nb_pt); */
/*   fprintf(fileVtk,"SCALARS data float 1\n"); */
/*   fprintf(fileVtk,"LOOKUP_TABLE default\n"); */
/*   for (int j = 1; j<=nb_pt; j++) */
/*   { */
/*     fprintf(fileVtk, "%i\n",(int) (j == connectionList[j])); */
/*   } */

/*   fprintf(fileVtk,"\n"); */

/*   fclose(fileVtk); */
/* } */


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
 * \param q pointer toward the global octree
 * \param nv number of vertices in the cell subtree
 *
 * Initialisation of the octree cell.
 *
 */
void _MMG3D_initOctree(_MMG3D_pOctree q, int nv)
{
  q->nv = nv;

  _MMG5_SAFE_MALLOC(q->q0,sizeof(_MMG3D_octree_s), _MMG3D_octree_s);

  _MMG3D_initOctree_s(q->q0);
}

/**
 * \param q pointer toward the octree cell
 * \param nv number of vertices in the cell subtree
 *
 * Free the octree cell.
 *
 */
void _MMG3D_freeOctree_s(_MMG3D_octree_s* q, int nv, int dim)
{
  int nbBitsInt,depthMax;

  nbBitsInt = sizeof(int)*8;
  depthMax  = nbBitsInt/dim - 1;

  if (q->nbVer>nv && q->depth < depthMax )
  {
    for (int i = 0; i<(1<<dim); i++)
      _MMG3D_freeOctree_s(&(q->branches[i]), nv, dim);
    _MMG5_SAFE_FREE(q->branches);
  }
  if (q->nbVer>0 && (q->nbVer<= nv || q->depth == depthMax))
  {
    _MMG5_SAFE_FREE(q->v);
  }
}

/**
 * \param q pointer toward the global octree.
 *
 * Free the global octree structure.
 *
 */
void _MMG3D_freeOctree(_MMG3D_pOctree q, int dim)
{
  _MMG3D_freeOctree_s(q->q0, q->nv, dim);
  _MMG5_SAFE_FREE(q->q0);
}


/**
 * \param q pointer toward the octree cell.
 * \param centre coordinates of the centre of the current subtree.
 * \param rect rectangle that we want to intersect with the subtree. We define
 * it given: the coordinates of one corner of the rectange and the length of
 * the rectangle in each dimension.
 * \param nv number of vertices in the subtree.
 * \param dim dimension work.
 * \param index number of octree cells that intersect \a rect
 *
 * Count the number of octree cells that intersect the rectangle \a rect. Result
 * is put in the \a index variable.
 *
 */
void _MMG3D_countListSquareRec(_MMG3D_octree_s* q, double* center, double rect[6],
                               int nv, int dim, int* index)
{
  double *recttemp;
  double *centertemp;

  double l = 1./(1<<(q->depth+1));

  _MMG5_SAFE_MALLOC(recttemp,2*dim,double);
  _MMG5_SAFE_MALLOC(centertemp,dim,double);

  if (q->branches==NULL)
  {
    (*index)++;
    //~ fprintf(stdout,"Indexe incrémenté %i\n", *index);
  }else
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

        _MMG3D_countListSquareRec(&(q->branches[0]), centertemp, recttemp, nv, dim, index);
      }
      if (rect[2]+rect[5] > center[2]) // branch 4
      {
        recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
        recttemp[5] = rect[2]+ rect[5] - recttemp[2];

        centertemp[2] = center[2]+l/2;

        _MMG3D_countListSquareRec(&(q->branches[4]), centertemp, recttemp, nv, dim, index);
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

        _MMG3D_countListSquareRec(&(q->branches[1]), centertemp, recttemp, nv, dim, index);
      }
      if (rect[2]+rect[5] > center[2]) // branch 5
      {
        recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
        recttemp[5] = rect[2]+ rect[5] - recttemp[2];

        centertemp[2] = center[2]+l/2;

        _MMG3D_countListSquareRec(&(q->branches[5]), centertemp, recttemp, nv, dim, index);
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

        _MMG3D_countListSquareRec(&(q->branches[2]), centertemp, recttemp, nv, dim, index);
      }
      if (rect[2]+rect[5] >= center[2]) // branch 6
      {
        recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
        recttemp[5] = rect[2]+ rect[5] - recttemp[2];

        centertemp[2] = center[2]+l/2;

        _MMG3D_countListSquareRec(&(q->branches[6]), centertemp, recttemp, nv, dim, index);
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

        _MMG3D_countListSquareRec(&(q->branches[3]), centertemp, recttemp, nv, dim, index);
      }
      if (rect[2]+rect[5] > center[2]) // branch 7
      {
        recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
        recttemp[5] = rect[2]+ rect[5] - recttemp[2];

        centertemp[2] = center[2]+l/2;

        _MMG3D_countListSquareRec(&(q->branches[7]), centertemp, recttemp, nv, dim, index);
      }
    }
  }
  _MMG5_SAFE_FREE(recttemp);
  _MMG5_SAFE_FREE(centertemp);
}

#warning remove one pointer on qlist +  try to mutualize the code with the previous one.
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
 * Count and list the number of octree cells that intersect the rectangle \a rect.
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
int _MMG3D_getListSquare(_MMG3D_pOctree q, double* rect,
                         _MMG3D_octree_s*** qlist, int dim)
{
  double *rect2, *center;
  int    i,index;

  index = 0;

  _MMG5_SAFE_MALLOC(rect2,2*dim,double);
  _MMG5_SAFE_MALLOC(center,dim,double);

  memcpy(rect2, rect, sizeof(double)*dim*2);

  for (i = 0; i < dim; ++i)
    center[i] = 0.5;

  //~ fprintf(stdout,"q nulle ? %i\n", q->branches == NULL);
  _MMG3D_countListSquareRec(q->q0,center, rect2, q->nv, dim, &index);
  //~ fprintf(stdout, "index : %i\n", index);

  _MMG5_SAFE_MALLOC(qlist,index,_MMG3D_octree_s**);

  for (i = 0; i<index; i++)
    (*qlist)[i] = NULL;

  //~ if (qlist == NULL)
  //~ fprintf(stdout, "niuratisecuinserirt");
  index = 0;
  memcpy(rect2, rect, sizeof(double)*dim*2);
  for (i = 0; i < dim; i++)
    center[i] = 0.5;

  _MMG3D_getListSquareRec(q->q0, center, rect2, qlist, q->nv, dim, &index);

  _MMG5_SAFE_FREE(rect2);
  _MMG5_SAFE_FREE(center);

  return index;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward an octree cell.
 * \param ver vertex coordinates scaled such that the quadrant is [0;1]x[0;1]x[0;1]
 * \param no vertex index in the mesh.
 * \param nv subtree number of points.
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
  int      quadrant;

  nbBitsInt = sizeof(int)*8;
  dim       = mesh->dim;
  depthMax  = nbBitsInt/dim - 1;

  _MMG5_SAFE_MALLOC(pt,dim,double);

  if ( q->depth < depthMax )
  {
    if (q->nbVer < nv)
    {
      //~ fprintf(stdout, "ok1 %f %f \n", ver[0], ver[1]);
      if (q->v == NULL)
      {
        //~ fprintf(stdout, "ok1bis\n");
        _MMG5_SAFE_MALLOC(q->v,nv,int);
      }

        q->v[q->nbVer] = no; // on peut faire mieux en triant les vertices ajoutés, ça permet ensuite de faire directement descendre les points à leur place lorsque l'on dépasse le nombre de points dans la feuille et de rechercher plus rapidement dans la feuille
      // intéressant seulement si nv est grand
        q->nbVer++;
      return;
    }else if (q->nbVer == nv && q->branches==NULL) //creation of sub-branch and relocation of vertices in the sub-branches
    {
      //~ fprintf(stdout, "ok2\n");
      _MMG5_SAFE_MALLOC(q->branches,(1<<dim),_MMG3D_octree_s);

      for ( i = 0; i<(1<<dim); i++)
      {
        _MMG3D_initOctree_s(&(q->branches[i]));
        q->branches[i].depth = q->depth+1;
      }
      //~ unsigned char quadrant = ver[0]>0.5+2*(ver[1]>0.5);
      //~ if ((q->branches[quadrant].v = (int *) malloc(nv*sizeof(int)))==NULL)
      //~ fprintf(stdout, "Probleme allocation vertex\n");
      //~ q->branches[quadrant].v[0] = no;
      //~ q->branches[quadrant].nbVer++;
      q->nbVer++;
      for (i = 0; i<nv; i++)
      {

        memcpy(pt,&mesh->point[q->v[i]] ,dim*sizeof(double));
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
      _MMG5_SAFE_FREE(q->v);

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
      _MMG5_SAFE_REALLOC(q->v,nv*(q->nbVer/nv+1),int,"octree");

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

  memcpy(pt, &mesh->point[no] ,dim*sizeof(double));
  //~ fprintf(stdout, "ajout du point %f %f %f\n", pt[0], pt[1], pt[2]);
  _MMG3D_addOctreeRec(mesh, q->q0, pt , no, q->nv);

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
  if ( q->depth < depth && q->nbVer > nv)
  {

    for (int i = 0; i < (1<<dim); i++)
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

  for (int i = 0; i<sizeof(int)*8/dim; i++)
  {
    fprintf(stdout,"\n profondeur %i \n", i);
    _MMG3D_printArbreDepth(q->q0, i, q->nv, dim);

  }
  fprintf(stdout,"\n fin \n");
}


/**
* \param q pointer toward an octree cell
* \param nv number of vertices in the subtree
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
      sizeBranches += _MMG3D_sizeArbreLinkRec(&(q->branches[i]), nv, dim)+sizeof(_MMG3D_octree_s)+(1<<dim)*sizeof(_MMG3D_octree_s*);
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
  for(int i = 1; i <= nb_pt; i++)
  {
    norm = 0;
    if(i != no)
    {
      for (int j = 0; j<dim; j++)
      {

        x[j] =  mesh->point[i].c[j] - mesh->point[no].c[j];
      }
      for (int j = 0; j<dim; j++)
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

  ppt = &mesh->point[no];
  rect[0] = ppt->c[0]-l;
  rect[1] = ppt->c[1]-l;
  rect[2] = ppt->c[2]-l;
  rect[3] = 2*l;
  rect[4] = 2*l;
  rect[5] = 2*l;
  //~ fprintf(stdout,"q nulle ? %i\n", q->branches == NULL);
  ns = _MMG3D_getListSquare(q, rect, &qlist, mesh->dim);

  //~ fprintf(stdout, "OK %i\n",ns);
  //~ if (ns>0)
  //~ {
  //~ fprintf(stdout, "OK %i\n",ns);
  //~ fprintf(stdout,"qlist %i\n", qlist[0]->nbVer);
  //~ fprintf(stdout, "OK %i\n",ns);
  //~
  //~ }
  for (int i = 0; i < ns; i++)
  {
    for (int j = 0; j<qlist[i]->nbVer; j++)
    {
      nver = qlist[i]->v[j];
      //~ fprintf(stdout,"nver %i\n",nver);
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
  //~ for (int i = 0; i < ns; i++)
  //~ {
  //~ for (int j = i+1; j < ns; j++)
  //~ {
  //~ if (qlist[i] == qlist[j])
  //~ fprintf(stdout,"Problème\n");
  //~ }
  //~ }

  _MMG5_SAFE_FREE(qlist);

  if (sqrt(lmin)<l)
    return nmin;
  else
    return -1;
}


//~
//~ void plotOctree( MMG5_pMesh mesh, _MMG3D_pOctree q, char* name)
//~ {
//~ //Declare postscript
//~ MMG2D_postscript ps;
//~ //Open postscript file
//~ MMG2D_openPostscript(&ps,name);
//~ double ver[2];
//~ int NN;
//~ double color[3];
//~ color[0] = 0;
//~ color[1] = 0;
//~ color[2] = 0;
//~ double lmoy = 0.3*sqrt(2./(3*sqrt(3)*q->q0->nbVer));
//~ lmoy = lmoy<0.05 ? lmoy:0.05;
//~
//~ MMG2D_sethueColor(&ps, 0.);
//~ fprintf(ps.file,"/Helvetica findfont\n");
//~ fprintf(ps.file,"%i scalefont setfont\n", 8);
//~ plotOctreePointRec(&ps, mesh, q->q0, lmoy, 0, 0);
//~ MMG2D_sethsbColor(&ps, color);
//~ MMG2D_setLineWidth(&ps, 0.0001);
//~ plotOctreeSquareRec(&ps, mesh, q->q0, lmoy, 0, 0);
//~ double x,y,r;
//~ MMG2D_setLineWidth(&ps, 0.3);
//~ int NNBF;
//~
//~
//~ struct timeval tm1;
//~ start(&tm1);
//~ for (int i = 1; i<= mesh->np; i++)
//~ NNBF = NearNeighborBrutForce(mesh, i, 0.01, q->dim, mesh->np);
//~ fprintf(stdout,"temps brutForce : ");
//~ stop(&tm1);
//~ start(&tm1);
//~ for (int i = 1; i<= mesh->np; i++)
//~ NN = NearNeighborSquare(mesh,q,i,0.01);
//~ fprintf(stdout,"temps _MMG3D_octree_s : ");
//~ stop(&tm1);
//~ for (int i = 1; i<= mesh->np; i++)
//~ {
//~ memcpy(ver,mesh->point[i].c,sizeof(double)*2);
//NN = getNearestNeighbor( mesh, q, ver, i);
//int NNBF = NearestNeighborBrutForce(mesh, i);
//~
//~ NN = NearNeighborSquare(mesh,q,i,0.001);
//~ NNBF = NearNeighborBrutForce(mesh, i, 0.001, q->dim);
//~
//~ if (NN != NNBF)
//~ {
//~ x = mesh->point[i].c[0] - mesh->point[NN].c[0];
//~ y = mesh->point[i].c[1] - mesh->point[NN].c[1];
//~ r = x*x+y*y;
//~ x = mesh->point[i].c[0] - mesh->point[NNBF].c[0];
//~ y = mesh->point[i].c[1] - mesh->point[NNBF].c[1];
//~ x = y*y+x*x;
//~ fprintf(stdout," erreur %i trouvé %i (%f) au lieu de %i (%f)\n", i, NN, r, NNBF, x);
//~ MMG2D_setLineWidth(&ps, 1.8);
//~ }
//~
//~ if (NN > 0)
//~ MMG2D_drawLine( mesh, &ps, i, NN);
//~ else
//fprintf(stdout, "Noeud non relié %i\n", i);
//~
//~ if (NN != NNBF)
//~ {
//~ MMG2D_setLineWidth(&ps, 0.3);
//~ }
//~ }

//~ MMG2D_closePostscript(&ps);
//~ }

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

  ncells = _MMG3D_getListSquare(octree, methalo, &lococ,3);

  /* Check the octree cells */
  for ( i=0; i<ncells; ++i ) {
    for (j=0; j<lococ[i]->nbVer; ++j)
    {
      ip1  = lococ[i]->v[j];
      pp1  = &mesh->point[ip1];


#warning to remove when delOctree will be implemented
      if ( !MG_VOK(pp1) ) continue;

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

  return(1);
}
