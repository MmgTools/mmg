#include "mmg3d.h"

double distMet(double* x, double* metric, int dim)
{
  double dist = 0;
  int ind = 0;
  for(int i = 0; i< dim; i++)
  {
    for (int j = i; j<dim; j++)
    {
      dist += metric[ind]*x[i]*x[j];
      ind++;
    }
  }
  return dist;
}

typedef struct preQuadtree
{
  int* v;
  unsigned int nbVer; // attention il faut mettre long sur proc 32 bits sinon le max est 65 535 noeuds
  unsigned char depth;
  struct preQuadtree* branches; // voir type recursif, si impossible : pointer void et typage à l'allocation.
} preQuadtree;

typedef struct quadtree
{
  int nv;
  int dim;
  preQuadtree* q0;
} quadtree;


void sorties_vtk(MMG5_pMesh mesh, int* connectionList, int nb_pt, char* fich)
{
  MMG5_pPoint ppt;
  int imax = nb_pt, jmax = nb_pt;

  FILE* fileVtk;

  if (!(fileVtk = fopen(fich,"wb")))
    fprintf(stdout,"sorties_vtk: unable to open fije %s\n", fich);

  // ----------------------Header-------------------------------------
  fprintf(fileVtk,"# vtk DataFile Version 3.1\n");
  fprintf(fileVtk,"Mesh nodes\n");
  fprintf(fileVtk,"ASCII\n");
  fprintf(fileVtk,"DATASET UNSTRUCTURED_GRID\n\n");

  //----------------------Nodes position------------------------------
  fprintf(fileVtk,"POINTS    %i    float\n", nb_pt);
  for (int i = 1; i<=nb_pt; i++)
  {
    ppt = &mesh->point[i];
    fprintf(fileVtk,"%g    %g    %g\n",ppt->c[0],ppt->c[1],ppt->c[2]);
  }

  fprintf(fileVtk,"\n");

  //---------------------Elements-------------------------------------
  fprintf(fileVtk,"CELLS    %i    %i\n", nb_pt, 3*nb_pt);
  for (int j = 1; j<=nb_pt; j++)
  {
    fprintf(fileVtk, "2    %i    %i\n", j-1, connectionList[j]-1);
  }

  fprintf(fileVtk,"\n");

  fprintf(fileVtk,"CELL_TYPES    %i\n",nb_pt);
  for (int j = 1; j<=nb_pt; j++)
  {
    fprintf(fileVtk,"3 ");
  }
  fprintf(fileVtk, "\n\n");

  //---------------------Data-----------------------------------------
  fprintf(fileVtk,"POINT_DATA    %i\n", nb_pt);
  fprintf(fileVtk,"SCALARS data float 1\n");
  fprintf(fileVtk,"LOOKUP_TABLE default\n");
  for (int j = 1; j<=nb_pt; j++)
  {
    fprintf(fileVtk, "%i\n",(int) (j == connectionList[j]));
  }

  fprintf(fileVtk,"\n");

  fclose(fileVtk);
}


void initPreQuadtree( preQuadtree* q)
{
  q->nbVer = 0;
  q->depth = 0;
  q->v = NULL;
  q->branches = NULL;
}

void initQuadtree(quadtree* q, int nv, int dim)
{
  q->nv = nv;
  q->dim = dim;
  q->q0 = (preQuadtree*) malloc(sizeof(preQuadtree));
  if (q->q0 == NULL)
    fprintf(stdout,"initQuadtree: allocation failed");
  initPreQuadtree(q->q0);
}

void freePreQuadtree(preQuadtree* q, int nv, int dim)
{
  if (q->nbVer>nv && q->depth < sizeof(int)*8/dim-1)
  {
    for (int i = 0; i<(1<<dim); i++)
      freePreQuadtree(&(q->branches[i]), nv, dim);
    free(q->branches);
  }
  if (q->nbVer>0 && (q->nbVer<= nv || q->depth == sizeof(int)*8/dim-1))
  {
    free(q->v);
  }
}

void freeQuadtree(quadtree* q)
{
  freePreQuadtree(q->q0, q->nv, q->dim);
  free(q->q0);
}


void start(struct timeval* tm1)
{
  gettimeofday(tm1, NULL);
}

void stop(struct timeval* tm1)
{
  struct timeval tm2;
  gettimeofday(&tm2, NULL);

  unsigned long long t = 1000 * (tm2.tv_sec - (*tm1).tv_sec) + (tm2.tv_usec - (*tm1).tv_usec) / 1000;
  printf("%llu ms\n", t);
}

// rect : x,y,z,l,h,p Very ugly way to do it but I couldn't find better
void countListSquareRec(preQuadtree* q, double* center, double* rect, int nv, int dim, int* index)
{
  if (q->branches==NULL)
  {
    (*index)++;
    //~ fprintf(stdout,"Indexe incrémenté %i\n", *index);
  }else
  {
    double recttemp[2*dim];
    double centertemp[dim];
    double l = 1./(1<<(q->depth+1));

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

        countListSquareRec(&(q->branches[0]), centertemp, recttemp, nv, dim, index);
      }
      if (rect[0]+rect[2] >center[0] && rect[1]<center[1]) // branch 1
      {
        recttemp[0] = rect[0]<center[0] ? center[0]:rect[0];
        recttemp[1] = rect[1];
        recttemp[2] = rect[0]+rect[2]-recttemp[0];
        recttemp[3] = rect[1]+rect[3] < center[1] ? rect[3]:center[1]-rect[1];

        centertemp[0] = center[0]+l/2;
        centertemp[1] = center[1]-l/2;

        countListSquareRec(&(q->branches[1]), centertemp, recttemp, nv, dim, index);
      }
      if (rect[0]<center[0] && rect[1]+rect[3]>center[1]) // branch 2
      {
        recttemp[0] = rect[0];
        recttemp[1] = rect[1]>center[1] ? rect[1]:center[1];
        recttemp[2] = rect[0]+rect[2] < center[0] ? rect[2]:center[0]-rect[0];
        recttemp[3] = rect[1] + rect[3]- recttemp[1];

        centertemp[0] = center[0]-l/2;
        centertemp[1] = center[1]+l/2;

        countListSquareRec(&(q->branches[2]), centertemp, recttemp, nv, dim, index);
      }
      if (rect[0]+rect[2] >center[0] && rect[1]+rect[3]>center[1]) // branch 3
      {
        recttemp[0] = rect[0]>center[0] ? rect[0]:center[0];
        recttemp[1] = rect[1]>center[1] ? rect[1]:center[1];
        recttemp[2] = rect[0]+rect[2]-recttemp[0];
        recttemp[3] = rect[1]+rect[3]-recttemp[1];

        centertemp[0] = center[0]+l/2;
        centertemp[1] = center[1]+l/2;

        countListSquareRec(&(q->branches[3]), centertemp, recttemp, nv, dim, index);
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

          countListSquareRec(&(q->branches[0]), centertemp, recttemp, nv, dim, index);
        }
        if (rect[2]+rect[5] > center[2]) // branch 4
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          countListSquareRec(&(q->branches[4]), centertemp, recttemp, nv, dim, index);
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

          countListSquareRec(&(q->branches[1]), centertemp, recttemp, nv, dim, index);
        }
        if (rect[2]+rect[5] > center[2]) // branch 5
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          countListSquareRec(&(q->branches[5]), centertemp, recttemp, nv, dim, index);
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

          countListSquareRec(&(q->branches[2]), centertemp, recttemp, nv, dim, index);
        }
        if (rect[2]+rect[5] >= center[2]) // branch 6
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          countListSquareRec(&(q->branches[6]), centertemp, recttemp, nv, dim, index);
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

          countListSquareRec(&(q->branches[3]), centertemp, recttemp, nv, dim, index);
        }
        if (rect[2]+rect[5] > center[2]) // branch 7
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          countListSquareRec(&(q->branches[7]), centertemp, recttemp, nv, dim, index);
        }
      }
    }
  }
}

// rect : x,y,z,l,h,p Very ugly way to do it but I couldn't find better
void getListSquareRec(preQuadtree* q, double* center, double* rect, preQuadtree*** qlist, int nv, int dim, int* index)
{
  if (q->branches==NULL)
  {
    (*qlist)[*index] = q;
    (*index)++;
    //~ fprintf(stdout,"Indexe incrémenté %i\n", *index);
  }else
  {
    double recttemp[2*dim];
    double centertemp[dim];
    double l = 1./(1<<(q->depth+1));

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

        getListSquareRec(&(q->branches[0]), centertemp, recttemp, qlist, nv, dim, index);
      }
      if (rect[0]+rect[2] >center[0] && rect[1]<center[1]) // branch 1
      {
        recttemp[0] = rect[0]<center[0] ? center[0]:rect[0];
        recttemp[1] = rect[1];
        recttemp[2] = rect[0]+rect[2]-recttemp[0];
        recttemp[3] = rect[1]+rect[3] < center[1] ? rect[3]:center[1]-rect[1];

        centertemp[0] = center[0]+l/2;
        centertemp[1] = center[1]-l/2;

        getListSquareRec(&(q->branches[1]), centertemp, recttemp, qlist, nv, dim, index);
      }
      if (rect[0]<center[0] && rect[1]+rect[3]>center[1]) // branch 2
      {
        recttemp[0] = rect[0];
        recttemp[1] = rect[1]>center[1] ? rect[1]:center[1];
        recttemp[2] = rect[0]+rect[2] < center[0] ? rect[2]:center[0]-rect[0];
        recttemp[3] = rect[1] + rect[3]- recttemp[1];

        centertemp[0] = center[0]-l/2;
        centertemp[1] = center[1]+l/2;

        getListSquareRec(&(q->branches[2]), centertemp, recttemp, qlist, nv, dim, index);
      }
      if (rect[0]+rect[2] >center[0] && rect[1]+rect[3]>center[1]) // branch 3
      {
        recttemp[0] = rect[0]>center[0] ? rect[0]:center[0];
        recttemp[1] = rect[1]>center[1] ? rect[1]:center[1];
        recttemp[2] = rect[0]+rect[2]-recttemp[0];
        recttemp[3] = rect[1]+rect[3]-recttemp[1];

        centertemp[0] = center[0]+l/2;
        centertemp[1] = center[1]+l/2;

        getListSquareRec(&(q->branches[3]), centertemp, recttemp, qlist, nv, dim, index);
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

          getListSquareRec(&(q->branches[0]), centertemp, recttemp, qlist, nv, dim, index);
        }
        if (rect[2]+rect[5] > center[2]) // branch 4
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          getListSquareRec(&(q->branches[4]), centertemp, recttemp, qlist, nv, dim, index);
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

          getListSquareRec(&(q->branches[1]), centertemp, recttemp, qlist, nv, dim, index);
        }
        if (rect[2]+rect[5] > center[2]) // branch 5
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          getListSquareRec(&(q->branches[5]), centertemp, recttemp, qlist, nv, dim, index);
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

          getListSquareRec(&(q->branches[2]), centertemp, recttemp, qlist, nv, dim, index);
        }
        if (rect[2]+rect[5] >= center[2]) // branch 6
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          getListSquareRec(&(q->branches[6]), centertemp, recttemp, qlist, nv, dim, index);
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

          getListSquareRec(&(q->branches[3]), centertemp, recttemp, qlist, nv, dim, index);
        }
        if (rect[2]+rect[5] > center[2]) // branch 7
        {
          recttemp[2] = rect[2]<center[2] ? center[2]:rect[2];
          recttemp[5] = rect[2]+ rect[5] - recttemp[2];

          centertemp[2] = center[2]+l/2;

          getListSquareRec(&(q->branches[7]), centertemp, recttemp, qlist, nv, dim, index);
        }
      }
    }
  }
}

int getListSquare(quadtree* q, double* rect, preQuadtree*** qlist)
{
  int dim = q->dim;
  double rect2[dim*2];
  memcpy(rect2, rect, sizeof(double)*dim*2);
  double center[dim];
  for (int i = 0; i < dim; i++)
    center[i] = 0.5;
  int index=0;
  //~ fprintf(stdout,"q nulle ? %i\n", q->branches == NULL);
  countListSquareRec(q->q0,center, rect2, q->nv, q->dim, &index);
  //~ fprintf(stdout, "index : %i\n", index);
  *qlist = (preQuadtree**) malloc(sizeof(preQuadtree*)*index);
  for (int i = 0; i<index; i++)
    (*qlist)[i] = NULL;
  //~ if (qlist == NULL)
  //~ fprintf(stdout, "niuratisecuinserirt");
  index = 0;
  memcpy(rect2, rect, sizeof(double)*dim*2);
  for (int i = 0; i < dim; i++)
    center[i] = 0.5;

  getListSquareRec(q->q0, center, rect2, qlist, q->nv, q->dim, &index);

  return index;
}

void addVertexRec3d(MMG5_pMesh mesh, preQuadtree* q, double* ver, const int no, int nv, int dim)
{

  if ( q->depth < sizeof(int)*8/dim -1 ) // size of int in bits divided by dimension
  {
    if (q->nbVer < nv)
    {
      //~ fprintf(stdout, "ok1 %f %f \n", ver[0], ver[1]);
      if (q->v == NULL)
      {
        //~ fprintf(stdout, "ok1bis\n");
        q->v = (int *) malloc(nv*sizeof(int));
        if (q->v ==NULL)
          fprintf(stdout, "addVertexRec3d: allocation branches failed at 1\n");
      }
      q->v[q->nbVer] = no; // on peut faire mieux en triant les vertices ajoutés, ça permet ensuite de faire directement descendre les points à leur place lorsque l'on dépasse le nombre de points dans la feuille et de rechercher plus rapidement dans la feuille
      // intéressant seulement si nv est grand
      q->nbVer++;
      return;
    }else if (q->nbVer == nv && q->branches==NULL) //creation of sub-branch and relocation of vertices in the sub-branches
    {
      //~ fprintf(stdout, "ok2\n");
      q->branches = (preQuadtree*) malloc((1<<dim)*sizeof(preQuadtree));
      if (q->branches==NULL)
        fprintf(stdout, "addVertexRec3d: allocation branches failed at 2\n");
      for (int i = 0; i<(1<<dim); i++)
      {
        initPreQuadtree(&(q->branches[i]));
        q->branches[i].depth = q->depth+1;
      }
      //~ unsigned char quadrant = ver[0]>0.5+2*(ver[1]>0.5);
      //~ if ((q->branches[quadrant].v = (int *) malloc(nv*sizeof(int)))==NULL)
      //~ fprintf(stdout, "Probleme allocation vertex\n");
      //~ q->branches[quadrant].v[0] = no;
      //~ q->branches[quadrant].nbVer++;
      q->nbVer++;
      for (int i = 0; i<nv; i++)
      {
        double pt[dim];

        memcpy(pt,&mesh->point[q->v[i]] ,dim*sizeof(double));
        for (int j =0; j < q->depth; j++)
        {
          for (int k = 0; k<dim; k++)
          {
            pt[k] -= ((float) (pt[k]>0.5))*0.5;
            pt[k] *= 2;
          }
        }
        //~ addVertexRec(mesh, q, mesh->point[q->v[i]].c, q->v[i]);
        addVertexRec3d(mesh, q, pt, q->v[i],nv,dim);
        q->nbVer--;
      }
      addVertexRec3d(mesh, q, ver, no, nv, dim);
      q->nbVer--;
      free(q->v);
      q->v = NULL;
    }else
    {
      int quadrant = 0;

      for (int i = 0; i<dim; i++)
      {
        quadrant += ((float) (ver[i]>0.5))*(1<<i);
        ver[i] -= ((float) (ver[i]>0.5))*0.5;
        ver[i] *= 2;
      }

      q->nbVer++;
      addVertexRec3d(mesh, &(q->branches[quadrant]), ver, no, nv, dim);
    }
  }else
  {

    if (q->nbVer%nv == 0)
      q->v = (int*) realloc(q->v, sizeof(int)*nv*(q->nbVer/nv+1));
    if (q->v == NULL)
    {
      fprintf(stdout,"addVertex : Reallocation failed\n");
      return;
    }
    q->v[q->nbVer] = no;
    q->nbVer++;
  }
}

void addVertex3d(MMG5_pMesh mesh, quadtree* q, const int no)
{
  double pt[q->dim];

  memcpy(pt, &mesh->point[no] ,q->dim*sizeof(double));
  //~ fprintf(stdout, "ajout du point %f %f %f\n", pt[0], pt[1], pt[2]);
  addVertexRec3d(mesh, q->q0, pt , no, q->nv, q->dim);
}

void printArbreDepth(preQuadtree* q, int depth, int nv, int dim)
{
  if ( q->depth < depth && q->nbVer > nv)
  {

    for (int i = 0; i < (1<<dim); i++)
    {
      printArbreDepth(&(q->branches[i]),depth, nv, dim);
    }
    //~ fprintf(stdout," | ");
  }else if (q->depth == depth)
  {
    fprintf(stdout,"%i ",q->nbVer);
  }
}
void printArbre(quadtree* q)
{
  for (int i = 0; i<sizeof(int)*8/q->dim; i++)
  {
    fprintf(stdout,"\n profondeur %i \n", i);
    printArbreDepth(q->q0, i, q->nv, q-> dim);

  }
  fprintf(stdout,"\n fin \n");
}

int sizeArbreRec(preQuadtree* q, int nv, int dim)
{

  if (q->branches != NULL)
  {
    int sizeBranches = 0;
    for (int i= 0; i <4; i++)
    {
      sizeBranches += sizeArbreRec(&(q->branches[i]),nv,dim)+sizeof(preQuadtree)+(1<<dim)*sizeof(preQuadtree*);
    }
    return sizeBranches;
  }else if(q->v != NULL)
  {
    return nv*sizeof(int)+sizeof(preQuadtree);
  }else
  {
    return sizeof(preQuadtree);
  }
}

int sizeArbre(quadtree* q)
{
  return sizeArbreRec(q->q0, q->nv, q->dim);
}

int sizeArbreLinkRec(preQuadtree* q, int nv, int dim)
{
  if (q->branches != NULL)
  {
    int sizeBranches = 0;
    for (int i= 0; i <4; i++)
    {
      sizeBranches += sizeArbreLinkRec(&(q->branches[i]), nv, dim)+sizeof(preQuadtree)+(1<<dim)*sizeof(preQuadtree*);
    }
    return sizeBranches;
  }else if(q->v != NULL)
  {
    return sizeof(int)+sizeof(preQuadtree);
  }else
  {
    return sizeof(preQuadtree);
  }
}

int sizeArbreLink(quadtree* q)
{
  return sizeArbreLinkRec(q->q0, q->nv, q->dim)+q->q0->nbVer*sizeof(int);
}
//~
//~ void drawSquare(MMG2D_postscript* ps, double x, double y, double l)
//~ {
//~ double x1 = x*ps->box_width+ps->point0[0];
//~ double y1 = y*ps->box_width+ps->point0[1];
//~ double l1 = l*ps->box_width;
//~ fprintf(ps->file,"%.6g %.6g %.6g %.6g rectstroke\n", x1, y1, l1, l1);
//~ }

//~ void plotQuadtreePointRec(MMG2D_postscript* ps, MMG5_pMesh mesh, preQuadtree* q, double lmoy, double x, double y)
//~ {
//~ double dl = 1./(1<<(q->depth+1));
//~ double color[3];
//~ color[0] = 0;
//~ color[1] = 0;
//~ color[2] = 0;
//~ double size;
//~ if (q->branches != NULL)
//~ {
//~ for (int i= 0; i <4; i++)
//~ {
//~ plotQuadtreePointRec(ps, mesh, &(q->branches[i]), lmoy, x+(i%2)*dl, y+(i/2)*dl);
//~ }
//~ }else if(q->v != NULL)
//~ {
//~ size = lmoy/3<dl/q->nbVer ? lmoy/3:dl/q->nbVer;
//~ MMG2D_sethueColor(ps, 0.);
//~ for (int i = 0; i<q->nbVer; i++)
//~ {
//~ if (q->v[i] == 25)
//~ {
//~ MMG2D_sethueColor(ps, 0.6);
//~ MMG2D_drawPoint(mesh, ps, q->v[i], size*2);
//~ MMG2D_sethueColor(ps, 0.);
//~ }else
//~ MMG2D_drawPoint(mesh, ps, q->v[i], size);
//~ }
//~ MMG2D_sethsbColor(ps, color);
//~ for (int i = 0; i<q->nbVer; i++)
//~ {
//~ fprintf(ps->file,"%.6g %.6g %i %i (%i) outputtext\n",  mesh->point[q->v[i]].c[0]*ps->box_width+ps->point0[0], mesh->point[q->v[i]].c[1]*ps->box_width+ps->point0[1], 2, 0, q->v[i]);
//~ }
//~ }
//~ }


//~ void plotQuadtreeSquareRec(MMG2D_postscript* ps, MMG5_pMesh mesh, preQuadtree* q, double lmoy, double x, double y)
//~ {
//~ double dl = 1./(1<<(q->depth+1));
//~ if (q->branches != NULL)
//~ {
//~ for (int i= 0; i <4; i++)
//~ {
//~ plotQuadtreeSquareRec(ps, mesh, &(q->branches[i]), lmoy, x+(i%2)*dl, y+(i/2)*dl);
//~ }
//~ }else
//~ {
//~ drawSquare(ps, x, y, 2*dl);
//~ }
//~ }

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

int NearNeighborSquare(MMG5_pMesh mesh, quadtree* q, int no, double l)
{
  MMG5_pPoint   ppt,ppt1;
  preQuadtree** qlist;
  int ns,nver;
  double rect[2*q->dim];
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
  ns = getListSquare(q, rect, &qlist);

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

  free(qlist);

  if (sqrt(lmin)<l)
    return nmin;
  else
    return -1;
}


//~
//~ void plotQuadtree( MMG5_pMesh mesh, quadtree* q, char* name)
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
//~ plotQuadtreePointRec(&ps, mesh, q->q0, lmoy, 0, 0);
//~ MMG2D_sethsbColor(&ps, color);
//~ MMG2D_setLineWidth(&ps, 0.0001);
//~ plotQuadtreeSquareRec(&ps, mesh, q->q0, lmoy, 0, 0);
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
//~ fprintf(stdout,"temps preQuadtree : ");
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

