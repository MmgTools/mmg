/**
 * \file common/MMG5_doSol_template.h
 * \brief common template for doSol functions
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \date 12 2021
 * \copyright GNU Lesser General Public License.
 **/

#define MMG_PASTER(x,y) x ## _ ## y
#define MMG_EVALUATOR(x,y)  MMG_PASTER(x,y)
#define MMG_FUNCTION_NAME(a,b) MMG_EVALUATOR(a, b)

static inline
int MMG_FUNCTION_NAME(MMG5_doSol_iso, MMG5_ELEMENT_TYPE)
  (MMG5_pMesh mesh,MMG5_pSol met,
   int (*Set_solSize)(MMG5_pMesh,MMG5_pSol,int,int,int),
   int ne,MMG5_ELEMENT_TYPE *elt,int8_t na,const uint8_t iare[][2]) {

  MMG5_ELEMENT_TYPE *pt;
  MMG5_pPoint  p1,p2;
  double       ux,uy,uz,dd;
  int          i,k,ia,ib,ipa,ipb,type;
  int         *mark;

  MMG5_SAFE_CALLOC(mark,mesh->np+1,int,return 0);

  /* Memory alloc */
  if ( met->size!=1 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,met->size);
    return 0;
  }

  type=1;
  if ( !Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
    return 0;

  /* Travel the triangles edges and add the edge contribution to edges
   * extermities */
  for (k=1; k<=ne; k++) {
    pt = &elt[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<na; i++) {
      ia  = iare[i][0];
      ib  = iare[i][1];
      ipb = pt->v[ib];
      p1  = &mesh->point[ipa];
      p2  = &mesh->point[ipb];

      ux  = p1->c[0] - p2->c[0];
      uy  = p1->c[1] - p2->c[1];
      uz  = p1->c[2] - p2->c[2];
      dd  = sqrt(ux*ux + uy*uy + uz*uz);

      met->m[ipa] += dd;
      mark[ipa]++;
      met->m[ipb] += dd;
      mark[ipb]++;
    }
  }

  /* if hmax is not specified, compute it from the metric */
  if ( mesh->info.hmax < 0. ) {
    dd = 0.;
    for (k=1; k<=mesh->np; k++) {
      if ( !mark[k] ) continue;
      dd = MG_MAX(dd,met->m[k]);
    }
    assert ( dd );
    mesh->info.hmax = 10.*dd;
  }

  /* vertex size */
  for (k=1; k<=mesh->np; k++) {
    if ( !mark[k] ) {
      met->m[k] = mesh->info.hmax;
      continue;
    }
    else
      met->m[k] = met->m[k] / (double)mark[k];
  }

  MMG5_SAFE_FREE(mark);
  return 1;
}

static inline
int MMG_FUNCTION_NAME(MMG5_doSol_ani, MMG5_ELEMENT_TYPE)
  (MMG5_pMesh mesh,MMG5_pSol met,
   int (*Set_solSize)(MMG5_pMesh,MMG5_pSol,int,int,int),
   int ne,MMG5_ELEMENT_TYPE *elt,int8_t na,const uint8_t iare[][2]) {

  MMG5_ELEMENT_TYPE *pt;
  MMG5_pPoint  p1,p2;
  double       u[3],dd,tensordot[6];
  int          i,j,k,iadr,ipa,ipb,type;
  int         *mark;

  MMG5_SAFE_CALLOC(mark,mesh->np+1,int,return 0);

  /* Memory alloc */
  if ( met->size!=6 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,met->size);
    return 0;
  }

  type = 3;
  if ( !Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
    return 0;

  /* Travel the triangles edges and add the edge contribution to edges
   * extermities */
  for (k=1; k<=ne; k++) {
    pt = &elt[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<na; i++) {
      ipa = pt->v[iare[i][0]];
      ipb = pt->v[iare[i][1]];
      p1  = &mesh->point[ipa];
      p2  = &mesh->point[ipb];

      u[0]  = p1->c[0] - p2->c[0];
      u[1]  = p1->c[1] - p2->c[1];
      u[2]  = p1->c[2] - p2->c[2];

      tensordot[0] = u[0]*u[0];
      tensordot[1] = u[0]*u[1];
      tensordot[2] = u[0]*u[2];
      tensordot[3] = u[1]*u[1];
      tensordot[4] = u[1]*u[2];
      tensordot[5] = u[2]*u[2];

      iadr = 6*ipa;
      for ( j=0; j<6; ++j ) {
        met->m[iadr+j]   += tensordot[j];
      }
      mark[ipa]++;

      iadr = 6*ipb;
      for ( j=0; j<6; ++j ) {
        met->m[iadr+j]   += tensordot[j];
      }
      mark[ipb]++;
    }
  }

  /* if hmax is not specified, compute it from the metric */
  double hmax = FLT_MAX;
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !mark[k] ) {
      continue;
    }

    /* Metric = nedges/dim * inv (sum(tensor_dot(edges,edges))).
     * sum(tensor_dot) is stored in sol->m so reuse tensordot to
     * compute M.  */
    iadr = 6*k;
    if ( !MMG5_invmat(met->m+iadr,tensordot) ) {
      /* Non invertible matrix: impose FLT_MIN, it will be truncated by hmax
       * later */
      fprintf(stdout, " ## Warning: %s: %d: non invertible matrix."
             " Impose hmax size at point\n",__func__,__LINE__);
      met->m[iadr+0] = FLT_MIN;
      met->m[iadr+1] = 0;
      met->m[iadr+2] = 0;
      met->m[iadr+3] = FLT_MIN;
      met->m[iadr+4] = 0;
      met->m[iadr+5] = FLT_MIN;
      continue;
    }

    dd = (double)mark[k]*0.5;

    for ( j=0; j<6; ++j ) {
      met->m[iadr+j] = dd*tensordot[j];
    }

    /* Check metric */
    double lambda[3],vp[3][3];
    if (!MMG5_eigenv(1,met->m+iadr,lambda,vp) ) {
      fprintf(stdout, " ## Warning: %s: %d: non diagonalizable metric."
              " Impose hmax size at point\n",__func__,__LINE__);
      met->m[iadr+0] = FLT_MIN;
      met->m[iadr+1] = 0;
      met->m[iadr+2] = 0;
      met->m[iadr+3] = FLT_MIN;
      met->m[iadr+4] = 0;
      met->m[iadr+5] = FLT_MIN;
      continue;
    }

    assert ( lambda[0] > 0. && lambda[1] > 0.  && lambda[2] > 0.
            && "Negative eigenvalue");

    /* If the vertex belongs to only colinear edges one of the eigenvalue is
     * infinite: do not take it into account, it will be truncated by hmax
     * later */
    for ( j=0; j<3; ++j ) {
      if ( isfinite(lambda[j]) ) {
        hmax = MG_MIN(hmax,lambda[0]);
      }
    }
  }
  if ( mesh->info.hmax < 0.) {
    assert ( hmax > 0. && hmax < FLT_MAX && "Wrong hmax value" );
    mesh->info.hmax = 10./sqrt(hmax);
  }

  /* vertex size at orphan points: impose hmax size */
  hmax = 1./(mesh->info.hmax*mesh->info.hmax);
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !mark[k] ) {
      iadr = 6*k;
      met->m[iadr]   = hmax;
      met->m[iadr+3] = met->m[iadr];
      met->m[iadr+5] = met->m[iadr];
    }
  }

  MMG5_SAFE_FREE(mark);
  return 1;
}
