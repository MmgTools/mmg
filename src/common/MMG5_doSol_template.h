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
