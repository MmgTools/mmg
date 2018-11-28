/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
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
 * \file common/isosiz.c
 * \brief Fonctions for isotropic size map computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the meric structure.
 * \param ptt pointer toward the triangle structure.
 * \return The computed area.
 *
 * Compute the area of the surface triangle \a ptt with respect to
 * the isotropic metric \a met.
 *
 */
double MMG5_surftri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) {
  double   *a,*b,*c,abx,aby,abz,acx,acy,acz,det,n[3];

  a = mesh->point[ptt->v[0]].c;
  b = mesh->point[ptt->v[1]].c;
  c = mesh->point[ptt->v[2]].c;

  /* area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];

  n[0] = aby*acz - abz*acy;
  n[1] = abz*acx - abx*acz;
  n[2] = abx*acy - aby*acx;
  det  = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];

  return  0.5*sqrt(det) ;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param funcname name of the calling function
 *
 * \return 1 if success, 0 if fail.
 *
 * Print that we enter in the defsiz function in high verbosity level and check
 * the hmax value.
 *
 */
int MMG5_defsiz_startingMessage (MMG5_pMesh mesh,MMG5_pSol met,const char * funcname ) {

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Defining %stropic map\n",(met->size==1)?"iso":"aniso");

  if ( mesh->info.hmax < 0.0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative hmax value.\n",funcname);
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ip0 index of the first edge extremity
 * \param ip1 index of the second edge extremity
 *
 * \return 1 if success, 0 if fail.
 *
 * Compute the euclidean length of the edge \a ip0 \a ip1,
 * add this length to the metric of the edge extremities and
 * increment the count of times we have processed this extremities.
 *
 */
int MMG5_sum_reqEdgeLengthsAtPoint ( MMG5_pMesh mesh,MMG5_pSol met,int ip0,int ip1 ) {
  MMG5_pPoint p0,p1;
  int         j;
  double      len,dist;

  /* Compute the euclidean edge length */
  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];

  len = 0.;
  for ( j=0; j<mesh->dim; ++j ) {
    dist = p1->c[j]-p0->c[j];
    len += dist*dist;
  }
  len = sqrt(len);

  /* Add the length to the point's metric and increment the number of
   * times the point has been seen */
  met->m[met->size*ip0] += len;
  met->m[met->size*ip1] += len;
  ++p0->s;
  ++p1->s;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Compute the mean metric at mesh points with a non-nul \a s field. At the
 * beginning, the metric of a given point contains the sum of n metrics and the
 * \a s field of the point the number of metrics summed in the point. Set the
 * flag of the processed points to 3.
 *
 */
int MMG5_compute_meanMetricAtMarkedPoints_iso ( MMG5_pMesh mesh,MMG5_pSol met ) {
  MMG5_pPoint p0;
  int         k;

  for ( k=1; k<=mesh->np; k++ ) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) )  continue;

    if ( !p0->s ) continue;

    met->m[k] /= p0->s;
    p0->flag = 3;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ismet 1 if the user has given an input metric (so we need to reset it)
 *
 * \return 1 if success, 0 if fail.
 *
 * For a triangle mesh, process the triangles and set to 0 the metrics at points
 * that are at the extremities of a required edge.
 *
 */
int MMG5_reset_metricAtReqEdges_surf ( MMG5_pMesh mesh,MMG5_pSol met,int8_t ismet ) {
  MMG5_pPoint p0;
  MMG5_pTria  pt;
  int         k,i,j,ip0,ip1;

  if ( ismet ) {
    for ( k=1; k<=mesh->nt; k++ ) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for ( i=0; i<3; ++i ) {
        if ( (pt->tag[i] & MG_REQ) || (pt->tag[i] & MG_NOSURF) ||
             (pt->tag[i] & MG_PARBDY) ) {

          ip0 = pt->v[MMG5_iprv2[i]];
          ip1 = pt->v[MMG5_inxt2[i]];

          for ( j=0; j<met->size; ++j ) {

            met->m[ met->size*ip0 + j ] = 0.;
            met->m[ met->size*ip1 + j ] = 0.;
          }
        }
      }
    }
  }

  return 1;
}
