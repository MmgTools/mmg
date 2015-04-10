/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
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
 * \file common/mmg_isosiz.c
 * \brief Fonctions for isotropic size map computation.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param ip1 index of edge's extremity.
 * \param ip2 index of edge's extremity.
 * \param isedg 1 if the edge is a ridge, 0 otherwise.
 * \return length of edge according to the prescribed metric.
 *
 * Compute length of edge \f$[i0;i1]\f$ according to the prescribed iso.
 * metric.
 *
 */
inline double _MMG5_lenedg_iso(MMG5_pMesh mesh,MMG5_pSol met,int ip1,int ip2, char isedg) {
  MMG5_pPoint   p1,p2;
  double   h1,h2,l,r,len;

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  h1 = met->m[ip1];
  h2 = met->m[ip2];
  l = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]) \
    + (p2->c[2]-p1->c[2])*(p2->c[2]-p1->c[2]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;
  len = fabs(r) < _MMG5_EPS ? l / h1 : l / (h2-h1) * log(r+1.0);

  return(len);
}

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
double _MMG5_surftri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) {
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

  return( 0.5*sqrt(det) );
}
