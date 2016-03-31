/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
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
 * \file mmg2d/intmet_2d.c
 * \brief Interpolation of metrics
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "mmg2d.h"

/* Interpolation of metric met along edge i of triangle k, according to parameter s; 
   ip = index of the new point */
int _MMG2_intmet_iso(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s) {
  MMG5_pTria  pt;
  int    ip1,ip2;
  char   i1,i2;
  
  pt  = &mesh->tria[k];
  i1  = _MMG5_inxt2[i];
  i2  = _MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  met->m[ip] = s * (met->m[ip1] + met->m[ip2]);
  
  return(1);
}
