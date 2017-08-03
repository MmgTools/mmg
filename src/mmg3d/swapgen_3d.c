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
 * \file mmg3d/swapgen_3d.c
 * \brief Functions for swapping process inside the mesh.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "inlined_functions_3d.h"

/**
 * \param mesh pointer toward the mesh structure
 * \param met pointer toward the metric structure.
 * \param start tetrahedra in which the swap should be performed
 * \param ia edge that we want to swap
 * \param ilist pointer to store the size of the shell of the edge
 * \param list pointer to store the shell of the edge
 * \param crit improvment coefficient
 * \return 0 if fail, the index of point corresponding to the swapped
 * configuration otherwise (\f$4*k+i\f$).
 * \param typchk type of checking permformed for edge length (hmin or LSHORT
 * criterion).
 *
 * Check whether swap of edge \a ia in \a start should be performed, and
 * return \f$4*k+i\f$ the index of point corresponding to the swapped
 * configuration. The shell of edge is built during the process.
 *
 */
int _MMG5_chkswpgen(MMG5_pMesh mesh,MMG5_pSol met,int start,int ia,
                    int *ilist,int *list,double crit,char typchk) {
  MMG5_pTetra    pt,pt0;
  MMG5_pPoint    p0;
  double    calold,calnew,caltmp;
  int       na,nb,np,adj,piv,npol,refdom,k,l,iel;
  int       *adja,pol[MMG3D_LMAX+2];
  char      i,ipa,ipb,ip,ier,ifac;

  pt  = &mesh->tetra[start];
  refdom = pt->ref;

  pt0 = &mesh->tetra[0];
  na  = pt->v[_MMG5_iare[ia][0]];
  nb  = pt->v[_MMG5_iare[ia][1]];
  calold = pt->qual;

  /* Store shell of ia in list, and associated pseudo polygon in pol */
  (*ilist) = 0;
  npol = 0;
  list[(*ilist)] = 6*start+ia;
  (*ilist)++;
  adja = &mesh->adja[4*(start-1)+1];
  adj  = adja[_MMG5_ifar[ia][0]];      // start travelling by face (ia,0)
  ifac = adj%4;
  piv  = pt->v[_MMG5_ifar[ia][1]];
  pol[npol] = 4*start + _MMG5_ifar[ia][1];
  npol++;

  /* Edge is on a boundary between two different domains */
  if ( mesh->info.opnbdy )
    if ( pt->xt && (mesh->xtetra[pt->xt].ftag[_MMG5_ifar[ia][1]] & MG_BDY) )
      return 0;

  while ( adj ) {
    adj /= 4;
    if ( adj ==start ) break;

    pt = &mesh->tetra[adj];
    if ( pt->tag & MG_REQ ) return(0);

    /* Edge is on a boundary between two different domains */
    if ( pt->ref != refdom )  return(0);
    else if ( mesh->info.opnbdy ) {
      if ( pt->xt && (mesh->xtetra[pt->xt].ftag[ifac] & MG_BDY) ) return 0;
    }

    calold = MG_MIN(calold, pt->qual);
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = _MMG5_iare[i][0];
      ipb = _MMG5_iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    list[(*ilist)] = 6*adj +i;
    (*ilist)++;
    /* overflow */
    if ( (*ilist) > MMG3D_LMAX-3 )  return(0);

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
      pol[npol] = 4*adj + _MMG5_ifar[i][1];
      npol++;
      adj = adja[ _MMG5_ifar[i][0] ];
      piv = pt->v[ _MMG5_ifar[i][1] ];
    }
    else {
      assert(pt->v[ _MMG5_ifar[i][1] ] == piv);
      pol[npol] = 4*adj + _MMG5_ifar[i][0];
      npol++;
      adj = adja[ _MMG5_ifar[i][1] ];
      piv = pt->v[ _MMG5_ifar[i][0] ];
    }
    ifac = adj%4;
  }
  //CECILE : je vois pas pourquoi ca ameliore de faire ce test
  //plus rapide mais du coup on elimine des swap...
  //4/01/14 commentaire
  // if ( calold*_MMG3D_ALPHAD > 0.5 )  return(0);

  /* Prevent swap of an external boundary edge */
#warning Profiling : why we don t check the edge tag instead of computing the shell
   if ( !adj ) return(0);

  assert(npol == (*ilist)); // du coup, apres on pourra virer npol

  /* Find a configuration that enhances the worst quality within the shell */
  for (k=0; k<npol; k++) {
    iel = pol[k] / 4;
    ip  = pol[k] % 4;
    np  = mesh->tetra[iel].v[ip];
    calnew = 1.0;
    ier = 1;

    if ( mesh->info.fem ) {
      p0 = &mesh->point[np];
      if ( p0->tag & MG_BDY ) {
        for (l=0; l<npol;l++) {
          if ( k < npol-1 ) {
            if ( l == k || l == k+1 )  continue;
          }
          else {
            if ( l == npol-1 || l == 0 )  continue;
          }
          iel = pol[l] / 4;
          ip  = pol[l] % 4;
          pt = &mesh->tetra[iel];
          p0 = &mesh->point[pt->v[ip]];
          if ( p0->tag & MG_BDY ) {
            ier = 0;
            break;
          }
        }
      }
      if ( !ier )  continue;
      ier = 1;
    }

    for (l=0; l<(*ilist); l++) {
      /* Do not consider tets of the shell of collapsed edge */
      if ( k < npol-1 ) {
        if ( l == k || l == k+1 )  continue;
      }
      else {
        if ( l == npol-1 || l == 0 )  continue;
      }
      iel = list[l] / 6;
      i   = list[l] % 6;
      pt  = &mesh->tetra[iel];

      /* Prevent from creating a tetra with 4 bdy vertices */
      if ( mesh->point[np].tag & MG_BDY ) {
        if ( ( mesh->point[pt->v[_MMG5_ifar[i][0]]].tag & MG_BDY ) &&
             ( mesh->point[pt->v[_MMG5_ifar[i][1]]].tag & MG_BDY ) ) {
          if ( ( mesh->point[pt->v[_MMG5_iare[i][0]]].tag & MG_BDY ) ||
               ( mesh->point[pt->v[_MMG5_iare[i][1]]].tag & MG_BDY ) ) {
            ier = 0;
            break;
          }
        }
      }

      /* First tetra obtained from iel */
      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      pt0->v[_MMG5_iare[i][0]] = np;


      if ( typchk==1 && met->size > 1 && met->m )
        caltmp = _MMG5_caltet33_ani(mesh,met,pt0);
      else
        caltmp = _MMG5_orcal(mesh,met,0);

      calnew = MG_MIN(calnew,caltmp);

      ier = (calnew > crit*calold);
      if ( !ier )  break;

      /* Second tetra obtained from iel */
      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      pt0->v[_MMG5_iare[i][1]] = np;

      if ( typchk==1 && met->size > 1 && met->m )
        caltmp = _MMG5_caltet33_ani(mesh,met,pt0);
      else
        caltmp = _MMG5_orcal(mesh,met,0);

      calnew = MG_MIN(calnew,caltmp);

      ier = (calnew > crit*calold);
      if ( !ier )  break;
    }
    if ( ier )  return(pol[k]);
  }
  return(0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param nconf configuration.
 * \param ilist number of tetrahedra in the shell of the edge that we want
 *  to swap.
 * \param list pointer toward the shell of the edge that we want to swap.
 * \param octree pointer toward the octree structure in Delaunay mode,
 * NULL pointer in pattern mode.
 * \param typchk type of checking permformed for edge length (hmin or LSHORT
 * criterion).
 * \return -1 if lack of memory, 0 if fail to swap, 1 otherwise.
 *
 * Perform swap of edge whose shell is passed according to configuration nconf.
 *
 */
int _MMG5_swpgen(MMG5_pMesh mesh,MMG5_pSol met,int nconf,int ilist,int *list,
                 _MMG3D_pOctree octree, char typchk) {
  MMG5_pTetra    pt;
  MMG5_pPoint    p0,p1;
  int       iel,na,nb,np,nball,ret,start;
  double    m[3];
  char      ia,ip,iq;
  int       ier;

  iel = list[0] / 6;
  ia  = list[0] % 6;

  pt = &mesh->tetra[iel];
  na = pt->v[_MMG5_iare[ia][0]];
  nb = pt->v[_MMG5_iare[ia][1]];
  p0 = &mesh->point[na];
  p1 = &mesh->point[nb];

  /* Temporarily create midpoint at swapped edge */
  m[0] = 0.5*(p0->c[0] + p1->c[0]);
  m[1] = 0.5*(p0->c[1] + p1->c[1]);
  m[2] = 0.5*(p0->c[2] + p1->c[2]);

  np  = _MMG3D_newPt(mesh,m,0);
  if(!np){
    _MMG5_POINT_REALLOC(mesh,met,np,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new point\n",__func__);
                        _MMG5_INCREASE_MEM_MESSAGE();
                        return(-1)
                        ,m,0,-1);
  }
  if ( met->m ) {
    if ( typchk == 1 && (met->size>1) ) {
      if ( _MMG3D_intmet33_ani(mesh,met,iel,ia,np,0.5)<=0 )  return(0);
    }
    else {
      if ( _MMG5_intmet(mesh,met,iel,ia,np,0.5)<=0 ) return(0);
    }
  }

  /** First step : split of edge (na,nb) */
  ret = 2*ilist + 0;
  ier = _MMG5_split1b(mesh,met,list,ret,np,0,typchk-1);
  /* pointer adress may change if we need to realloc memory during split */
  pt = &mesh->tetra[iel];

  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Warning: %s: unable to swap internal edge.\n",
      __func__);
    return(-1);
  }
  else if ( !ier )  {
    _MMG3D_delPt(mesh,np);
    return(0);
  }

  /** Second step : collapse of np towards enhancing configuration */
  start = nconf / 4;
  iq = nconf % 4;

  pt = &mesh->tetra[start];
  for (ip=0; ip<4; ip++) {
    if ( pt->v[ip] == np )  break;
  }
  assert(ip<4);

  memset(list,0,(MMG3D_LMAX+2)*sizeof(int));
  nball = _MMG5_boulevolp(mesh,start,ip,list);

  ier = _MMG5_colver(mesh,met,list,nball,iq,typchk);
  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Warning: %s: unable to swap internal edge.\n",
      __func__);
    return(-1);
  }
  else if ( ier ) _MMG3D_delPt(mesh,ier);

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param k index of the tetrahedron with multiple boundary faces (to be swapped).
 * \param metRidTyp metric storage (classic or special)
 * \return -1 if lack of memory, 0 if fail to swap, 1 otherwise.
 *
 * Search an adjacent to the tetra \a k and perform swap 2->3 (the common face
 * of the 2 tetra is destroyed and replaced by a common edge used by the three
 * new elts).
 *
 * \remark used in anatet4 to remove the tetra with multiple boundary faces.
 *
 */
int MMG3D_swap23(MMG5_pMesh mesh,MMG5_pSol met,int k,char metRidTyp) {
  MMG5_pTetra          pt0,pt1,ptnew;
  MMG5_xTetra          xt[3];
  MMG5_pxTetra         pxt0,pxt1;
  double               vold0,vold,vnew;
  int                  k1,conf0,conf1,*adja,iel,np,xt1;
  int                  adj0_2,adj0_3,adj1_1,adj1_2,adj1_3;
  char                 j0,j1,i,isxt[3];
  unsigned char        tau0[4],tau1[4];
  const unsigned char *taued0,*taued1;

  pt0   = &mesh->tetra[k];
  vold0 = _MMG5_orvol(mesh->point,pt0->v);

  assert ( pt0->xt );

  pxt0 = &mesh->xtetra[pt0->xt];
  for (j0=0; j0<4; j0++) {
    if ( pxt0->ftag[j0] & MG_BDY ) continue;

    /** Neighbouring element with which we will try to swap */
    adja = &mesh->adja[4*(k-1)+1];
    k1   = adja[j0]/4;
    j1   = adja[j0]%4;

    assert(k1);

    /* Search in which configurations are the tetrahedra (default is case 0-0)
     *
     *           3                    2------------- 0
     *         ,/|`\                  |`\          /|
     *       ,/  |  `\                |  `\       /.|
     *     ,/    '.   `\              '.   `\    / |
     *   ,/       |     `\             |     `\ / .|
     * ,/         |       `\           |       /\.|
     * 0-----------'.--------2          '     /  3
     * `\.         |      ,/            |    / ,/
     *    `\.      |    ,/              |   /,/
     *       `\.   '. ,/                '. ,/
     *          `\. |/                   |/
     *             `1                    1
     */

    /* k may be in configuration 0, 3, 6 or 9. Default is case 0 */
    conf0 = 3*j0;

    switch(conf0) {
    case 0:
      tau0[0] = 0; tau0[1] = 1; tau0[2] = 2; tau0[3] = 3;
      taued0 = &MMG5_permedge[0][0];
      break;
    case 3:
      tau0[0] = 1; tau0[1] = 0; tau0[2] = 3; tau0[3] = 2;
      taued0 = &MMG5_permedge[3][0];
      break;
    case 6:
      tau0[0] = 2; tau0[1] = 0; tau0[2] = 1; tau0[3] = 3;
      taued0 = &MMG5_permedge[6][0];
      break;
    case 9:
      tau0[0] = 3; tau0[1] = 0; tau0[2] = 2; tau0[3] = 1;
      taued0 = &MMG5_permedge[9][0];
      break;
    }

    /* k1 may be in configuration j1, j1+1, j1+2 */
    pt1 = &mesh->tetra[k1];

    if ( pt1->tag & MG_REQ ) continue;

    assert(pt0->ref == pt1->ref);
    for ( i=0; i<3; ++i )
      if ( pt0->v[_MMG5_idir[j0][0]] == pt1->v[_MMG5_idir[j1][i]] ) break;

    assert(i<3);
    conf1 = 3*j1+i;

    switch(conf1) {
    case 0:
      tau1[0] = 0; tau1[1] = 1; tau1[2] = 2; tau1[3] = 3;
      taued1 = &MMG5_permedge[0][0];
      break;
    case 1:
      tau1[0] = 0; tau1[1] = 2; tau1[2] = 3; tau1[3] = 1;
      taued1 = &MMG5_permedge[1][0];
      break;
    case 2:
      tau1[0] = 0; tau1[1] = 3; tau1[2] = 1; tau1[3] = 2;
      taued1 = &MMG5_permedge[2][0];
      break;
    case 3:
      tau1[0] = 1; tau1[1] = 0; tau1[2] = 3; tau1[3] = 2;
      taued1 = &MMG5_permedge[3][0];
      break;
    case 4:
      tau1[0] = 1; tau1[1] = 3; tau1[2] = 2; tau1[3] = 0;
      taued1 = &MMG5_permedge[5][0];
      break;
    case 5:
      tau1[0] = 1; tau1[1] = 2; tau1[2] = 0; tau1[3] = 3;
      taued1 = &MMG5_permedge[4][0];
      break;
    case 6:
      tau1[0] = 2; tau1[1] = 0; tau1[2] = 1; tau1[3] = 3;
      taued1 = &MMG5_permedge[6][0];
      break;
    case 7:
      tau1[0] = 2; tau1[1] = 1; tau1[2] = 3; tau1[3] = 0;
      taued1 = &MMG5_permedge[7][0];
      break;
    case 8:
      tau1[0] = 2; tau1[1] = 3; tau1[2] = 0; tau1[3] = 1;
      taued1 = &MMG5_permedge[8][0];
      break;
    case 9:
      tau1[0] = 3; tau1[1] = 0; tau1[2] = 2; tau1[3] = 1;
      taued1 = &MMG5_permedge[9][0];
      break;
    case 10:
      tau1[0] = 3; tau1[1] = 2; tau1[2] = 1; tau1[3] = 0;
      taued1 = &MMG5_permedge[11][0];
      break;
    case 11:
      tau1[0] = 3; tau1[1] = 1; tau1[2] = 0; tau1[3] = 2;
      taued1 = &MMG5_permedge[10][0];
      break;
    }

    /* Test volume of the 3 created tets */
    vold = MG_MIN(vold0,_MMG5_orvol(mesh->point,pt1->v));

    ptnew = &mesh->tetra[0];
    memcpy(ptnew,pt0,sizeof(MMG5_Tetra));
    np    = pt1->v[tau1[0]];

#warning which threshold to refuse the swap?
    ptnew->v[tau0[1]] = np;
    vnew = _MMG5_orvol(mesh->point,ptnew->v);
    if ( vnew < _MMG5_EPSD2 ) continue;
    else if ( vold > _MMG5_NULKAL && vnew < _MMG5_NULKAL ) continue;

    ptnew->v[tau0[1]] = pt0->v[tau0[1]];
    ptnew->v[tau0[2]] = np;
    vnew = _MMG5_orvol(mesh->point,ptnew->v);
    if ( vnew < _MMG5_EPSD2 ) continue;
    else if ( vold > _MMG5_NULKAL && vnew < _MMG5_NULKAL ) continue;

    ptnew->v[tau0[2]] = pt0->v[tau0[2]];
    ptnew->v[tau0[3]] = np;
    vnew = _MMG5_orvol(mesh->point,ptnew->v);
    if ( vnew < _MMG5_EPSD2 ) continue;
    else if ( vold > _MMG5_NULKAL && vnew < _MMG5_NULKAL ) continue;

    /** Swap */
    xt1 = pt1->xt;
    memcpy(pt1,pt0,sizeof(MMG5_Tetra));

    iel = _MMG3D_newElt(mesh);
    if ( !iel ) {
      _MMG5_TETRA_REALLOC(mesh,iel,mesh->gap,
                          fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                  " a new element.\n",__func__);
                          _MMG5_INCREASE_MEM_MESSAGE();
                          fprintf(stderr,"  Exit program.\n");
                          return -1,-1);
    }
    ptnew = &mesh->tetra[iel];
    memcpy(ptnew,pt0,sizeof(MMG5_Tetra));

    /* First tetra: k */
    pt0->v[tau0[1]] = np;

    /* Second tetra: k1 */
    pt1->v[tau0[2]] = np;

    /* Third tetra: iel */
    ptnew->v[tau0[3]] = np;

    /* xtetra and adjacency update */
    pxt0 = &mesh->xtetra[pt0->xt];
    memcpy(&xt[0],pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt[1],pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt[2],pxt0,sizeof(MMG5_xTetra));

    /* Store the old adja */
    adja = &mesh->adja[4*(k-1) +1];
    adj0_2 = adja[tau0[2]];
    adj0_3 = adja[tau0[3]];

    adja = &mesh->adja[4*(k1-1) +1];
    adj1_1 = adja[tau1[1]];
    adj1_2 = adja[tau1[2]];
    adj1_3 = adja[tau1[3]];

    /* New adja for the new tets */
    adja = &mesh->adja[4*(k-1) +1];
    adja[tau0[0]] = adj1_1;
    adja[tau0[2]] = 4*k1  + tau0[1] ;
    adja[tau0[3]] = 4*iel + tau0[1] ;
    if ( adj1_1 )
      mesh->adja[4*(adj1_1/4-1) + 1 + adj1_1%4] = 4*k + tau0[0];

    adja = &mesh->adja[4*(k1-1) +1];
    adja[tau0[0]] = adj1_3;
    adja[tau0[1]] = 4*k   + tau0[2] ;
    adja[tau0[2]] = adj0_2;
    adja[tau0[3]] = 4*iel + tau0[2] ;
    if ( adj1_3 )
      mesh->adja[4*(adj1_3/4-1) + 1 + adj1_3%4] = 4*k1 + tau0[0];
    if ( adj0_2 )
      mesh->adja[4*(adj0_2/4-1) + 1 + adj0_2%4] = 4*k1 + tau0[2];

    adja = &mesh->adja[4*(iel-1) +1];
    adja[tau0[0]] = adj1_2;
    adja[tau0[1]] = 4*k   + tau0[3] ;
    adja[tau0[2]] = 4*k1  + tau0[3] ;
    adja[tau0[3]] = adj0_3;
    if ( adj1_2 )
      mesh->adja[4*(adj1_2/4-1) + 1 + adj1_2%4] = 4*iel + tau0[0];
    if ( adj0_3 )
      mesh->adja[4*(adj0_3/4-1) + 1 + adj0_3%4] = 4*iel + tau0[3];

    if ( !pt1->xt ) {
      /* Assignation of the xt fields to the appropriate tets */
      /* xt[0] */
      xt[0].tag[taued0[0]] = 0;
      xt[0].tag[taued0[3]] = 0;
      xt[0].tag[taued0[4]] = 0;

      xt[0].edg[taued0[0]] = 0;
      xt[0].edg[taued0[3]] = 0;
      xt[0].edg[taued0[4]] = 0;

      xt[0].ref[ tau0[0]] = 0;
      xt[0].ref[ tau0[2]] = 0;
      xt[0].ref[ tau0[3]] = 0;
      xt[0].ftag[tau0[0]] = 0;
      xt[0].ftag[tau0[2]] = 0;
      xt[0].ftag[tau0[3]] = 0;

      MG_SET(xt[0].ori, tau0[0]);
      MG_SET(xt[0].ori, tau0[2]);
      MG_SET(xt[0].ori, tau0[3]);

      /* xt[1] */
      xt[1].tag[taued0[1]] = 0;
      xt[1].tag[taued0[3]] = 0;
      xt[1].tag[taued0[5]] = 0;

      xt[1].edg[taued0[1]] = 0;
      xt[1].edg[taued0[3]] = 0;
      xt[1].edg[taued0[5]] = 0;

      xt[1].ref[ tau0[0]] = 0;
      xt[1].ref[ tau0[1]] = 0;
      xt[1].ref[ tau0[3]] = 0;
      xt[1].ftag[tau0[0]] = 0;
      xt[1].ftag[tau0[1]] = 0;
      xt[1].ftag[tau0[3]] = 0;

      MG_SET(xt[1].ori, tau0[0]);
      MG_SET(xt[1].ori, tau0[1]);
      MG_SET(xt[1].ori, tau0[3]);

      /* xt[2] */
      xt[1].tag[taued0[2]] = 0;
      xt[1].tag[taued0[4]] = 0;
      xt[1].tag[taued0[5]] = 0;

      xt[1].edg[taued0[2]] = 0;
      xt[1].edg[taued0[4]] = 0;
      xt[1].edg[taued0[5]] = 0;

      xt[1].ref[ tau0[0]] = 0;
      xt[1].ref[ tau0[1]] = 0;
      xt[1].ref[ tau0[2]] = 0;
      xt[1].ftag[tau0[0]] = 0;
      xt[1].ftag[tau0[1]] = 0;
      xt[1].ftag[tau0[2]] = 0;

      MG_SET(xt[1].ori, tau0[0]);
      MG_SET(xt[1].ori, tau0[1]);
      MG_SET(xt[1].ori, tau0[2]);

     }
    else {
      pxt1 = &mesh->xtetra[xt1];

      /* Assignation of the xt fields to the appropriate tets */
      /* xt[0] */
      xt[0].tag[taued0[0]] = 0;
      xt[0].tag[taued0[3]] = pxt1->tag[taued1[1]];
      xt[0].tag[taued0[4]] = pxt1->tag[taued1[2]];

      xt[0].edg[taued0[0]] = 0;
      xt[0].edg[taued0[3]] =  pxt1->edg[taued1[1]];
      xt[0].edg[taued0[4]] =  pxt1->edg[taued1[2]];

      xt[0].ref[ tau0[0]] = pxt1->ref[tau1[1]];
      xt[0].ref[ tau0[2]] = 0;
      xt[0].ref[ tau0[3]] = 0;
      xt[0].ftag[tau0[0]] = pxt1->ftag[tau1[1]];
      xt[0].ftag[tau0[2]] = 0;
      xt[0].ftag[tau0[3]] = 0;

      if ( MG_GET(pxt1->ori,tau1[1]) ) MG_SET(xt[0].ori, tau0[0]);
      MG_SET(xt[0].ori, tau0[2]);
      MG_SET(xt[0].ori, tau0[3]);

      /* xt[1] */
      xt[1].tag[taued0[1]] = 0;
      xt[1].tag[taued0[3]] = pxt1->tag[taued1[0]];
      xt[1].tag[taued0[5]] = pxt1->tag[taued1[2]];

      xt[1].edg[taued0[1]] = 0;
      xt[1].edg[taued0[3]] =  pxt1->edg[taued1[0]];
      xt[1].edg[taued0[5]] =  pxt1->edg[taued1[2]];

      xt[1].ref[ tau0[0]] = pxt1->ref[tau1[3]];
      xt[1].ref[ tau0[1]] = 0;
      xt[1].ref[ tau0[3]] = 0;
      xt[1].ftag[tau0[0]] = pxt1->ftag[tau1[3]];
      xt[1].ftag[tau0[1]] = 0;
      xt[1].ftag[tau0[3]] = 0;

      if ( MG_GET(pxt1->ori,tau1[3]) ) MG_SET(xt[1].ori, tau0[0]);
      MG_SET(xt[1].ori, tau0[1]);
      MG_SET(xt[1].ori, tau0[3]);

      /* xt[2] */
      xt[2].tag[taued0[2]] = 0;
      xt[2].tag[taued0[4]] = pxt1->tag[taued1[0]];
      xt[2].tag[taued0[5]] = pxt1->tag[taued1[1]];

      xt[2].edg[taued0[2]] = 0;
      xt[2].edg[taued0[4]] = pxt1->edg[taued1[0]];
      xt[2].edg[taued0[5]] = pxt1->edg[taued1[1]];

      xt[2].ref[ tau0[0]] = pxt1->ref[tau1[2]];
      xt[2].ref[ tau0[1]] = 0;
      xt[2].ref[ tau0[2]] = 0;
      xt[2].ftag[tau0[0]] = pxt1->ftag[tau1[2]];
      xt[2].ftag[tau0[1]] = 0;
      xt[2].ftag[tau0[2]] = 0;

      if ( MG_GET(pxt1->ori,tau1[2]) ) MG_SET(xt[2].ori, tau0[0]);
      MG_SET(xt[2].ori, tau0[1]);
      MG_SET(xt[2].ori, tau0[2]);
    }

    /* Assignation of the xt fields to the appropriate tets */
    isxt[0] = isxt[1] = isxt[2] = 0;
    for (i=0; i<4; i++) {
      if ( xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
      if ( xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
      if ( xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
    }

    assert ( (isxt[0] && isxt[1]) || (isxt[1] && isxt[2]) || (isxt[2] && isxt[0]) );

    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));

      if ( xt1 ) {
        if ( isxt[1] ) {
          pt1->xt = xt1;
          memcpy(pxt1,&xt[1],sizeof(MMG5_xTetra));
          if ( isxt[2] ) {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              /* realloc of xtetras table */
              _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");
                                 return -1,-1);
            }
            ptnew->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[2],sizeof(MMG5_xTetra));
          }
          else ptnew->xt = 0;
        }
        else {
          pt1->xt   = 0;
          ptnew->xt = xt1;
          memcpy(pxt1,&xt[2],sizeof(MMG5_xTetra));
        }
      }
      else {
        if ( isxt[1] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return -1,-1);
          }
          pt1->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[1],sizeof(MMG5_xTetra));
        }
        else pt1->xt = 0;

        if ( isxt[2] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return -1,-1);
          }
          ptnew->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[2],sizeof(MMG5_xTetra));
        }
        else ptnew->xt = 0;
      }
    }
    else {
      pt0->xt = 0;
      memcpy(pxt0 ,&xt[2],sizeof(MMG5_xTetra));

      if ( xt1 ) {
        pt1->xt = xt1;
        memcpy(pxt1,&xt[1],sizeof(MMG5_xTetra));
      }
      else {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             fprintf(stderr,"  Exit program.\n");
                             return -1,-1);
        }
        pt1->xt = mesh->xt;
        pxt0 = &mesh->xtetra[mesh->xt];
        memcpy(pxt0,&xt[1],sizeof(MMG5_xTetra));
      }
    }

    /** Quality Update */
    if ( (!metRidTyp) && met->m && met->size>1 ) {
      pt0->qual   = _MMG5_caltet33_ani(mesh,met,pt0);
      pt1->qual   = _MMG5_caltet33_ani(mesh,met,pt1);
      ptnew->qual = _MMG5_caltet33_ani(mesh,met,ptnew);
    }
    else
    {
      pt0->qual   = _MMG5_orcal(mesh,met,k);
      pt1->qual   = _MMG5_orcal(mesh,met,k1);
      ptnew->qual = _MMG5_orcal(mesh,met,iel);
    }
    pt0->mark   = mesh->mark;
    pt1->mark   = mesh->mark;
    ptnew->mark = mesh->mark;

    return 1;
  }
  return 0;
}
