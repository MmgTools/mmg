/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
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
 * \file mmg3d/quality.c
 * \brief Functions to compute elements quality and edge lengths.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmg3d.h"

extern char ddb;

/** compute tetra oriented quality of iel (return 0.0 when element is inverted) */
inline double _MMG5_orcal(MMG5_pMesh mesh,MMG5_pSol met,int iel) {
  MMG5_pTetra     pt;

  pt = &mesh->tetra[iel];

  if ( met->m )
    return(_MMG5_caltet(mesh,met,pt->v[0],pt->v[1],pt->v[2],pt->v[3]));
  else // with -A option we are in aniso but without metric.
    return(_MMG5_caltet_iso(mesh,met,pt->v[0],pt->v[1],pt->v[2],pt->v[3]));
}


/** compute tetra quality iso */
inline double _MMG5_caltet_iso(MMG5_pMesh mesh,MMG5_pSol met,int ia,int ib,int ic,int id) {
  double     abx,aby,abz,acx,acy,acz,adx,ady,adz,bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     vol,v1,v2,v3,rap;
  double    *a,*b,*c,*d;

  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;
  d = mesh->point[id].c;

  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  rap = abx*abx + aby*aby + abz*abz;

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  rap += acx*acx + acy*acy + acz*acz;

  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];
  rap += adx*adx + ady*ady + adz*adz;

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;
  if ( vol < _MMG5_EPSD2 )  return(0.0);

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];
  rap += bcx*bcx + bcy*bcy + bcz*bcz;

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];
  rap += bdx*bdx + bdy*bdy + bdz*bdz;

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];
  rap += cdx*cdx + cdy*cdy + cdz*cdz;
  if ( rap < _MMG5_EPSD2 )  return(0.0);

  /* quality = vol / len^3/2 */
  rap = rap * sqrt(rap);
  return(vol / rap);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the meric structure.
 * \param ia, ib, ic, id point index
 * \return The anisotropic quality of the tet (ia,ib,ic,id).
 *
 * Compute the quality of the tet K (\a ia,\a ib,\a ic,\a id) with respect to
 * the anisotropic metric \a met.
 *    Q = V_met(K) / (sum(len(edge_K)^2)^(3/2)
 *
 * \todo test with the square of this measure
 */
inline double _MMG5_caltet_ani(MMG5_pMesh mesh,MMG5_pSol met,int ia,int ib,int ic,int id) {
  double     cal,abx,aby,abz,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     h1,h2,h3,h4,h5,h6,det,vol,rap,v1,v2,v3,num;
  double    *a,*b,*c,*d;
  double     *ma,*mb,*mc,*md,mm[6];
  int        j,iadr;

  cal = _MMG5_NULKAL;

  /* average metric */
  memset(mm,0,6*sizeof(double));
  iadr = (ia)*met->size;
  ma   = &met->m[iadr];
  iadr = (ib)*met->size;
  mb   = &met->m[iadr];
  iadr = (ic)*met->size;
  mc   = &met->m[iadr];
  iadr = (id)*met->size;
  md   = &met->m[iadr];
  for (j=0; j<6; j++)
    mm[j] = 0.25 * (ma[j]+mb[j]+mc[j]+md[j]);
  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;
  d = mesh->point[id].c;

  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;            
  if ( vol <= 0. )  return(cal);
  det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \
      - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \
      + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]);   
  if ( det < _MMG5_EPSOK )   {
    //printf("--- INVALID METRIC : DET (%d) %e\n",iel,det);
    return(cal);
  }
  det = sqrt(det) * vol;
  /* edge lengths */
  h1 =      mm[0]*abx*abx + mm[3]*aby*aby + mm[5]*abz*abz \
     + 2.0*(mm[1]*abx*aby + mm[2]*abx*abz + mm[4]*aby*abz);
  h2 =      mm[0]*acx*acx + mm[3]*acy*acy + mm[5]*acz*acz \
     + 2.0*(mm[1]*acx*acy + mm[2]*acx*acz + mm[4]*acy*acz);
  h3 =      mm[0]*adx*adx + mm[3]*ady*ady + mm[5]*adz*adz \
     + 2.0*(mm[1]*adx*ady + mm[2]*adx*adz + mm[4]*ady*adz);

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  h4 =      mm[0]*bdx*bdx + mm[3]*bdy*bdy + mm[5]*bdz*bdz \
     + 2.0*(mm[1]*bdx*bdy + mm[2]*bdx*bdz + mm[4]*bdy*bdz);
  h5 =      mm[0]*cdx*cdx + mm[3]*cdy*cdy + mm[5]*cdz*cdz \
     + 2.0*(mm[1]*cdx*cdy + mm[2]*cdx*cdz + mm[4]*cdy*cdz);
  h6 =      mm[0]*bcx*bcx + mm[3]*bcy*bcy + mm[5]*bcz*bcz \
     + 2.0*(mm[1]*bcx*bcy + mm[2]*bcx*bcz + mm[4]*bcy*bcz);

  /* quality */
  rap = h1 + h2 + h3 + h4 + h5 + h6;
  num = sqrt(rap) * rap;  

  cal = det / num;  
  if(cal <= _MMG5_NULKAL) {
    printf(" TOO BAD QUALITY %e %e %e %e\n",cal,num,det,vol);  
    return(_MMG5_NULKAL);
  }
  //printf("cal %e %e %e\n",cal,num,det);
  assert(cal > _MMG5_NULKAL);
  return(cal); 
}

/* identify type of element :
   ityp= 0: 4 faces bonnes          (elt ok)
   1: 4 faces bonnes, vol nul (sliver) ou "quasi sliver" ie 4 faces ok
   2: 4 faces ok, vol nul+sommet proche face   (chapeau)
   3: 3 faces bonnes, 1 obtuse    (aileron)
   4: 2 faces bonnes, 2 faces aigu => 1 petite arete
   5: 1 face bonne, 3 petites aretes
   6: 2 faces grandes aretes, 2 faces petites _MMG5_iaretes
   7: 4 faces grandes aretes
   8: 2 faces obtus, 1 faces aigu et une face OK
   item: bad entity
*/


/* nb face obtuse :    nb faces aigu :
   ityp :  0: 0        0
   1: 0        0
   2: 0        0
   3: 1        0
   4: 0        2
   5: 0        3
   6: 2        2
   7: 0        4
*/
/* nb gde arete :    nb petite arete :
   ityp :  0: 0        0
   1: 0        0
   2: 0        0
   3: 1        0
   4: 0        1
   5: 0        3
   6: 1        1
   7: 0        2
*/
#define _MMG5_EPSVOL 0.001
#define RAPMAX    0.4//0.3//0.25
unsigned char inxt[7]    = { 1,2,0,1,2,0,1 };

/**
 * \warning Not used.
 */
int _MMG5_typelt(MMG5_pMesh mesh,int iel,int *item) {
  MMG5_pTetra    pt;
  MMG5_pPoint    pa,pb,pc,pd;
  double    abx,aby,abz,acx,acy,acz,adx,ady,adz,v1,v2,v3,vol;
  double    bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz,h[6],volchk,ssmall;
  double    s[4],dd,rapmin,rapmax,surmin,surmax;
  int       i,k,ia,ib,ic,id,isur,isurmax,isurmin,iarmax,iarmin;
  int       nobtus,naigu,aigu;
  short     i0,i1,i2;
  double lmoy;

  pt = &mesh->tetra[iel];
  if ( !pt->v[0] )  return(-1);

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];
  pa = &mesh->point[ia];
  pb = &mesh->point[ib];
  pc = &mesh->point[ic];
  pd = &mesh->point[id];

  /* volume */
  abx = pb->c[0] - pa->c[0];
  aby = pb->c[1] - pa->c[1];
  abz = pb->c[2] - pa->c[2];

  acx = pc->c[0] - pa->c[0];
  acy = pc->c[1] - pa->c[1];
  acz = pc->c[2] - pa->c[2];

  adx = pd->c[0] - pa->c[0];
  ady = pd->c[1] - pa->c[1];
  adz = pd->c[2] - pa->c[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;

  /* max edge */
  h[0] = abx*abx + aby*aby + abz*abz;
  h[1] = acx*acx + acy*acy + acz*acz;
  h[2] = adx*adx + ady*ady + adz*adz;

  bcx = pc->c[0] - pb->c[0];
  bcy = pc->c[1] - pb->c[1];
  bcz = pc->c[2] - pb->c[2];

  bdx = pd->c[0] - pb->c[0];
  bdy = pd->c[1] - pb->c[1];
  bdz = pd->c[2] - pb->c[2];

  cdx = pd->c[0] - pc->c[0];
  cdy = pd->c[1] - pc->c[1];
  cdz = pd->c[2] - pc->c[2];

  h[3] = bcx*bcx + bcy*bcy + bcz*bcz;
  h[4] = bdx*bdx + bdy*bdy + bdz*bdz;
  h[5] = cdx*cdx + cdy*cdy + cdz*cdz;

  /* face areas */
  dd = cdy*bdz - cdz*bdy;
  s[0] = dd * dd;
  dd = cdz*bdx - cdx*bdz;
  s[0] = s[0] + dd * dd;
  dd = cdx*bdy - cdy*bdx;
  s[0] = s[0] + dd * dd;
  s[0] = sqrt(s[0]);

  s[1] = sqrt(v1*v1 + v2*v2 + v3*v3);

  dd = bdy*adz - bdz*ady;
  s[2] = dd * dd;
  dd = bdz*adx - bdx*adz;
  s[2] = s[2] + dd * dd;
  dd = bdx*ady - bdy*adx;
  s[2] = s[2] + dd * dd;
  s[2] = sqrt(s[2]);

  dd = aby*acz - abz*acy;
  s[3] = dd * dd;
  dd = abz*acx - abx*acz;
  s[3] = s[3] + dd * dd;
  dd = abx*acy - aby*acx;
  s[3] = s[3] + dd * dd;
  s[3] = sqrt(s[3]);

  /* classification */
  rapmin = h[0];
  rapmax = h[0];
  iarmin = 0;
  iarmax = 0;
  for (i=1; i<6; i++) {
    if ( h[i] < rapmin ) {
      rapmin = h[i];
      iarmin = i;
    }
    else if ( h[i] > rapmax ) {
      rapmax = h[i];
      iarmax = i;
    }
  }
  rapmin = sqrt(rapmin);
  rapmax = sqrt(rapmax);
  volchk = _MMG5_EPSVOL * rapmin*rapmin*rapmin;
  //printf("$$$$$$$$$$$$$$$$$$$$$ iel %d %e %e\n",iel,volchk,vol);
  /* small volume: types 1,2,3,4 */
  if ( vol < volchk ) {
    puts("volume nul : type 1,2,3,4");

    ssmall = 0.4 * (s[0]+s[1]+s[2]+s[3]);
    isur   = 0;
    for (i=0; i<4; i++)
      isur += s[i] > ssmall;

    /* types 2,3 */
    item[0] = iarmax;
    item[1] = _MMG5_isar[iarmax][0];
    if ( isur == 1 ) {
      surmin   = s[0];
      isurmin = 0;
      surmax   = s[0];
      isurmax = 0;
      for (i=1; i<4; i++) {
        if ( s[i] < surmin ) {
          surmin  = s[i];
          isurmin = i;
        }
        else if ( s[i] > surmax ) {
          surmax  = s[i];
          isurmax = i;
        }
      }
      dd = surmin / surmax;
      if ( dd < RAPMAX ) {
        item[1] = _MMG5_isar[iarmax][0];
        return(3);
      }
      else {
        item[0] = isurmax;
        item[1] = isurmin;
        return(2);
      }
    }

    /* types 1 */
    isur = 0;
    if ( s[0]+s[1] > ssmall )  isur = 1;
    if ( s[0]+s[2] > ssmall )  isur++;
    if ( s[0]+s[3] > ssmall )  isur++;

    if ( isur > 2 ) {
      dd = rapmin / rapmax;
      item[0] = iarmin;
      item[1] = _MMG5_idir[iarmin][0];
      if ( dd < 0.01 )  return(4);
      if ( s[0]+s[1] > ssmall ) {
        item[0] = 0;
        return(1);
      }
      if ( s[0]+s[2] > ssmall ) {
        item[0] = 1;
        return(1);
      }
      if ( s[0]+s[3] > ssmall ) {
        item[0] = 2;
        return(1);
      }
    }

    //puts("default");
    item[0] = 0;
    return(1);
  }/*end chkvol*/

  dd = rapmin / rapmax;
  // printf("dd %e %e %e\n",dd,RAPMAX,0.7*RAPMAX);
  /* types 3,6,7 */
  if ( dd < RAPMAX ) { /*ie une arete 3 fois plus gde qu'une autre*/
    lmoy = 0;
    for (i=0; i<6; i++)  {h[i] = sqrt(h[i]); lmoy+=h[i];}
    lmoy *= 1./6.;

    nobtus = 0;
    naigu = 0;
    for (k=0; k<4; k++) {
      aigu = 0;
      for(i=0 ; i<3 ; i++) {
        i0 = _MMG5_idir[k][i];
        i1 = _MMG5_idir[k][inxt[i]];
        i2 = _MMG5_idir[k][inxt[i+1]];
        if((h[i0]>h[i1] && h[i0]>h[i2]) && (h[i0]>2.5*h[i1] || h[i0]>2.5*h[i2])) {//obtu ? > /*130*/ 150
          //be carefull isocele triangle --> opposite edge only
          if(!( fabs(1-h[i0]/h[i1]) < 0.1 || fabs(1-h[i0]/h[i2])< 0.1 ) ) {
            if(h[i0] > 0.9659258/*0.9063078*/*(h[i1]+h[i2])) break;
          }
        }
        if((h[i0]<h[i1] && h[i0]<h[i2]) && (2.5*h[i0]<h[i1] || 2.5*h[i0]<h[i2])) {//aigu ? <20
          if(h[i0] < 0.17*(h[i1]+h[i2])) aigu++;
        }
      }
      if(i<3) nobtus++;
      else if(aigu) naigu++;
    }
    //printf("%d on trouve %d aigu et %d obtu\n",iel,naigu,nobtus);
    switch(nobtus) {
    case(3):case(4):
      return(8);
    case(2):
      if(naigu==2) {
        return(6);
      } else {
        return(8);
      }
    case(1):
      return(3);
    case(0):
      if(naigu==4) {
        return(7);
      } else if(naigu==3) {
        return(5);
      } else if(naigu==2) {
        return(4);
      } else {
        //printf("%d on trouve %d aigu et %d obtu\n",iel,naigu,nobtus);
        return(9); //toutes les faces sont ok
      }
    }
  }


  item[0] = 0;
  // printf("edge %e %e -- vol %e\n",rapmin,rapmax,vol);
  //puts("default");
  return(9);
}

/**
 * \warning Not used.
 */
int _MMG5_badelt(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra   pt;
  double   kal;
  int      k,it,maxit,nd/*,item[2],typ*/;
  /*int      ntyp[10];*/
  int      list[_MMG5_LMAX+2],i,ilist,nconf,ns;

  it = 0;
  maxit = 3;
  do {
    nd = 0;
    ns = 0;
    //for (k=1; k<10; k++) ntyp[k]=0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      kal = /*_MMG5_ALPHAD **/ pt->qual;
      if ( kal > 0.0096225 /*_MMG5_BADKAL/_MMG5_ALPHAD*/ )  continue;
      //typ =  _MMG5_typelt(mesh,k,item);
      //ntyp[typ]++;
      nd++;
      /*treat bad elt*/
      /*1) try to swp one edge*/
      for(i=0 ; i<6 ; i++) {
        nconf = _MMG5_chkswpgen(mesh,met,k,i,&ilist,list,1.01);
        if ( nconf ) {
          ns++;
          if(!_MMG5_swpgen(mesh,met,nconf,ilist,list,NULL)) return(-1);
          break;
        }
      }
    }
    /*printf("on trouve %d bad elt\n",nd);
      for (k=0; k<=9; k++)
      if ( ntyp[k] )
      printf("  optim [%d]      = %5d  %6.2f %%\n",k,ntyp[k],100.0*ntyp[k]/nd);
    */if ( ns > 0 )
      fprintf(stdout,"     %8d edge swapped\n",ns);
  }
  while ( ++it < maxit && nd > 0 );
  return(nd);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Compute sizes of edges of the mesh, and displays histo.
 *
 */
int _MMG5_prilen(MMG5_pMesh mesh, MMG5_pSol met) {
  MMG5_pTetra     pt;
  MMG5_pxTetra    pxt;
  _MMG5_Hash      hash;
  double          len,avlen,lmin,lmax;
  int             k,np,nq,amin,bmin,amax,bmax,ned,hl[9];
  char            ia,i0,i1,ier,i;
  static double   bd[9]= {0.0, 0.3, 0.6, 0.7071, 0.9, 1.3, 1.4142, 2.0, 5.0};
  //{0.0, 0.2, 0.5, 0.7071, 0.9, 1.111, 1.4142, 2.0, 5.0};

  memset(hl,0,9*sizeof(int));
  ned = 0;
  avlen = 0.0;
  lmax = 0.0;
  lmin = 1.e30;
  amin = amax = bmin = bmax = 0;

  /* Hash all edges in the mesh */
  if ( !_MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return(0);

  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = _MMG5_iare[ia][0];
      i1 = _MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      if(!_MMG5_hashEdge(mesh,&hash,np,nq,0)){
        fprintf(stdout,"%s:%d: Error: function _MMG5_hashEdge return 0\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
      }
    }
  }


  /* Pop edges from hash table, and analyze their length */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    for(ia=0; ia<6; ia++) {
      i0 = _MMG5_iare[ia][0];
      i1 = _MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      /* Remove edge from hash ; ier = 1 if edge has been found */
      ier = _MMG5_hashPop(&hash,np,nq);
      if( ier ) {
        ned ++;
        if ( pt->xt )
          len = _MMG5_lenedg(mesh,met,np,nq,(pxt->tag[ia] & MG_GEO));
        else
          len = _MMG5_lenedg(mesh,met,np,nq,0);

        avlen += len;

        if( len < lmin ) {
          lmin = len;
          amin = np;
          bmin = nq;
        }

        if ( len > lmax ) {
          lmax = len;
          amax = np;
          bmax = nq;
        }

        /* Locate size of edge among given table */
        for(i=0; i<8; i++) {
          if ( bd[i] <= len && len < bd[i+1] ) {
            hl[i]++;
            break;
          }
        }
        if( i == 8 ) hl[8]++;
      }
    }
  }

  /* Display histogram */
  _MMG5_displayHisto(mesh, ned, &avlen, amin, bmin, lmin,
                     amax, bmax, lmax, &bd[0], &hl[0]);

  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 *
 * Print histogram of mesh qualities.
 *
 */
void _MMG5_outqua(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra    pt;
  double   rap,rapmin,rapmax,rapavg,med,good;
  int      i,k,iel,ok,ir,imax,nex,his[5];

  /*compute tet quality*/
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
     if( !MG_EOK(pt) )   continue;
     pt->qual = _MMG5_orcal(mesh,met,k);
  }
  if ( abs(mesh->info.imprim) <= 0 ) return;

  rapmin  = 2.0;
  rapmax  = 0.0;
  rapavg  = med = good = 0.0;
  iel     = 0;

  for (k=0; k<5; k++)  his[k] = 0;

  nex = ok = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) ) {
      nex++;
      continue;
    }
    ok++;
    if ( _MMG5_orvol(mesh->point,pt->v) < 0.0 ) {
      fprintf(stdout,"dans quality vol negatif\n");
    }
    rap = _MMG5_ALPHAD * pt->qual;
    if ( rap < rapmin ) {
      rapmin = rap;
      iel    = ok;
    }
    if ( rap > 0.5 )  med++;
    if ( rap > 0.12 ) good++;
    if ( rap < _MMG5_BADKAL )  mesh->info.badkal = 1;
    rapavg += rap;
    rapmax  = MG_MAX(rapmax,rap);
    ir = MG_MIN(4,(int)(5.0*rap));
    his[ir] += 1;
  }

#ifndef DEBUG
  fprintf(stdout,"\n  -- MESH QUALITY   %d\n",mesh->ne - nex);
  fprintf(stdout,"     BEST   %8.6f  AVRG.   %8.6f  WRST.   %8.6f (%d)\n",
          rapmax,rapavg / (mesh->ne-nex),rapmin,iel);
#else
  fprintf(stdout,"     BEST   %e  AVRG.   %e  WRST.   %e (%d)\n => %d %d %d %d\n",
          rapmax,rapavg / (mesh->ne-nex),rapmin,iel,
          _MMG5_indPt(mesh,mesh->tetra[iel].v[0]),_MMG5_indPt(mesh,mesh->tetra[iel].v[1]),
          _MMG5_indPt(mesh,mesh->tetra[iel].v[2]),_MMG5_indPt(mesh,mesh->tetra[iel].v[3]));
#endif
  if ( abs(mesh->info.imprim) < 3 ){
    if (rapmin == 0){
      fprintf(stdout,"  ## WARNING: TOO BAD QUALITY FOR THE WORST ELEMENT\n");
      _MMG5_unscaleMesh(mesh,met);
      MMG5_saveMesh(mesh);
      MMG5_saveMet(mesh,met);
      exit(EXIT_FAILURE);
    }
    return;
  }

  /* print histo */
  fprintf(stdout,"     HISTOGRAMM:");
  fprintf(stdout,"  %6.2f %% > 0.12\n",100.0*(good/(float)(mesh->ne-nex)));
  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"                  %6.2f %% >  0.5\n",100.0*( med/(float)(mesh->ne-nex)));
    imax = MG_MIN(4,(int)(5.*rapmax));
    for (i=imax; i>=(int)(5*rapmin); i--) {
      fprintf(stdout,"     %5.1f < Q < %5.1f   %7d   %6.2f %%\n",
              i/5.,i/5.+0.2,his[i],100.*(his[i]/(float)(mesh->ne-nex)));
    }
  }
  if (rapmin == 0){
    fprintf(stdout,"  ## WARNING: TOO BAD QUALITY FOR THE WORST ELEMENT\n");
    _MMG5_unscaleMesh(mesh,met);
    MMG5_saveMesh(mesh);
    MMG5_saveMet(mesh,met);
    exit(EXIT_FAILURE);
  }
}

/**
 *
 * Approximation of the final number of vertex.
 *
 * \warning Not used.
 */
int _MMG5_countelt(MMG5_pMesh mesh,MMG5_pSol sol, double *weightelt, long *npcible) {
  MMG5_pTetra pt;
  MMG5_pxTetra pxt;
  double      len;
  int         k,ia,ipa,ipb,lon,l;
  //int   npbdry;
  int    *pdel,lenint,loc,nedel,longen;
  int      isbdry;
  double   dned,dnface,dnint/*,dnins*/,w,lenavg,lent[6];
  double   dnpdel,dnadd,leninv,dnaddloc,dnpdelloc;
  int      list[_MMG5_LMAX],ddebug=0,ib;
  long     nptot;
  //FILE *inm;

  pdel = (int*) calloc(mesh->np+1,sizeof(int));
  nptot = (long) mesh->np;

  // svg des poids
  // npbdry = 0;
  // inm = fopen("poid.sol","w");
  // fprintf(inm,"MeshVersionFormatted 2\n Dimension 3 \n SolAtTetrahedra \n %d\n 1 1 \n",mesh->ne);

  // substraction of the half of the number of bdry vertex to avoid the surestimation due of the interface
  // for (k=1; k<=mesh->np; k++) {
  //   if(mesh->point[k].tag & MG_BDY) npbdry++;
  // }
  // nptot -= 0.5*npbdry;

  dnadd = dnpdel = 0;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /*longueur moyenne*/
    lenavg = 0;
    for(ib=0 ; ib<6 ; ib++) {
      ipa = _MMG5_iare[ib][0];
      ipb = _MMG5_iare[ib][1];
      if ( pt->xt )
        lent[ib] = _MMG5_lenedg(mesh,sol,pt->v[ipa],pt->v[ipb],
                                (pxt->tag[ib] & MG_GEO ));
      else
        lent[ib] = _MMG5_lenedg(mesh,sol,pt->v[ipa],pt->v[ipb],0);
      lenavg+=lent[ib];
    }
    lenavg /= 6.;

    w = 0;
    if(weightelt)
      weightelt[k] = 0;
    nedel = 0;

    for (ia=0; ia<6; ia++) {
      //lon = MMG5_coquil(mesh,k,ia,&list);
      longen = _MMG5_coquil(mesh,k,ia,list);
      lon = longen/2;
      isbdry = 0;//longen%2;
      if(!lon) continue;
      /* if ( isbdry )  { */
      /*    assert(longen%2); */
      /*    //printf("_MMG5_coquil %d\n",longen/2); */
      /*    continue; */
      /* } */
      //assert(!(longen%2));
      for (l=1; l<lon; l++)
        if ( list[l] < 6*k )  break;

      if ( l < lon )  {
        loc = 1;
        //continue;
      } else {
        loc = 0;
      }

      dnaddloc = 0;
      dnpdelloc = 0;

      len = lent[ia];

      if(ddebug) printf("len %e\n",len);
      if(len > 3) {
        loc = 0;
        len = lenavg;
        lenint = ((int) len);
        if(fabs(lenint -len) > 0.5) lenint++;
        //POURQUOI SURESTIMER ???lenint++;
        //nb de point a inserer sur l'arete : longueur - 1
        dned = lenint - 1;
        //nb de point a l'interieur de la face si toutes les aretes sont coupees le meme nb de fois
        dnface = (lenint+2)*(lenint+1) / 2. - 3 - 3*dned;
        //nb de point a l'interieur du tetra si toutes les aretes sont coupees le meme nb de fois
        dnint = (lenint+3)*(lenint+2)*(lenint+1) / 6. - 4 - 4*dnface - 6*dned;
        //nb de point a inserer pour cette arete de ce tetra : on divise par lon
        // dnins = dned*(1./lon) + (dnface/3. + dnint/6.);//(dnface/12. + dnint/6.);
        if(!isbdry) {
          //nb points sur l'arete +
          //lon*(2/3 nb point sur la face (ie 1/3 de face et 2 faces adj a l'arete) + 1/6 nb de point interne)
          dnaddloc = dned + lon*(2*dnface/3. + dnint/6.);
        } else {
          dnaddloc = 0.5*(dned + lon*(2*dnface/3. + dnint/6.));
        }
        dnaddloc *= 1./lon;
        if(!loc) {
          if(_MMG5_ALPHAD * pt->qual >= 0.5) /*on ne compte les points internes que pour les (tres) bons tetras*/
            dnaddloc = dnaddloc;
          else if(_MMG5_ALPHAD * pt->qual >= 1./5.)
            dnaddloc = dned / lon + 2*dnface/3.;
          else
            dnaddloc = dned / lon ;
          //rajout de 30% parce que 1) on vise des longueurs de 0.7 et
          //2) on ne tient pas compte du fait qu'on divise tjs par 2 dans la generation
          if( (_MMG5_ALPHAD*pt->qual <= 1./50.) )
            dnaddloc = 0;
          else  if((_MMG5_ALPHAD*pt->qual <= 1./10.) )
            dnaddloc =  0.2*dnaddloc; //CEDRIC : essayer 0.3 0.4
          else if((len > 10) && (_MMG5_ALPHAD*pt->qual >= 1./1.5) ) //on sous-estime uniquement pour les tres bons
            dnaddloc = dnaddloc*0.3 + dnaddloc; //CEDRIC : essayer 0.3 ?
          else if(len < 6 && len>3) //CEDRIC : essayer len < 3,4, 6,7 mais aussi en commentant le test sur len puis pour la qual > 3, 5,8
            dnaddloc = 0.7*dnaddloc; //CEDRIC : essayer 0.9 0.7 0.6


          dnadd += dnaddloc;
        }
      } else if(len > 2.8) {
        if(!isbdry) {
          dnaddloc = 2.;
        } else {
          dnaddloc = 1;
        }
        if(!loc){
          if(!isbdry) {
            dnadd += 2.;
          } else {
            dnadd++;
          }
        }
        // dnins = 2;
      } else if(len > 1.41) {
        if(!isbdry)
          dnaddloc = 1;
        if(!loc) {
          if(!isbdry) dnadd += 1.;
        }
        // dnins = 1;
      } else if(len < 0.6) {
        nedel = 1;

        leninv = 1./len;
        if(pt->v[ipa]<pt->v[ipb]) {
          if(!pdel[pt->v[ipa]]) {
            if(!isbdry) {
              dnpdelloc = (leninv - 1.)/leninv;
            } else {
              dnpdelloc = 0.5*(leninv - 1.)/leninv;
            }
            if(!loc) {
              dnpdel+=dnpdelloc;
              pdel[pt->v[ipa]]=1;
            }
          } else if(!pdel[pt->v[ipb]]) {
            if(!isbdry) {
              dnpdelloc = (leninv - 1.)/leninv;
            } else {
              dnpdelloc = 0.5*(leninv - 1.)/leninv;
            }
            if(!loc) {
              dnpdel +=dnpdelloc;
              pdel[pt->v[ipb]]=1;
            }
          }
        } else {
          if(!pdel[pt->v[ipb]]) {
            if(!isbdry) {
              dnpdelloc = (leninv - 1.)/leninv;
            } else {
              dnpdelloc = 0.5*(leninv - 1.)/leninv;
            }
            if(!loc) {
              dnpdel+=dnpdelloc;
              pdel[pt->v[ipb]]=1;
            }
          } else if(!pdel[pt->v[ipa]]) {
            if(!isbdry) {
              dnpdelloc = (leninv - 1.)/leninv;
            } else {
              dnpdelloc = 0.5*(leninv - 1.)/leninv;
            }
            if(!loc) {
              dnpdel+=dnpdelloc;
              pdel[pt->v[ipa]]=1;
            }
          }
        }
        //ndel++;
      }

      //pour cette arete de ce tetra :
      //PHASE 1 = dnaddloc + nedel (on compte un si arete trop petite)
      //PHASE 2 = dnaddloc
      if(ddebug) printf("on ajoute %e\n",dnaddloc);
      w += (2*dnaddloc);//1./lon*(2*dnaddloc + dnpdelloc);

    }/*for ia*/
    if(ddebug) printf("on soustrait %d\n",nedel);

    w += 0.5*nedel;

    //si l'elt ne doit pas etre ni splitte ni collapse est ce qu'on le compte ???
    //if(w==0) w+=1;

    //fprintf(inm,"%e\n",w);
    if(weightelt)
      weightelt[k] = 10*w;
  } /*For k*/


  nptot += (long) dnadd - (long) dnpdel;
  *npcible = nptot;
  fprintf(stdout,"ESTIMATION OF THE FINAL NUMBER OF NODES : %ld  _MMG5_ADD %f  _MMG5_DEL %f\n",nptot,dnadd,dnpdel);

  free(pdel);

  //fclose(inm);
  return(1);
}
