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
 * \file mmg3d/opttyp.c
 * \brief Functions for the optimization of very bad elements.
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Algiane Froehly
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmg3d.h"

/* identify type of element :
   ityp= 0: 4 faces bonnes          (elt ok)
   1: 4 faces bonnes, vol nul (sliver)
   2: 4 faces ok, vol nul+sommet proche face   (chapeau)
   3: 3 faces bonnes, 1 obtuse    (aileron)
   4: 2 faces bonnes, 2 faces aigu => 1 petite arete
   5: 1 face bonne, 3 petites aretes
   6: 2 faces grandes aretes, 2 faces petites iaretes
   7: 4 faces grandes aretes
   item: bad entity
*/


/* nb face obtuse :    nb faces aigu :
   ityp :  0: 0             0
   1: 0             0
   2: 0             0
   3: 1             0
   4: 0             2
   5: 0             3
   6: 2             2
   7: 0             4
*/
/* nb gde arete :    nb petite arete :
   ityp :  0: 0             0
   1: 0             0
   2: 0             0
   3: 1             0
   4: 0             1
   5: 0             3
   6: 1             1
   7: 0             2
*/

int MMG3D_typelt(MMG5_pMesh mesh,int iel,int *item) {
  MMG5_pTetra    pt;
  MMG5_pPoint    pa,pb,pc,pd;
  double    abx,aby,abz,acx,acy,acz,adx,ady,adz,v1,v2,v3,vol;
  double    bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz,h[6],volchk,ssmall;
  double    s[4],dd,rapmin,rapmax,surmin,surmax;
  int       i,k,ia,ib,ic,id,ityp,isur,isurmax,isurmin,iarmax,iarmin;
  int       nobtus,naigu;
  short     i0,i1,i2;
  double    EPSVOL = 0.001;
  double    RAPMAX = 0.25;

  ityp = 0;
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
  volchk = EPSVOL * rapmin*rapmin*rapmin;

  /* small volume: types 1,2,3,4 */
  if ( vol < volchk ) {
    //puts("volume nul : type 1,2,3,4");
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
  /* types 3,6,7 */
  if ( dd < RAPMAX ) { /*ie une arete 4 fois plus gde qu'une autre*/

    for (i=0; i<6; i++)  h[i] = sqrt(h[i]);

    nobtus = 0;
    for (k=0; k<4; k++) {
      for (i=0; i<3; i++) {
        i0 = _MMG5_idir[k][i];
        i1 = _MMG5_idir[k][_MMG5_inxt3[i]];
        i2 = _MMG5_idir[k][_MMG5_inxt3[i+1]];
        if ( h[i0]+h[i1] < 1.2*h[i2] ) {/*1.4 ie une face obtus*/
          nobtus++;
          item[0] = i2;
          item[1] = _MMG5_idir[k][_MMG5_inxt3[i+1]];
        }
      }
    }

    switch(nobtus){
    case 0 :
      break;
    case 1:
      item[0] = iarmax;
      item[1] = _MMG5_isar[iarmax][0];
      return(3);
    case 2:
      item[0] = iarmin;
      item[1] = iarmax;
      return(6);
    default:
      item[0] = iarmin;
      item[1] = iarmax;
      //printf("default obtus %d\n",nobtus);
      return(7);
    }
  }

  /* type 4,5,7 */
  else if ( dd < 0.7*RAPMAX ) {
    naigu = 0;
    for (k=0; k<4; k++) {
      for (i=0; i<3; i++) {
        i0 = _MMG5_idir[k][i];
        i1 = _MMG5_idir[k][_MMG5_inxt3[i]];
        i2 = _MMG5_idir[k][_MMG5_inxt3[i+1]];
        if ( h[i0]+h[i1] > 1.5*h[i2] )  naigu++;/*1.5*/
      }
    }
    switch(naigu){
    case 0 :
      break;
    case 1:
      break;
    case 2:
      item[0] = iarmin;
      return(4);
    case 3:
      /*#warning definir item*/
      return(5);
    default:
      item[0] = iarmin;
      item[1] = iarmax;
      //printf("default aigu\n");
      return(7);
    }
  }
  item[0] = 0;
  return(1);
}
static inline int _MMG3D_swpItem(MMG5_pMesh mesh,  MMG5_pSol met,_MMG5_pBucket bucket,int k,int iar) {
  MMG5_pTetra   pt,pt1;
  MMG5_pxTetra  pxt;
  int           l,list[MMG3D_LMAX+2],lon,iel,nconf,ier;
  double        crit;
  double        OCRIT = 1.01;
  int           nf = 0,j;

  ier = 0;
  pt = &mesh->tetra[k];
  lon = _MMG5_coquil(mesh,k,iar,&list[0]);
  if(lon%2) return(0);
  lon = lon/2;
  if ( lon > 2 ) {
    crit = pt->qual;
    for (l=0; l<lon; l++) {
      iel = list[l] / 6;
      pt1 = &mesh->tetra[iel];
      if(pt1->tag & MG_REQ) break;
      if ( pt1->qual < crit )  crit = pt1->qual;
      if ( pt1->xt ) {
        pxt = &mesh->xtetra[pt1->xt];
        for (j=0; j<4; j++)
          if ( pxt->ftag[j] & MG_BDY )  nf++;
      }
    }
    if(nf > 1) return(0);//printf("on risque de crreerrrr %d\n",nf);
    if(l<lon)  {
      ier = 0;
    } else {
      crit *= OCRIT;
      /* Prevent swap of a ref or tagged edge */
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        if ( pxt->edg[iar] || pxt->tag[iar] ) return(0);
      }
      
      nconf = _MMG5_chkswpgen(mesh,met,k,iar,&lon,list,OCRIT,2);
      if ( nconf ) {
        ier = _MMG5_swpgen(mesh,met,nconf,lon,list,bucket,2);
        if ( ier < 0 ) return(-1);
        else 
          return(ier);
      }           
    }
  }

  return(ier);
}
static inline int _MMG3D_swpalmostall(MMG5_pMesh mesh,  MMG5_pSol met,_MMG5_pBucket bucket,int k,int iar) {
  int           i,ier;
  
  ier = 0;
  for(i=0 ; i<6 ; i++) {
    if(i==iar) continue;
    ier = _MMG3D_swpItem(mesh,met,bucket,k,i);
    if ( ier < 0 ) return(-1);
    else if(ier) 
      return(ier);
  }
  return(ier);
}

static inline int _MMG3D_splitItem(MMG5_pMesh mesh,  MMG5_pSol met,_MMG5_pBucket bucket,
                                   int k,int iar,double OCRIT) {
  MMG5_pTetra   pt;
  double        len;
  double        LLONG2 = 0.1;
  int           ier,j,ip;

  ier = 0;
  pt = &mesh->tetra[k];
  len = _MMG5_lenedg(mesh,met,iar,pt);
  if(len > LLONG2) {
    ier = _MMG5_splitedg(mesh,met,k,iar,OCRIT);
  }

  if(ier) {
    for(j=0 ; j<4 ; j++) {
      if(pt->v[j] == ier) break;	
    }
    assert(j<4);
    ier = _MMG3D_movv_ani(mesh,met,k,j);
  }
  return(ier);
}

static inline int _MMG3D_splitalmostall(MMG5_pMesh mesh,  MMG5_pSol met,_MMG5_pBucket bucket,
                                        int k,int iar) {
  MMG5_pTetra   pt,pt1;
  int           i,l,list[MMG3D_LMAX+2],lon,iel,nconf,ier;
  double        len;
  double        OCRIT=1.01;

  ier = 0;
  
  pt = &mesh->tetra[k];
  /* if(_MMG5_orvolnorm(mesh,k) < 5.e-9) { */
  /*   OCRIT *= 0.5; */
  /* } else */
  /*   OCRIT *= 0.75; */
  
  for(i=0 ; i<6 ; i++) {
    if(i==iar) continue;
    ier = _MMG3D_splitItem(mesh,met,bucket,k,i,OCRIT);
    if(ier) return(ier);
  }           

  return(ier);
}

int MMG3D_opttyp(MMG5_pMesh mesh, MMG5_pSol met,_MMG5_pBucket bucket) {
  MMG5_pTetra    pt,pt1;
  MMG5_pxTetra   pxt;
  double         crit,critswp,len;
  int            k,j,ityp,cs[10],ds[10],item[2],lon,*adja,iadr;
  int            list[MMG3D_LMAX+2],l,iel,ier,i,nd,nconf,ne,npeau;
  int            it,maxit,ntot;
  double         OCRIT = 1.01;
  double         LLONG2 = 0.1;
  int ddebug ;

  //printf("on lance anatet4\n");
  ier = 0;//_MMG5_anatet4(mesh,met,1);
  if(ier) {
    printf("on a splitte %d\n",ier);
 /* memory free */
    _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));
     if ( !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stdout,"  ## Hashing problem. Exit program.\n");
      return(0);
    }
  }

  ntot = 0;
  crit = 0.2 / _MMG5_ALPHAD;

  it = 0;
  maxit = 10;
  do {
    ne = mesh->ne;
    nd = 0;
    memset(cs,0,10*sizeof(int)); 
    memset(ds,0,10*sizeof(int));
    
    for (k=1 ; k<=ne ; k++) {
      pt = &mesh->tetra[k];
      if(!pt->v[0]) continue;
      /* if(pt->qual <= 3.117138e-07) {printf("k %d ityp %d %e npeau %d\n",k,ityp,pt->qual,npeau);
        ddebug = 1;}
        else*/
        ddebug = 0;

      if(pt->qual > crit) continue;
    
      ityp = MMG3D_typelt(mesh,k,item);
      cs[ityp]++;

      /*tet with bdry faces*/
      iadr = 4*(k-1) + 1;
      adja = &mesh->adja[iadr];
      /*optim bdry tetra*/
      npeau = 0;
      for(i=0 ; i<4 ; i++) {
        if(!adja[i]) npeau++;     
      }
      if(npeau>1 && !mesh->info.noswap) {
        if(1 || mesh->info.imprim<-4) printf("%d faces de peau!!!! %d (typ %d) %e\n",npeau,k,ityp,pt->qual);
        /* if ( !_MMG5_split4bar(mesh,met,k,1) ) return(-1); */
        /* else { */
        /*   nd++;  */
        /*   ds[ityp]++; */
        /*   printf("on a splitte\n"); */
        /*   continue; */
        /* } */
      }  
      if(npeau) {       
        ier = 0;//MMG3D_optbdry(mesh,sol,k);
        if(ier) { 
          nd++;
          ds[ityp]++;
          continue;
        }
      }

      switch(ityp) {
      
      case 1:  /* sliver */
      case 3:  /* aileron*/
      case 6:  /* O good face: move away closest vertices */
      case 7:
      default:
        if(mesh->info.noswap) break;
        ier = _MMG3D_swpItem(mesh,met,bucket,k,item[0]);
        if(ddebug) printf("on swp %d ?\n",ier);
        if(ier > 0) {
          nd++;
          ds[ityp]++;
          break;
        } else if(!ier) {
          /*second try to split the biggest edge*/
          if(!mesh->info.noinsert) {
            /* if(_MMG5_orvolnorm(mesh,k) < 5.e-9) { */
            /*   OCRIT *= 0.5; */
            /* } else */
            /*   OCRIT *= 0.75; */
            ier = 0;//_MMG3D_splitItem(mesh,met,bucket,k,item[0],OCRIT);
            if(ddebug) printf("on split %d ?\n",ier);

            if(ier) {
              nd++;
              ds[ityp]++;
              break;
            }
          } /*end noinsert*/

          ier = _MMG3D_swpalmostall(mesh,met,bucket,k,item[0]);
          if(ddebug) printf("on swp2 %d ?\n",ier);

          if(ier > 0) { 
            nd++;
            ds[ityp]++;
            break;
          }
          ier = 0;//_MMG3D_splitalmostall(mesh,met,bucket,k,item[0]);
          if(ddebug) printf("on split2 %d ?\n",ier);

          if(ier > 0) { 
            nd++;
            ds[ityp]++;
            break;
          }
        }
        break;
      case 2: /*chapeau*/
        if((met->size!=1) && _MMG3D_movv_ani(mesh,met,k,item[0])) {
          nd++;
          ds[ityp]++;
        } else {
          for(i=0 ; i<4 ; i++) {
            if(item[0]==i) continue;
            if((met->size!=1) && _MMG3D_movv_ani(mesh,met,k,i)) {
              nd++;
              ds[ityp]++;
              break;
            }
          }
        }
        break;
      }/*end switch*/
    }/*end for k*/
    /* for (k=0; k<=7; k++)
      if ( cs[k] )
        printf("  optim [%d]      = %5d   %5d  %6.2f %%\n",k,cs[k],ds[k],100.0*ds[k]/cs[k]);
    */
    ntot += nd;
  } while (nd && it++<maxit);
          
  
    return(ntot);
  }
