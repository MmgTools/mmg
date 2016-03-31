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
 * \file mmg2d/unused_2d.c
 * \brief A bunch of unused functions (for now)
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "mmg2d.h"

#define EPSD 1.e-10
#define LLONG1 1.9

/*************************************************************************************/
/********************************** Begin solmap_d2.c *************************************/
/*************************************************************************************/

/**
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward the sol structure
 * \return 1 if success
 *
 * Compute isotropic size map according to the mean of the length of the edges
 * passing through a point.
 *
 */
int MMG2_doSol(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria      ptt,pt;
  MMG5_pPoint     p1,p2;
  double          ux,uy,dd;
  int             i,k,ib,ipa,ipb;
  int             MMG_inxtt[5] = {0,1,2,0,1};
  
  sol->np = mesh->np;
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    p1->tagdel = 0;
  }
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !ptt->v[0] )  continue;
    
    for (i=0; i<3; i++) {
      ib  = MMG_inxtt[i+1];
      ipa = ptt->v[i];
      ipb = ptt->v[ib];
      p1  = &mesh->point[ipa];
      p2  = &mesh->point[ipb];
      
      ux  = p1->c[0] - p2->c[0];
      uy  = p1->c[1] - p2->c[1];
      dd  = sqrt(ux*ux + uy*uy);
      
      sol->m[ipa] += dd;
      p1->tagdel++;
      sol->m[ipb] += dd;
      p2->tagdel++;
    }
  }
  
  /* vertex size */
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !p1->tagdel )  {
      sol->m[k] = mesh->info.hmax;
      continue;
    }
    
    sol->m[k] = MG_MIN(mesh->info.hmax,MG_MAX(mesh->info.hmin,sol->m[k] / (double)p1->tagdel));
    p1->tagdel = 0;
  }
  
  /* compute quality */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    pt->qual = MMG2_caltri_in(mesh,sol,pt);
  }
  
  if ( mesh->info.imprim < -4 )
    fprintf(stdout,"     HMIN %f   HMAX %f\n",mesh->info.hmin,mesh->info.hmax);
  return(1);
}

/*************************************************************************************/
/********************************** End solmap_d2.c *************************************/
/*************************************************************************************/

/*************************************************************************************/
/********************************** Begin evalgeom_d2.c *************************************/
/*************************************************************************************/

/* read mesh data */
int MMG2_evalgeom(MMG5_pMesh mesh) {
  MMG5_pTria     pt,pt1;
  MMG5_pPoint    ppt,ppa,ppb;
  double    capx,capy,cbpx,cbpy,alpha,rbound;
  int       *list,k,j,lon,iadr,*adja,nbdry,ibdry[2],ip,i,iel;
  int       ref,nc;
  
  nc = 0;
  rbound = mesh->info.dhd*M_PI/180.;
  /*corners detection*/
  _MMG5_SAFE_MALLOC(list,MMG2D_LMAX,int);
  
  ++mesh->base;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if(!M_EOK(pt)) continue;
    for(j=0 ; j<3 ; j++) {
      ppt = &mesh->point[pt->v[j]];
      if(!(ppt->tag & M_BDRY)) continue;
      if(ppt->tagdel == mesh->base) continue;
      
      lon = MMG2_boulep(mesh,k,j,list);
      assert(lon);
      /*bdry triangles*/
      nbdry = 0;
      ref = mesh->tria[list[1]/3].ref;
      for(i=1 ; i<=lon ; i++) {
        iel = list[i]/3;
        ip  = list[i]%3;
        pt1 = &mesh->tria[iel];
        if(pt1->ref!=ref) continue;
        iadr = 3*(iel-1) + 1;
        adja = &mesh->adja[iadr];
        if(!adja[MMG2_iopp[ip][0]] || (mesh->tria[adja[MMG2_iopp[ip][0]]/3].ref != ref)) {
          
          if(nbdry>=2) fprintf(stdout,"NON MANIFOLD MESH\n");
          if(MMG2_iare[MMG2_iopp[ip][0]][0]==ip)
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][1]];
          else {
            assert(MMG2_iare[MMG2_iopp[ip][0]][1]==ip) ;
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][0]];
          }
          //printf("1) tr %d : %d -- edge %d = %d %d\n",k,iel,MMG2_iopp[ip][0],MMG2_iare[MMG2_iopp[ip][0]][0],
          //                  MMG2_iare[MMG2_iopp[ip][0]][1]);
        }
        if(!adja[MMG2_iopp[ip][1]] || (mesh->tria[adja[MMG2_iopp[ip][1]]/3].ref != ref)) {
          if(nbdry>=2) fprintf(stdout,"NON MANIFOLD MESH\n");
          if(MMG2_iare[MMG2_iopp[ip][1]][0]==ip)
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][1]];
          else {
            assert(MMG2_iare[MMG2_iopp[ip][1]][1]==ip) ;
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][0]];
          }
          //printf("2) tr %d : %d -- edge %d = %d %d\n",k,iel);
        }
      }
      if(nbdry!=2) {
        fprintf(stdout,"NON MANIFOLD DOMAIN (NO CORNERS DETECTION) %d -- vertex %d\n",nbdry,k);
        continue;
      }
      //calcul de l'angle forme par les 3 points
      ppa  = &mesh->point[ibdry[0]];
      ppb  = &mesh->point[ibdry[1]];
      capx = ppt->c[0] - ppa->c[0];
      capy = ppt->c[1] - ppa->c[1];
      cbpx = ppt->c[0] - ppb->c[0];
      cbpy = ppt->c[1] - ppb->c[1];
      alpha = capx*cbpx + capy*cbpy;
      alpha /= sqrt(capx*capx+capy*capy)*sqrt(cbpx*cbpx+cbpy*cbpy);
      alpha = acos(alpha);
      //printf("point %d : %e (= %e)-- %e %e\n",pt->v[j],alpha,alpha*180./M_PI,capx,capy);
      if(alpha < rbound && (!(ppt->tag & M_CORNER)) ) {
        ++nc;
        ppt->tag |= M_CORNER;
      }
      
      ppt->tagdel++;
    }
  }
  
  _MMG5_SAFE_FREE(list);
  
  if ( abs(mesh->info.imprim) > 3 && nc )
    fprintf(stdout,"     %d corners detected\n",nc);
  return(1);
}

/*************************************************************************************/
/********************************** End evalgeom_d2.c *************************************/
/*************************************************************************************/

/*************************************************************************************/
/********************************** Begin cendel_d2.c *************************************/
/*************************************************************************************/


int MMG2_cendel(MMG5_pMesh mesh,MMG5_pSol sol,double declic,int base) {
  MMG5_pTria      pt,pt1;
  pQueue     queue;
  double      crit;
  int       *adja,*list,adj,iadr,i,k,ns,np;
  
  /* queue on quality */
  queue = MMG2_kiuini(mesh,mesh->nt,declic,-1);
  assert(queue);
  list  = (int*)malloc(MMG2D_LMAX*sizeof(int));
  assert(list);
  ns = 0;
  np = 0;
  do {
    k = MMG2_kiupop(queue);
    if ( !k )  break;
    np++;
    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;
    
    /* base internal edges */
    iadr  = 3*(k-1) + 1;
    adja  = &mesh->adja[iadr];
    for (i=0; i<3; i++) {
      adj = adja[i] / 3;
      if ( !adj || pt->ref != mesh->tria[adj].ref )  continue;
      //check required
      if((mesh->point[pt->v[MMG2_iare[i][0]]].tag & M_REQUIRED)
         && (mesh->point[pt->v[MMG2_iare[i][1]]].tag & M_REQUIRED)) {
        continue;
      }
      
      pt1  = &mesh->tria[adj];
      crit = 0.99 * M_MAX(pt->qual,pt1->qual);
      if ( _MMG2_swapdelone(mesh,sol,k,i,crit,list) ) {
        ns++;
        break;
      }
      
    }
  }
  while ( k );
  
  if ( mesh->info.imprim < - 4 )
    fprintf(stdout,"     %7d PROPOSED  %7d SWAPPED\n",np,ns);
  
  MMG2_kiufree(queue);
  free(list);
  return(ns);
}

/*************************************************************************************/
/********************************** End cendel_d2.c *************************************/
/*************************************************************************************/


/*************************************************************************************/
/********************************** Begin split_d2.c *************************************/
/*************************************************************************************/


#define QSEUIL 1e4
#define LSHORT1 0.65

/*insert ip on edge between k1 and adj1/3 */
int MMG2_split(MMG5_pMesh mesh,MMG5_pSol sol,int ip,int k1,int adj1) {
  MMG5_pTria     pt1,pt2,pt3,pt4,ptmp;
  MMG5_pEdge     ped,ped1;
  int       k2,adj2,jel,kel,voy1,voy2,iar1,iar2,iara1,iara2;
  int       *adja,*adja1,*adja2,tmp1,tmp2,piar1,piar2,pvoy1,piara1,piara2,pvoy2;
  int       iadr,tmp,voy,num,newed,num1,num2;
  double    air,cal1,cal2,cal3,cal4,coe,*ca,*cb,*ma,*mb,len;
  
  coe = QSEUIL/ALPHA;
  
  voy2  = adj1%3;
  k2    = adj1/3;
  adja2 = &mesh->adja[3*(k2-1) + 1];
  adj2  = adja2[voy2];
  assert(adj2/3==k1);
  voy1  = adj2%3;
  
  iar1  = MMG2_iare[voy1][0];
  iar2  = MMG2_iare[voy1][1];
  iara1 = MMG2_iare[voy2][0];
  iara2 = MMG2_iare[voy2][1];
  
  pt1  = &mesh->tria[k1];
  pt2  = &mesh->tria[k2];
  
  assert(pt2->v[iara1]==pt1->v[iar2]);
  assert(pt2->v[iara2]==pt1->v[iar1]);
  
  piar1  = pt1->v[iar1];
  piar2  = pt1->v[iar2];
  pvoy1  = pt1->v[voy1];
  piara1 = pt2->v[iara1];
  piara2 = pt2->v[iara2];
  pvoy2  = pt2->v[voy2];
  
  /*test split ok*/
  /*test area > 0*/
  air = MMG2_quickarea(mesh->point[piar2].c,mesh->point[pvoy1].c,mesh->point[ip].c);
  if(air < EPSA) return(0);
  air = MMG2_quickarea(mesh->point[piara1].c,mesh->point[ip].c,mesh->point[pvoy2].c);
  if(air < EPSA) return(0);
  air = MMG2_quickarea(mesh->point[ip].c,mesh->point[pvoy1].c,mesh->point[piar1].c);
  if(air < EPSA) return(0);
  air = MMG2_quickarea(mesh->point[ip].c,mesh->point[piara2].c,mesh->point[pvoy2].c);
  if(air < EPSA) return(0);
  
  /*test qual*/
  ptmp = &mesh->tria[0];
  ptmp->v[0] = piar2;
  ptmp->v[1] = pvoy1;
  ptmp->v[2] = ip;
  cal1       = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal1 > coe) return(0);
  
  ptmp->v[0] = piara1;
  ptmp->v[1] = ip;
  ptmp->v[2] = pvoy2;
  cal2       = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal2 > coe) return(0);
  
  ptmp->v[0] = ip;
  ptmp->v[1] = pvoy1;
  ptmp->v[2] = piar1;
  cal3       = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal3 > coe) return(0);
  
  ptmp->v[0] = ip;
  ptmp->v[1] = piara2;
  ptmp->v[2] = pvoy2;
  cal4       = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal4 > coe) return(0);
  
  /*test length : pvoy1-ip and pvoy2-ip*/
  ca   = &mesh->point[ip].c[0];
  iadr = (ip-1)*sol->size + 1;
  ma   = &sol->m[iadr];
  
  cb  = &mesh->point[pvoy1].c[0];
  iadr = (pvoy1-1)*sol->size + 1;
  mb   = &sol->m[iadr];
  
  len = MMG2_length(ca,cb,ma,mb);
  //printf("edg %d %d : %e\n",pvoy1,ip,len);
  if(len < LSHORT1) return(0);
  
  cb  = &mesh->point[pvoy2].c[0];
  iadr = (pvoy2-1)*sol->size + 1;
  mb   = &sol->m[iadr];
  
  len = MMG2_length(ca,cb,ma,mb);
  //printf("edg %d %d : %e\n",pvoy2,ip,len);
  if(len < LSHORT1) return(0);
  
  /*check*/
  cb  = &mesh->point[piara2].c[0];
  iadr = (piara2-1)*sol->size + 1;
  mb   = &sol->m[iadr];
  
  len = MMG2_length(ca,cb,ma,mb);
  //printf("edg %d (%d)  %d : %e\n",piar1,piara2,ip,len);
  if(len < LSHORT1) return(0);
  
  cb  = &mesh->point[piar2].c[0];
  iadr = (piar2-1)*sol->size + 1;
  mb   = &sol->m[iadr];
  
  len = MMG2_length(ca,cb,ma,mb);
  //printf("edg %d %d : %e\n",piar2,ip,len);
  if(len < LSHORT1) return(0);
  
  pt1->v[0] = piar2;
  pt1->v[1] = pvoy1;
  pt1->v[2] = ip;
  pt1->qual = cal1;
  
  pt2->v[0] = piara1;
  pt2->v[1] = ip;
  pt2->v[2] = pvoy2;
  pt2->qual = cal2;
  
  
  jel  = _MMG2D_newElt(mesh);
  if ( !jel ) {
    _MMG5_TRIA_REALLOC(mesh,jel,mesh->gap,
                       printf("  ## Error: unable to allocate a new element.\n");
                       _MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");
                       exit(EXIT_FAILURE));
    pt1  = &mesh->tria[k1];
    pt2  = &mesh->tria[k2];
    ptmp = &mesh->tria[0];
    adja2 =  &mesh->adja[3*(k2-1) + 1];
  }
  kel  = _MMG2D_newElt(mesh);
  if ( !kel ) {
    _MMG5_TRIA_REALLOC(mesh,kel,mesh->gap,
                       printf("  ## Error: unable to allocate a new element.\n");
                       _MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");
                       exit(EXIT_FAILURE));
    pt1  = &mesh->tria[k1];
    pt2  = &mesh->tria[k2];
    ptmp = &mesh->tria[0];
    adja2 =  &mesh->adja[3*(k2-1) + 1];
  }
  pt3  = &mesh->tria[jel];
  pt3->v[0] = ip;
  pt3->v[1] = pvoy1;
  pt3->v[2] = piar1;
  pt3->ref  = pt1->ref;
  pt3->qual = cal3;
  
  pt4  = &mesh->tria[kel];
  pt4->v[0] = ip;
  pt4->v[1] = piara2;
  pt4->v[2] = pvoy2;
  pt4->ref  = pt2->ref;
  pt4->qual = cal4;
  
  /*adj*/
  adja1 = &mesh->adja[3*(k1-1) + 1];
  /*printf("adj1 : %d %d %d\n",adja1[0]/3,adja1[1]/3,adja1[2]/3);
   printf("adj2 : %d %d %d\n",adja2[0]/3,adja2[1]/3,adja2[2]/3);
   */
  adja = &mesh->adja[3*(jel-1) + 1];
  adja[0] = adja1[iar2];
  pt3->edg[0] = pt1->edg[iar2];
  tmp = 3*(adja1[iar2]/3 - 1) + 1;
  voy = adja1[iar2]%3;
  if(adja1[iar2]) (&mesh->adja[tmp])[voy] = 3*jel + 0;
  adja[1] = 3*kel + 2;
  num = 0;
  if(pt1->edg[voy1]) {
    /*edge creation*/
    /*split edge piar1 piar2 */
    // #warning same tangent
    num = pt1->edg[voy1];
    assert(num);
    ped = &mesh->edge[num];
    newed = _MMG5_newEdge(mesh);
    if ( !newed ) {
      _MMG5_EDGE_REALLOC(mesh,newed,mesh->gap,
                         printf("  ## Error: unable to allocate a new edge.\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         printf("  Exit program.\n");
                         exit(EXIT_FAILURE));
      ped = &mesh->edge[num];
    }
    
    ped1 = &mesh->edge[newed];
    memcpy(ped1,ped,sizeof(MMG5_Edge));
    ped1->a = ip;
    
    ped->b = ip;
    
    if(ped->a==piar1) {
      pt3->edg[1] = num;
      pt4->edg[2] = num;
      // printf("on met num pt3  %d %d\n",pt3->v[MMG2_iare[1][0]],pt3->v[MMG2_iare[1][1]]);
      //printf("on met num pt4 %d %d\n",pt4->v[MMG2_iare[2][0]],pt4->v[MMG2_iare[2][1]]);
      
    } else {
      pt3->edg[1] = newed;
      pt4->edg[2] = newed;
    }
  }
  adja[2] = 3*k1  + 0;
  pt3->edg[2] = 0;
  
  if(kel) {
    adja = &mesh->adja[3*(kel-1) + 1];
    adja[0] = adja2[iara1];
  }
  pt4->edg[0] = pt2->edg[iara1];
  
  
  if(adja2[iara1]) (&mesh->adja[3*(adja2[iara1]/3-1)+1])[adja2[iara1]%3] = 3*kel + 0;
  adja[1] = 3*k2  + 0;
  pt4->edg[1] = 0;
  adja[2] = 3*jel + 1;
  if(pt1->edg[voy1]) assert(pt4->edg[2]);
  
  tmp1 = adja1[iar1];
  num1 = pt1->edg[iar1];
  tmp2 = adja2[iara2];
  num2 = pt2->edg[iara2];
  adja1[0] = 3*jel + 2;
  pt1->edg[0] = 0;
  adja1[1] = 3*k2  + 2;
  
  adja1[2] = tmp1;
  pt1->edg[2] = num1;
  if(tmp1)
    (&mesh->adja[3*(tmp1/3-1)+1])[tmp1%3] = 3*k1 + 2;
  
  
  adja2[0] = 3*kel + 1;
  pt2->edg[0] = 0;
  adja2[1] = tmp2;
  pt2->edg[1] = num2;
  if(tmp2)
    (&mesh->adja[3*(tmp2/3-1)+1])[tmp2%3] = 3*k2 + 1;
  adja2[2] = 3*k1 + 1;
  if(num) {
    if(ped->a==piar1) {
      pt1->edg[1] = newed;
      pt2->edg[2] = newed;
    } else {
      pt1->edg[1] = num;
      pt2->edg[2] = num;
      //printf("on met num k1 %d %d\n",pt1->v[MMG2_iare[1][0]],pt1->v[MMG2_iare[1][1]]);
      //printf("on met num k2 %d %d\n",pt2->v[MMG2_iare[2][0]],pt2->v[MMG2_iare[2][1]]);
      
    }
  } else {
    pt1->edg[1] = 0;
    pt2->edg[2] = 0;
  }
  
  
  
  if(MMG2D_callbackinsert)
    MMG2D_callbackinsert((int) ip,(int) k1,(int) k2,(int)jel,(int) kel);
  
  return(1);
}

/*insert ip on edge in k1 */
int MMG2_splitbdry(MMG5_pMesh mesh,MMG5_pSol sol,int ip,int k1,int voy1,double *tang) {
  MMG5_pTria     pt1,pt3,ptmp;
  MMG5_pEdge     ped,ped1;
  MMG5_pPoint    ppt;
  int       jel,iar1,iar2,i,num,newed,num1,num2;
  int       *adja,*adja1,tmp1,piar1,piar2,pvoy1,ref1,ref2;
  double    air,cal1,cal2,coe;
  
  coe = QSEUIL/ALPHA;
  
  iar1  = MMG2_iare[voy1][0];
  iar2  = MMG2_iare[voy1][1];
  
  pt1  = &mesh->tria[k1];
  //printf("tr1 %d %d %d -- voy %d : %d %d \n",pt1->v[0],pt1->v[1],pt1->v[2],voy1,iar1,iar2);
  
  piar1  = pt1->v[iar1];
  piar2  = pt1->v[iar2];
  pvoy1  = pt1->v[voy1];
  
  /*test split ok*/
  air = MMG2_quickarea(mesh->point[piar2].c,mesh->point[pvoy1].c,mesh->point[ip].c);
  if(air < EPSA) return(0);
  air = MMG2_quickarea(mesh->point[ip].c,mesh->point[pvoy1].c,mesh->point[piar1].c);
  if(air < EPSA) return(0);
  
  ptmp = &mesh->tria[0];
  ptmp->v[0] = piar2;
  ptmp->v[1] = pvoy1;
  ptmp->v[2] = ip;
  cal1 = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal1 > coe) return(0);
  
  ptmp->v[0] = ip;
  ptmp->v[1] = pvoy1;
  ptmp->v[2] = piar1;
  cal2 = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal2 > coe) return(0);
  
  
  pt1->v[0] = piar2;
  pt1->v[1] = pvoy1;
  pt1->v[2] = ip;
  pt1->qual = cal1;
  
  jel  = _MMG2D_newElt(mesh);
  if ( !jel ) {
    _MMG5_TRIA_REALLOC(mesh,jel,mesh->gap,
                       printf("  ## Error: unable to allocate a new element.\n");
                       _MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");
                       exit(EXIT_FAILURE));
    pt1  = &mesh->tria[k1];
    ptmp = &mesh->tria[0];
    
  }
  pt3  = &mesh->tria[jel];
  pt3->v[0] = ip;
  pt3->v[1] = pvoy1;
  pt3->v[2] = piar1;
  pt3->ref  = pt1->ref;
  pt3->qual = cal2;
  
  /*adj*/
  adja1 = &mesh->adja[3*(k1-1) + 1];
  
  adja = &mesh->adja[3*(jel-1) + 1];
  adja[0] = adja1[iar2];
  if(adja1[iar2])
    (&mesh->adja[3*(adja1[iar2]/3-1)+1])[adja1[iar2]%3] = 3*jel + 0;
  adja[1] = 0;
  adja[2] = 3*k1  + 0;
  
  tmp1 = adja1[iar1];
  adja1[0] = 3*jel + 2;
  adja1[1] = 0;
  adja1[2] = tmp1;
  if(tmp1)
    (&mesh->adja[3*(tmp1/3-1)+1])[tmp1%3] = 3*k1 + 2;
  
  /*si dep alors on met la moy des dep dans ip*/
  if(mesh->info.lag >=0) {
    printf("  ## Error: option not implemented: merge option 9\n");
    exit(EXIT_FAILURE);
    /* mesh->disp.mv[2*(ip-1) + 1 + 0] = 0.5*(mesh->disp.mv[2*(piar1-1) + 1 + 0] + mesh->disp.mv[2*(piar2-1) + 1 + 0]); */
    /* mesh->disp.mv[2*(ip-1) + 1 + 1] = 0.5*(mesh->disp.mv[2*(piar1-1) + 1 + 1] + mesh->disp.mv[2*(piar2-1) + 1 + 1]);     */
    /* d1 = mesh->disp.mv[2*(ip-1) + 1 + 0]*mesh->disp.mv[2*(ip-1) + 1 + 0] */
    /*   + mesh->disp.mv[2*(ip-1) + 1 + 1]*mesh->disp.mv[2*(ip-1) + 1 + 1]; */
    /* if ( d1 > 1e-24 )  mesh->point[ip].tag  |= M_MOVE; */
  }
  
  /*propagation des ref de peau*/
  ref1 = mesh->point[piar1].ref;
  ref2 = mesh->point[piar2].ref;
  if( ref1 || ref2 ) mesh->point[ip].ref = ref1 > ref2 ? ref1 : ref2;
  
  /*split edge piar1 piar2 if exist */
  num = pt1->edg[voy1];
  assert(num);
  ped = &mesh->edge[num];
  
  newed = _MMG5_newEdge(mesh);
  if ( !newed ) {
    _MMG5_EDGE_REALLOC(mesh,newed,mesh->gap,
                       printf("  ## Error: unable to allocate a new edge.\n");
                       _MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");
                       exit(EXIT_FAILURE));
    ped = &mesh->edge[num];
  }
  
  ped1 = &mesh->edge[newed];
  memcpy(ped1,ped,sizeof(MMG5_Edge));
  ped1->a = ip;
  ppt = &mesh->point[ip];
  for(i=0 ; i<2 ; i++)
    ppt->n[i] = tang[i];
  
  ped->b = ip;
  
  num1 = pt1->edg[iar1];
  num2 = pt1->edg[iar2];
  pt1->edg[2] = num1;
  pt3->edg[0] = num2;
  
  if(ped->a==piar1) {
    pt3->edg[1] = num;
    pt1->edg[1] = newed;
  } else {
    pt3->edg[1] = newed;
    pt1->edg[1] = num;
  }
  
  pt1->edg[0]  = 0;
  //end add edge
  
  return(1);
}

/*************************************************************************************/
/********************************** End split_d2.c *************************************/
/*************************************************************************************/

/*************************************************************************************/
/********************************** Begin colpoi_d2.c *************************************/
/*************************************************************************************/


/*collapse edge ppb-->ppa*/
int MMG2_colpoi(MMG5_pMesh mesh, MMG5_pSol sol,int iel,int iar,int ia,int ib,double coe) {
  MMG5_pTria     pt,pt1;
  MMG5_pPoint    ppa,ppb,pp1,pp2,ppa1,ppb1;
  MMG5_pEdge     ped;
  int       pib,pia,jel,iadr,a1,v1,*adja,iaa,voy,a,a2,v2,adj,adj1,adjj1;
  int     *list,lon,i,kel,num,ed,i1,i2;
  double    declic,*cal,air,coor[2],solu[3],*c1,*c2,*m1,*m2,len;
  double    cbound,capx,capy,cbpx,cbpy,alpha;
  int       nbdry,ref,ip,iadri,*adjai,ibdry[2],it1;
  
  pt  = &mesh->tria[iel];
  pib = pt->v[ib];
  pia = pt->v[ia];
  ppa = &mesh->point[pia];
  ppb = &mesh->point[pib];
  
  if(ppb->tag & M_BDRY) return(0);
  
  iadr = 3*(iel-1) + 1;
  adja = &mesh->adja[iadr];
  jel  = adja[iar]/3;
  
  list  = (int*)malloc(MMG2D_LMAX*sizeof(int));
  assert(list);
  
  lon = MMG2_boulep(mesh,iel,ib,list);
  if(!lon) {
    free(list);
    return(0);
  }
  
  //if vertex between two SD, check geometry criterion
  if(ppb->tag & M_SD) {
    if(!(ppa->tag & M_SD)) {
      free(list);
      return(0);
    }
    /*check geom*/
    cbound = 178.*M_PI/180.;
    nbdry = 0;
    ref = mesh->tria[list[lon]/3].ref;
    for(i=1 ; i<=lon ; i++) {
      kel = list[i]/3;
      ip  = list[i]%3;
      pt1 = &mesh->tria[kel];
      iadri = 3*(kel-1) + 1;
      adjai =  &mesh->adja[iadri];
      if(pt1->ref!=ref) {
        it1 = (i==1)?lon:i-1;
        if(adjai[MMG2_iopp[ip][0]]/3 == list[it1]/3 ) {
          //assert(nbdry<2);
          if(MMG2_iare[MMG2_iopp[ip][0]][0]==ip)
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][1]];
          else {
            assert(MMG2_iare[MMG2_iopp[ip][0]][1]==ip) ;
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][0]];
          }
        } else {
          assert(adjai[MMG2_iopp[ip][1]]/3 == list[it1]/3);
          if(MMG2_iare[MMG2_iopp[ip][1]][0]==ip)
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][1]];
          else {
            assert(MMG2_iare[MMG2_iopp[ip][1]][1]==ip) ;
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][0]];
          }
        }
      }
      ref = pt1->ref;
    }
    if(nbdry!=2) return(0);
    //assert(nbdry==2); //sinon non manifold
    //calcul de l'angle forme par les 3 points
    ppa1  = &mesh->point[ibdry[0]];
    ppb1  = &mesh->point[ibdry[1]];
    capx = ppb->c[0] - ppa1->c[0];
    capy = ppb->c[1] - ppa1->c[1];
    cbpx = ppb->c[0] - ppb1->c[0];
    cbpy = ppb->c[1] - ppb1->c[1];
    alpha = capx*cbpx + capy*cbpy;
    alpha /= sqrt(capx*capx+capy*capy)*sqrt(cbpx*cbpx+cbpy*cbpy);
    alpha = acos(alpha);
    
    if(alpha < cbound ) {
      free(list);
      return(0);
    }
  }
  
  cal  = (double*)malloc((lon+1)*sizeof(double));
  assert(cal);
  
  /*simu colps ppb-->ppa*/
  memcpy(coor, ppb->c,2*sizeof(double));
  memcpy(ppb->c,ppa->c,2*sizeof(double));
  //memcpy(solu,&sol->m[3*(pib-1) + 1],sol->size*sizeof(double));
  for(i=0 ; i<sol->size ; i++) {
    solu[i] = sol->m[sol->size*(pib-1) + 1 + i];
    sol->m[sol->size*(pib-1) + 1 + i] = sol->m[sol->size*(pia-1) + 1 + i];
  }
  
  //memcpy(&sol->m[3*(pib-1) + 1],&sol->m[3*(pia-1) + 1],sol->size*sizeof(double));
  
  
  /*check config*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    if(kel==jel) continue;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    air  = MMG2_quickarea(mesh->point[pt1->v[0]].c,mesh->point[pt1->v[1]].c,
                          mesh->point[pt1->v[2]].c);
    if(air < EPSA) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*(pib-1) + 1],solu,sol->size*sizeof(double));
      free(cal);
      free(list);
      return(0);
    }
    declic = coe*pt1->qual;
    cal[i] = MMG2_caltri_in(mesh,sol,pt1);
    if (cal[i] > declic) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*(pib-1) + 1],solu,sol->size*sizeof(double));
      free(cal);
      free(list);
      return(0);
    }
  }
  
  /*check lengths*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    if(kel==jel) continue;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    /*check first edge containing ib*/
    ed   = (voy+2)%3;
    i1   = pt1->v[MMG2_iare[ed][0]];
    i2   = pt1->v[MMG2_iare[ed][1]];
    pp1  = &mesh->point[i1];
    pp2  = &mesh->point[i2];
    c1   = &pp1->c[0];
    c2   = &pp2->c[0];
    iadr = (i1-1)*sol->size + 1;
    m1   = &sol->m[iadr];
    iadr = (i2-1)*sol->size + 1;
    m2   = &sol->m[iadr];
    
    len = MMG2_length(c1,c2,m1,m2);
    if (len > LLONG1) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*(pib-1) + 1],solu,sol->size*sizeof(double));
      free(cal);
      free(list);
      return(0);
    }
  }
  
  /*update tria*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    if(kel==jel) continue;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    
    assert(pt1->v[voy]==pib);
    pt1->v[voy] = pia;
    
    /*edge*/
    num  = pt1->edg[MMG2_iare[voy][0]];
    if(num) {
      ped = &mesh->edge[num];
      if(ped->a==pib) ped->a = pia;
      if(ped->b==pib) ped->b = pia;
    }
    num = pt1->edg[ MMG2_iare[voy][1]];
    if(num) {
      ped = &mesh->edge[num];
      if(ped->a==pib) ped->a = pia;
      if(ped->b==pib) ped->b = pia;
      
    }
    pt1->qual = cal[i];
  }
  
  /*adj of iel*/
  adj  = adja[ib];
  a1   = adja[ia]/3;
  v1   = adja[ia]%3;
  
  iaa  = MMG2_iare[adja[iar]%3][1];
  
  adjj1 = MMG2_iare[adja[iar]%3][0];
  
  /*adj of jel*/
  iadr = 3*(jel-1) + 1;
  adja = &mesh->adja[iadr];
  adj1 = adja[adjj1];
  a2   = adja[iaa]/3;
  v2   = adja[iaa]%3;
  
  
  /*adja*/
  iadr = 3*(a1-1) + 1;
  a    = adj/3;
  voy  = adj%3;
  adja = &mesh->adja[iadr];
  adja[v1] = 3*a + voy;
  
  if(a) {
    iadr = 3*(a-1) + 1;
    adja = &mesh->adja[iadr];
    adja[voy] = 3*a1 + v1;
  }
  /*update edge*/
  num = pt->edg[ib];
  if(num) {
    mesh->tria[a1].edg[v1] = num;
  }
  
  iadr = 3*(a2-1) + 1;
  a    = adj1/3;
  voy  = adj1%3;
  if(a2) {
    adja = &mesh->adja[iadr];
    adja[v2] = 3*a + voy;
  }
  if(a) {
    iadr = 3*(a-1) + 1;
    adja = &mesh->adja[iadr];
    adja[voy] = 3*a2 + v2;
  }
  /*update edge*/
  pt1  = &mesh->tria[jel];
  num = pt1->edg[adjj1];
  if(a2 && num) {
    if(!((mesh->edge[num].a==mesh->tria[a2].v[MMG2_iare[v2][0]] || mesh->edge[num].a==mesh->tria[a2].v[MMG2_iare[v2][1]])
         && (mesh->edge[num].b==mesh->tria[a2].v[MMG2_iare[v2][0]] ||mesh->edge[num].b==mesh->tria[a2].v[MMG2_iare[v2][1]]))) {
      /* printf("on a un soucis 1\n"); */
      /* printf("pnum %d %d dans %d %d %d\n",mesh->edge[num].a,mesh->edge[num].b,pt1->v[0],pt1->v[1],pt1->v[2]); */
      /* printf("edgea %d %d\n",mesh->tria[a2].v[MMG2_iare[v2][0]],mesh->tria[a2].v[MMG2_iare[v2][1]]); */
      if(mesh->info.imprim > 6)
        fprintf(stdout," ## Warning: bad configuration for collapse\n");
    }
    mesh->tria[a2].edg[v2] = num;
  }
  num = pt1->edg[iaa];
  if(a && num) {
    /* printf("tr %d : %d %d %d et num %d\n",a,mesh->tria[a].edg[0],mesh->tria[a].edg[1],mesh->tria[a].edg[2],num); */
    /* printf("pnum %d %d dans %d %d %d\n",mesh->edge[num].a,mesh->edge[num].b,pt1->v[0],pt1->v[1],pt1->v[2]); */
    /* printf("edgea %d %d\n",mesh->tria[a].v[MMG2_iare[voy][0]],mesh->tria[a].v[MMG2_iare[voy][1]]); */
    if(!((mesh->edge[num].a==mesh->tria[a].v[MMG2_iare[voy][0]] || mesh->edge[num].a==mesh->tria[a].v[MMG2_iare[voy][1]])
         && (mesh->edge[num].b==mesh->tria[a].v[MMG2_iare[voy][0]] || mesh->edge[num].b==mesh->tria[a].v[MMG2_iare[voy][1]]))) {
      /* printf("on a un soucis\n"); */
      /* printf("pnum %d %d dans %d %d %d\n",mesh->edge[num].a,mesh->edge[num].b,pt1->v[0],pt1->v[1],pt1->v[2]); */
      /* printf("edgea %d %d\n",mesh->tria[a].v[MMG2_iare[voy][0]],mesh->tria[a].v[MMG2_iare[voy][1]]); */
      if(mesh->info.imprim > 6) fprintf(stdout," ## Warning: bad configuration for collapse\n");
    }
    mesh->tria[a].edg[voy] = num;
  }
  
  num = pt->edg[iar];
  if(num) {
    _MMG5_delEdge(mesh,num);
  }
  
  
  _MMG2D_delElt(mesh,iel);
  _MMG2D_delElt(mesh,jel);
  memcpy(ppb->c,coor,3*sizeof(double));
  memcpy(&sol->m[sol->size*(pib-1) + 1],solu,sol->size*sizeof(double));
  
  free(list);
  free(cal);
  return(1);
}

/**
 * \brief return 1 if the edge does not verify hausd criterion
 */
int MMG2_chkedg(MMG5_pMesh mesh, MMG5_pPoint ppa,MMG5_pPoint ppb) {
  double t0[2],t1[2];
  double l,ux,uy,cosn,ps;
  int    i;
  
  ux = ppa->c[0] - ppb->c[0];
  uy = ppa->c[1] - ppb->c[1];
  l = ux*ux + uy*uy;
  l = sqrt(l);
  
  //with tangent, compute control point
  for(i=0 ; i<2 ; i++) {
    t0[i] = l*ppa->n[i];
    t1[i] = l*ppb->n[i];
  }
  if(ppa->tag & M_CORNER) {
    for(i=0 ; i<2 ; i++)
      t0[i] = ppa->c[i] - ppb->c[i];
  }
  if(ppb->tag & M_CORNER) {
    for(i=0 ; i<2 ; i++)
      t1[i] = ppa->c[i] - ppb->c[i];
  }
  /*check if t0 has the same sens of vect(P0P1)*/
  if(t0[0]/(ppb->c[0]-ppa->c[0]) < 0 || t0[1]/(ppb->c[1]-ppa->c[1])<0) {
    //printf("t0/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
    for(i=0 ; i<2 ; i++) {
      t0[i] *= -1;
    }
  }
  /*check if t1 has the opposite sens of vect(P0P1)*/
  if(t1[0]/(ppb->c[0]-ppa->c[0]) > 0 || t1[1]/(ppa->c[1]-ppb->c[1])>0) {
    //printf("t1/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
    for(i=0 ; i<2 ; i++) {
      t1[i] *= -1;
    }
  }
  //compute the distance between mid point and curve (with angle between edge and tang)
  ps = t0[0]*ux + t0[1]*uy;
  ps /= l*l;
  cosn = ps*ps ;
  cosn *= fabs(1.0-cosn);
  
  cosn *= (0.25*l);
  if(cosn > mesh->info.hausd*mesh->info.hausd) return(1);
  //idem for ppb
  ps = -(t1[0]*ux + t1[1]*uy);
  ps /= l*l;
  cosn = ps*ps ;
  cosn *= fabs(1.0-cosn);
  cosn *= (0.25*l);
  if(cosn > mesh->info.hausd*mesh->info.hausd) return(1);
  
  
  
  return(0);
}

/**
 * \return -1 if fail, 0 if we can't collapse, 1 otherwise.
 *
 * collapse edge ppb-->ppa and ppbppa is boundary edge
 *
 */
int MMG2_colpoibdry(MMG5_pMesh mesh, MMG5_pSol sol,int iel,int iar,int ia,int ib,double coe) {
  MMG5_pTria     pt,pt1;
  MMG5_pEdge     ped;
  MMG5_pPoint    ppa,ppb,pp1,pp2,ppa1,ppb1;
  int       pib,pia,jel,iadr,a1,v1,*adja,voy,a,adj;
  int     *list,lon,i,kel,num,i1,i2,ed;
  double    declic,*cal,air,coor[2],solu[3],*c1,*c2,*m1,*m2,len;
  //double    capx,capy,cbpx,cbpy,alpha,cbound;
  int       iadri,*adjai,nbdry,ibdry[2],ip;
  pt  = &mesh->tria[iel];
  pib = pt->v[ib];
  pia = pt->v[ia];
  ppa = &mesh->point[pia];
  ppb = &mesh->point[pib];
  
  if (ppb->tag & M_CORNER) return(0);
  iadr = 3*(iel-1) + 1;
  adja = &mesh->adja[iadr];
  jel  = adja[iar]/3;
  assert(!jel);
  _MMG5_SAFE_MALLOC(list,MMG2D_LMAX,int);
  
  lon = MMG2_boulep(mesh,iel,ib,list);
  if(!lon) {
    _MMG5_SAFE_FREE(list);
    return(0);
  }
  
  /*check geom*/
  nbdry = 0;
  for(i=1 ; i<=lon ; i++) {
    kel = list[i]/3;
    ip  = list[i]%3;
    pt1 = &mesh->tria[kel];
    iadri = 3*(kel-1) + 1;
    adjai = &mesh->adja[iadri];
    if(!adjai[MMG2_iopp[ip][0]]) {
      assert(nbdry<2);
      if(MMG2_iare[MMG2_iopp[ip][0]][0]==ip)
        ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][1]];
      else {
        assert(MMG2_iare[MMG2_iopp[ip][0]][1]==ip) ;
        ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][0]];
      }
    }
    if(!adjai[MMG2_iopp[ip][1]]) {
      assert(nbdry<2); //sinon non manifold
      if(MMG2_iare[MMG2_iopp[ip][1]][0]==ip)
        ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][1]];
      else {
        assert(MMG2_iare[MMG2_iopp[ip][1]][1]==ip) ;
        ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][0]];
      }
    }
  }
  assert(nbdry==2); //sinon non manifold
  /*first check that the two edges verify the hausd criterion*/
  ppa1  = &mesh->point[ibdry[0]];
  ppb1  = &mesh->point[ibdry[1]];
  
  if(MMG2_chkedg(mesh,ppb,ppa1))   {
    _MMG5_SAFE_FREE(list);
    return(0);
  }
  if(MMG2_chkedg(mesh,ppb,ppb1))  {
    _MMG5_SAFE_FREE(list);
    return(0);
  }
  
  /*second check that the new edge verify the hausd criteron*/
  if(MMG2_chkedg(mesh,ppb1,ppa1))  {
    _MMG5_SAFE_FREE(list);
    return(0);
  }
  
  /* //comment from here */
  /*   //calcul de l'angle forme par les 3 points  */
  /*   capx = ppb->c[0] - ppa1->c[0];  */
  /*   capy = ppb->c[1] - ppa1->c[1];  */
  /*   cbpx = ppb->c[0] - ppb1->c[0];  */
  /*   cbpy = ppb->c[1] - ppb1->c[1];  */
  /*   alpha = capx*cbpx + capy*cbpy; */
  /*   alpha /= sqrt(capx*capx+capy*capy)*sqrt(cbpx*cbpx+cbpy*cbpy);   */
  /*   alpha = acos(alpha); */
  /*   //printf("point %d : %e (= %e)-- %e %e\n",pt->v[j],alpha,alpha*180./M_PI,capx,capy); */
  
  /*   if(alpha < cbound ) { */
  /*     free(list); */
  /*     return(0); */
  /*     //to here */
  /*   } else */
  if(lon > 100) {
    _MMG5_SAFE_FREE(list);
    return(0);
  }
  _MMG5_SAFE_CALLOC(cal,lon+1,double);
  
  /*simu colps ppb-->ppa*/
  memcpy(coor, ppb->c,2*sizeof(double));
  memcpy(ppb->c,ppa->c,2*sizeof(double));
  memcpy(solu,&sol->m[sol->size*(pib-1) + 1],sol->size*sizeof(double));
  memcpy(&sol->m[sol->size*(pib-1) + 1],&sol->m[sol->size*(pia-1) + 1],sol->size*sizeof(double));
  
  /*check config*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    air  = MMG2_quickarea(mesh->point[pt1->v[0]].c,mesh->point[pt1->v[1]].c,
                          mesh->point[pt1->v[2]].c);
    if(air < EPSA) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*(pib-1) + 1],solu,sol->size*sizeof(double));
      _MMG5_SAFE_FREE(cal);
      _MMG5_SAFE_FREE(list);
      return(0);
    }
    declic = coe*pt1->qual;
    cal[i] = MMG2_caltri_in(mesh,sol,pt1);
    if (cal[i] > declic) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*(pib-1) + 1],solu,sol->size*sizeof(double));
      _MMG5_SAFE_FREE(cal);
      _MMG5_SAFE_FREE(list);
      return(0);
    }
  }
  
  /*check lengths*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    if(kel==jel) continue;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    /*check second edge containing ib*/
    ed   = (voy+2)%3;
    i1   = pt1->v[MMG2_iare[ed][0]];
    i2   = pt1->v[MMG2_iare[ed][1]];
    pp1  = &mesh->point[i1];
    pp2  = &mesh->point[i2];
    c1   = &pp1->c[0];
    c2   = &pp2->c[0];
    iadr = (i1-1)*sol->size + 1;
    m1   = &sol->m[iadr];
    iadr = (i2-1)*sol->size + 1;
    m2   = &sol->m[iadr];
    
    len = MMG2_length(c1,c2,m1,m2);
    if (len > LLONG1) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*(pib-1) + 1],solu,sol->size*sizeof(double));
      _MMG5_SAFE_FREE(cal);
      _MMG5_SAFE_FREE(list);
      return(0);
    }
  }
  /*update tria*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    assert(pt1->v[voy]==pib);
    pt1->v[voy] = pia;
    
    pt1->qual = cal[i];
    
    /*edge*/
    num = pt1->edg[ MMG2_iare[voy][0]];
    if(num) {
      ped = &mesh->edge[num];
      if(ped->a==pib) ped->a = pia;
      if(ped->b==pib) ped->b = pia;
    }
    num = pt1->edg[ MMG2_iare[voy][1]];
    if(num) {
      ped = &mesh->edge[num];
      if(ped->a==pib) ped->a = pia;
      if(ped->b==pib) ped->b = pia;
    }
  }
  
  adj  = adja[ib];
  a1   = adja[ia]/3;
  v1   = adja[ia]%3;
  pt1  = &mesh->tria[a1];
  
  /*adja*/
  a    = adj/3;
  voy  = adj%3;
  if(a1) {
    iadr = 3*(a1-1) + 1;
    adja = &mesh->adja[iadr];
    adja[v1] = 3*a + voy;
  }
  if(a) {
    iadr = 3*(a-1) + 1;
    adja = &mesh->adja[iadr];
    adja[voy] = 3*a1 + v1;
  }
  
  /*del edge pib pia if exist*/
  num = pt->edg[iar];
  if(!num) {
    /* printf("la edge %d %d iar %d tr %d\n",pia,pib,iar,iel); */
    /* printf("pt %d %d %d\n",pt->v[0],pt->v[1],pt->v[2]); */
    /* printf("pt->ned %d %d %d\n",pt->edg[0],pt->edg[1],pt->edg[2]); */
    fprintf(stdout," ## Error: problem with an edge."
            " Check your data and/or report the bug\n");
    return(-1);
  }
  assert(num);
  _MMG5_delEdge(mesh,num);
  
  /*check if tr iel has other edge*/
  if(pt->edg[ia]) {
    assert(a);
    ped = &mesh->edge[pt->edg[ia]];
    if(ped->a==pib)
      ped->a = pia;
    else {
      assert(ped->b==pib);
      ped->b=pia;
    }
    mesh->tria[a].edg[voy]=pt->edg[ia];
  }
  if(pt->edg[ib]) {
    assert(a1);
    ped = &mesh->edge[pt->edg[ib]];
    mesh->tria[a1].edg[v1]=pt->edg[ib];
  }
  _MMG2D_delElt(mesh,iel);
  memcpy(ppb->c,coor,3*sizeof(double));
  memcpy(&sol->m[sol->size*(pib-1) + 1],solu,sol->size*sizeof(double));
  
  _MMG5_SAFE_FREE(cal);
  _MMG5_SAFE_FREE(list);
  
  // if ( !MG2_chkmsh(mesh,0) ) exit(EXIT_FAILURE);
  
  return(1);
}

/*************************************************************************************/
/********************************** End colpoi_d2.c *************************************/
/*************************************************************************************/

/*************************************************************************************/
/********************************** Begin hash_d2.c *************************************/
/*************************************************************************************/



/**
 * \param mesh pointer toward the mesh structure.
 * \param ip1,ip2,ip3 integer
 * \param t  : computed tangent at vertex ip2
 *
 * Compute the tangent at vertex ip2
 *
 */
int MMG2_computetangent(MMG5_pMesh mesh,int ip1,int ip2,int ip3,double *t) {
  double c2[2],c1[2],c3[2],dd1,dd2;
  double dd,n1[2],n2[2],n[2];
  //  double t1[2],t2[2],theta,kappa;
  int    i;
  
  for(i=0 ; i<2 ; i++) {
    c1[i] =  mesh->point[ip1].c[i];
    c2[i] =  mesh->point[ip2].c[i];
    c3[i] =  mesh->point[ip3].c[i];
  }
  if(mesh->point[ip2].tag & M_CORNER) {
    t[0] = c2[0]-c1[0];
    t[1] = c2[1]-c1[1];
  }
  
  
  /*mean normal*/
  n1[0] = - (c2[1] - c1[1]);
  n1[1] =  (c2[0] - c1[0]);
  dd1 = 1;//sqrt(n1[0]*n1[0]+n1[1]*n1[1]);
  n2[0] = - (c3[1] - c2[1]);
  n2[1] =  (c3[0] - c2[0]);
  dd2 = 1;//sqrt(n2[0]*n2[0]+n2[1]*n2[1]);
  
  dd = 0;
  for(i=0 ; i<2 ; i++) {
    n[i] = 0.5*(n1[i]/dd1+n2[i]/dd2);
    dd += n[i]*n[i];
  }
  dd = sqrt(dd);
  n[0] /= dd;
  n[1]/=dd;
  t[0] = n[1];
  t[1] = -n[0];
  /* t1[0] = c2[0] - c1[0]; */
  /* t1[1] = c2[1] - c1[1]; */
  /* t2[0] = c3[0] - c2[0]; */
  /* t2[1] = c3[1] - c2[1]; */
  /* theta = (t1[0]*t2[0]+t1[1]*t2[1])/(sqrt(t1[0]*t1[0]+t1[1]*t1[1])*sqrt(t2[0]*t2[0]+t2[1]*t2[1])); */
  /* //printf("theta %e %e == %e\n",theta,acos(theta),M_Pi/2.); */
  /* theta = acos(theta); */
  /* //sign */
  /* dd = t1[0]*t2[1]-t1[1]*t2[0]; */
  /* if(dd<0) {/\*puts("eheheheheheheheh 222");*\/theta *= -1.;} */
  /* dd = sqrt(t1[0]*t1[0]+t1[1]*t1[1]) +sqrt(t2[0]*t2[0]+t2[1]*t2[1]); */
  /* //printf("arc length %e\n",dd); */
  /* kappa = theta / (0.5*dd); */
  /* printf("pour %d on trouve kappa %e\n",ip2,kappa); */
  
  /* if(fabs(kappa)<1e-8){ */
  /*   printf("STRAIGHT EDGE\n"); */
  /*   return(1); */
  /* }  */
  
  /* r = 1./kappa; */
  /* if(r>0) { */
  /*   n[0] = - (c2[1] - c1[1]); */
  /*   n[1] =  (c2[0] - c1[0]); */
  /* } else { */
  /*   n[0] =  (c2[1] - c1[1]); */
  /*   n[1] = - (c2[0] - c1[0]); */
  /*   r *= -1; */
  /* } */
  /* dd = sqrt(n[0]*n[0]+n[1]*n[1]); */
  /* n[0] /= dd; */
  /* n[1] /= dd; */
  /* //can be < 0 !! ==> straight edge */
  /* if((r*r - 0.25*dd*dd) < 0) { */
  /*   printf("STRAIGHT EDGE 2\n"); */
  /*   return(1); */
  /* } */
  /* center[0] = 0.5*(c1[0]+c2[0]) + (sqrt(r*r - 0.25*dd*dd))*(n[0]); */
  /* center[1] = 0.5*(c1[1]+c2[1]) + (sqrt(r*r - 0.25*dd*dd))*(n[1]); */
  
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param ped pointer toward an Edge.
 * \param k integer
 * \param ik
 * \param i2
 * \param nv
 * \return 1
 *
 * Compute the tangent.
 *
 */
int MMG2_tangent(MMG5_pMesh mesh,MMG5_pEdge ped,int k,int ik,int i2,int nv) {
  MMG5_pTria pt,ptiel;
  MMG5_pEdge pediel;
  MMG5_pPoint pta,ptb,ptaiel,ptbiel;
  
  int iadj,ip,i,iadr,*adja,iel,voyiel,vi1,vi2,num,ipnv;
  
  pt = &mesh->tria[k];
  iadr = 3*(k-1)+1;
  adja = &mesh->adja[iadr];
  iadj = i2;
  ip = pt->v[ik];
  pta = &mesh->point[ped->a];
  ptb = &mesh->point[ped->b];
  
  if(!adja[iadj]) {
    if(nv) {
      MMG2_computetangent(mesh,pt->v[i2],ped->b,ip,ptb->n);
    }
    else {
      MMG2_computetangent(mesh,pt->v[i2],ped->a,ip,pta->n);
    }
    num = pt->edg[i2];
    assert(num);
    pediel = &mesh->edge[num];
    ptaiel = &mesh->point[pediel->a];
    ptbiel = &mesh->point[pediel->b];
    ipnv = (nv==0)?ped->a:ped->b;
    if(pediel->a==ipnv){
      for(i=0 ; i<2 ; i++) {
        if(nv)
          ptaiel->n[i] = ptb->n[i];
        else
          ptaiel->n[i] = pta->n[i];
      }
    } else {
      assert(pediel->b==ipnv);
      for(i=0 ; i<2 ; i++) {
        if(nv)
          ptbiel->n[i] = ptb->n[i];
        else
          ptbiel->n[i] = pta->n[i];
      }
    }
    return(1);
    
  }
  do {
    iel = adja[iadj]/3;
    ptiel = &mesh->tria[iel];
    voyiel = adja[iadj]%3;
    /*find ptiel*/
    iadr = 3*(iel-1) + 1;
    adja = &mesh->adja[iadr];
    
    vi1 = MMG2_idir[voyiel+1];
    vi2 = MMG2_idir[voyiel+2];
    if(ptiel->v[vi2]==ip) {
      num = ptiel->edg[vi2];
      iadj = vi2;
    } else {
      num = ptiel->edg[vi1];
      iadj = vi1;
    }
    ip = ptiel->v[voyiel];
  } while (!num && iel);
  if(nv) {
    MMG2_computetangent(mesh,pt->v[i2],ped->b,ptiel->v[voyiel],ptb->n);
  }
  else {
    MMG2_computetangent(mesh,pt->v[i2],ped->a,ptiel->v[voyiel],pta->n);
  }
  
  return(1);
}


/*************************************************************************************/
/********************************** End hash_d2.c *************************************/
/*************************************************************************************/

/*************************************************************************************/
/********************************** Begin mmg2d2.c *************************************/
/*************************************************************************************/


/**
 *
 * \return <0 value if triangle outside ; > 0 if triangle inside
 *
 * Seek if triangle is inside, outside or not determinated (thus, the edge must
 * be enforced or the case is non convexe.
 *
 */
int MMG2_findpos(MMG5_pMesh mesh,MMG5_pTria pt,int ip1,int ip2,int ip3,int ip4,int base) {
  MMG5_pEdge     ped;
  int            i,j,nb,ped0,ped1;
  unsigned int   MMG2_iare2[3][2] = { {0,1}, {0,2}, {1,2} };
  
  /*au moins 1 point de la BB ??*/
  nb = 0;
  for(i=0 ; i<3 ; i++) {
    if(pt->v[i]==ip1 || pt->v[i]==ip2 || pt->v[i]==ip3 || pt->v[i]==ip4) {
      nb++;
    }
  }
  if(nb==3 || nb==2) {
    pt->base = -base;
    pt->ref  = 3;
    return(-base);
  } else if(!nb){
    pt->base = base;
    //pt->ref  = base;
    return(base);
  }
  /*1 pt de la BB => on ne sait pas si le tr est dedans ou dehors*/
  /*contient un edge ?*/
  for(i=0 ; i<3 ; i++) {
    ped0 = pt->v[MMG2_iare2[i][0]];
    ped1 = pt->v[MMG2_iare2[i][1]];
    for(j=1 ; j<=mesh->na ; j++) {
      ped = &mesh->edge[j];
      if((ped->a == ped0 && ped->b==ped1) || (ped->b == ped0 && ped->a==ped1)){
        ped->base = -1;
        break;
      }
    }
    if(j<=mesh->na) break;
  }
  if(i<3) { /*le tr est dehors*/
    //printf("on a trouve edge : %d %d\n",ped0,ped1);
    pt->base = -base;
    pt->ref  = 3;
    return(base);
  }
  
  /*le tr est indetermine*/
  pt->base = 0;
  
  return(0);
}

/* base boundary vertices and compute tangents*/
int MMG2_baseBdry(MMG5_pMesh mesh) {
  MMG5_pTria     pt,pt1;
  MMG5_pPoint    ppt;
  MMG5_pEdge     ped;
  int      *adja,adj,iadr,k,i,ip,ned,num,i1,i2;
  HashTable edgeT;
  
  ned = 0;
  if(mesh->na) {
    ned = mesh->na;
    /*edge treatment*/
    edgeT.size  = mesh->namax;
    edgeT.nxtmax = 3*mesh->namax+1;
    edgeT.hnxt  = mesh->namax;
    _MMG5_SAFE_CALLOC(edgeT.item,edgeT.nxtmax,Hedge);
    
    memset(edgeT.item,0,edgeT.nxtmax*sizeof(Hedge));
    for (k=edgeT.size; k<edgeT.nxtmax; k++)
      edgeT.item[k].nxt = k+1;
    for(k=1 ; k<=mesh->na ; k++) {
      ped = &mesh->edge[k];
      if(!ped->a) continue;
      MMG2_hashEdge(&edgeT,k,ped->a,ped->b);
    }
  }
  
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    
    iadr = 3*(k-1) + 1;
    adja = &mesh->adja[iadr];
    
    for (i=0; i<3; i++) {
      adj = adja[i]/3;
      pt1 = &mesh->tria[adj];
      pt->edg[i] = 0;
      if ( !adj  ) {
        num = 0;
        if(ned) {
          num = MMG2_hashEdge(&edgeT,k,pt->v[MMG2_iopp[i][0]],pt->v[MMG2_iopp[i][1]]);
        }
        ip  = pt->v[MMG2_iopp[i][0]];
        mesh->point[ip].tag |= M_BDRY;
        if ( mesh->info.nosurf && ( !( mesh->point[ip].tag & M_REQUIRED) ) ) {
          mesh->point[ip].tag |= M_REQUIRED;
          mesh->point[ip].tag |= M_NOSURF;
        }
        
        ip  = pt->v[MMG2_iopp[i][1]];
        mesh->point[ip].tag |= M_BDRY;
        if ( mesh->info.nosurf && ( !( mesh->point[ip].tag & M_REQUIRED) ) ) {
          mesh->point[ip].tag |= M_REQUIRED;
          mesh->point[ip].tag |= M_NOSURF;
        }
        
        if(num) {
          pt->edg[i] = num;
        } else {
          num = _MMG5_newEdge(mesh);
          if ( !num ) {
            _MMG5_EDGE_REALLOC(mesh,num,mesh->gap,
                               printf("  ## Error: unable to allocate a new edge.\n");
                               _MMG5_INCREASE_MEM_MESSAGE();
                               printf("  Exit program.\n");
                               exit(EXIT_FAILURE));
          }
          pt->edg[i] = num;
          ped = &mesh->edge[num];
          ped->a = pt->v[MMG2_iopp[i][0]];
          ped->b = pt->v[MMG2_iopp[i][1]];
        }
      } else if(pt->ref != pt1->ref) {
        num = 0;
        if(ned) {
          num = MMG2_hashEdge(&edgeT,k,pt->v[MMG2_iopp[i][0]],pt->v[MMG2_iopp[i][1]]);
        }
        ip  = pt->v[MMG2_iopp[i][0]];
        mesh->point[ip].tag |= M_SD;
        if ( mesh->info.nosurf && ( !( mesh->point[ip].tag & M_REQUIRED) ) ) {
          mesh->point[ip].tag |= M_REQUIRED;
          mesh->point[ip].tag |= M_NOSURF;
        }
        
        ip  = pt->v[MMG2_iopp[i][1]];
        mesh->point[ip].tag |= M_SD;
        if ( mesh->info.nosurf  && ( !( mesh->point[ip].tag & M_REQUIRED) )  ) {
          mesh->point[ip].tag |= M_REQUIRED;
          mesh->point[ip].tag |= M_NOSURF;
        }
        
        if(num) {
          pt->edg[i] = num;
        } else {
          num = _MMG5_newEdge(mesh);
          if ( !num ) {
            _MMG5_EDGE_REALLOC(mesh,num,mesh->gap,
                               printf("  ## Error: unable to allocate a new edge.\n");
                               _MMG5_INCREASE_MEM_MESSAGE();
                               printf("  Exit program.\n");
                               exit(EXIT_FAILURE));
          }
          pt->edg[i] = num;
          ped = &mesh->edge[num];
          ped->a = pt->v[MMG2_iopp[i][0]];
          ped->b = pt->v[MMG2_iopp[i][1]];
        }
      }
    }
  }
  
  if(ned && edgeT.item)  _MMG5_SAFE_FREE(edgeT.item);
  
  /*compute tangents*/
  mesh->base++;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;
    iadr = 3*(k-1) + 1;
    adja = &mesh->adja[iadr];
    for (i=0; i<3; i++) {
      num = pt->edg[i];
      if(!num) continue;
      ped = &mesh->edge[num];
      
      i1 = MMG2_idir[i+1];
      i2 = MMG2_idir[i+2];
      ppt = &mesh->point[ped->a];
      /*tangent on ped->a*/
      if(ppt->flag != mesh->base) {
        if(pt->v[i1]==ped->a) {
          assert(pt->v[i2]==ped->b);
          MMG2_tangent(mesh,ped,k,i,i2,0);
          ppt->flag = mesh->base;
        } else {
          assert(pt->v[i2]==ped->a);
          MMG2_tangent(mesh,ped,k,i,i1,0);
          ppt->flag = mesh->base;
        }
        
      } else {
        //printf("on a deja calcule pour %d %e %e\n",ped->a,ped->t0[0],ped->t0[1]);
      }
      
      
      
      /* ppt = &mesh->point[ped->b]; */
      /* /\*tangent on ped->b*\/ */
      /* if(ppt->flag != mesh->flag) { */
      /*  if(pt->v[i1]==ped->b) { */
      /*    MMG2_tangent(mesh,ped,k,i,i2,1); */
      /*  } else { */
      /*    MMG2_tangent(mesh,ped,k,i,i1,1);     */
      /*  } */
      /*  ppt->flag = mesh->flag; */
      
      /* } else { */
      /*  //printf("on a deja calcule pour %d %e %e\n",ped->b,ped->t1[0],ped->t1[1]); */
      /* } */
      
    }
  }
  return(1);
}


/*************************************************************************************/
/********************************** End mmg2d2.c *************************************/
/*************************************************************************************/

/*************************************************************************************/
/********************************** Begin mmg2d0.c *************************************/
/*************************************************************************************/

int MMG2_mmg2d0(MMG5_pMesh mesh,MMG5_pSol sol) {
  int       ns,nm,nsiter,nmiter,nmbar,it,maxtou;
  double    declic;
  
  /*optim*/
  ns     = 0;
  nm     = 0;
  nmiter = 0;
  nsiter = 0;
  it     = 0;
  maxtou = 100;
  do {
    /*edge flip*/
    if(!mesh->info.noswap) {
      declic = 1.5 / ALPHA;
      nsiter = MMG2_cendel(mesh,sol,declic,-1);
      if ( nsiter && mesh->info.imprim < 0)
        fprintf(stdout,"     %7d SWAPPED\n",nsiter);
      
      ns+=nsiter;
    }
    /*point relocation*/
    if(!mesh->info.nomove) {
      declic = 1.5 / ALPHA;
      nmiter = MMG2_optlen(mesh,sol,declic,-1);
      if(sol->size==1) nmbar =  optlen_iso_bar(mesh,sol,declic,-1);
      else nmbar=0;
      nm += nmiter+nmbar;
      if ( mesh->info.imprim < 0)
        fprintf(stdout,"     %7d + %7d MOVED \n",nmiter,nmbar);
    }
    
    
  } while((nmiter+nsiter > 0) && (++it <= maxtou));
  if ( mesh->info.imprim )
    fprintf(stdout,"     %7d SWAPPED %7d MOVED\n",ns,nm);
  
  return(1);
}

/*************************************************************************************/
/********************************** End mmg2d0.c *************************************/
/*************************************************************************************/

/*************************************************************************************/
/********************************** Begin mmg2d1.c *************************************/
/*************************************************************************************/

#define BUCKSIZ    64
#define M_LONG     1.4//1.85//1.4//1.421
#define M_SHORT    0.65//0.8//0.65//0.707

int MMG2_invmat(double *m,double *minv) {
  double        det;
  
  if(fabs(m[1]) < EPSD) { /*mat diago*/
    minv[0] = 1./m[0];
    minv[1] = 0;
    minv[2] = 1./m[2];
  } else {
    det = m[0]*m[2] - m[1]*m[1];
    det = 1. / det;
    minv[0] = det * m[2];
    minv[1] = - det * m[1];
    minv[2] = det * m[0];
  }
  return(1);
}

int interp_ani(double *ma,double *mb,double *mp,double t) {
  double  dma[3],dmb[3],mai[3],mbi[3],mi[3];
  int   i;
  
  for (i=0; i<3; i++) {
    dma[i] = ma[i];
    dmb[i] = mb[i];
  }
  
  if ( !MMG2_invmat(dma,mai) || !MMG2_invmat(dmb,mbi) ) {
    fprintf(stderr,"  ## Error: unable to interpole the metric.\n");
    return(0);
  }
  
  for (i=0; i<3; i++)
    mi[i] = (1.0-t)*mai[i] + t*mbi[i];
  
  if ( !MMG2_invmat(mi,mai) ) {
    fprintf(stderr,"  ## Error: invalid metric.\n");
    return(0);
  }
  
  for (i=0; i<3; i++)  mp[i] = mai[i];
  return(1);
}

int interp_iso(double *ma,double *mb,double *mp,double t) {
  
  *mp = (1.0-t)*(*ma) + t*(*mb);
  return(1);
}

static int cassar(MMG5_pMesh mesh,MMG5_pSol sol,int ia,int ib,double t) {
  MMG5_pPoint   p1,p2;
  //Displ      pd;
  double   c[2],t1,*ma,*mb,*mp;
  int      ip,iadr,memlack;
  
  memlack = 0;
  
  p1 = &mesh->point[ia];
  p2 = &mesh->point[ib];
  t1 = 1.0 - t;
  
  c[0] = t1*p1->c[0] +  t*p2->c[0];
  c[1] = t1*p1->c[1] +  t*p2->c[1];
  ip   = _MMG2D_newPt(mesh,c,0);
  if ( !ip ) {
    /* reallocation of point table */
    
    _MMG2D_POINT_REALLOC(mesh,sol,ip,mesh->gap,
                         printf("  ## Error: unable to allocate a new point\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         memlack=1;
                         return(-1)
                         ,c,0);
    p1  = &mesh->point[ia];
    p2  = &mesh->point[ib];
  }
  
  
  /*interpol metric*/
  iadr = (ia-1)*sol->size + 1;
  ma  = &sol->m[iadr];
  
  iadr = (ib-1)*sol->size + 1;
  mb  = &sol->m[iadr];
  
  iadr = (ip-1)*sol->size + 1;
  mp  = &sol->m[iadr];
  
  if ( sol->size==1 ) {
    if(!interp_iso(ma,mb,mp,t1) ) return(-1);
  }
  else {
    if(!interp_ani(ma,mb,mp,t1) ) return(-1);
  }
  
  /*interpol dep si option 9*/
  if( mesh->info.lag >=0) {
    printf(" ## Error: option not available:"
           " comment because of merge needs mmg2d1 option 9\n");
    exit(EXIT_FAILURE);
    // pd = mesh->disp;
    //pd.mv[2*(ip-1) + 1 + 0] = t1*pd.mv[2*(ia-1) + 1 + 0] + t*pd.mv[2*(ib-1) + 1 + 0];
    //pd.mv[2*(ip-1) + 1 + 1] = t1*pd.mv[2*(ia-1) + 1 + 1] + t*pd.mv[2*(ib-1) + 1 + 1];
  }
  if ( memlack )  return(-1);
  return(ip);
}

/*compute new vertex on bdry*/
static int cassarbdry(MMG5_pMesh mesh,MMG5_pSol sol,int ied,int ia,int ib,double t,double *tang) {
  MMG5_pPoint   p0,p1,ppt;
  MMG5_pEdge    ped;
  // Displ      pd;
  double   c[2],pc1[2],pc2[2],t0[2],t1[2],t_1,*ma,*mb,*mp;//,dx,dy;
  double   l;
  int      ip,iadr,i,inv,memlack;
  inv = 0;
  memlack = 0;
  p0 = &mesh->point[ia];
  p1 = &mesh->point[ib];
  
  l = (p0->c[0] - p1->c[0])*(p0->c[0] - p1->c[0]) +
  (p0->c[1] - p1->c[1])*(p0->c[1] - p1->c[1]);
  l = sqrt(l);
  ped = &mesh->edge[ied];
  if(ia != ped->a) {
    assert(ib == ped->a);
    for(i=0 ; i<2 ; i++) {
      t0[i] = l*p0->n[i];
      t1[i] = l*p1->n[i];
    }
    /*if corner, recompute the tangent*/
    if(p0->tag & M_CORNER) {
      for(i=0 ; i<2 ; i++)
        t0[i] = p0->c[i] - p1->c[i];
    }
    if(p1->tag & M_CORNER) {
      for(i=0 ; i<2 ; i++)
        t1[i] = p0->c[i] - p1->c[i];
    }
  } else {
    assert(ib == ped->b);
    for(i=0 ; i<2 ; i++) {
      t0[i] = l*p0->n[i];
      t1[i] = l*p1->n[i];
    }
    /*if corner, recompute the tangent*/
    if(p0->tag & M_CORNER) {
      for(i=0 ; i<2 ; i++)
        t0[i] = p1->c[i] - p0->c[i];
    }
    if(p1->tag & M_CORNER) {
      for(i=0 ; i<2 ; i++)
        t1[i] = p1->c[i] - p0->c[i];
    }
    
  }
  
  
  /*check if t0 has the same sens of vect(P0P1)*/
  if(t0[0]/(p1->c[0]-p0->c[0]) < 0 || t0[1]/(p1->c[1]-p0->c[1])<0) {
    //printf("t0/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
    for(i=0 ; i<2 ; i++) {
      t0[i] *= -1;
    }
    inv = 1;
  }
  /*check if t1 has the opposite sens of vect(P0P1)*/
  if(t1[0]/(p1->c[0]-p0->c[0]) > 0 || t1[1]/(p1->c[1]-p0->c[1])>0) {
    //printf("t1/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
    for(i=0 ; i<2 ; i++) {
      t1[i] *= -1;
    }
  }
  /*control points*/
  for(i=0 ; i<2 ; i++) {
    pc1[i] = (t0[i]+3*p0->c[i])/3.;
    pc2[i] = (t1[i]+3*p1->c[i])/3.;
  }
  
  /*coor new point*/
  t_1 = 1.0 - t;
  for(i=0 ; i<2 ; i++) {
    c[i] = t_1*t_1*t_1*p0->c[i] + 3*t*t_1*t_1*pc1[i] + 3*t*t*t_1*pc2[i] +  t*t*t*p1->c[i];
  }
  // printf("c %e %e -- mid %e %e\n",c[0],c[1],0.5*(p0->c[0]+p1->c[0]),0.5*(p0->c[1]+p1->c[1]));
  ip   = _MMG2D_newPt(mesh,c,0);
  if ( !ip ) {
    /* reallocation of point table */
    
    _MMG2D_POINT_REALLOC(mesh,sol,ip,mesh->gap,
                         printf("  ## Error: unable to allocate a new bdry point\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         memlack=1;
                         return(-1);
                         ,c,0);
    p0 = &mesh->point[ia];
    p1 = &mesh->point[ib];
  }
  
  /*tangent new point*/
  ppt = &mesh->point[ip];
  for(i=0 ; i<2 ; i++) {
    // tang[i] = -(-3*t_1*t_1*p0->c[i] +(3*t_1*(1-3*t))*pc1[i]
    //    + (3*t*(2-3*t))*pc2[i] +  3*t*t*p1->c[i]);
    // printf("tang %e %e diff %e\n",tang[i],(3./8.)*(p1->c[i]+pc2[i]-pc1[i]-p0->c[i]),
    //     fabs(tang[i]-(3./8.)*(p1->c[i]+pc2[i]-pc1[i]-p0->c[i])));
    tang[i] = (3./8.)*(p1->c[i]+pc2[i]-pc1[i]-p0->c[i]);
  }
  l = tang[0]*tang[0] + tang[1]*tang[1];
  l = 1./sqrt(l);
  tang[0] *= l;
  tang[1] *= l;
  /*check if tang has the same sens than P0P*/
  if(inv) {
    if ((tang[0]/(ppt->c[0]-p0->c[0]) > 0 || tang[1]/(ppt->c[1]-p0->c[1])>0)) {
      //printf("tang/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
      for(i=0 ; i<2 ; i++) {
        tang[i] *= -1;
      }
    }
  } else if((tang[0]/(ppt->c[0]-p0->c[0]) < 0 || tang[1]/(ppt->c[1]-p0->c[1])<0)) {
    //printf("tang/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
    for(i=0 ; i<2 ; i++) {
      tang[i] *= -1;
    }
  }
  
  /*   /\*change the two other tangents*\/ */
  /*   if(ia != ped->a) { */
  /*     assert(ib == ped->a); */
  /*     for(i=0 ; i<2 ; i++) { */
  /*       /\*t0*\/ */
  /*       p0->n[i] = -(3./2.)*(pc1[i]-p0->c[i]); */
  /*       /\*t1*\/ */
  /*       p1->n[i] = -(3./2.)*(p1->c[i]-pc2[i]); */
  /*     } */
  /*     /\*check tangent orientation*\/ */
  /* #warning remove ? check orientation tangent */
  /*     dx = t1[0]/p0->n[0]; */
  /*     dy = t1[1]/p0->n[1]; */
  /*     if(dy < 0 || dx <0) { */
  /*       //printf("1) pbs de colinearite %e\n",fabs(dx-dy)); */
  /*     } */
  /*     dx = t0[0]/p1->n[0]; */
  /*     dy = t0[1]/p1->n[1]; */
  /*     if(dx < 0 || dy <0) { */
  /*       //printf("3) pbs de colinearite %e %e %e\n",fabs(dx-dy),dx,dy); */
  /*     } */
  /*   } else { */
  /*     assert(ib == ped->b); */
  /*     for(i=0 ; i<2 ; i++) { */
  /*       /\*t0*\/ */
  /*       p0->n[i] = -(3./2.)*(pc1[i]-p0->c[i]); */
  /*       /\*t1*\/ */
  /*       p1->n[i] = -(3./2.)*(p1->c[i]-pc2[i]); */
  /*     } */
  /*     /\*check tangent orientation*\/ */
  /* #warning remove ? check orientation tangent */
  /*     dx = t1[0]/p1->n[0]; */
  /*     dy = t1[1]/p1->n[1]; */
  /*     if(dx < 0 || dy <0) { */
  /*       //printf("2) pbs de colinearite %e %e %e\n",fabs(dx-dy),dx,dy); */
  /*     } */
  /*     dx = t0[0]/p0->n[0]; */
  /*     dy = t0[1]/p0->n[1]; */
  /*     if(dx < 0 || dy <0) { */
  /*       //printf("4) pbs de colinearite pts %d %e %e %e\n",ia,fabs(dx-dy),dx,dy); */
  /*       for(i=0 ; i<2 ; i++) { */
  /*         p0->n[i] *= -1; */
  /*       } */
  
  /*     } */
  /*   } */
  
  /*interpol metric*/
  iadr = (ia-1)*sol->size + 1;
  ma  = &sol->m[iadr];
  
  iadr = (ib-1)*sol->size + 1;
  mb  = &sol->m[iadr];
  
  iadr = (ip-1)*sol->size + 1;
  mp  = &sol->m[iadr];
  
  if ( sol->size==1 ) {
    if(!interp_iso(ma,mb,mp,t_1) ) return(-1);
  }
  else {
    if(!interp_ani(ma,mb,mp,t_1) ) return(-1);
  }
  
  /*interpol dep si option 9*/
  if( mesh->info.lag >= 0) {
    //#warning option 9
    printf(" ## Error: option not available:"
           " comment because of merge needs mmg2d1 option 9\n");
    exit(EXIT_FAILURE);
    /* pd = mesh->disp; */
    /* pd.mv[2*(ip-1) + 1 + 0] = t_1*pd.mv[2*(ia-1) + 1 + 0] + t*pd.mv[2*(ib-1) + 1 + 0]; */
    /* pd.mv[2*(ip-1) + 1 + 1] = t_1*pd.mv[2*(ia-1) + 1 + 1] + t*pd.mv[2*(ib-1) + 1 + 1]; */
  }
  if ( memlack )  return(-1);
  return(ip);
}

/**
 * \param mesh poitner toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param bucket pointer toward the bucket structure.
 * \param declic quality threshold.
 * \param alert if 1, we are unable to create a new vertex.
 * \param ni number of inserted points.
 * \param nc nuber of collapsed points.
 * \return 0 if fail, 1 otherwise.
 *
 * Analyse the edges, split the longer and collapse the shorter one.
 *
 */
static int analar(MMG5_pMesh mesh,MMG5_pSol sol,pBucket bucket,
                  double declic,int *alert, int *ni, int *nc) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppa,ppb;
  double  *ca,*cb,*ma,*mb,tail,t,tang[2];
  int     *adja,voi[3],k,iadr,adj,/*base,*/nbp,npp,ip;
  int     nt,ier;
  int     i,i1,i2;
  int     ins,i0,ii0;
  
  //  base  = ++mesh->base;
  (*ni)  = 0;
  (*nc)  = 0;
  nt  = mesh->nt;
  npp = 0.0;
  
  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;
    //else if ( /*pt->flag = base-1 || pt->qual < declic*/ )  continue;
    
    /* base internal edges */
    iadr  = 3*(k-1) + 1;
    adja  = &mesh->adja[iadr];
    voi[0] = adja[0];
    voi[1] = adja[1];
    voi[2] = adja[2];
    i0 = 0;
    if (!voi[1]) i0 = 1;
    if (!voi[2]) i0 = 2;
    for (ii0=i0; ii0<i0+3; ii0++) {
      i = ii0%3;
      adj = voi[i] / 3;
      
      i1   = pt->v[MMG2_idir[i+1]];
      i2   = pt->v[MMG2_idir[i+2]];
      
      ppa  = &mesh->point[i1];
      ppb  = &mesh->point[i2];
      //#warning bad test for edge required
      if((ppa->tag & M_REQUIRED) && (ppb->tag & M_REQUIRED)) {
        //printf("edge required %d %d\n",i1,i2);
        continue;
      }
      ca   = &ppa->c[0];
      cb   = &ppb->c[0];
      iadr = (i1-1)*sol->size + 1;
      ma   = &sol->m[iadr];
      iadr = (i2-1)*sol->size + 1;
      mb   = &sol->m[iadr];
      tail = MMG2_length(ca,cb,ma,mb);
      
      if ( tail > M_LONG && *alert <= 1 ) {
        npp++;
        nbp = tail + 0.5;
        if ( nbp*(nbp+1) < 0.99*tail*tail )  nbp++;
        t = 1.0 / (float)nbp;
        if ( nbp < 3 || nbp > 15 )  t = 0.5;
        if ( !adj || pt->ref != mesh->tria[adj].ref )  {
          /*add bdry*/
          if(!pt->edg[i])  {
            /* if(mesh->info.ddebug) { */
            /*   printf("tr %d : %d %d %d mais %d\n",k,pt->edg[0],pt->edg[1],pt->edg[2],i); */
            /*   printf("%d %d %d\n",pt->v[0],pt->v[1],pt->v[2]); */
            /* } */
            assert(mesh->tria[adj].ref!=pt->ref);
            assert((mesh->tria[adj]).edg[voi[i]%3]);
            //#warning find why we have to do that
            pt->edg[i] = (mesh->tria[adj]).edg[voi[i]%3];
          }
          assert(pt->edg[i]);
          ip = cassarbdry(mesh,sol,pt->edg[i],i1,i2,0.5,tang);
        } else {
          ip = cassar(mesh,sol,i1,i2,t);
        }
        if(ip < 0) {
          if(mesh->info.imprim > 6)
            printf("  ## Warning: impossible to create new vertex\n");
          //return(0);
          //printf("ahhhhhhhhhhhhhhhh\n");
          *alert = 2;
        } else {
          if ( !adj || pt->ref != mesh->tria[adj].ref )  {
            /*boundary edge*/
            if(!adj) {
              ins = MMG2_splitbdry(mesh,sol,ip,k,i,tang);
              if(!ins) {
                _MMG2D_delPt(mesh,ip);
                continue;
              }
              mesh->point[ip].tag |= M_BDRY;
              (*ni) += 1;
              break;
            } else {
              mesh->point[ip].tag |= M_SD;
              ins = MMG2_split(mesh,sol,ip,k,voi[i]);
              if(!ins) {
                _MMG2D_delPt(mesh,ip);
                continue;
              }
              (*ni) += 1;
              break;
            }
            continue;
          } else {
            ins = MMG2_split(mesh,sol,ip,k,voi[i]);
            if(!ins) {
              _MMG2D_delPt(mesh,ip);
              continue;
            }
            (*ni) += 1;
            break;
          }
        }
      }
      
      else if ( tail < M_SHORT ) {
        if ( !adj || pt->ref != mesh->tria[adj].ref )  {
          if(!adj) {
            
            ier = MMG2_colpoibdry(mesh,sol,k,i,MMG2_iare[i][0],
                                  MMG2_iare[i][1],2.75);
            if ( ier ==-1 ) return(0);
            
            else if ( !ier ){
              ier = MMG2_colpoibdry(mesh,sol,k,i,MMG2_iare[i][1],
                                    MMG2_iare[i][0],2.75);
              if ( ier==-1 ) return(0);
              else if ( !ier ){
                continue;
              } else {
                (*nc)++;
                _MMG2D_delPt(mesh,i1);
                break;
              }
            }
            (*nc)++;
            _MMG2D_delPt(mesh,i2);
            break;
          } else {
            if(!MMG2_colpoi(mesh,sol,k,i,MMG2_iare[i][0],MMG2_iare[i][1],2.75)) {
              if(!MMG2_colpoi(mesh,sol,k,i,MMG2_iare[i][1],MMG2_iare[i][0],2.75)) {
                continue;
              } else {
                (*nc)++;
                _MMG2D_delPt(mesh,i1);
                break;
              }
            }
            _MMG2D_delPt(mesh,i2);
            (*nc)++;
            break;
          }
        } else {
          if(!MMG2_colpoi(mesh,sol,k,i,MMG2_iare[i][0],MMG2_iare[i][1],2.75)) {
            if(!MMG2_colpoi(mesh,sol,k,i,MMG2_iare[i][1],MMG2_iare[i][0],2.75)) {
              continue;
            } else {
              (*nc)++;
              _MMG2D_delPt(mesh,i1);
              break;
              
            }
          }
          (*nc)++;
          _MMG2D_delPt(mesh,i2);
          break;
        }
      }
    }
  }
  if ( mesh->info.imprim > 5 ) {
    fprintf(stdout,"    %8d INSERTED %8d COLLAPSED\n",*ni,*nc);
  }
  return(1);
}
static int analargeom(MMG5_pMesh mesh,MMG5_pSol sol,int *alert) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppa,ppb;
  double  tail,tang[2];
  int     *adja,voi[3],k,iadr,adj,/* base ,*/npp,ip;
  int     nt;
  int     i,i1,i2,ni,maxtou,it;
  int     ins,i0,ii0,nitot;
  maxtou = 30;
  it     = 0;
  nitot = 0;
  do {
    ni  = 0;
    // base  = ++mesh->base;
    nt    = mesh->nt;
    npp   = 0;
    for (k=1; k<=nt; k++) {
      pt = &mesh->tria[k];
      if ( !M_EOK(pt) )  continue;
      //else if ( /*pt->flag = base-1 || pt->qual < declic*/ )  continue;
      
      /* base internal edges */
      iadr  = 3*(k-1) + 1;
      adja  = &mesh->adja[iadr];
      voi[0] = adja[0];
      voi[1] = adja[1];
      voi[2] = adja[2];
      i0 = 0;
      if (!voi[1]) i0 = 1;
      if (!voi[2]) i0 = 2;
      for (ii0=i0; ii0<i0+3; ii0++) {
        i = ii0%3;
        adj = voi[i] / 3;
        
        i1   = pt->v[MMG2_idir[i+1]];
        i2   = pt->v[MMG2_idir[i+2]];
        
        ppa  = &mesh->point[i1];
        ppb  = &mesh->point[i2];
        //#warning bad test for edge required
        if((ppa->tag & M_REQUIRED) && (ppb->tag & M_REQUIRED)) {
          //printf("edge required %d %d\n",i1,i2);
          continue;
        }
        if ( adj )  continue;
        
        tail = MMG2_chkedg(mesh,ppa,ppb);
        if(!tail) continue;
        if ( *alert <= 1 ) {
          npp++;
          assert(pt->edg[i]);
          ip = cassarbdry(mesh,sol,pt->edg[i],i1,i2,0.5,tang);
          
          if(ip < 0) {
            if(mesh->info.imprim > 6) printf("  ## Warning: impossible to create new vertex\n");
            *alert = 2;
          }
          ins = MMG2_splitbdry(mesh,sol,ip,k,i,tang);
          if(!ins) {
            _MMG2D_delPt(mesh,ip);
            continue;
          }
          mesh->point[ip].tag |= M_BDRY;
          ni += 1;
          break;
        }
      }
    }
    if ( mesh->info.imprim > 5 ) {
      fprintf(stdout,"    %8d INSERTED \n",ni);
    }
    nitot +=ni;
  } while(ni > 0 && it++ < maxtou);
  return(nitot);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return 1 if success, 0 if strongly fail.
 *
 * Mesh adaptation.
 *
 **/
int MMG2_mmg2d1(MMG5_pMesh mesh,MMG5_pSol sol) {
  pBucket  bucket;
  double   declic;
  int      ns,base,alert,it,maxtou;
  int      nadd,ndel,nswp,ngeom,ni,nc;
  
  nadd = ndel = nswp = 0;
  if ( mesh->info.imprim < 0 ) {
    MMG2_outqua(mesh,sol);
    MMG2_prilen(mesh,sol);
  }
  
  /* 1. Delaunization */
  if ( mesh->info.imprim < -4 )
    fprintf(stdout,"  -- DELAUNIZATION\n");
  
  /* edge flip */
  if(!mesh->info.noswap) {
    declic = 1.1 / ALPHA;
    base   = mesh->base;
    ns = MMG2_cendel(mesh,sol,declic,base);
    nswp += ns;
    if ( mesh->info.imprim > 5 )
      fprintf(stdout,"  -- %8d SWAPPED\n",ns);
  }
  alert  = 0;
  
  /* 1. Geometric mesh */
  if ( mesh->info.imprim > 3 )
    fprintf(stdout,"  -- GEOMETRIC MESH\n");
  ngeom = analargeom(mesh,sol,&alert);
  if ( mesh->info.imprim && (abs(mesh->info.imprim) < 6) )
    fprintf(stdout,"     %8d splitted\n",ngeom);
  
  /* 2. field points */
  bucket = MMG2_newBucket(mesh,M_MAX(mesh->info.bucket,BUCKSIZ));
  assert(bucket);
  declic = 1.5 / ALPHA;
  maxtou = 30;
  it     = 0;
  do {
    ni = 0;
    nc = 0;
    if ( !analar(mesh,sol,bucket,declic,&alert,&ni,&nc) ) {
      return(0);
    }
    nadd += ni;
    ndel += nc;
    if(!mesh->info.noswap) {
      ns = MMG2_cendel(mesh,sol,declic,mesh->base);
      nswp += ns;
      if ( mesh->info.imprim > 5 )
        fprintf(stdout,"  -- %8d SWAPPED\n",ns);
    }
  }
  while ( ++it < maxtou && (ni+nc > 0));
  //while ( ++it < maxtou && (ni+nc > 0));//> 0.05*mesh->np));
  
  if ( mesh->info.imprim && (abs(mesh->info.imprim) < 6) && ( nadd || ndel ) ) {
    fprintf(stdout,"     %8d splitted, %8d collapsed,"
            " %8d swapped.\n",nadd,ndel,nswp);
  }
  
  if ( mesh->info.imprim < 0 ) {
    MMG2_outqua(mesh,sol);
    MMG2_prilen(mesh,sol);
  }
  
  /* free memory */
  _MMG5_SAFE_FREE(bucket->head);
  _MMG5_SAFE_FREE(bucket->link);
  _MMG5_SAFE_FREE(bucket);
  
  return(1);
}

/*************************************************************************************/
/********************************** End mmg2d1.c *************************************/
/*************************************************************************************/

/** UNUSED FUNCTION */
/*insertion of the list of points inside the mesh*/
/*return 0 if pbs occur*/
int MMG2_insertpoint(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria   pt,pt1;
  MMG5_pPoint  ppt;
  double  declic;
  int     k,nsiter,lel,mel,nel,ia,ib,ic,aext0,aext1,aext2;
  int     iadr,*adja,nflat,ie1,ie2,ie3,voy,text,atext1,atext2;
  double  EPSDD = 1e+6;
  
  for(k=1 ; k<=mesh->np - 4 ; k++) {
    if(1){
      if(mesh->info.ddebug) printf("------------------------------------try swap to increase mesh quality\n");
      declic = 1.1 / ALPHA;
      nsiter = 1;
      while (nsiter) {
        nsiter = MMG2_cendel(mesh,sol,declic,-1);
        if ( nsiter && mesh->info.imprim < 0)
          fprintf(stdout,"     %7d SWAPPED\n",nsiter);
      }
    }
    ppt = &mesh->point[k];
    if(ppt->tmp==1) continue;
    /*recherche du triangle contenant le point : lel*/
    lel = MMG2_findTria(mesh,k);
    assert(lel);
    
    nflat = 0; // to avoid bad triangle  1, 2, 4, 3, 5, 7
    pt  = &mesh->tria[lel];
    iadr = 3*(lel-1) + 1;
    adja = &mesh->adja[iadr];
    ia  = pt->v[0];
    ib  = pt->v[1];
    ic  = pt->v[2];
    aext0 = adja[0];
    aext1 = adja[1];
    aext2 = adja[2];
    
    /*creation de trois triangles*/
    mel = _MMG2D_newElt(mesh);
    if ( !mel ) {
      _MMG5_TRIA_REALLOC(mesh,mel,mesh->gap,
                         printf("  ## Error: unable to allocate a new element.\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         printf("  Exit program.\n");
                         exit(EXIT_FAILURE));
      pt  = &mesh->tria[lel];
    }
    nel = _MMG2D_newElt(mesh);
    if ( !nel ) {
      _MMG5_TRIA_REALLOC(mesh,mel,mesh->gap,
                         printf("  ## Error: unable to allocate a new element.\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         printf("  Exit program.\n");
                         exit(EXIT_FAILURE));
      pt  = &mesh->tria[lel];
    }
    pt->v[2] = k;
    pt->qual = MMG2_caltri_in(mesh,sol,pt);
    if(pt->qual > EPSDD) {
      nflat = 1;
    } else {
      adja[0] = 3*nel + 1;
      adja[1] = 3*mel + 0;
      if(aext2)
        (&mesh->adja[3*(aext2/3-1) + 1])[aext2%3] = 3*lel + 2;
      //printf("adj of %d : %d %d %d -- %d(%d) : %d\n",lel,nel,mel,aext2/3,aext2/3,aext2%2,lel);
    }
    pt  = &mesh->tria[nel];
    pt->base = mesh->base;
    iadr = 3*(nel-1) + 1;
    adja = &mesh->adja[iadr];
    pt->v[0] = ib;
    pt->v[1] = ic;
    pt->v[2] = k;
    pt->qual = MMG2_caltri_in(mesh,sol,pt);
    if(pt->qual > EPSDD) {
      nflat += 2;
    } else {
      adja[0] = 3*mel + 1;
      adja[1] = 3*lel + 0;
      adja[2] = aext0;
      if(aext0)
        (&mesh->adja[3*(aext0/3-1) + 1])[aext0%3] = 3*nel + 2;
      //printf("adj of %d : %d %d %d -- %d(%d) : %d\n",nel,mel,lel,aext0/3,aext0/3,aext0%3,nel);
    }
    pt  = &mesh->tria[mel];
    pt->base = mesh->base;
    iadr = 3*(mel-1) + 1;
    adja = &mesh->adja[iadr];
    pt->v[0] = ic;
    pt->v[1] = ia;
    pt->v[2] = k;
    pt->qual = MMG2_caltri_in(mesh,sol,pt);
    if(pt->qual > EPSDD) {
      nflat += 4;
      /*on coupe le tr aext1 en 2*/
      text = aext1/3;
      voy  = aext1%3;
      pt1  = &mesh->tria[text];
      ie1  = pt1->v[voy];
      ie2  = pt1->v[MMG2_inxt[voy]];
      ie3  = pt1->v[MMG2_inxt[voy+1]];
      pt1->v[0] = ie1;
      pt1->v[1] = ie2;
      pt1->v[2] = k;
      pt1->qual = MMG2_caltri_in(mesh,sol,pt1);
      pt->v[0] = ie3;
      pt->v[1] = ie1;
      pt->v[2] = k;
      pt->qual = MMG2_caltri_in(mesh,sol,pt);
      /*maj des adj*/
      /*text*/
      iadr = 3*(text-1) + 1;
      adja = &mesh->adja[iadr];
      atext1 = adja[MMG2_inxt[voy]];
      atext2 = adja[MMG2_inxt[voy+1]];
      adja[0] = 3*lel + 1;
      adja[1] = 3*mel + 0;
      adja[2] = atext2;
      if(atext2)
        (&mesh->adja[3*(atext2/3-1) + 1])[atext2%3] = 3*text + 2;
      //if(ddebug) printf("adj of %d : %d %d %d \n",text,adja[0]/3,adja[1]/3,adja[2]/3);
      /*mel*/
      iadr = 3*(mel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[0] = 3*text + 1;
      adja[1] = 3*nel + 0;
      adja[2] = atext1;
      if(atext1)
        (&mesh->adja[3*(atext1/3-1) + 1])[atext1%3] = 3*mel + 2;
      if(nflat==6 /*|| (*adj)==7*/) {//CECILE 14/04/14:*adj ne veut rien dire!!!
        adja[1] = 3*(aext0/3);
        //pas de suite sinon on perd de l'info (&mesh->adja[3*(aext0/3-1) + 1])[0] = 3*mel + 1;
        
      } else {
        /*nel*/
        iadr = 3*(nel-1) + 1;
        adja = &mesh->adja[iadr];
        adja[0] = 3*mel + 1;
      }
      
      /*lel*/
      iadr = 3*(lel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[1] = 3*text + 0;
      
    } else {
      adja[0] = 3*lel + 1;
      adja[1] = 3*nel + 0;
      adja[2] = aext1;
      if(aext1)
        (&mesh->adja[3*(aext1/3-1) + 1])[aext1%3] = 3*mel + 2;
      //printf("adj of %d : %d %d %d -- %d(%d) : %d\n",mel,lel,nel,aext1/3,aext1/3,aext1%3,mel);
      
    }
    if(nflat==1 || nflat==3 || nflat==7 || nflat==5){
      /*on coupe le tr aext2 en 2*/
      text = aext2/3;
      voy  = aext2%3;
      pt1  = &mesh->tria[text];
      ie1  = pt1->v[voy];
      ie2  = pt1->v[MMG2_inxt[voy]];
      ie3  = pt1->v[MMG2_inxt[voy+1]];
      pt1->v[0] = ie1;
      pt1->v[1] = ie2;
      pt1->v[2] = k;
      pt1->qual = MMG2_caltri_in(mesh,sol,pt1);
      pt = &mesh->tria[lel];
      pt->v[0] = ie3;
      pt->v[1] = ie1;
      pt->v[2] = k;
      pt->qual = MMG2_caltri_in(mesh,sol,pt);
      /*maj des adj*/
      /*text*/
      iadr = 3*(text-1) + 1;
      adja = &mesh->adja[iadr];
      atext1 = adja[MMG2_inxt[voy]];
      atext2 = adja[MMG2_inxt[voy+1]];
      adja[0] = 3*nel + 1;
      adja[1] = 3*lel + 0;
      adja[2] = atext2;
      if(atext2)
        (&mesh->adja[3*(atext2/3-1) + 1])[atext2%3] = 3*text + 2;
      
      /*lel*/
      iadr = 3*(lel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[0] = 3*text + 1;
      adja[1] = 3*mel + 0;
      adja[2] = atext1;
      if(atext1)
        (&mesh->adja[3*(atext1/3-1) + 1])[atext1%3] = 3*lel + 2;
      if(nflat==5 || nflat==7) {
        adja[1] = 3*(aext1/3);
        if(aext1)
          (&mesh->adja[3*(aext1/3-1) + 1])[0] = 3*lel + 1;
      }
      
      /*nel*/
      iadr = 3*(nel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[1] = 3*text + 0;
      
      /*mel*/
      iadr = 3*(mel-1) + 1;
      adja = &mesh->adja[iadr];
      if(!(nflat==5 || nflat==7))
        adja[0] = 3*lel + 1;
      
    }
    if(nflat==2 || nflat==3 || nflat==7 || nflat == 6){
      /*on coupe le tr aext0 en 2*/
      text = aext0/3;
      voy  = aext0%3;
      pt1  = &mesh->tria[text];
      ie1  = pt1->v[voy];
      ie2  = pt1->v[MMG2_inxt[voy]];
      ie3  = pt1->v[MMG2_inxt[voy+1]];
      pt1->v[0] = ie1;
      pt1->v[1] = ie2;
      pt1->v[2] = k;
      pt1->qual = MMG2_caltri_in(mesh,sol,pt1);
      pt = &mesh->tria[nel];
      pt->v[0] = ie3;
      pt->v[1] = ie1;
      pt->v[2] = k;
      pt->qual = MMG2_caltri_in(mesh,sol,pt);
      /*maj des adj*/
      /*text*/
      iadr = 3*(text-1) + 1;
      adja = &mesh->adja[iadr];
      atext1 = adja[MMG2_inxt[voy]];
      atext2 = adja[MMG2_inxt[voy+1]];
      adja[0] = 3*mel + 1;
      adja[1] = 3*nel + 0;
      adja[2] = atext2;
      if(atext2)
        (&mesh->adja[3*(atext2/3-1) + 1])[atext2%3] = 3*text + 2;
      /*nel*/
      iadr = 3*(nel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[0] = 3*text + 1;
      adja[1] = 3*lel + 0;
      adja[2] = atext1;
      if(atext1)
        (&mesh->adja[3*(atext1/3-1) + 1])[atext1%3] = 3*nel + 2;
      if (nflat==3 || nflat==7) {
        adja[1] = 3*(aext2/3);
        //(&mesh->adja[3*(aext2/3-1) + 1])[aext2%3] = 3*nel + 1;
      }
      
      /*mel*/
      iadr = 3*(mel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[1] = 3*text + 0;
    }
    if(mesh->info.ddebug) {
      if ( !_MMG5_chkmsh(mesh,0,0) ) exit(EXIT_FAILURE);
    }
    
  }
  return(1);
}

