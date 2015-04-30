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
 * \file mmg3d/chkmsh.c
 * \brief Check the input mesh validity.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmg3d.h"

#define  _MMG5_EPSLOC   1.00005
#define  IEDG(a,b) (((a) > 0) && ((b) > 0)) ? ((a)+(b)) : (((a)+(b))-(1))

extern char ddb;

/**
 *
 * \warning Not used.
 */
void _MMG5_chkvol(MMG5_pMesh mesh) {
  MMG5_pTetra    pt;
  int       k;
#ifdef DEBUG
  int       ier=1;
#endif

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( _MMG5_orvol(mesh->point,pt->v) < _MMG5_NULKAL ) {
      printf("  tetra %d  volume %e\n",k,_MMG5_orvol(mesh->point,pt->v));
#ifdef DEBUG
      ier = 0;
#endif
    }
  }
#ifdef DEBUG
  assert(ier);
#endif
}

/**
 *
 * \warning Not used.
 */
int _MMG5_chkmshsurf(MMG5_pMesh mesh){
  MMG5_pTria      pt;
  int        k,k1;
  int        *adja,*adja1;
  char       i,voy;

  for (k=1; k<=mesh->nt; k++) {
    pt   = &mesh->tria[k];
    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_NOM )  continue;
      k1  = adja[i] / 3;
      voy = adja[i] % 3;

      if(!k1) continue;
      adja1 = &mesh->adjt[3*(k1-1)+1];

      if(adja1[voy] / 3 != k){
        printf("Wrong adjacency relation for triangles : %d %d \n",k,k1);
        exit(EXIT_FAILURE);
      }
    }
  }
  return(1);
}

int _MMG5_mmg3dChkmsh(MMG5_pMesh mesh,int severe,int base) {
  MMG5_pTetra    pt,pt1,pt2;
  MMG5_pxTetra   pxt;
  int            *adja,*adja1,adj,adj1,k,i,iadr;
  int            iel,a0,a1,a2,b0,b1,b2;
  unsigned char  voy,voy1;
  /* commentated part variables
     MMG5_pTetra        pt0;
     MMG5_xTetra       pxt0,pxt1,pxt2;
     int           ilists,ilistv,lists[_MMG5_LMAX+2],listv[_MMG5_LMAX+2];
     int           ielprv,ielnxt,l,nump,np,nq;
     unsigned char j,iface,ifaceprv,ifacenxt,indp,indpprv,indpnxt,tag0,tag1,tag2,ia;
  */

  for (k=1; k<=mesh->ne; k++) {
    pt1 = &mesh->tetra[k];
    if ( !MG_EOK(pt1) || pt1->ref < 0 )   continue;
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];

    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      voy = adja[i] % 4;
      if ( !adj )  continue;

      if ( adj == k ) {
        fprintf(stdout,"  1. Wrong adjacency %d %d\n",k,adj);
        printf("k %d: %d %d %d %d\n",k,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3]);
        printf("adj (%d): %d %d %d %d\n",k,adja[0]/4,adja[1]/4,adja[2]/4,adja[3]/4);
        exit(EXIT_FAILURE);
      }
      pt2 = &mesh->tetra[adj];
      if ( !MG_EOK(pt2) || pt2->ref < 0 ){
        fprintf(stdout,"  4. Invalid adjacent %d %d\n",adj,k);
        printf("sommets k   %d: %d %d %d %d\n",k,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3]);
        printf("sommets adj %d: %d %d %d %d\n",adj,pt2->v[0],pt2->v[1],pt2->v[2],pt2->v[3]);
        printf("numeros adj %d: %d %d %d %d\n",k,adja[0]/4,adja[1]/4,adja[2]/4,adja[3]/4);
        exit(EXIT_FAILURE);
      }
      iadr  = (adj-1)*4 + 1;
      adja1 = &mesh->adja[iadr];
      adj1  = adja1[voy] / 4;
      voy1  = adja1[voy] % 4;
      if ( adj1 != k || voy1 != i ) {
        fprintf(stdout,"  2. Wrong adjacency %d %d\n",k,adj1);
        printf("k %d: %d %d %d %d\n",k,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3]);
        printf("a %d: %d %d %d %d\n",adj,pt2->v[0],pt2->v[1],pt2->v[2],pt2->v[3]);
        printf("adj(%d): %d %d %d %d\n",k,adja[0]/4,adja[1]/4,adja[2]/4,adja[3]/4);
        printf("adj(%d): %d %d %d %d\n",adj,adja1[0]/4,adja1[1]/4,adja1[2]/4,adja1[3]/4);
        exit(EXIT_FAILURE);
      }

      a0 = pt1->v[_MMG5_idir[i][0]];
      a1 = pt1->v[_MMG5_idir[i][1]];
      a2 = pt1->v[_MMG5_idir[i][2]];

      b0 = pt2->v[_MMG5_idir[voy][0]];
      b1 = pt2->v[_MMG5_idir[voy][1]];
      b2 = pt2->v[_MMG5_idir[voy][2]];

      if(!(((a0 == b0)&&(a1 == b1)&&(a2 ==b2))||((a0 == b0)&&(a1 == b2)&&(a2 ==b1))\
           || ((a0 == b1)&&(a1 == b0)&&(a2 ==b2)) || ((a0 == b1)&&(a1 == b2)&&(a2 ==b0))\
           || ((a0 == b2)&&(a1 == b0)&&(a2 ==b1)) || ((a0 == b2)&&(a1 == b1)&&(a2 ==b0)) )){
        printf("Inconsistent faces : tetra %d face %d ; tetra %d face %i \n",k,i,adj,voy);
        printf("Tet 1 : %d %d %d \n",a0,a1,a2);
        printf("Tet 2 : %d %d %d \n",b0,b1,b2);
        exit(EXIT_FAILURE);
      }
    }
  }

  /* This test may have to be disactivated : check wether boundary faces (i.e. no neighbour)
     arise only when a BDY face is hit */
  for(k=1;k<=mesh->ne;k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )   continue;
    adja = &mesh->adja[4*(k-1)+1];
    for(i=0;i<4;i++){
      if(!adja[i]){
        if(!pt->xt){
          printf("Tetra %d : boundary face not tagged : %d \n",k,i);
          MMG5_saveMesh(mesh);
          exit(EXIT_FAILURE);
        }
        else{
          pxt = &mesh->xtetra[pt->xt];
          if(!(pxt->ftag[i] & MG_BDY)){
            printf("Tetra %d : boundary face not tagged : %d \n",k,i);
            MMG5_saveMesh(mesh);
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  }

  /* Case of an implicit surface embedded in mesh : check whether a face separating two
     different subdomains is tagged bdy */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )   continue;

    adja = &mesh->adja[4*(k-1)+1];
    for(i=0; i<4; i++){
      if(!adja[i]) continue;
      iel = adja[i] / 4;
      pt1 = &mesh->tetra[iel];

      if(pt->ref != pt1->ref){
        if(!pt->xt){
          printf("Tetra %d face %d : common face is a limit of two subdomains and has not xt : %d %d %d  \n",k,i,pt->v[_MMG5_idir[i][0]],pt->v[_MMG5_idir[i][1]],pt->v[_MMG5_idir[i][2]]);
          MMG5_saveMesh(mesh);
          exit(EXIT_FAILURE);
        }
        else{
          pxt = &mesh->xtetra[pt->xt];
          if(!(pxt->ftag[i] & MG_BDY)){
            printf("Tetra %d %d : common face is a limit of two subdomains and is not tagged %d %d %d -->%d\n",k,i,pt->v[_MMG5_idir[i][0]],pt->v[_MMG5_idir[i][1]],pt->v[_MMG5_idir[i][2]], pxt->ftag[i]);
            MMG5_saveMesh(mesh);
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  }

  /*for(k=1;k<=mesh->ne;k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )   continue;
    if(!pt->xt) continue;

    pxt = &mesh->xtetra[pt->xt];

    for(ia = 0;ia<6;ia++){
    np = pt->v[_MMG5_iare[ia][0]];
    nq = pt->v[_MMG5_iare[ia][1]];

    if(!(((np == 6)&&(nq == 3204))||((nq == 6)&&(np == 3204))))
    continue;

    if(!(pxt->tag[ia] & MG_REF)){
    MMG5_saveMesh(mesh);
    exit(EXIT_FAILURE);
    }
    }
    }*/

  /* Test Boundary references : do they match, from one face to its neighbour ? */
  /*for(k=1;k<=mesh->ne;k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )   continue;
    if(!pt->xt) continue;

    pxt = &mesh->xtetra[pt->xt];

    for(i=0;i<4;i++){
    if(!(pxt->ftag[i] & MG_BDY)) continue;

    for(j=0;j<3;j++){
    ip = _MMG5_idir[i][j];
    nump = pt->v[ip];

    if(!_MMG5_boulesurfvolp(mesh,k,ip,i,listv,&ilistv,lists,&ilists)){
    printf("%s:%d: Error: function _MMG5_boulesurfvolp return 0\n",__FILE__,__LINE__)
    exit(EXIT_FAILURE);
    }

    for(l=0;l<ilists;l++){
    iel = lists[l]/4;
    iface = lists[l]%4;

    if(l==0){
    ielprv = lists[ilists-1]/4;
    ifaceprv = lists[ilists-1]%4;
    }
    else{
    ielprv = lists[l-1]/4;
    ifaceprv = lists[l-1]%4;
    }

    if(l==ilists-1){
    ielnxt = lists[0]/4;
    ifacenxt = lists[0]%4;
    }
    else{
    ielnxt = lists[l+1]/4;
    ifacenxt = lists[l+1]%4;
    }

    pt0 = &mesh->tetra[iel];
    pt1 = &mesh->tetra[ielprv];
    pt2 = &mesh->tetra[ielnxt];
    assert(pt0->xt && pt1->xt && pt2->xt);
    pxt0 = &mesh->xtetra[pt0->xt];
    pxt1 = &mesh->xtetra[pt1->xt];
    pxt2 = &mesh->xtetra[pt2->xt];

    for(indp = 0;indp<3;indp++){
    if(pt0->v[_MMG5_idir[iface][indp]] == nump){
    break;
    }
    }
    assert(indp < 3);

    for(indpprv = 0;indpprv<3;indpprv++){
    if(pt1->v[_MMG5_idir[ifaceprv][indpprv]] == nump){
    break;
    }
    }
    assert(indpprv < 3);

    for(indpnxt = 0;indpnxt<3;indpnxt++){
    if(pt2->v[_MMG5_idir[ifacenxt][indpnxt]] == nump){
    break;
    }
    }
    assert(indpnxt < 3);

    assert(pxt0->ftag[iface] & MG_BDY);
    assert(pxt1->ftag[ifaceprv] & MG_BDY);
    assert(pxt2->ftag[ifacenxt] & MG_BDY);*/

  /* Cannot rely on tag boundary of edges */
  /*tag0 = (pxt0->tag[_MMG5_iarf[iface][_MMG5_iprv2[indp]]]) & (~MG_BDY);
    tag1 = (pxt1->tag[_MMG5_iarf[ifaceprv][inxt2[indpprv]]]) & (~MG_BDY);

    if(tag0 != tag1){
    printf("Unconsistent tag of edge : tetra %d %d pour le point %d\n",iel,ielprv,nump);
    printf("tags : %d %d \n",tag0,tag1);

    MMG5_saveMesh(mesh);
    exit(EXIT_FAILURE);
    }

    tag0 = (pxt0->tag[_MMG5_iarf[iface][inxt2[indp]]]) & (~MG_BDY);
    tag2 = (pxt2->tag[_MMG5_iarf[ifacenxt][_MMG5_iprv2[indpnxt]]]) & (~MG_BDY);

    if(tag0 != tag2){
    printf("Unconsistent tag of edge : tetra %d %d pour le point %d\n",iel,ielnxt,nump);
    printf("tags : %d %d \n",tag0,tag2);
    MMG5_saveMesh(mesh);
    exit(EXIT_FAILURE);
    }
    }
    }
    }

    }*/

  /* Delaunay criterion */
  /*
    for (k=1; k<=mesh->ne; k++) {
    pt1 = &mesh->tetra[k];
    if ( !pt1->v[0] )  continue;
    iadr = (k-1)*4 + 1;
    adja = &mesh->adja[iadr];
    if ( !cenrad(mesh,k,c,&ray) )  continue;

    for (i=0; i<4; i++) {
    if ( !adja[i] )  continue;
    adj = adja[i] / 4;
    voy = adja[i] % 4;
    pt2 = &mesh->tetra[adj];

    ppt = &mesh->point[ pt2->v[voy] ];
    dd = (ppt->c[0] - c[0]) * (ppt->c[0] - c[0]) \
    + (ppt->c[1] - c[1]) * (ppt->c[1] - c[1]) \
    + (ppt->c[2] - c[2]) * (ppt->c[2] - c[2]);
    if ( _MMG5_EPSLOC*_MMG5_EPSLOC*dd < ray ) {
    fprintf(stdout,"  ## Non-Delaunay mesh:  %.14f < %.14f\n",dd,ray);
    exit(EXIT_FAILURE);
    }
    }
    }
  */

  /*  if ( !severe )  return(1);

      for (k=1; k<=mesh->ne; k++) {
      pt1 = &mesh->tetra[k];
      if ( !pt1->v[0] )  continue;
      else if (pt1->flag < base )  continue;
      iadr = 4*(k-1) + 1;
      adja = &mesh->adja[iadr];

      for (i=0; i<4; i++) {
      adj = (adja[i]-1) / 4 + 1;
      voy = (adja[i]-1) % 4;
      if ( !adj )  continue;

      ip  = pt1->v[i];
      ppt = &mesh->point[ip];
      if ( ppt->tag & M_UNUSED ) {
      fprintf(stdout,"  6. Unused vertex %d  %d\n",k,ip);
      printf("%d %d %d %d\n",pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3]);
      exit(EXIT_FAILURE);
      }
      lon = MMG_boulep(mesh,k,i,&list);
      for (l=1; l<=lon; l++) {
      kk  = list.tetra[l] / 4;
      nk  = list.tetra[l] % 4;
      pt2 = &mesh->tetra[kk];
      if ( pt2->v[nk] != ip ) {
      fprintf(stdout,"  5. Wrong ball %d, %d\n",ip,pt2->v[nk]);
      exit(EXIT_FAILURE);
      }
      }
      if ( lon < 1 )  continue;
      len = 0;
      for (kk=1; kk<=mesh->ne; kk++) {
      pt2 = &mesh->tetra[kk];
      if ( !pt2->v[0] )  continue;
      for (j=0; j<4; j++)
      if ( pt2->v[j] == ip ) {
      len++;
      break;
      }
      }
      if ( len != lon ) {
      fprintf(stdout,"  7. Incorrect ball %d: %d %d\n",pt1->v[i],lon,len);
      exit(EXIT_FAILURE);
      }
      }
      }*/

  //fprintf(stdout,"  ** MESH STRUCTURE IS OK\n");
  return(1);
}

/**
 * Search boundary faces containing point np.
 *
 * \warning Not used.
 **/
int _MMG5_chkptonbdy(MMG5_pMesh mesh,int np){
  MMG5_pTetra      pt;
  MMG5_pxTetra     pxt;
  MMG5_pPoint      p0;
  int         k;
  char        i,j,ip;

  for(k=1;k<=mesh->np;k++)
    mesh->point[k].flag = 0;

  /* Put flag = 1 at each point belonging to a boundary face */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if(!MG_EOK(pt)) continue;
    if(!pt->xt) continue;
    pxt = &mesh->xtetra[pt->xt];
    for(i=0; i<4; i++){
      if(!(pxt->ftag[i] & MG_BDY)) continue;
      for(j=0; j<3; j++){
        ip = _MMG5_idir[i][j];
        if(pt->v[ip] == np) printf("Le pt : %d sur la face %d du tetra %d : %d %d %d %d \n",pt->v[ip],i,k,pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
        p0 = &mesh->point[pt->v[ip]];
        p0->flag = 1;
      }
    }
  }

  /* Make sure that all the remaining points are not tagged BDY */
  for(k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    if(!MG_VOK(p0)) continue;
    if(p0->flag) continue;
    if(p0->tag & MG_BDY){
      printf("      Fct. chkptonbdy : point %d tagged bdy while belonging to no BDY face\n",k);
      exit(EXIT_FAILURE);
    }
  }

  return(1);
}

/**
 *
 * Count how many boundary faces share point nump.
 *
 * \warning Not used.
 */
int _MMG5_cntbdypt(MMG5_pMesh mesh, int nump){
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  int k,nf;
  char i,j,ip;

  nf = 0;

  for(k=1; k<=mesh->ne;k++){
    pt = &mesh->tetra[k];
    if(!MG_EOK(pt)) continue;
    if(!pt->xt) continue;
    pxt = &mesh->xtetra[pt->xt];
    for(i=0; i<4; i++){
      if(!(pxt->ftag[i] & MG_BDY)) continue;
      for(j=0; j<3; j++){
        ip = _MMG5_idir[i][j];
        if(pt->v[ip] == nump){
          printf("La face : %d %d %d \n dans le tetra : %d %d %d %d \n",pt->v[_MMG5_idir[i][0]],pt->v[_MMG5_idir[i][1]],pt->v[_MMG5_idir[i][2]],pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
          nf++;
        }
      }
    }
  }
  return(nf);
}

/** Count the number of tetras that have several boundary faces, as well as the number of internal
    edges connecting points of the boundary */
int _MMG5_chkfemtopo(MMG5_pMesh mesh) {
  MMG5_pTetra      pt,pt1;
  MMG5_pxTetra     pxt;
  MMG5_pPoint      p0,p1;
  int         k,nf,ntet,ned,np,ischk,ilist,list[_MMG5_LMAX+2],l,np1,npchk,iel;
  char        i0,j,i,i1,ia;

  ntet = ned = 0;
  for(k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Count elements with at least two boundary faces */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    else if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    nf = 0;
    for (i=0; i<4; i++) {
      if ( pxt->ftag[i] & MG_BDY )  nf++;
    }
    if ( nf >= 2 )  ntet++;
  }
  if ( ntet )  printf("  *** %d tetras with at least 2 boundary faces.\n",ntet);

  /* Count internal edges connecting two points of the boundary */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<4; i++) {
      np = pt->v[i];
      p0 = &mesh->point[np];
      if ( !(p0->tag & MG_BDY) )  continue;

      ischk = p0->flag % 2;
      if ( ischk )  continue;
      p0->flag += 1;

      ilist = _MMG5_boulevolp(mesh,k,i,list);
      for (l=0; l<ilist; l++) {
        iel = list[l] / 4;
        i0  = list[l] % 4;
        i1  = i0;

        pt1 = &mesh->tetra[iel];
        for (j=0; j<3; j++) {
          i1  = _MMG5_inxt3[i1];
          np1 = pt1->v[i1];
          if ( np1 < np )  continue;
          p1 = &mesh->point[np1];
          if ( !(p1->tag & MG_BDY) )  continue;

          ischk = p1->flag % 2;
          npchk = p1->flag / 2;
          if ( npchk == np )  continue;

          ia = IEDG(i0,i1);
          p1->flag = 2*np + ischk;
          if ( !_MMG5_srcbdy(mesh,iel,ia) )  ned++;
        }
      }
    }
  }
  if ( ned )  printf("  *** %d internal edges connecting boundary points.\n",ned);
  return(1);
}

/**
 *
 * Search face n0,n1,n2 in mesh, and get the support tetras, with the
 * corresponding refs.
 *
 * \warning Not used.
 */
int srcface(MMG5_pMesh mesh,int n0,int n1,int n2) {
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  int       k,ip0,ip1,ip2,minn,maxn,sn,mins,maxs,sum,ref;
  char      i,tag;

  minn = MG_MIN(n0,MG_MIN(n1,n2));
  maxn = MG_MAX(n0,MG_MAX(n1,n2));
  sn   = n0 + n1 + n2;
  pxt = 0;

  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) ) continue;

    if( pt->xt ) pxt = &mesh->xtetra[pt->xt];
    for(i=0; i<4; i++) {
      ip0 = pt->v[_MMG5_idir[i][0]];
      ip1 = pt->v[_MMG5_idir[i][1]];
      ip2 = pt->v[_MMG5_idir[i][2]];

      mins = MG_MIN(ip0,MG_MIN(ip1,ip2));
      maxs = MG_MAX(ip0,MG_MAX(ip1,ip2));
      sum  = ip0 + ip1 + ip2;
      tag  = pt->xt ? pxt->ftag[i] : 0;
      ref  = pt->xt ? pxt->ref[i] : 0;

      if( mins == minn && maxs == maxn && sum == sn ) {
        printf("Face %d in tetra %d with ref %d : corresponding ref %d , tag : %d\n",i,k,pt->ref,ref,tag);
      }
    }
  }


  return(1);
}
