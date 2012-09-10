#include "mmg3d.h"

#define  EPSLOC   1.00005

extern char ddb;

int chkmsh(pMesh mesh,int severe,int base) {
  pTetra	pt,pt0,pt1,pt2;
  pxTetra   pxt,pxt0,pxt1,pxt2;
  int		*adja,*adja1,adj,adj1,k,i,iadr,ilists,ilistv,lists[LMAX+2],listv[LMAX+2],ip;
  int       iel,ielprv,ielnxt,l,nump,np,nq;
  int       a0,a1,a2,b0,b1,b2;
  unsigned char	voy,voy1,j,iface,ifaceprv,ifacenxt,indp,indpprv,indpnxt,tag0,tag1,tag2,ia;

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
	      exit(1);
      }
	    pt2 = &mesh->tetra[adj];
	    if ( !MG_EOK(pt2) || pt2->ref < 0 ){
	      fprintf(stdout,"  4. Invalid adjacent %d %d\n",adj,k);
	      printf("sommets k   %d: %d %d %d %d\n",k,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3]);
	      printf("sommets adj %d: %d %d %d %d\n",adj,pt2->v[0],pt2->v[1],pt2->v[2],pt2->v[3]);
	      printf("numeros adj %d: %d %d %d %d\n",k,adja[0]/4,adja[1]/4,adja[2]/4,adja[3]/4);
	      exit(1);
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
	      exit(1);
	    }
      
      a0 = pt1->v[idir[i][0]];
      a1 = pt1->v[idir[i][1]];      
      a2 = pt1->v[idir[i][2]]; 
      
      b0 = pt2->v[idir[voy][0]];
      b1 = pt2->v[idir[voy][1]];      
      b2 = pt2->v[idir[voy][2]]; 
      
      if(!(((a0 == b0)&&(a1 == b1)&&(a2 ==b2))||((a0 == b0)&&(a1 == b2)&&(a2 ==b1))\
         || ((a0 == b1)&&(a1 == b0)&&(a2 ==b2)) || ((a0 == b1)&&(a1 == b2)&&(a2 ==b0))\
         || ((a0 == b2)&&(a1 == b0)&&(a2 ==b1)) || ((a0 == b2)&&(a1 == b1)&&(a2 ==b0)) )){
        printf("Inconsistent faces : tetra %d face %d ; tetra %d face %i \n",k,i,adj,voy);
        printf("Tet 1 : %d %d %d \n",a0,a1,a2);
        printf("Tet 2 : %d %d %d \n",b0,b1,b2);
        exit(0);
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
		      saveMesh(mesh);
		      exit(0);
		    }
		    else{
		      pxt = &mesh->xtetra[pt->xt];
		      if(!(pxt->ftag[i] & MG_BDY)){
			      printf("Tetra %d : boundary face not tagged : %d \n",k,i);
			      saveMesh(mesh);
			      exit(0);
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
	  np = pt->v[iare[ia][0]];
	  nq = pt->v[iare[ia][1]];
	  
	  if(!(((np == 6)&&(nq == 3204))||((nq == 6)&&(np == 3204))))
	    continue;
	  	  
	  if(!(pxt->tag[ia] & MG_REF)){
		saveMesh(mesh);
	    exit(0);
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
	    ip = idir[i][j];
		nump = pt->v[ip];
		
	    assert(boulesurfvolp(mesh,k,ip,i,listv,&ilistv,lists,&ilists));
		
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
		    if(pt0->v[idir[iface][indp]] == nump){
			  break;
			}
		  }
		  assert(indp < 3);
		  
		  for(indpprv = 0;indpprv<3;indpprv++){
		    if(pt1->v[idir[ifaceprv][indpprv]] == nump){
			  break;
			}
		  }
		  assert(indpprv < 3);
		  
		  for(indpnxt = 0;indpnxt<3;indpnxt++){
		    if(pt2->v[idir[ifacenxt][indpnxt]] == nump){
			  break;
			}
		  }
		  assert(indpnxt < 3);
		  
		  assert(pxt0->ftag[iface] & MG_BDY);
		  assert(pxt1->ftag[ifaceprv] & MG_BDY);
		  assert(pxt2->ftag[ifacenxt] & MG_BDY);*/
		  
		  /* Cannot rely on tag boundary of edges */ 
		  /*tag0 = (pxt0->tag[iarf[iface][iprv2[indp]]]) & (~MG_BDY);
		  tag1 = (pxt1->tag[iarf[ifaceprv][inxt2[indpprv]]]) & (~MG_BDY);
		  
		  if(tag0 != tag1){
		    printf("Unconsistent tag of edge : tetra %d %d pour le point %d\n",iel,ielprv,nump);
			printf("tags : %d %d \n",tag0,tag1);

			saveMesh(mesh);
			exit(0);
		  }
		  
		  tag0 = (pxt0->tag[iarf[iface][inxt2[indp]]]) & (~MG_BDY);
		  tag2 = (pxt2->tag[iarf[ifacenxt][iprv2[indpnxt]]]) & (~MG_BDY);
		  
		  if(tag0 != tag2){
			printf("Unconsistent tag of edge : tetra %d %d pour le point %d\n",iel,ielnxt,nump);
			printf("tags : %d %d \n",tag0,tag2);
			saveMesh(mesh);
			exit(0);
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
      if ( EPSLOC*EPSLOC*dd < ray ) {
        fprintf(stdout,"  ## Non-Delaunay mesh:  %.14f < %.14f\n",dd,ray);
	exit(1);
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
	exit(1);
      }
      lon = MMG_boulep(mesh,k,i,&list);
      for (l=1; l<=lon; l++) {
        kk  = list.tetra[l] / 4;
	nk  = list.tetra[l] % 4;
	pt2 = &mesh->tetra[kk];
	if ( pt2->v[nk] != ip ) {
	  fprintf(stdout,"  5. Wrong ball %d, %d\n",ip,pt2->v[nk]);
	  exit(1);
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
        exit(1);
      }
    }
  }*/

  //fprintf(stdout,"  ** MESH STRUCTURE IS OK\n");
  return(1);
}
