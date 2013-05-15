#include "mmg3d.h"

extern Info  info;
extern char  ddb;
#ifdef DEBUG
double lmintmp,lmaxtmp,lentmp;
extern double      tabtmp[12][7];
int tabeltmp[4][7],iatmp,ibtmp,ictmp;
#endif
/** Table that associates to each (even) permutation of the 4 vertices of a tetrahedron
 *  the corresponding permutation of its edges. Labels :
 *  0  : [0,1,2,3]
 *  1  : [0,2,3,1]
 *  2  : [0,3,1,2]
 *  3  : [1,0,3,2]
 *  4  : [1,2,0,3]
 *  5  : [1,3,2,0]
 *  6  : [2,0,1,3]
 *  7  : [2,1,3,0]
 *  8  : [2,3,0,1]
 *  9  : [3,0,2,1]
 *  10 : [3,1,0,2]
 *  11 : [3,2,1,0]
 *  The edge 0 of the config 1 become the edge 1 of the reference config so permedge[1][0]=1 ...
 */

unsigned char permedge[12][6] = {
  {0,1,2,3,4,5}, {1,2,0,5,3,4}, {2,0,1,4,5,3}, {0,4,3,2,1,5},
  {3,0,4,1,5,2}, {4,3,0,5,2,1}, {1,3,5,0,2,4}, {3,5,1,4,0,2},
  {5,1,3,2,4,0}, {2,5,4,1,0,3}, {4,2,5,0,3,1}, {5,4,2,3,1,0} };

/** simulate split 1 edge of tetra : return 0 if split leads to invalid situation, else 1 */
int split1_sim(pMesh mesh,pSol met,int k,int vx[6]) {
  pTetra   pt,pt0;
  double   vold,vnew;
  unsigned char tau[4],*taued;

  /* tau = sigma^-1 = permutation that sends the reference config (edge 01 split) to the current */
  pt = &mesh->tetra[k];
  vold = orvol(mesh->point,pt->v);
  pt0 = &mesh->tetra[0];

  /* default is case 1 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &permedge[0][0];
  switch(pt->flag) {
  case 2:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &permedge[6][0];
    break;
  case 4:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &permedge[2][0];
    break;
  case 8:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &permedge[4][0];
    break;
  case 16:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &permedge[10][0];
    break;
  case 32:
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    taued = &permedge[11][0];
    break;
  }

  /* Test volume of the two created tets */
  memcpy(pt0,pt,sizeof(Tetra));
  pt0->v[tau[1]] = vx[taued[0]];
  vnew = orvol(mesh->point,pt0->v);
  if ( vnew < EPSD2 )  return(0);
  else if ( vold > NULKAL && vnew < NULKAL )  return(0);

  memcpy(pt0,pt,sizeof(Tetra));
  pt0->v[tau[0]] = vx[taued[0]];
  vnew = orvol(mesh->point,pt0->v);
  if ( vnew < EPSD2 )  return(0);
  else if ( vold > NULKAL && vnew < NULKAL )  return(0);

  return(1);
}

/** split 1 edge of tetra */
void split1(pMesh mesh,pSol met,int k,int vx[6]) {
  pTetra   pt,pt1;
  xTetra   xt,xt1;
  pxTetra  pxt0;
  int      iel,ref;
  char     i,tag,isxt,isxt1;
  unsigned char tau[4],*taued;

  /* create a new tetra */
  pt  = &mesh->tetra[k];
#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt->v[iare[i][0]],pt->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        //       tabtmp[0][1]++;
        lmintmp=lentmp;
      }else if(lentmp>lmaxtmp) {
        //tabtmp[0][2]++;
        lmaxtmp=lentmp;
      }
    }
#endif
  iel = newElt(mesh);
  assert(iel);
  pt1 = &mesh->tetra[iel];
  pt1 = memcpy(pt1,pt,sizeof(Tetra));
  pxt0 = 0;
  if ( pt->xt ) {
    pxt0 = &mesh->xtetra[pt->xt];
    memcpy(&xt,pxt0,sizeof(xTetra));
    memcpy(&xt1,pxt0,sizeof(xTetra));
  }
  else {
    memset(&xt,0,sizeof(xTetra));
    memset(&xt1,0,sizeof(xTetra));
  }

  /* default is case 1 */
  tau[0] = 0; tau[1] = 1; tau[2] = 2; tau[3] = 3;
  taued = &permedge[0][0];
  switch(pt->flag) {
  case 2:
    tau[0] = 2; tau[1] = 0; tau[2] = 1; tau[3] = 3;
    taued = &permedge[6][0];
    break;
  case 4:
    tau[0] = 0; tau[1] = 3; tau[2] = 1; tau[3] = 2;
    taued = &permedge[2][0];
    break;
  case 8:
    tau[0] = 1; tau[1] = 2; tau[2] = 0; tau[3] = 3;
    taued = &permedge[4][0];
    break;
  case 16:
    tau[0] = 3; tau[1] = 1; tau[2] = 0; tau[3] = 2;
    taued = &permedge[10][0];
    break;
  case 32:
    tau[0] = 3; tau[1] = 2; tau[2] = 1; tau[3] = 0;
    taued = &permedge[11][0];
    break;
  }
  /* delete old edge */
  if ( hPop(&mesh->htab,pt->v[tau[0]],pt->v[tau[1]],&ref,&tag) ) {
    /* Generic formulation of split of 1 edge */
    pt->v[tau[1]] = pt1->v[tau[0]] = vx[taued[0]];
    hEdge(&mesh->htab,pt->v[tau[0]],pt->v[tau[1]],ref,tag);
    hEdge(&mesh->htab,pt1->v[tau[0]],pt1->v[tau[1]],ref,tag);
  }
  else
    pt->v[tau[1]] = pt1->v[tau[0]] = vx[taued[0]];

  if ( pt->xt ) {
    /* Reset edge tag */
    xt.tag[ taued[3]] = xt.ftag[ tau[3]];
    xt.tag[ taued[4]] = xt.ftag[ tau[2]];
    xt1.tag[taued[1]] = xt1.ftag[tau[3]];
    xt1.tag[taued[2]] = xt1.ftag[tau[2]];
    xt.ref[ tau[0]] = 0;  xt.ftag[ tau[0]] = 0;
    xt1.ref[tau[1]] = 0;  xt1.ftag[tau[1]] = 0;
  }

  pt->flag = pt1->flag = 0;
  isxt = 0 ;
  isxt1 = 0;
  for (i=0; i<4; i++) {
    if ( xt.ref[i]  || xt.ftag[i] )  isxt = 1;
    if ( xt1.ref[i] || xt1.ftag[i] ) isxt1 = 1;
    if ( isxt && isxt1 )  break;
  }

  if ( pt->xt ) {
    if ( isxt && !isxt1 ) {
      pt1->xt = 0;
      memcpy(pxt0,&xt,sizeof(xTetra));
    }
    else if ( !isxt && isxt1 ) {
      pt1->xt = pt->xt;
      pt->xt  = 0;
      pxt0 = &mesh->xtetra[pt1->xt];
      memcpy(pxt0,&xt1,sizeof(xTetra));
    }
    else if ( isxt && isxt1 ) {
      mesh->xt++;
      assert(mesh->xt < mesh->xtmax);
      pt1->xt = mesh->xt;
      memcpy(pxt0,&xt,sizeof(xTetra));
      pxt0 = &mesh->xtetra[mesh->xt];
      memcpy(pxt0,&xt1,sizeof(xTetra));
    }
    else {
      pt->xt = 0;
      pt1->xt = 0;
    }
  }
  /* Quality update */
  pt->qual=orcal(mesh,k);
  pt1->qual=orcal(mesh,iel);
#ifdef DEBUG
  tabtmp[0][3]=lmintmp;  tabtmp[0][4]=lmaxtmp;
  tabtmp[0][5]=DBL_MAX;  tabtmp[0][6]=DBL_MIN;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt->v[iare[i][0]],pt->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[0][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[0][2]++;
      }
      lentmp=lenedg(mesh,met,pt1->v[iare[i][0]],pt1->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[0][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[0][2]++;
      }
    }
#endif
}

/** Simulate at the same time creation and bulging of one point, with new position o,
    to be inserted at an edge, whose shell is passed :
    return 0 if final position is invalid, 1 if all checks are ok */
int simbulgept(pMesh mesh,int *list,int ret,double o[3]) {
  pTetra    pt,pt0;
  pPoint    ppt0;
  double    calold,calnew,caltmp;
  int       k,iel,ilist;
  char      ie,ia,ib;

  ilist = ret / 2;
  pt0  = &mesh->tetra[0];
  ppt0 = &mesh->point[0];
  ppt0->c[0] = o[0];
  ppt0->c[1] = o[1];
  ppt0->c[2] = o[2];

  calold = calnew = DBL_MAX;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 6;
    ie  = list[k] % 6;
    ia = iare[ie][0];
    ib = iare[ie][1];

    pt = &mesh->tetra[iel];
    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[ia] = 0;
    calold = MG_MIN(calold,pt->qual);
    caltmp = orcal(mesh,0);
    if ( caltmp < EPSD )  return(0);
    calnew = MG_MIN(calnew,caltmp);

    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[ib] = 0;
    caltmp = orcal(mesh,0);
    if ( caltmp < EPSD )  return(0);
    calnew = MG_MIN(calnew,caltmp);
  }
  /*if ( calold < NULKAL && calnew <= calold )  return(0);
    else if ( calnew < 0.3*calold )  return(0);*/


  return(1);
}

/** Split edge list[0]%6, whose shell list is passed, introducing point ip
    Beware : shell has to be enumerated in ONLY ONE TRAVEL (always same sense) */
int split1b(pMesh mesh, pSol met,int *list, int ret, int ip,int cas){
  pTetra    pt,pt1;
  xTetra    xt,xt1;
  pxTetra   pxt0;
  int       ilist,k,open,iel,jel,*newtet,nump,numq,*adja;
  int       *adjan,nei2,nei3,mel,ref;
  char      ie,tau[4],*taued,isxt,isxt1,i,j,voy,tag;
  double    lmin,lmax,len;

  ilist = ret / 2;
  open  = ret % 2;

  if(met->m){
    lmin = 0.6;
    lmax = 1.3;
    for (j=0; j<ilist; j++) {
      for(i=0;i<6;i++) {
        len = lenedg(mesh,met, mesh->tetra[list[j]/6].v[iare[i][0]],
                     mesh->tetra[list[j]/6].v[iare[i][1]]);
        if ( len < lmin) {
          lmin = len;
        }
        else if ( len > lmax) {
          lmax = len;
        }
      }
    }
    /* cree-t-on une trop petite arete ? (voir le bug de BUG_Split1b_SpereIso_0.125h_met) */
    if ( cas ) {
      //lmin=0.6;
      for (j=0; j<ilist; j++) {
        iel = list[j] / 6;
        pt  = &mesh->tetra[iel];
        ie  = list[j] % 6;
        len = lenedg(mesh,met, pt->v[isar[ie][0]],ip);
        if ( len < lmin )  break;
        len = lenedg(mesh,met, pt->v[isar[ie][1]],ip);
        if ( len < lmin )  break;
      }
      if ( j < ilist )  return(0);
    }
  }


  newtet = (int*)calloc(ilist,sizeof(int));
  assert(newtet);

  iel = list[0] / 6;
  ie  = list[0] % 6;
  pt  = &mesh->tetra[iel];

  nump = pt->v[iare[ie][0]];
  numq = pt->v[iare[ie][1]];

  /* If need be, update edge references */
  hPop(&mesh->htab,nump,numq,&ref,&tag);
  if ( ref || tag ) {
    hEdge(&mesh->htab,nump,ip,ref,tag);
    hEdge(&mesh->htab,numq,ip,ref,tag);
  }

  /* Fill list newtet[k] = +_created tetra for list[k]/6 : + if kept tetra (= one associated to
     pt->v[tau[0]]) is associated with nump, - if with numq */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 6;
    ie  = list[k] % 6;
    pt  = &mesh->tetra[iel];
    /* identity : case 0 */
    tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
    switch(ie) {
    case 1:
      tau[0] = 2; tau[1] = 0; tau[2] = 1; tau[3] = 3;
      break;
    case 2:
      tau[0] = 0; tau[1] = 3; tau[2] = 1; tau[3] = 2;
      break;
    case 3:
      tau[0] = 1; tau[1] = 2; tau[2] = 0; tau[3] = 3;
      break;
    case 4:
      tau[0] = 3; tau[1] = 1; tau[2] = 0; tau[3] = 2;
      break;
    case 5:
      tau[0] = 3; tau[1] = 2; tau[2] = 1; tau[3] = 0;
      break;
    }
    jel = newElt(mesh);
    if(!jel){
      fprintf(stdout,"%s:%d: Error: unable to allocate a new element\n"
              ,__FILE__,__LINE__);
      for(;k>0;--k){
        delElt(mesh,newtet[k]);
      }
      return(-1);
    }
    pt1 = &mesh->tetra[jel];
    memcpy(pt1,pt,sizeof(Tetra));

    if ( pt->v[tau[0]] == nump )
      newtet[k] = jel;
    else
      newtet[k] = -jel;
  }

  /* Special case : only one element in the shell */
  if ( ilist == 1 ) {
    assert(open);
    iel = list[0] / 6;
    ie  = list[0] % 6;
    pt = &mesh->tetra[iel];
    jel = fabs(newtet[0]);
    pt1 = &mesh->tetra[jel];

    pxt0 = 0;
    if ( pt->xt ) {
      pxt0 = &mesh->xtetra[pt->xt];
      memcpy(&xt,pxt0,sizeof(xTetra));
      memcpy(&xt1,pxt0,sizeof(xTetra));
    }
    else {
      memset(&xt,0, sizeof(xTetra));
      memset(&xt1,0, sizeof(xTetra));
    }

    /* tau = sigma^-1 = permutation that sends the reference config (edge 01 split) to current */
    tau[0] = 0; tau[1] = 1; tau[2] = 2; tau[3] = 3;
    taued = &permedge[0][0];
    switch(ie){
    case 1:
      tau[0] = 2; tau[1] = 0; tau[2] = 1; tau[3] = 3;
      taued = &permedge[6][0];
      break;
    case 2:
      tau[0] = 0; tau[1] = 3; tau[2] = 1; tau[3] = 2;
      taued = &permedge[2][0];
      break;
    case 3:
      tau[0] = 1; tau[1] = 2; tau[2] = 0; tau[3] = 3;
      taued = &permedge[4][0];
      break;
    case 4:
      tau[0] = 3; tau[1] = 1; tau[2] = 0; tau[3] = 2;
      taued = &permedge[10][0];
      break;
    case 5:
      tau[0] = 3; tau[1] = 2; tau[2] = 1; tau[3] = 0;
      taued = &permedge[11][0];
      break;
    }

    /* Generic formulation of split of 1 edge */
    pt->v[tau[1]] = pt1->v[tau[0]] = ip;
    if(pt->xt){
      /* Reset edge tag */
      xt.tag[ taued[3]] = xt.ftag[ tau[3]];
      xt.tag[ taued[4]] = xt.ftag[ tau[2]];
      xt1.tag[taued[1]] = xt1.ftag[tau[3]];
      xt1.tag[taued[2]] = xt1.ftag[tau[2]];
      xt.ref[tau[0]]   = 0;  xt.ftag[tau[0]]  = 0;
      xt1.ref[tau[1]]  = 0;  xt1.ftag[tau[1]] = 0;
    }

    pt->flag = pt1->flag = 0;

    isxt = 0 ;
    isxt1 = 0;

    for(i=0;i<4;i++){
      if(xt.ref[i] || xt.ftag[i]) isxt = 1;
      if(xt1.ref[i] || xt1.ftag[i]) isxt1 = 1;
    }

    if(pt->xt){
      if((isxt)&&(!isxt1)){
        pt1->xt = 0;
        pxt0 = &mesh->xtetra[pt->xt];
        memcpy(pxt0,&xt,sizeof(xTetra));
      }
      else if ((!isxt)&&(isxt1)){
        pt1->xt = pt->xt;
        pt->xt = 0;
        pxt0 = &mesh->xtetra[pt1->xt];
        memcpy(pxt0,&xt1,sizeof(xTetra));
      }
      else if (isxt && isxt1){
        mesh->xt++;
        assert(mesh->xt < mesh->xtmax);
        pt1->xt = mesh->xt;
        pxt0 = &mesh->xtetra[pt->xt];
        memcpy(pxt0,&xt,sizeof(xTetra));
        pxt0 = &mesh->xtetra[pt1->xt];
        memcpy(pxt0,&xt1,sizeof(xTetra));
      }
      else{
        pt->xt = 0;
        pt1->xt = 0;
      }
    }

    /* Update of adjacency relations */
    adja = &mesh->adja[4*(iel-1)+1];
    adjan = &mesh->adja[4*(jel-1)+1];

    adja[tau[2]] = 0;
    adja[tau[3]] = 0;
    adjan[tau[2]] = 0;
    adjan[tau[3]] = 0;

    mel = adja[tau[0]] / 4;
    voy = adja[tau[0]] % 4;
    adja[tau[0]] = 4*jel + tau[1];
    adjan[tau[0]] = 4*mel + voy;
    adjan[tau[1]] = 4*iel + tau[0];

    if(mel){
      adjan = &mesh->adja[4*(mel -1) +1];
      adjan[voy] = 4*jel + tau[0];
    }
    /* Quality update */
    pt->qual=orcal(mesh,iel);
    pt1->qual=orcal(mesh,jel);

#ifdef DEBUG
    for(j=0;j<ilist;j++){
      for(i=0;i<6;i++)
        {
          len=lenedg(mesh,met,
                     mesh->tetra[list[j]/6].v[iare[i][0]],
                     mesh->tetra[list[j]/6].v[iare[i][1]]);
          if(len<lmin) {
            tabtmp[0][1]+=0.5;
          }else if(len>lmax) {
            tabtmp[0][2]+=0.5;
          }
        }
    }
    tabtmp[0][0]+=ilist;
#endif
    free(newtet);
    newtet=NULL;
    return(1);
  }

  /*if(ddb) printf("le tet 0 : %d son new : %d, le tet 1 %d et son new : %d \n",list[0]/6,newtet[0],list[1]/6,newtet[1]);

    if(ddb){
    printf("Ecrit la coquille de %d %d ouverte ? %d: \n",nump,numq,open);
    for(k=0;k<ilist;k++){
    printf("coquille ; %d arete %d et le nouveau %d \n",list[k]/6,list[k]%6,newtet[k]);
    }
    }*/

  /* General case : update each element of the shell */
  for(k=0; k<ilist; k++){
    iel = list[k] / 6;
    ie  = list[k] % 6;
    pt = &mesh->tetra[iel];
    jel = fabs(newtet[k]);
    pt1 = &mesh->tetra[jel];

    pxt0 = 0;
    if ( pt->xt ) {
      pxt0 = &mesh->xtetra[pt->xt];
      memcpy(&xt,pxt0,sizeof(xTetra));
      memcpy(&xt1,pxt0,sizeof(xTetra));
    }
    else {
      memset(&xt,0, sizeof(xTetra));
      memset(&xt1,0, sizeof(xTetra));
    }

    /* tau = sigma^-1 = permutation that sends the reference config (edge 01 split) to current */
    tau[0] = 0; tau[1] = 1; tau[2] = 2; tau[3] = 3;
    taued = &permedge[0][0];
    switch(ie){
    case 1:
      tau[0] = 2; tau[1] = 0; tau[2] = 1; tau[3] = 3;
      taued = &permedge[6][0];
      break;
    case 2:
      tau[0] = 0; tau[1] = 3; tau[2] = 1; tau[3] = 2;
      taued = &permedge[2][0];
      break;
    case 3:
      tau[0] = 1; tau[1] = 2; tau[2] = 0; tau[3] = 3;
      taued = &permedge[4][0];
      break;
    case 4:
      tau[0] = 3; tau[1] = 1; tau[2] = 0; tau[3] = 2;
      taued = &permedge[10][0];
      break;
    case 5:
      tau[0] = 3; tau[1] = 2; tau[2] = 1; tau[3] = 0;
      taued = &permedge[11][0];
     break;
    }

    /* Generic formulation of split of 1 edge */
    pt->v[tau[1]] = pt1->v[tau[0]] = ip;
    if(pt->xt){
      /* Reset edge tag */
      xt.tag[ taued[3]] = xt.ftag[ tau[3]];
      xt.tag[ taued[4]] = xt.ftag[ tau[2]];
      xt1.tag[taued[1]] = xt1.ftag[tau[3]];
      xt1.tag[taued[2]] = xt1.ftag[tau[2]];
      xt.ref[tau[0]]   = 0;  xt.ftag[tau[0]]  = 0;
      xt1.ref[tau[1]]  = 0;  xt1.ftag[tau[1]] = 0;
    }

    pt->flag = pt1->flag = 0;

    isxt = 0 ;
    isxt1 = 0;

    for(i=0;i<4;i++){
      if(xt.ref[i] || xt.ftag[i]) isxt = 1;
      if(xt1.ref[i] || xt1.ftag[i]) isxt1 = 1;
    }

    if(pt->xt){
      if((isxt)&&(!isxt1)){
        pt1->xt = 0;
        pxt0 = &mesh->xtetra[pt->xt];
        memcpy(pxt0,&xt,sizeof(xTetra));
      }
      else if ((!isxt)&&(isxt1)){
        pt1->xt = pt->xt;
        pt->xt = 0;
        pxt0 = &mesh->xtetra[pt1->xt];
        memcpy(pxt0,&xt1,sizeof(xTetra));
      }
      else if (isxt && isxt1){
        mesh->xt++;
        assert(mesh->xt < mesh->xtmax);
        pt1->xt = mesh->xt;
        pxt0 = &mesh->xtetra[pt->xt];
        memcpy(pxt0,&xt,sizeof(xTetra));
        pxt0 = &mesh->xtetra[pt1->xt];
        memcpy(pxt0,&xt1,sizeof(xTetra));
      }
      else{
        pt->xt = 0;
        pt1->xt = 0;
      }
    }

    /* Update of adjacency relations */
    adja = &mesh->adja[4*(iel-1)+1];
    adjan = &mesh->adja[4*(jel-1)+1];

    nei2 = adja[tau[2]];
    nei3 = adja[tau[3]];

    /* Adjacency relations through both splitted faces */
    if(k == 0){
      if((list[1] / 6) == (nei2 / 4)){
        if( MG_SMSGN(newtet[0],newtet[1])){  //new elt of list[0] goes with new elt of list[1]
          adja[tau[2]] = nei2;
          adjan[tau[2]] = 4*fabs(newtet[1])+(nei2 %4);
        }
        else{
          adja[tau[2]] = 4*fabs(newtet[1])+(nei2 %4);
          adjan[tau[2]] = nei2;
        }

        if(open){
          adja[tau[3]] = 0;
          adjan[tau[3]] = 0;
        }

        else{
          assert((list[ilist-1] / 6) == (nei3 / 4));
          if( MG_SMSGN(newtet[0],newtet[ilist-1])){
            adja[tau[3]] = nei3;
            adjan[tau[3]] = 4*fabs(newtet[ilist-1])+(nei3 %4);
          }
          else{
            adja[tau[3]] = 4*fabs(newtet[ilist-1])+(nei3 %4);
            adjan[tau[3]] = nei3;
          }
        }
      }

      else{
        assert((list[1] / 6) == (nei3 / 4));
        if( MG_SMSGN(newtet[0],newtet[1])){
          adja[tau[3]] = nei3;
          adjan[tau[3]] = 4*fabs(newtet[1])+(nei3 %4);
        }
        else{
          adja[tau[3]] = 4*fabs(newtet[1])+(nei3 %4);
          adjan[tau[3]] = nei3;
        }

        if(open){
          adja[tau[2]] = 0;
          adjan[tau[2]] = 0;
        }

        else{
          assert((list[ilist-1]) / 6 == (nei2 / 4));
          if( MG_SMSGN(newtet[0],newtet[ilist-1]) ){
            adja[tau[2]] = nei2;
            adjan[tau[2]] = 4*fabs(newtet[ilist-1])+(nei2 %4);
          }
          else{
            adja[tau[2]] = 4*fabs(newtet[ilist-1])+(nei2 %4);
            adjan[tau[2]] = nei2;
          }
        }
      }
    }

    else if(k==ilist-1){
      if((list[ilist-2] / 6) == (nei2 / 4)){
        if( MG_SMSGN(newtet[ilist-1],newtet[ilist-2]) ){
          adja[tau[2]] = nei2;
          adjan[tau[2]] = 4*fabs(newtet[ilist-2])+(nei2 %4);
        }
        else{
          adja[tau[2]] = 4*fabs(newtet[ilist-2])+(nei2 %4);
          adjan[tau[2]] = nei2;
        }

        if(open){
          adja[tau[3]] = 0;
          adjan[tau[3]] = 0;
        }

        else{
          assert((list[0]) / 6 == (nei3 / 4));
          if( MG_SMSGN(newtet[ilist-1],newtet[0]) ){
            adja[tau[3]] = nei3;
            adjan[tau[3]] = 4*fabs(newtet[0])+(nei3 %4);
          }
          else{
            adja[tau[3]] = 4*fabs(newtet[0])+(nei3 %4);
            adjan[tau[3]] = nei3;
          }
        }
      }

      else{
        assert((list[ilist-2] / 6) == (nei3 / 4));
        if( MG_SMSGN(newtet[ilist-1],newtet[ilist-2]) ){
          adja[tau[3]] = nei3;
          adjan[tau[3]] = 4*fabs(newtet[ilist-2])+(nei3 %4);
        }
        else{
          adja[tau[3]] = 4*fabs(newtet[ilist-2])+(nei3 %4);
          adjan[tau[3]] = nei3;
        }

        if(open){
          adja[tau[2]] = 0;
          adjan[tau[2]] = 0;
        }

        else{
          assert((list[0]) / 6 == (nei2 / 4));
          if( MG_SMSGN(newtet[ilist-1],newtet[0]) ){
            adja[tau[2]] = nei2;
            adjan[tau[2]] = 4*fabs(newtet[0])+(nei2 %4);
          }
          else{
            adja[tau[2]] = 4*fabs(newtet[0])+(nei2 %4);
            adjan[tau[2]] = nei2;
          }
        }
      }
    }

    else{
      if((list[k-1] / 6) == (nei2 / 4)){
        if( MG_SMSGN(newtet[k],newtet[k-1]) ){
          adja[tau[2]] = nei2;
          adjan[tau[2]] = 4*fabs(newtet[k-1])+(nei2 %4);
        }
        else{
          adja[tau[2]] = 4*fabs(newtet[k-1])+(nei2 %4);
          adjan[tau[2]] = nei2;
        }

        assert((list[k+1]) / 6 == (nei3 / 4));
        if( MG_SMSGN(newtet[k],newtet[k+1]) ){
          adja[tau[3]] = nei3;
          adjan[tau[3]] = 4*fabs(newtet[k+1])+(nei3 %4);
        }
        else{
          adja[tau[3]] = 4*fabs(newtet[k+1])+(nei3 %4);
          adjan[tau[3]] = nei3;
        }
      }

      else{
        assert((list[k-1] / 6) == (nei3 / 4));
        if( MG_SMSGN(newtet[k],newtet[k-1]) ){
          adja[tau[3]] = nei3;
          adjan[tau[3]] = 4*fabs(newtet[k-1])+(nei3 %4);
        }
        else{
          adja[tau[3]] = 4*fabs(newtet[k-1])+(nei3 %4);
          adjan[tau[3]] = nei3;
        }

        assert((list[k+1]) / 6 == (nei2 / 4));
        if( MG_SMSGN(newtet[k],newtet[k+1]) ){
          adja[tau[2]] = nei2;
          adjan[tau[2]] = 4*fabs(newtet[k+1])+(nei2 %4);
        }
        else{
          adja[tau[2]] = 4*fabs(newtet[k+1])+(nei2 %4);
          adjan[tau[2]] = nei2;
        }
      }
    }

    /* Internal adjacency relations update */
    mel = adja[tau[0]] / 4;
    voy = adja[tau[0]] % 4;
    adja[tau[0]] = 4*jel + tau[1];
    adjan[tau[0]] = 4*mel + voy;
    adjan[tau[1]] = 4*iel + tau[0];

    if(mel){
      adjan = &mesh->adja[4*(mel -1) +1];
      adjan[voy] = 4*jel + tau[0];
    }
    /* Quality update */
    pt->qual=orcal(mesh,iel);
    pt1->qual=orcal(mesh,jel);
  }

#ifdef DEBUG
  tabtmp[0][3]=lmin;     tabtmp[0][4]=lmax;
  tabtmp[0][5]=DBL_MAX;  tabtmp[0][6]=DBL_MIN;
  for(j=0;j<ilist;j++){
    for(i=0;i<6;i++){
      len=lenedg(mesh,met,
                 mesh->tetra[list[j]/6].v[iare[i][0]],
                 mesh->tetra[list[j]/6].v[iare[i][1]]);
      if(len<lmin) {
        tabtmp[0][1]+=0.5;
      }else if(len>lmax) {
        tabtmp[0][2]+=0.5;
      }
    }
  }
  tabtmp[0][0]+=ilist;
#endif
  free(newtet);
  newtet=NULL;
  return(1);
}

/** Simulate split of two edges that belong to a common face */
int split2sf_sim(pMesh mesh,pSol met,int k,int vx[6]){
  pTetra        pt,pt0;
  double   vold,vnew;
  unsigned char tau[4],*taued,imin;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = orvol(mesh->point,pt->v);

  /* identity is case 48 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &permedge[0][0];
  switch(pt->flag){
  case 24 :
    tau[0] = 0 ; tau[1] = 2 ; tau[2] = 3 ; tau[3] = 1;
    taued = &permedge[1][0];
    break;
  case 40 :
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &permedge[2][0];
    break;
  case 6 :
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &permedge[5][0];
    break;
  case 34 :
    tau[0] = 1 ; tau[1] = 0 ; tau[2] = 3 ; tau[3] = 2;
    taued = &permedge[3][0];
    break;
  case 36 :
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &permedge[4][0];
    break;
  case 20 :
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &permedge[6][0];
    break;
  case 5 :
    tau[0] = 2 ; tau[1] = 1 ; tau[2] = 3 ; tau[3] = 0;
    taued = &permedge[7][0];
    break;
  case 17 :
    tau[0] = 2 ; tau[1] = 3 ; tau[2] = 0 ; tau[3] = 1;
    taued = &permedge[8][0];
    break;
  case 9 :
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    taued = &permedge[9][0];
    break;
  case 3 :
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    taued = &permedge[11][0];
    break;
  case 10 :
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &permedge[10][0];
    break;
  }

  /* Test orientation of the three tets to be created */
  imin = (pt->v[tau[1]] < pt->v[tau[2]]) ? tau[1] : tau[2] ;

  memcpy(pt0,pt,sizeof(Tetra));
  pt0->v[tau[1]] = vx[taued[4]];
  pt0->v[tau[2]] = vx[taued[5]];
  vnew = orvol(mesh->point,pt0->v);
  if ( vnew < EPSD2 )  return(0);
  else if ( vold > NULKAL && vnew < NULKAL )  return(0);

  if ( imin == tau[1] ) {
    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[tau[2]] = vx[taued[5]];
    pt0->v[tau[3]] = vx[taued[4]];
    vnew = orvol(mesh->point,pt0->v);
    if ( vnew < EPSD2 )  return(0);
    else if ( vold > NULKAL && vnew < NULKAL )  return(0);

    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[tau[3]] = vx[taued[5]];
    vnew = orvol(mesh->point,pt0->v);
    if ( vnew < EPSD2 )  return(0);
    else if ( vold > NULKAL && vnew < NULKAL )  return(0);
  }
  else {
    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[tau[3]] = vx[taued[4]];
    vnew = orvol(mesh->point,pt0->v);
    if ( vnew < EPSD2 )  return(0);
    else if ( vold > NULKAL && vnew < NULKAL )  return(0);

    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[tau[1]] = vx[taued[4]];
    pt0->v[tau[3]] = vx[taued[5]];
    vnew = orvol(mesh->point,pt0->v);
    if ( vnew < EPSD2 )  return(0);
    else if ( vold > NULKAL && vnew < NULKAL )  return(0);
  }
  return(1);
}

/** Split of two edges that belong to a common face : 1 tetra becomes 3 */
void split2sf(pMesh mesh,pSol met,int k,int vx[6]){
  pTetra        pt[3];
  xTetra        xt[3];
  pxTetra       pxt0;
  int           iel,i,ref4,ref5;
  int           newtet[3];
  char          flg,imin,firstxt,isxt[3],tag4,tag5;
  unsigned char tau[4],*taued;

  pt[0] = &mesh->tetra[k];
#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        lmintmp=lentmp;
      }else if(lentmp>lmaxtmp) {
        lmaxtmp=lentmp;
      }
    }
#endif
  flg   = pt[0]->flag;
  pt[0]->flag = 0;
  newtet[0]=k;

  iel = newElt(mesh);
  assert(iel);
  pt[1] = &mesh->tetra[iel];
  memcpy(pt[1],pt[0],sizeof(Tetra));
  newtet[1]=iel;

  iel = newElt(mesh);
  assert(iel);
  pt[2] = &mesh->tetra[iel];
  memcpy(pt[2],pt[0],sizeof(Tetra));
  newtet[2]=iel;

  if ( pt[0]->xt ) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    memcpy(&xt[0],pxt0,sizeof(xTetra));
    memcpy(&xt[1],pxt0,sizeof(xTetra));
    memcpy(&xt[2],pxt0,sizeof(xTetra));
  }
  else {
    pxt0 = 0;
    memset(&xt[0],0,sizeof(xTetra));
    memset(&xt[1],0,sizeof(xTetra));
    memset(&xt[2],0,sizeof(xTetra));
  }
  /* identity is case 48 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &permedge[0][0];
  switch(flg){
  case 24 :
    tau[0] = 0 ; tau[1] = 2 ; tau[2] = 3 ; tau[3] = 1;
    taued = &permedge[1][0];
    break;
  case 40 :
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &permedge[2][0];
    break;
  case 6 :
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &permedge[5][0];
    break;
  case 34 :
    tau[0] = 1 ; tau[1] = 0 ; tau[2] = 3 ; tau[3] = 2;
    taued = &permedge[3][0];
    break;
  case 36 :
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &permedge[4][0];
    break;
  case 20 :
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &permedge[6][0];
    break;
  case 5 :
    tau[0] = 2 ; tau[1] = 1 ; tau[2] = 3 ; tau[3] = 0;
    taued = &permedge[7][0];
    break;
  case 17 :
    tau[0] = 2 ; tau[1] = 3 ; tau[2] = 0 ; tau[3] = 1;
    taued = &permedge[8][0];
    break;
  case 9 :
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    taued = &permedge[9][0];
    break;
  case 3 :
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    taued = &permedge[11][0];
    break;
  case 10 :
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &permedge[10][0];
    break;
  }

  hPop(&mesh->htab,pt[0]->v[tau[1]],pt[0]->v[tau[3]],&ref4,&tag4);
  if ( ref4 || tag4 ) {
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[4]],ref4,tag4);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[4]],ref4,tag4);
  }
  hPop(&mesh->htab,pt[0]->v[tau[2]],pt[0]->v[tau[3]],&ref5,&tag5);
  if ( ref5 || tag5 ) {
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[5]],ref5,tag5);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[5]],ref5,tag5);
  }

  /* Generic formulation for the split of 2 edges belonging to a common face */
  imin = (pt[0]->v[tau[1]] < pt[0]->v[tau[2]]) ? tau[1] : tau[2] ;
  pt[0]->v[tau[1]]  = vx[taued[4]] ;  pt[0]->v[tau[2]] = vx[taued[5]];
  xt[0].tag[taued[0]] = xt[0].ftag[tau[2]];
  xt[0].tag[taued[1]] = xt[0].ftag[tau[1]];
  xt[0].tag[taued[3]] = xt[0].ftag[tau[0]];
  xt[0].ref[tau[3]] = 0;  xt[0].ftag[tau[3]] = 0;

  if ( imin == tau[1] ) {
    pt[1]->v[tau[2]] = vx[taued[5]];  pt[1]->v[tau[3]] = vx[taued[4]];
    pt[2]->v[tau[3]] = vx[taued[5]];

    xt[1].tag[taued[1]] = xt[1].ftag[tau[1]];
    xt[1].tag[taued[2]] = xt[1].ftag[tau[2]];
    xt[1].tag[taued[3]] = xt[1].ftag[tau[0]];
    xt[1].tag[taued[5]] = xt[1].ftag[tau[0]];
    xt[1].ref[  tau[1]] = 0;  xt[1].ref[  tau[3]] = 0;
    xt[1].ftag[ tau[1]] = 0;  xt[1].ftag[ tau[3]] = 0;

    xt[2].tag[taued[2]] = xt[2].ftag[tau[1]];
    xt[2].tag[taued[4]] = xt[2].ftag[tau[0]];
    xt[2].ref[ tau[2]]  = 0;  xt[2].ftag[ tau[2]] = 0;
  }
  else {
    pt[1]->v[tau[3]] = vx[taued[4]];
    pt[2]->v[tau[1]] = vx[taued[4]];  pt[2]->v[tau[3]] = vx[taued[5]];

    xt[1].tag[taued[2]] = xt[1].ftag[tau[2]];
    xt[1].tag[taued[5]] = xt[1].ftag[tau[0]];
    xt[1].ref[  tau[1]] = 0;  xt[1].ftag[ tau[1]] = 0;

    xt[2].tag[taued[0]] = xt[2].ftag[tau[2]];
    xt[2].tag[taued[2]] = xt[2].ftag[tau[1]];
    xt[2].tag[taued[3]] = xt[2].ftag[tau[0]];
    xt[2].tag[taued[4]] = xt[2].ftag[tau[0]];
    xt[2].ref[  tau[2]] = 0;  xt[2].ref[  tau[3]] = 0;
    xt[2].ftag[ tau[2]] = 0;  xt[2].ftag[ tau[3]] = 0;
  }

  /* Assignation of the xt fields to the appropriate tets */
  isxt[0] = isxt[1] = isxt[2] = 0;
  for (i=0; i<4; i++) {
    if ( xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
    if ( xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
    if ( xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
  }

  if ( pt[0]->xt ) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(xTetra));
      pt[1]->xt = pt[2]->xt = 0;
      for (i=1; i<3; i++) {
        if ( isxt[i] ) {
          assert(mesh->xt < mesh->xtmax-1);
          pt[i]->xt = ++mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&(xt[i]),sizeof(xTetra));
        }
      }
    }
    else {
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = 0;
      for (i=1; i<3; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[pt[i]->xt];
            memcpy(pxt0,&(xt[i]),sizeof(xTetra));
          }
          else {
            assert(mesh->xt < mesh->xtmax-1);
            pt[i]->xt = ++mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
        }
      }
      pt[0]->xt = 0;
    }
  }
  /* Quality update */
  pt[0]->qual=orcal(mesh,newtet[0]);
  pt[1]->qual=orcal(mesh,newtet[1]);
  pt[2]->qual=orcal(mesh,newtet[2]);

#ifdef DEBUG
  tabtmp[1][3]=lmintmp;  tabtmp[1][4]=lmaxtmp;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[1][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[1][2]++;
      }
      lentmp=lenedg(mesh,met,pt[1]->v[iare[i][0]],pt[1]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[1][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[1][2]++;
      }
      lentmp=lenedg(mesh,met,pt[2]->v[iare[i][0]],pt[2]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[1][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[1][2]++;
      }
    }
#endif
}

/** Split of two OPPOSITE edges */
void split2(pMesh mesh,pSol met,int k,int vx[6]) {
  pTetra   pt[4];
  xTetra   xt[4];
  pxTetra  pxt0;
  int      i,iel,ref0,ref5;
  int      newtet[4];
  char     flg,firstxt,isxt[4],tag0,tag5;
  unsigned char tau[4],*taued;

  pt[0] = &mesh->tetra[k];
#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        lmintmp=lentmp;
      }else if(lentmp>lmaxtmp) {
        lmaxtmp=lentmp;
      }
    }
  iatmp=0;ibtmp=0;
  for(i=0;i<6;i++){
    if(iatmp==0 && vx[i]>0){
      iatmp=i;
    }else if(iatmp>0 && vx[i]>0){
      ibtmp=i;
    }
  }
#endif

  flg   = pt[0]->flag;
  pt[0]->flag = 0;
  newtet[0]=k;

  iel = newElt(mesh);
  assert(iel);
  pt[1] = &mesh->tetra[iel];
  memcpy(pt[1],pt[0],sizeof(Tetra));
  newtet[1]=iel;

  iel = newElt(mesh);
  assert(iel);
  pt[2] = &mesh->tetra[iel];
  memcpy(pt[2],pt[0],sizeof(Tetra));
  newtet[2]=iel;

  iel = newElt(mesh);
  assert(iel);
  pt[3] = &mesh->tetra[iel];
  memcpy(pt[3],pt[0],sizeof(Tetra));
  newtet[3]=iel;

  pxt0 = 0;
  if ( pt[0]->xt) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    memcpy(&xt[0],pxt0,sizeof(xTetra));
    memcpy(&xt[1],pxt0,sizeof(xTetra));
    memcpy(&xt[2],pxt0,sizeof(xTetra));
    memcpy(&xt[3],pxt0,sizeof(xTetra));
  }
  else {
    memset(&xt[0],0,sizeof(xTetra));
    memset(&xt[1],0,sizeof(xTetra));
    memset(&xt[2],0,sizeof(xTetra));
    memset(&xt[3],0,sizeof(xTetra));
  }
  /* identity : case 33 */
  tau[0] = 0;  tau[1] = 1;  tau[2] = 2;  tau[3] = 3;
  taued = &permedge[0][0];
  switch(flg){
  case 18:
    tau[0] = 3;  tau[1] = 1;  tau[2] = 0;  tau[3] = 2;
    taued = &permedge[10][0];
    break;
  case 12:
    tau[0] = 0;  tau[1] = 3;  tau[2] = 1;  tau[3] = 2;
    taued = &permedge[2][0];
    break;
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[1]],&ref0,&tag0);
  if( ref0 || tag0 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[0]],ref0,tag0);
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[0]],ref0,tag0);
  }

  hPop(&mesh->htab,pt[0]->v[tau[2]],pt[0]->v[tau[3]],&ref5,&tag5);
  if( ref5 || tag5 ){
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[5]],ref5,tag5);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[5]],ref5,tag5);
  }

  /* Generic formulation for the split of 2 opposite edges */
  pt[0]->v[tau[1]] = vx[taued[0]];  pt[0]->v[tau[2]] = vx[taued[5]];
  pt[1]->v[tau[1]] = vx[taued[0]];  pt[1]->v[tau[3]] = vx[taued[5]];
  pt[2]->v[tau[0]] = vx[taued[0]];  pt[2]->v[tau[2]] = vx[taued[5]];
  pt[3]->v[tau[0]] = vx[taued[0]];  pt[3]->v[tau[3]] = vx[taued[5]];

  xt[0].ref[tau[0]]  = 0;  xt[0].ref[tau[3]]  = 0;
  xt[0].ftag[tau[0]] = 0;  xt[0].ftag[tau[3]] = 0;

  xt[1].ref[tau[0]]  = 0;  xt[1].ref[tau[2]]  = 0;
  xt[1].ftag[tau[0]] = 0;  xt[1].ftag[tau[2]] = 0;

  xt[2].ref[tau[1]]  = 0;  xt[2].ref[tau[3]]  = 0;
  xt[2].ftag[tau[1]] = 0;  xt[2].ftag[tau[3]] = 0;

  xt[3].ref[tau[1]]  = 0;  xt[3].ref[tau[2]]  = 0;
  xt[3].ftag[tau[1]] = 0;  xt[3].ftag[tau[2]] = 0;

  /* Assignation of the xt fields to the appropriate tets */
  memset(isxt,0,4*sizeof(char));
  for (i=0; i<4; i++) {
    if ( xt[0].ref[i] || xt[0].ftag[i] )  isxt[0] = 1;
    if ( xt[1].ref[i] || xt[1].ftag[i] )  isxt[1] = 1;
    if ( xt[2].ref[i] || xt[2].ftag[i] )  isxt[2] = 1;
    if ( xt[3].ref[i] || xt[3].ftag[i] )  isxt[3] = 1;
  }
  if ( pt[0]->xt) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(xTetra));
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          assert(mesh->xt < mesh->xtmax-1);
          pt[i]->xt = ++mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(xTetra));
        }
        else {
          pt[i]->xt = 0;
        }
      }
    }
    else {
      firstxt = 1;
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[pt[i]->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
          else {
            assert(mesh->xt < mesh->xtmax-1);
            pt[i]->xt = ++mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
        }
        else {
          pt[i]->xt = 0;
        }
      }
      pt[0]->xt = 0;
    }
  }
  /* Quality update */
    pt[0]->qual=orcal(mesh,newtet[0]);
    pt[1]->qual=orcal(mesh,newtet[1]);
    pt[2]->qual=orcal(mesh,newtet[2]);
    pt[3]->qual=orcal(mesh,newtet[3]);

#if DEBUG
  tabtmp[2][3]=lmintmp;  tabtmp[2][4]=lmaxtmp;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        if((pt[0]->v[iare[i][0]]==iatmp && pt[0]->v[iare[i][1]]==ibtmp)||
           (pt[0]->v[iare[i][0]]==ibtmp && pt[0]->v[iare[i][1]]==iatmp)) {
          tabtmp[2][1]= tabtmp[2][1]+0.25;
        }else{
          tabtmp[2][1]= tabtmp[2][1]+0.5;
        }
      }else if(lentmp>lmaxtmp) {
        if((pt[0]->v[iare[i][0]]==iatmp && pt[0]->v[iare[i][1]]==ibtmp)||
           (pt[0]->v[iare[i][0]]==ibtmp && pt[0]->v[iare[i][1]]==iatmp)) {
          tabtmp[2][2]= tabtmp[2][2]+0.25;
        }else{
          tabtmp[2][2]= tabtmp[2][2]+0.5;
        }
      }
      lentmp=lenedg(mesh,met,pt[1]->v[iare[i][0]],pt[1]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        if((pt[1]->v[iare[i][0]]==iatmp && pt[1]->v[iare[i][1]]==ibtmp)||
           (pt[1]->v[iare[i][0]]==ibtmp && pt[1]->v[iare[i][1]]==iatmp)) {
          tabtmp[2][1]= tabtmp[2][1]+0.25;
        }else{
          tabtmp[2][1]= tabtmp[2][1]+0.5;
        }
      }else if(lentmp>lmaxtmp) {
        if((pt[1]->v[iare[i][0]]==iatmp && pt[1]->v[iare[i][1]]==ibtmp)||
           (pt[1]->v[iare[i][0]]==ibtmp && pt[1]->v[iare[i][1]]==iatmp)) {
          tabtmp[2][2]= tabtmp[2][2]+0.25;
        }else{
          tabtmp[2][2]= tabtmp[2][2]+0.5;
        }

      }
      lentmp=lenedg(mesh,met,pt[2]->v[iare[i][0]],pt[2]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        if((pt[2]->v[iare[i][0]]==iatmp && pt[2]->v[iare[i][1]]==ibtmp)||
           (pt[2]->v[iare[i][0]]==ibtmp && pt[2]->v[iare[i][1]]==iatmp)) {
          tabtmp[2][1]= tabtmp[2][1]+0.25;
        }else{
          tabtmp[2][1]= tabtmp[2][1]+0.5;
        }
      }else if(lentmp>lmaxtmp) {
        if((pt[2]->v[iare[i][0]]==iatmp && pt[2]->v[iare[i][1]]==ibtmp)||
           (pt[2]->v[iare[i][0]]==ibtmp && pt[2]->v[iare[i][1]]==iatmp)) {
          tabtmp[2][2]= tabtmp[2][2]+0.25;
        }else{
          tabtmp[2][2]= tabtmp[2][2]+0.5;
        }
      }
      lentmp=lenedg(mesh,met,pt[3]->v[iare[i][0]],pt[3]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        if((pt[3]->v[iare[i][0]]==iatmp && pt[3]->v[iare[i][1]]==ibtmp)||
           (pt[3]->v[iare[i][0]]==ibtmp && pt[3]->v[iare[i][1]]==iatmp)) {
          tabtmp[2][1]= tabtmp[2][1]+0.25;
        }else{
          tabtmp[2][1]= tabtmp[2][1]+0.5;
        }

        tabtmp[2][1]++;
      }else if(lentmp>lmaxtmp) {
        if((pt[3]->v[iare[i][0]]==iatmp && pt[3]->v[iare[i][1]]==ibtmp)||
           (pt[3]->v[iare[i][0]]==ibtmp && pt[3]->v[iare[i][1]]==iatmp)) {
          tabtmp[2][2]= tabtmp[2][2]+0.25;
        }else{
          tabtmp[2][2]= tabtmp[2][2]+0.5;
        }
      }
    }
#endif
}

/** Simulate split of 1 face (3 edges) */
int split3_sim(pMesh mesh,pSol met,int k,int vx[6]) {
  pTetra    pt,pt0;
  double    vold,vnew;
  unsigned char tau[4],*taued;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = orvol(mesh->point,pt->v);

  /* identity is case 11 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &permedge[0][0];
  switch(pt->flag) {
  case 21:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &permedge[2][0];
    break;
  case 38:
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    taued = &permedge[9][0];
    break;
  case 56:
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &permedge[5][0];
    break;
  }

  /* Check orientation of the 4 newly created tets */
  memcpy(pt0,pt,sizeof(Tetra));
  pt0->v[tau[1]] = vx[taued[0]];
  pt0->v[tau[2]] = vx[taued[1]];
  vnew = orvol(mesh->point,pt0->v);
  if ( vnew < EPSD2 )  return(0);
  else if ( vold > NULKAL && vnew < NULKAL )  return(0);

  memcpy(pt0,pt,sizeof(Tetra));
  pt0->v[tau[0]] = vx[taued[0]];
  pt0->v[tau[2]] = vx[taued[3]];
  vnew = orvol(mesh->point,pt0->v);
  if ( vnew < EPSD2 )  return(0);
  else if ( vold > NULKAL && vnew < NULKAL )  return(0);

  memcpy(pt0,pt,sizeof(Tetra));
  pt0->v[tau[0]] = vx[taued[1]];
  pt0->v[tau[1]] = vx[taued[3]];
  vnew = orvol(mesh->point,pt0->v);
  if ( vnew < EPSD2 )  return(0);
  else if ( vold > NULKAL && vnew < NULKAL )  return(0);

  memcpy(pt0,pt,sizeof(Tetra));
  pt0->v[tau[0]] = vx[taued[0]];
  pt0->v[tau[1]] = vx[taued[3]];
  pt0->v[tau[2]] = vx[taued[1]];
  vnew = orvol(mesh->point,pt0->v);
  if ( vnew < EPSD2 )  return(0);
  else if ( vold > NULKAL && vnew < NULKAL )  return(0);

  return(1);
}

/** 1 face (3 edges) subdivided */
void split3(pMesh mesh,pSol met,int k,int vx[6]) {
  pTetra    pt[4];
  xTetra    xt[4];
  pxTetra   pxt0;
  int       iel,i,ref0,ref1,ref3;
  int       newtet[4];
  char      flg,firstxt,isxt[4],tag0,tag1,tag3;
  unsigned char tau[4],*taued;

  pt[0] = &mesh->tetra[k];
  flg   = pt[0]->flag;
  pt[0]->flag  = 0;
  newtet[0]=k;
#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        //       tabtmp[0][1]++;
        lmintmp=lentmp;
      }else if(lentmp>lmaxtmp) {
        //tabtmp[0][2]++;
        lmaxtmp=lentmp;
      }
    }
#endif
  /* create 3 new tetras */
  iel = newElt(mesh);
  assert(iel);
  pt[1] = &mesh->tetra[iel];
  pt[1] = memcpy(pt[1],pt[0],sizeof(Tetra));
  newtet[1]=iel;

  iel = newElt(mesh);
  assert(iel);
  pt[2] = &mesh->tetra[iel];
  pt[2] = memcpy(pt[2],pt[0],sizeof(Tetra));
  newtet[2]=iel;

  iel = newElt(mesh);
  assert(iel);
  pt[3] = &mesh->tetra[iel];
  pt[3] = memcpy(pt[3],pt[0],sizeof(Tetra));
  newtet[3]=iel;

  pxt0 = 0;
  if ( pt[0]->xt ) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    memcpy(&xt[0],pxt0, sizeof(xTetra));
    memcpy(&xt[1],pxt0, sizeof(xTetra));
    memcpy(&xt[2],pxt0, sizeof(xTetra));
    memcpy(&xt[3],pxt0, sizeof(xTetra));
  }
  else {
    memset(&xt[0],0, sizeof(xTetra));
    memset(&xt[1],0, sizeof(xTetra));
    memset(&xt[2],0, sizeof(xTetra));
    memset(&xt[3],0, sizeof(xTetra));
  }

  /* update vertices, case 11 is default */
  tau[0] = 0; tau[1] = 1; tau[2] = 2; tau[3] = 3;
  taued = &permedge[0][0];
  switch(flg) {
  case 21:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &permedge[2][0];
    break;
  case 38:
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    taued = &permedge[9][0];
    break;
  case 56:
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &permedge[5][0];
    break;
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[1]],&ref0,&tag0);
  if ( ref0 || tag0 ) {
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[0]],ref0,tag0);
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[0]],ref0,tag0);
  }
  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[2]],&ref1,&tag1);
  if ( ref1 || tag1 ) {
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[1]],ref1,tag1);
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[1]],ref1,tag1);
  }
  hPop(&mesh->htab,pt[0]->v[tau[1]],pt[0]->v[tau[2]],&ref3,&tag3);
  if ( ref3 || tag3 ) {
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[3]],ref3,tag3);
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[3]],ref3,tag3);
  }

  /* Generic formulation of split of 3 edges */
  pt[0]->v[tau[1]] = vx[taued[0]];  pt[0]->v[tau[2]] = vx[taued[1]];
  pt[1]->v[tau[0]] = vx[taued[0]];  pt[1]->v[tau[2]] = vx[taued[3]];
  pt[2]->v[tau[0]] = vx[taued[1]];  pt[2]->v[tau[1]] = vx[taued[3]];
  pt[3]->v[tau[0]] = vx[taued[0]];  pt[3]->v[tau[1]] = vx[taued[3]];  pt[3]->v[tau[2]] = vx[taued[1]];

  xt[0].tag[taued[3]] = xt[0].ftag[tau[3]];
  xt[0].tag[taued[4]] = xt[0].ftag[tau[2]];
  xt[0].tag[taued[5]] = xt[0].ftag[tau[1]];
  xt[0].ref[  tau[0]] = 0;  xt[0].ftag[ tau[0]] = 0;

  xt[1].tag[taued[1]] = xt[1].ftag[tau[3]];
  xt[1].tag[taued[2]] = xt[1].ftag[tau[2]];
  xt[1].tag[taued[5]] = xt[1].ftag[tau[0]];
  xt[1].ref[  tau[1]] = 0;  xt[1].ftag[ tau[1]] = 0;

  xt[2].tag[taued[0]] = xt[2].ftag[tau[3]];
  xt[2].tag[taued[2]] = xt[2].ftag[tau[1]];
  xt[2].tag[taued[4]] = xt[2].ftag[tau[0]];
  xt[2].ref[  tau[2]] = 0;  xt[2].ftag[ tau[2]] = 0;

  xt[3].tag[taued[0]] = xt[3].ftag[tau[3]];
  xt[3].tag[taued[1]] = xt[3].ftag[tau[3]];
  xt[3].tag[taued[2]] = xt[3].ftag[tau[2]];
  xt[3].tag[taued[3]] = xt[3].ftag[tau[3]];
  xt[3].tag[taued[4]] = xt[3].ftag[tau[0]];
  xt[3].tag[taued[5]] = xt[3].ftag[tau[1]];
  xt[3].ref[  tau[0]] = 0;  xt[3].ref[  tau[1]] = 0;  xt[3].ref[  tau[2]] = 0;
  xt[3].ftag[ tau[0]] = 0;  xt[3].ftag[ tau[1]] = 0;  xt[3].ftag[ tau[2]] = 0;

  /* Assignation of the xt fields to the appropriate tets */
  memset(isxt,0,4*sizeof(char));
  for (i=0; i<4; i++) {
    if ( xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
    if ( xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
    if ( xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
    if ( xt[3].ref[i] || xt[3].ftag[i] ) isxt[3] = 1;
  }

  if ( pt[0]->xt ) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(xTetra));
      pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          assert(mesh->xt < mesh->xtmax-1);
          pt[i]->xt = ++mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(xTetra));
        }
      }
    }
    else {
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[(pt[i])->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
          else {
            assert(mesh->xt < mesh->xtmax-1);
            pt[i]->xt = ++mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&(xt[i]),sizeof(xTetra));
          }
        }
      }
      pt[0]->xt = 0;
    }
  }
  /* Quality update */
  pt[0]->qual=orcal(mesh,newtet[0]);
  pt[1]->qual=orcal(mesh,newtet[1]);
  pt[2]->qual=orcal(mesh,newtet[2]);
  pt[3]->qual=orcal(mesh,newtet[3]);

#ifdef DEBUG
  tabtmp[3][3]=lmintmp;  tabtmp[3][4]=lmaxtmp;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[3][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[3][2]++;
      }
      lentmp=lenedg(mesh,met,pt[1]->v[iare[i][0]],pt[1]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[3][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[3][2]++;
      }
      lentmp=lenedg(mesh,met,pt[2]->v[iare[i][0]],pt[2]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[3][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[3][2]++;
      }
      lentmp=lenedg(mesh,met,pt[3]->v[iare[i][0]],pt[3]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[3][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[3][2]++;
      }
    }
#endif
}

/** Split 3 edge in cone configuration */
void split3cone(pMesh mesh,pSol met,int k,int vx[6]) {
  pTetra    pt[4];
  xTetra    xt[4];
  pxTetra   pxt0;
  int       iel,i,ref0,ref1,ref2;
  int       newtet[4];
  char      flg,firstxt,isxt[4],ia,ib,tag0,tag1,tag2;//ic;
  unsigned char tau[4],*taued;

  pt[0]  = &mesh->tetra[k];
  flg = pt[0]->flag;
  pt[0]->flag  = 0;
  newtet[0]=k;

#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        lmintmp=lentmp;
      }else if(lentmp>lmaxtmp) {
        lmaxtmp=lentmp;
      }
    }
#endif
  /* create 3 new tetras */
  iel = newElt(mesh);
  assert(iel);
  pt[1] = &mesh->tetra[iel];
  memcpy(pt[1],pt[0],sizeof(Tetra));
  newtet[1]=iel;

  iel = newElt(mesh);
  assert(iel);
  pt[2] = &mesh->tetra[iel];
  memcpy(pt[2],pt[0],sizeof(Tetra));
  newtet[2]=iel;

  iel = newElt(mesh);
  assert(iel);
  pt[3] = &mesh->tetra[iel];
  memcpy(pt[3],pt[0],sizeof(Tetra));
  newtet[3]=iel;

  if(pt[0]->xt){
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    memcpy(&xt[0],pxt0, sizeof(xTetra));
    memcpy(&xt[1],pxt0, sizeof(xTetra));
    memcpy(&xt[2],pxt0, sizeof(xTetra));
    memcpy(&xt[3],pxt0, sizeof(xTetra));
  }
  else{
    pxt0 = 0;
    memset(&xt[0],0, sizeof(xTetra));
    memset(&xt[1],0, sizeof(xTetra));
    memset(&xt[2],0, sizeof(xTetra));
    memset(&xt[3],0, sizeof(xTetra));
  }

  /* Set permutation of vertices : reference configuration is 7 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &permedge[0][0];

  switch(flg) {
  case 25:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &permedge[4][0];
    break;

  case 42:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &permedge[6][0];
    break;

  case 52:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &permedge[10][0];
    break;
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[1]],&ref0,&tag0);
  if( ref0 || tag0 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[0]],ref0,tag0);
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[0]],ref0,tag0);
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[2]],&ref1,&tag1);
  if( ref1 || tag1 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[1]],ref1,tag1);
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[1]],ref1,tag1);
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[3]],&ref2,&tag2);
  if( ref2 || tag2 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[2]],ref2,tag2);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[2]],ref2,tag2);
  }

  /* Generic formulation of split of 3 edges in cone configuration (edges 0,1,2 splitted) */
  /* Fill ia,ib,ic so that pt->v[ia] < pt->v[ib] < pt->v[ic] */
  if((pt[0])->v[tau[1]] < (pt[0])->v[tau[2]]){
    ia = tau[1];
    ib = tau[2];
  }
  else{
    ia = tau[2];
    ib = tau[1];
  }

  if((pt[0])->v[tau[3]] < (pt[0])->v[ia]){
    //ic = ib;
    ib = ia;
    ia = tau[3];
  }
  else{
    if((pt[0])->v[tau[3]] < (pt[0])->v[ib]){
      //ic = ib;
      ib = tau[3];
    }
    else{
      //ic = tau[3];
    }
  }

  pt[0]->v[tau[1]] = vx[taued[0]] ; pt[0]->v[tau[2]] = vx[taued[1]] ; pt[0]->v[tau[3]] = vx[taued[2]];
  xt[0].tag[taued[3]] = xt[0].ftag[tau[3]];
  xt[0].tag[taued[4]] = xt[0].ftag[tau[2]];
  xt[0].tag[taued[5]] = xt[0].ftag[tau[1]];
  xt[0].ref[  tau[0]] = 0;
  xt[0].ftag[ tau[0]] = 0;

  if(ia == tau[3]){
    pt[1]->v[tau[0]] = vx[taued[2]] ; pt[1]->v[tau[1]] = vx[taued[0]] ; pt[1]->v[tau[2]] = vx[taued[1]];
    xt[1].tag[taued[0]] = xt[1].ftag[tau[2]];
    xt[1].tag[taued[1]] = xt[1].ftag[tau[1]];
    xt[1].tag[taued[3]] = xt[1].ftag[tau[3]];
    xt[1].tag[taued[4]] = xt[1].ftag[tau[2]];
    xt[1].tag[taued[5]] = xt[1].ftag[tau[1]];
    xt[1].ref[  tau[0]] = 0;  xt[1].ref[  tau[3]] = 0;
    xt[1].ftag[ tau[0]] = 0;  xt[1].ftag[ tau[3]] = 0;

    if(ib == tau[1]){
      pt[2]->v[tau[0]] = vx[taued[0]] ; pt[2]->v[tau[2]] = vx[taued[1]] ;
      xt[2].tag[taued[1]] = xt[2].ftag[tau[3]];
      xt[2].tag[taued[2]] = xt[2].ftag[tau[2]];
      xt[2].tag[taued[3]] = xt[2].ftag[tau[3]];
      xt[2].tag[taued[5]] = xt[2].ftag[tau[1]];
      xt[2].ref[  tau[0]] = 0;  xt[2].ref[  tau[1]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[1]] = 0;

      pt[3]->v[tau[0]] = vx[taued[1]] ;
      xt[3].tag[taued[0]] = xt[3].ftag[tau[3]];
      xt[3].tag[taued[2]] = xt[3].ftag[tau[1]];
      xt[3].ref[  tau[2]] = 0;
      xt[3].ftag[ tau[2]] = 0;
    }
    else{
      assert(ib == tau[2]);

      pt[2]->v[tau[0]] = vx[taued[1]] ; pt[2]->v[tau[1]] = vx[taued[0]] ;
      xt[2].tag[taued[0]] = xt[2].ftag[tau[3]];
      xt[2].tag[taued[2]] = xt[2].ftag[tau[1]];
      xt[2].tag[taued[3]] = xt[2].ftag[tau[3]];
      xt[2].tag[taued[4]] = xt[2].ftag[tau[2]];
      xt[2].ref[  tau[0]] = 0;  xt[2].ref[  tau[2]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[2]] = 0;

      pt[3]->v[tau[0]] = vx[taued[0]] ;
      xt[3].tag[taued[1]] = xt[3].ftag[tau[3]];
      xt[3].tag[taued[2]] = xt[3].ftag[tau[2]];
      xt[3].ref[  tau[1]] = 0;
      xt[3].ftag[ tau[1]] = 0;
    }
  }

  else if (ia == tau[2]){
    pt[1]->v[tau[0]] = vx[taued[1]] ; pt[1]->v[tau[1]] = vx[taued[0]] ; pt[1]->v[tau[3]] = vx[taued[2]];
    xt[1].tag[taued[0]] = xt[1].ftag[tau[3]];
    xt[1].tag[taued[2]] = xt[1].ftag[tau[1]];
    xt[1].tag[taued[3]] = xt[1].ftag[tau[3]];
    xt[1].tag[taued[4]] = xt[1].ftag[tau[2]];
    xt[1].tag[taued[5]] = xt[1].ftag[tau[1]];
    xt[1].ref[  tau[0]] = 0;  xt[1].ref[  tau[2]] = 0;
    xt[1].ftag[ tau[0]] = 0;  xt[1].ftag[ tau[2]] = 0;

    if(ib == tau[3]){
      pt[2]->v[tau[0]] = vx[taued[2]] ; pt[2]->v[tau[1]] = vx[taued[0]] ;
      xt[2].tag[taued[0]] = xt[2].ftag[tau[2]];
      xt[2].tag[taued[1]] = xt[2].ftag[tau[1]];
      xt[2].tag[taued[3]] = xt[2].ftag[tau[3]];
      xt[2].tag[taued[4]] = xt[2].ftag[tau[2]];
      xt[2].ref[  tau[0]] = 0;  xt[2].ref[  tau[3]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[3]] = 0;

      pt[3]->v[tau[0]] = vx[taued[0]] ;
      xt[3].tag[taued[1]] = xt[3].ftag[tau[3]];
      xt[3].tag[taued[2]] = xt[3].ftag[tau[2]];
      xt[3].ref[  tau[1]] = 0;
      xt[3].ftag[ tau[1]] = 0;
    }
    else{
      assert(ib == tau[1]);

      pt[2]->v[tau[0]] = vx[taued[0]] ; pt[2]->v[tau[3]] = vx[taued[2]] ;
      xt[2].tag[taued[1]] = xt[2].ftag[tau[3]];
      xt[2].tag[taued[2]] = xt[2].ftag[tau[2]];
      xt[2].tag[taued[4]] = xt[2].ftag[tau[2]];
      xt[2].tag[taued[5]] = xt[2].ftag[tau[1]];
      xt[2].ref[  tau[0]] = 0;  xt[2].ref[  tau[1]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[1]] = 0;

      pt[3]->v[tau[0]] = vx[taued[2]] ;
      xt[3].tag[taued[0]] = xt[3].ftag[tau[2]];
      xt[3].tag[taued[1]] = xt[3].ftag[tau[1]];
      xt[3].ref[  tau[3]] = 0;
      xt[3].ftag[ tau[3]] = 0;
    }
  }
  else{
    assert(ia == tau[1]);

    pt[1]->v[tau[0]] = vx[taued[0]] ; pt[1]->v[tau[2]] = vx[taued[1]] ; pt[1]->v[tau[3]] = vx[taued[2]];
    xt[1].tag[taued[1]] = xt[1].ftag[tau[3]];
    xt[1].tag[taued[2]] = xt[1].ftag[tau[2]];
    xt[1].tag[taued[3]] = xt[1].ftag[tau[3]];
    xt[1].tag[taued[4]] = xt[1].ftag[tau[2]];
    xt[1].tag[taued[5]] = xt[1].ftag[tau[1]];
    xt[1].ref[  tau[0]] = 0;  xt[1].ref[  tau[1]] = 0;
    xt[1].ftag[ tau[0]] = 0;  xt[1].ftag[ tau[1]] = 0;

    if(ib == tau[2]){
      pt[2]->v[tau[0]] = vx[taued[1]] ; pt[2]->v[tau[3]] = vx[taued[2]] ;
      xt[2].tag[taued[0]] = xt[2].ftag[tau[3]];
      xt[2].tag[taued[2]] = xt[2].ftag[tau[1]];
      xt[2].tag[taued[4]] = xt[2].ftag[tau[2]];
      xt[2].tag[taued[5]] = xt[2].ftag[tau[1]];
      xt[2].ref[  tau[0]] = 0;  xt[2].ref[  tau[2]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[2]] = 0;

      pt[3]->v[tau[0]] = vx[taued[2]] ;
      xt[3].tag[taued[0]] = xt[3].ftag[tau[2]];
      xt[3].tag[taued[1]] = xt[3].ftag[tau[1]];
      xt[3].ref[  tau[3]] = 0;
      xt[3].ftag[ tau[3]] = 0;
    }
    else{
      assert(ib == tau[3]);

      pt[2]->v[tau[0]] = vx[taued[2]] ; pt[2]->v[tau[2]] = vx[taued[1]] ;
      xt[2].tag[taued[0]] = xt[2].ftag[tau[2]];
      xt[2].tag[taued[1]] = xt[2].ftag[tau[1]];
      xt[2].tag[taued[3]] = xt[2].ftag[tau[3]];
      xt[2].tag[taued[5]] = xt[2].ftag[tau[1]];
      xt[2].ref[  tau[0]] = 0;  xt[2].ref[  tau[3]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[3]] = 0;

      pt[3]->v[tau[0]] = vx[taued[1]] ;
      xt[3].tag[taued[0]] = xt[3].ftag[tau[3]];
      xt[3].tag[taued[2]] = xt[3].ftag[tau[1]];
      xt[3].ref[  tau[2]] = 0;
      xt[3].ftag[ tau[2]] = 0;
    }
  }

  /* Assignation of the xt fields to the appropriate tets */
  isxt[0] = isxt[1] = isxt[2] = isxt[3] = 0;

  for(i=0;i<4;i++){
    if(xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
    if(xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
    if(xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
    if(xt[3].ref[i] || xt[3].ftag[i] ) isxt[3] = 1;
  }

  if((pt[0])->xt){
    if(isxt[0]){
      memcpy(pxt0,&xt[0],sizeof(xTetra));
      pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;
      for(i=1;i<4;i++){
        if(isxt[i]){
          mesh->xt++;
          assert(mesh->xt < mesh->xtmax);
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(xTetra));
        }
      }
    }
    else{
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;
      for(i=1;i<4;i++){
        if(isxt[i]){
          if(firstxt){
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[(pt[i])->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
          else{
            mesh->xt++;
            assert(mesh->xt < mesh->xtmax);
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
        }
      }
      (pt[0])->xt = 0;
    }
  }
  /* Quality update */
  pt[0]->qual=orcal(mesh,newtet[0]);
  pt[1]->qual=orcal(mesh,newtet[1]);
  pt[2]->qual=orcal(mesh,newtet[2]);
  pt[3]->qual=orcal(mesh,newtet[3]);

#ifdef DEBUG
  tabtmp[4][3]=lmintmp;  tabtmp[4][4]=lmaxtmp;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[4][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[4][2]++;
      }
      lentmp=lenedg(mesh,met,pt[1]->v[iare[i][0]],pt[1]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[4][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[4][2]++;
      }
      lentmp=lenedg(mesh,met,pt[2]->v[iare[i][0]],pt[2]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[4][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[4][2]++;
      }
      lentmp=lenedg(mesh,met,pt[3]->v[iare[i][0]],pt[3]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[4][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[4][2]++;
      }
    }
#endif
}

void split3op(pMesh mesh, pSol met, int k, int vx[6]){
  pTetra        pt[5];
  xTetra        xt[5];
  pxTetra       pxt0;
  char          flg,tag0,tag1,tag5;
  int           iel,ref0,ref1,ref5;
  int           newtet[5];
  unsigned char imin12,imin03,tau[4],*taued,sym[4],symed[6],ip0,ip1,ip2,ip3,ie0,ie1;
  unsigned char ie2,ie3,ie4,ie5,isxt[5],firstxt,i;

  pt[0]  = &mesh->tetra[k];
  flg = pt[0]->flag;
  pt[0]->flag  = 0;
  newtet[0]=k;

#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        lmintmp=lentmp;
      }else if(lentmp>lmaxtmp) {
        lmaxtmp=lentmp;
      }
    }
#endif

  /* Set permutation /symmetry of vertices : generic case : 35 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &permedge[0][0];

  sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
  symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
  symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;

  switch(flg) {
  case 19:
    tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
    taued = &permedge[0][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;

  case 13:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &permedge[2][0];

    sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
    symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;
    break;

  case 37:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &permedge[2][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;

  case 22:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &permedge[10][0];

    sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
    symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;
    break;

  case 28:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &permedge[10][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;

  case 26:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &permedge[6][0];

    sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
    symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;
    break;

  case 14:
    tau[0] = 0 ; tau[1] = 2 ; tau[2] = 3 ; tau[3] = 1;
    taued = &permedge[1][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;

  case 49:
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    taued = &permedge[11][0];

    sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
    symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;
    break;

  case 50:
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    taued = &permedge[11][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;

  case 44:
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    taued = &permedge[9][0];

    sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
    symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;
    break;

  case 41:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &permedge[4][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;
  }

  ip0 = tau[sym[0]];
  ip1 = tau[sym[1]];
  ip2 = tau[sym[2]];
  ip3 = tau[sym[3]];

  ie0 = taued[symed[0]];
  ie1 = taued[symed[1]];
  ie2 = taued[symed[2]];
  ie3 = taued[symed[3]];
  ie4 = taued[symed[4]];
  ie5 = taued[symed[5]];

  /* Test : to be removed eventually */
  assert(vx[ie0] > 0);
  assert(vx[ie1] > 0);
  assert(vx[ie5] > 0);
  assert(vx[ie2] <= 0);
  assert(vx[ie3] <= 0);
  assert(vx[ie4] <= 0);

  hPop(&mesh->htab,pt[0]->v[ip0],pt[0]->v[ip1],&ref0,&tag0);
  if( ref0 || tag0 ){
    hEdge(&mesh->htab,pt[0]->v[ip0],vx[ie0],ref0,tag0);
    hEdge(&mesh->htab,pt[0]->v[ip1],vx[ie0],ref0,tag0);
  }

  hPop(&mesh->htab,pt[0]->v[ip0],pt[0]->v[ip2],&ref1,&tag1);
  if( ref1 || tag1 ){
    hEdge(&mesh->htab,pt[0]->v[ip0],vx[ie1],ref1,tag1);
    hEdge(&mesh->htab,pt[0]->v[ip2],vx[ie1],ref1,tag1);
  }

  hPop(&mesh->htab,pt[0]->v[ip2],pt[0]->v[ip3],&ref5,&tag5);
  if( ref5 || tag5 ){
    hEdge(&mesh->htab,pt[0]->v[ip2],vx[ie5],ref5,tag5);
    hEdge(&mesh->htab,pt[0]->v[ip3],vx[ie5],ref5,tag5);
  }

  imin03 = ((pt[0])->v[ip0] < (pt[0])->v[ip3]) ? ip0 : ip3;
  imin12 = ((pt[0])->v[ip1] < (pt[0])->v[ip2]) ? ip1 : ip2;

  /* Create new elements according to the current configuration */
  iel = newElt(mesh);
  assert(iel);
  pt[1] = &mesh->tetra[iel];
  pt[1] = memcpy(pt[1],pt[0],sizeof(Tetra));
  newtet[1]=iel;

  iel = newElt(mesh);
  assert(iel);
  pt[2] = &mesh->tetra[iel];
  pt[2] = memcpy(pt[2],pt[0],sizeof(Tetra));
  newtet[2]=iel;

  iel = newElt(mesh);
  assert(iel);
  pt[3] = &mesh->tetra[iel];
  pt[3] = memcpy(pt[3],pt[0],sizeof(Tetra));
  newtet[3]=iel;

  if((pt[0])->xt){
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    memcpy(&xt[0],pxt0, sizeof(xTetra));
    memcpy(&xt[1],pxt0, sizeof(xTetra));
    memcpy(&xt[2],pxt0, sizeof(xTetra));
    memcpy(&xt[3],pxt0, sizeof(xTetra));
  }
  else{
    pxt0 = 0;
    memset(&xt[0],0, sizeof(xTetra));
    memset(&xt[1],0, sizeof(xTetra));
    memset(&xt[2],0, sizeof(xTetra));
    memset(&xt[3],0, sizeof(xTetra));
  }

  if(!((imin12 == ip1) && (imin03 == ip3))){
    iel = newElt(mesh);
    assert(iel);
    pt[4] = &mesh->tetra[iel];
    pt[4] = memcpy(pt[4],pt[0],sizeof(Tetra));
    newtet[4]=iel;

    if(pt[0]->xt){
      pxt0 = &mesh->xtetra[(pt[0])->xt];
      memcpy(&xt[4],pxt0, sizeof(xTetra));
    }

    else{
      pxt0 = 0;
      memset(&xt[4],0, sizeof(xTetra));
    }
  }

  /* Generic formulation of split of 3 edges in op configuration (edges 0,1,5 splitted) */
  if((imin12 == ip2) && (imin03 == ip0)){
    pt[0]->v[ip0] = vx[ie1] ;  pt[0]->v[ip1] = vx[ie0] ; pt[0]->v[ip3] = vx[ie5] ;
    xt[0].ref[ip0] = 0 ; xt[0].ref[ip2] = 0 ;
    xt[0].ftag[ip0] = 0 ; xt[0].ftag[ip2] = 0 ;

    pt[1]->v[ip0] = vx[ie0] ; pt[1]->v[ip3] = vx[ie5] ;
    xt[1].ref[ip1] = 0 ; xt[1] .ref[ip2] = 0 ;
    xt[1].ftag[ip1] = 0 ; xt[1].ftag[ip2] = 0 ;

    pt[2]->v[ip0] = vx[ie0] ; pt[2]->v[ip2] = vx[ie5] ;
    xt[2].ref[ip1] = 0 ; xt[2].ref[ip3] = 0 ;
    xt[2].ftag[ip1] = 0 ; xt[2].ftag[ip3] = 0 ;

    pt[3]->v[ip1] = vx[ie0] ; pt[3]->v[ip2] = vx[ie1] ; pt[3]->v[ip3] = vx[ie5] ;
    xt[3].ref[ip0] = 0 ; xt[3].ref[ip2] = 0 ;
    xt[3].ftag[ip0] = 0 ; xt[3].ftag[ip2] = 0 ;

    pt[4]->v[ip1] = vx[ie0] ; pt[4]->v[ip2] = vx[ie5];
    xt[4].ref[ip0] = 0 ; xt[4].ref[ip3] = 0 ;
    xt[4].ftag[ip0] = 0 ; xt[4].ftag[ip3] = 0 ;
  }

  else if((imin12 == ip1) && (imin03 == ip0)){
    pt[0]->v[ip0] = vx[ie1] ; pt[0]->v[ip3] = vx[ie5] ;
    xt[0].ref[ip2] = 0 ;
    xt[0].ftag[ip2] = 0 ;

    pt[1]->v[ip0] = vx[ie0] ; pt[1]->v[ip2] = vx[ie1] ; pt[1]->v[ip3] = vx[ie5];
    xt[1].ref[ip0] = 0 ; xt[1].ref[ip1] = 0 ; xt[1].ref[ip2] = 0 ;
    xt[1].ftag[ip0] = 0 ; xt[1].ftag[ip1] = 0 ; xt[1].ftag[ip2] = 0 ;

    pt[2]->v[ip0] = vx[ie0] ; pt[2]->v[ip2] = vx[ie5] ;
    xt[2].ref[ip1] = 0 ; xt[2].ref[ip3] = 0 ;
    xt[2].ftag[ip1] = 0 ; xt[2].ftag[ip3] = 0 ;

    pt[3]->v[ip1] = vx[ie0] ; pt[3]->v[ip2] = vx[ie5];
    xt[3].ref[ip0] = 0 ; xt[3].ref[ip3] = 0 ;
    xt[3].ftag[ip0] = 0 ; xt[3].ftag[ip3] = 0 ;

    pt[4]->v[ip1] = vx[ie0] ; pt[4]->v[ip2] = vx[ie1]; pt[4]->v[ip3] = vx[ie5];
    xt[4].ref[ip0] = 0 ; xt[4].ref[ip2] = 0 ;
    xt[4].ftag[ip0] = 0 ; xt[4].ftag[ip2] = 0 ;

  }

  else if((imin12 == ip2) && (imin03 == ip3)){
    pt[0]->v[ip1] = vx[ie0] ; pt[0]->v[ip2] = vx[ie1] ;
    xt[0].ref[ip0] = 0 ;
    xt[0].ftag[ip0] = 0 ;

    pt[1]->v[ip0] = vx[ie1] ; pt[1]->v[ip1] = vx[ie0] ; pt[1]->v[ip2] = vx[ie5];
    xt[1].ref[ip1] = 0 ; xt[1].ref[ip2] = 0 ;  xt[1].ref[ip3] = 0 ;
    xt[1].ftag[ip1] = 0 ; xt[1].ftag[ip2] = 0 ; xt[1].ftag[ip3] = 0 ;

    pt[2]->v[ip0] = vx[ie0] ; pt[2]->v[ip2] = vx[ie5] ;
    xt[2].ref[ip1] = 0 ; xt[2].ref[ip3] = 0 ;
    xt[2].ftag[ip1] = 0 ; xt[2].ftag[ip3] = 0 ;

    pt[3]->v[ip0] = vx[ie1] ; pt[3]->v[ip1] = vx[ie0]; pt[3]->v[ip3] = vx[ie5];
    xt[3].ref[ip0] = 0 ; xt[3].ref[ip2] = 0 ;
    xt[3].ftag[ip0] = 0 ; xt[3].ftag[ip2] = 0 ;

    pt[4]->v[ip0] = vx[ie0] ; pt[4]->v[ip3] = vx[ie5];
    xt[4].ref[ip1] = 0 ; xt[4].ref[ip2] = 0 ;
    xt[4].ftag[ip1] = 0 ; xt[4].ftag[ip2] = 0 ;
  }
  else{
    assert((imin12 == ip1) && (imin03 == ip3)) ;

    pt[0]->v[ip1] = vx[ie0] ; pt[0]->v[ip2] = vx[ie1] ;
    xt[0].ref[ip0] = 0 ;
    xt[0].ftag[ip0] = 0 ;

    pt[1]->v[ip0] = vx[ie1] ; pt[1]->v[ip3] = vx[ie5] ;
    xt[1].ref[ip2] = 0 ;
    xt[1].ftag[ip2] = 0 ;

    pt[2]->v[ip0] = vx[ie0] ; pt[2]->v[ip2] = vx[ie1] ;
    xt[2].ref[ip0] = 0 ; xt[2].ref[ip1] = 0 ;
    xt[2].ftag[ip0] = 0 ; xt[2].ftag[ip1] = 0 ;

    pt[3]->v[ip0] = vx[ie1] ; pt[3]->v[ip2] = vx[ie5] ;
    xt[3].ref[ip2] = 0 ; xt[3].ref[ip3] = 0 ;
    xt[3].ftag[ip2] = 0 ; xt[3].ftag[ip3] = 0 ;
  }

  /* Assignation of the xt fields to the appropriate tets */
  if((imin12 == ip1) && (imin03 == ip3)){
    isxt[0] = isxt[1] = isxt[2] = isxt[3] = 0;

    for(i=0;i<4;i++){
      if((xt[0]).ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
      if((xt[1]).ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
      if((xt[2]).ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
      if((xt[3]).ref[i] || xt[3].ftag[i] ) isxt[3] = 1;
    }

    if(pt[0]->xt){
      if(isxt[0]){
        memcpy(pxt0,&xt[0],sizeof(xTetra));
        pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;

        for(i=1;i<4;i++){
          if(isxt[i]){
            mesh->xt++;
            assert(mesh->xt < mesh->xtmax);
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
        }
      }
      else{
        firstxt = 1;
        pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;

        for(i=1;i<4;i++){
          if(isxt[i]){
            if(firstxt){
              firstxt = 0;
              pt[i]->xt = pt[0]->xt;
              pxt0 = &mesh->xtetra[(pt[i])->xt];
              memcpy(pxt0,&(xt[i]),sizeof(xTetra));
            }
            else{
              mesh->xt++;
              assert(mesh->xt < mesh->xtmax);
              pt[i] ->xt = mesh->xt;
              pxt0 = &mesh->xtetra[mesh->xt];
              memcpy(pxt0,&xt[i],sizeof(xTetra));
            }
          }
        }
        pt[0]->xt = 0;
      }
    }

  }
  else{
    isxt[0] = isxt[1] = isxt[2] = isxt[3] = isxt[4] = 0;

    for(i=0;i<4;i++){
      if((xt[0]).ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
      if((xt[1]).ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
      if((xt[2]).ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
      if((xt[3]).ref[i] || xt[3].ftag[i] ) isxt[3] = 1;
      if((xt[4]).ref[i] || xt[4].ftag[i] ) isxt[4] = 1;
    }

    if(pt[0]->xt){
      if(isxt[0]){
        memcpy(pxt0,&(xt[0]),sizeof(xTetra));
        pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = 0;

        for(i=1;i<5;i++){
          if(isxt[i]){
            mesh->xt++;
            assert(mesh->xt < mesh->xtmax);
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
        }
      }
      else{
        firstxt = 1;
        pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = 0;

        for(i=1; i<5; i++){
          if( isxt[i] ){
            if(firstxt){
              firstxt = 0;
              pt[i]->xt = pt[0]->xt;
              pxt0 = &mesh->xtetra[pt[i]->xt];
              memcpy(pxt0,&xt[i],sizeof(xTetra));
            }
            else{
              mesh->xt++;
              assert(mesh->xt < mesh->xtmax);
              pt[i]->xt = mesh->xt;
              pxt0 = &mesh->xtetra[mesh->xt];
              memcpy(pxt0,&xt[i],sizeof(xTetra));
            }
          }
        }
        pt[0]->xt = 0;
      }
    }
  }
  /* Quality update */
  pt[0]->qual=orcal(mesh,newtet[0]);
  pt[1]->qual=orcal(mesh,newtet[1]);
  pt[2]->qual=orcal(mesh,newtet[2]);
  pt[3]->qual=orcal(mesh,newtet[3]);
  if(!((imin12 == ip1) && (imin03 == ip3))){
    pt[4]->qual=orcal(mesh,newtet[4]);
  }

#ifdef DEBUG
  if(!((imin12 == ip1) && (imin03 == ip3))){
    tabtmp[5][0]--;
    tabtmp[6][0]++;
  }
  iatmp=vx[ie0];ibtmp=vx[ie1];ictmp=vx[ie5];

  tabtmp[5][3]=lmintmp;  tabtmp[5][4]=lmaxtmp;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        if(pt[0]->v[iare[i][0]]==iatmp ||pt[0]->v[iare[i][1]]==iatmp){
          if(pt[0]->v[iare[i][0]]==ibtmp || pt[0]->v[iare[i][1]]==ibtmp
             || pt[0]->v[iare[i][0]]==ictmp || pt[0]->v[iare[i][0]]==ictmp){
            tabtmp[5][1]=tabtmp[5][1]+1./3.;
          }else{
            tabtmp[5][1]=tabtmp[5][1]+1./2.;
          }}else tabtmp[5][1]=tabtmp[5][1]+1./2.;
      }else if(lentmp>lmaxtmp) {
        if(pt[0]->v[iare[i][0]]==iatmp ||pt[0]->v[iare[i][1]]==iatmp){
          if(pt[0]->v[iare[i][0]]==ibtmp || pt[0]->v[iare[i][1]]==ibtmp
             || pt[0]->v[iare[i][0]]==ictmp || pt[0]->v[iare[i][0]]==ictmp){
            tabtmp[5][2]=tabtmp[5][2]+1./3.;
          }else{
            tabtmp[5][2]=tabtmp[5][2]+1./2.;
          }}else tabtmp[5][2]=tabtmp[5][2]+1./2.;
      }
      lentmp=lenedg(mesh,met,pt[1]->v[iare[i][0]],pt[1]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        if(pt[1]->v[iare[i][0]]==iatmp ||pt[1]->v[iare[i][1]]==iatmp){
          if(pt[1]->v[iare[i][0]]==ibtmp || pt[1]->v[iare[i][1]]==ibtmp
             || pt[1]->v[iare[i][0]]==ictmp || pt[1]->v[iare[i][0]]==ictmp){
            tabtmp[5][1]=tabtmp[5][1]+1./3.;
          }else{
            tabtmp[5][1]=tabtmp[5][1]+1./2.;
          }}else tabtmp[5][1]=tabtmp[5][1]+1./2.;
      }else if(lentmp>lmaxtmp) {
        if(pt[1]->v[iare[i][0]]==iatmp ||pt[1]->v[iare[i][1]]==iatmp){
          if(pt[1]->v[iare[i][0]]==ibtmp || pt[1]->v[iare[i][1]]==ibtmp
             || pt[1]->v[iare[i][0]]==ictmp || pt[1]->v[iare[i][0]]==ictmp){
            tabtmp[5][2]=tabtmp[5][2]+1./3.;
          }else{
            tabtmp[5][2]=tabtmp[5][2]+1./2.;
          }}else tabtmp[5][2]=tabtmp[5][2]+1./2.;
      }
      lentmp=lenedg(mesh,met,pt[2]->v[iare[i][0]],pt[2]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        if(pt[2]->v[iare[i][0]]==iatmp ||pt[2]->v[iare[i][1]]==iatmp){
          if(pt[2]->v[iare[i][0]]==ibtmp || pt[2]->v[iare[i][1]]==ibtmp
             || pt[2]->v[iare[i][0]]==ictmp || pt[2]->v[iare[i][0]]==ictmp){
            tabtmp[5][1]=tabtmp[5][1]+1./3.;
          }else{
            tabtmp[5][1]=tabtmp[5][1]+1./2.;
          }}else tabtmp[5][1]=tabtmp[5][1]+1./2.;
      }else if(lentmp>lmaxtmp) {
        if(pt[2]->v[iare[i][0]]==iatmp ||pt[2]->v[iare[i][1]]==iatmp){
          if(pt[2]->v[iare[i][0]]==ibtmp || pt[2]->v[iare[i][1]]==ibtmp
             || pt[2]->v[iare[i][0]]==ictmp || pt[2]->v[iare[i][0]]==ictmp){
            tabtmp[5][2]=tabtmp[5][2]+1./3.;
          }else{
            tabtmp[5][2]=tabtmp[5][2]+1./2.;
          }}else tabtmp[5][2]=tabtmp[5][2]+1./2.;
      }
      lentmp=lenedg(mesh,met,pt[3]->v[iare[i][0]],pt[3]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        if(pt[3]->v[iare[i][0]]==iatmp ||pt[3]->v[iare[i][1]]==iatmp){
          if(pt[3]->v[iare[i][0]]==ibtmp || pt[3]->v[iare[i][1]]==ibtmp
             || pt[3]->v[iare[i][0]]==ictmp || pt[3]->v[iare[i][0]]==ictmp){
            tabtmp[5][1]=tabtmp[5][1]+1./3.;
          }else{
            tabtmp[5][1]=tabtmp[5][1]+1./2.;
          }}else tabtmp[5][1]=tabtmp[5][1]+1./2.;
      }else if(lentmp>lmaxtmp) {
        if(pt[3]->v[iare[i][0]]==iatmp ||pt[3]->v[iare[i][1]]==iatmp){
          if(pt[3]->v[iare[i][0]]==ibtmp || pt[3]->v[iare[i][1]]==ibtmp
             || pt[3]->v[iare[i][0]]==ictmp || pt[3]->v[iare[i][0]]==ictmp){
            tabtmp[5][2]=tabtmp[5][2]+1./3.;
          }else{
            tabtmp[5][2]=tabtmp[5][2]+1./2.;
          }}else tabtmp[5][2]=tabtmp[5][2]+1./2.;
      }
      if(!((imin12 == ip1) && (imin03 == ip3))){
        lentmp=lenedg(mesh,met,pt[4]->v[iare[i][0]],pt[4]->v[iare[i][1]]);
        if(lentmp<lmintmp) {
          if(pt[4]->v[iare[i][0]]==iatmp ||pt[4]->v[iare[i][1]]==iatmp){
            if(pt[4]->v[iare[i][0]]==ibtmp || pt[4]->v[iare[i][1]]==ibtmp
               || pt[4]->v[iare[i][0]]==ictmp || pt[4]->v[iare[i][0]]==ictmp){
            }else{
              // tabtmp[6][1]=tabtmp[5][1]+1./2.;
            }}else {/*tabtmp[6][1]=tabtmp[6][1]+1./2.;*/}
        }else if(lentmp>lmaxtmp) {
          if(pt[4]->v[iare[i][0]]==iatmp ||pt[4]->v[iare[i][1]]==iatmp){
            if(pt[4]->v[iare[i][0]]==ibtmp || pt[4]->v[iare[i][1]]==ibtmp
               || pt[4]->v[iare[i][0]]==ictmp || pt[4]->v[iare[i][0]]==ictmp){
              // tabtmp[6][2]++;
            }else{
              // tabtmp[6][2]=tabtmp[6][2]+1./2.;
            }}else{ /*tabtmp[6][2]=tabtmp[6][2]+1./2.;*/}
        }
      }
    }
#endif
}

/** Split a tetra in 4 tetras by introducing its barycenter
    FOR NOW : flags, that tell which edge should be split, are not updated (erased) : UPDATE NEEDED ?*/
int split4bar(pMesh mesh, pSol met, int k){
  pTetra   pt[4];
  pPoint   ppt;
  xTetra   xt[4];
  pxTetra  pxt0;
  double   o[3],hnew;
  int      i,ib,iel;
  int      newtet[4];
  unsigned char isxt[4],firstxt;

  pt[0] = &mesh->tetra[k];
  pt[0]->flag = 0;
  newtet[0]=k;

#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp = lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if ( lentmp < lmintmp ) {
        lmintmp = lentmp;
      }
			else if ( lentmp > lmaxtmp ) {
        lmaxtmp = lentmp;
      }
    }
#endif
  o[0] = o[1] = o[2] = 0.0;
  hnew = 0.0;
  for (i=0; i<4; i++) {
    ib    = pt[0]->v[i];
    ppt   = &mesh->point[ib];
    o[0] += ppt->c[0];
    o[1] += ppt->c[1];
    o[2] += ppt->c[2];
    if ( met->m )  hnew += met->m[ib];
  }
  o[0] *= 0.25;
  o[1] *= 0.25;
  o[2] *= 0.25;
  hnew *= 0.25;

  ib = newPt(mesh,o,0);
  if(!ib){
    fprintf(stdout,"%s:%d: Error: unable to allocate a new point\n"
            ,__FILE__,__LINE__);
    return(0);
  }
  if ( met->m )  met->m[ib] = hnew;

  /* create 3 new tetras */
  iel = newElt(mesh);
  if(!iel){
    fprintf(stdout,"%s:%d: Error: unable to allocate a new element\n"
            ,__FILE__,__LINE__);
    delPt(mesh,ib);
    return(0);
  }
  pt[1] = &mesh->tetra[iel];
  pt[1] = memcpy(pt[1],pt[0],sizeof(Tetra));
  newtet[1]=iel;

  iel = newElt(mesh);
  if(!iel){
    fprintf(stdout,"%s:%d: Error: unable to allocate a new element\n"
            ,__FILE__,__LINE__);
    delPt(mesh,ib);
    delElt(mesh,newtet[1]);
    return(0);
  }
  pt[2] = &mesh->tetra[iel];
  pt[2] = memcpy(pt[2],pt[0],sizeof(Tetra));
  newtet[2]=iel;

  iel = newElt(mesh);
  if(!iel){
    fprintf(stdout,"%s:%d: Error: unable to allocate a new element\n"
            ,__FILE__,__LINE__);
    delPt(mesh,ib);
    delElt(mesh,newtet[1]);
    delElt(mesh,newtet[2]);
    return(0);
  }
  pt[3] = &mesh->tetra[iel];
  pt[3] = memcpy(pt[3],pt[0],sizeof(Tetra));
  newtet[3]=iel;

  memset(&xt[0],0, sizeof(xTetra));
  memset(&xt[1],0, sizeof(xTetra));
  memset(&xt[2],0, sizeof(xTetra));
  memset(&xt[3],0, sizeof(xTetra));
  pxt0 = 0;
  if ( pt[0]->xt ) {
    pxt0 = &mesh->xtetra[pt[0]->xt];
    memcpy(&xt[0],pxt0,sizeof(xTetra));
    memcpy(&xt[1],pxt0,sizeof(xTetra));
    memcpy(&xt[2],pxt0,sizeof(xTetra));
    memcpy(&xt[3],pxt0,sizeof(xTetra));
  }

  /* Update vertices and xt fields */
  pt[0]->v[0] = pt[1]->v[1] = pt[2]->v[2] = pt[3]->v[3] = ib;

  xt[0].tag[0]  = 0;
  xt[0].tag[1]  = 0;
  xt[0].tag[2]  = 0;
  xt[0].ref[1]  = 0;  xt[0].ref[2]  = 0;  xt[0].ref[3]  = 0;
  xt[0].ftag[1] = 0;  xt[0].ftag[2] = 0;  xt[0].ftag[3] = 0;

  xt[1].tag[0]  = 0;
  xt[1].tag[3]  = 0;
  xt[1].tag[4]  = 0;
  xt[1].ref[0]  = 0;  xt[1].ref[2]  = 0;  xt[1].ref[3]  = 0;
  xt[1].ftag[0] = 0;  xt[1].ftag[2] = 0;  xt[1].ftag[3] = 0;

  xt[2].tag[1]  = 0;
  xt[2].tag[3]  = 0;
  xt[2].tag[5]  = 0;
  xt[2].ref[0]  = 0;  xt[2].ref[1]  = 0;  xt[2].ref[3]  = 0;
  xt[2].ftag[0] = 0;  xt[2].ftag[1] = 0;  xt[2].ftag[3] = 0;

  xt[3].tag[2]  = 0;
  xt[3].tag[4]  = 0;
  xt[3].tag[5]  = 0;
  xt[3].ref[0]  = 0;  xt[3].ref[1]  = 0;  xt[3].ref[2]  = 0;
  xt[3].ftag[0] = 0;  xt[3].ftag[1] = 0;  xt[3].ftag[2] = 0;

  /* Assignation of the xt fields to the appropriate tets */
  memset(isxt,0,4*sizeof(char));
  for (i=0; i<4; i++) {
    if ( xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
    if ( xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
    if ( xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
    if ( xt[3].ref[i] || xt[3].ftag[i] ) isxt[3] = 1;
  }

  if ( pt[0]->xt ) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(xTetra));
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          assert(mesh->xt < mesh->xtmax-1);
          pt[i]->xt = ++mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(xTetra));
        }
        else {
          pt[i]->xt = 0;
        }
      }
    }
    else {
      firstxt = 1;
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[(pt[i])->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
          else {
            assert(mesh->xt < mesh->xtmax-1);
            pt[i]->xt = ++mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
        }
        else {
          pt[i]->xt = 0;
        }
      }
      pt[0]->xt = 0;
    }
  }
  /* Quality update */
  pt[0]->qual=orcal(mesh,newtet[0]);
  pt[1]->qual=orcal(mesh,newtet[1]);
  pt[2]->qual=orcal(mesh,newtet[2]);
  pt[3]->qual=orcal(mesh,newtet[3]);

#ifdef DEBUG
  tabtmp[7][3]=lmintmp;  tabtmp[7][4]=lmaxtmp;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[7][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[7][2]++;
      }
      lentmp=lenedg(mesh,met,pt[1]->v[iare[i][0]],pt[1]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[7][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[7][2]++;
      }
      lentmp=lenedg(mesh,met,pt[2]->v[iare[i][0]],pt[2]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[7][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[7][2]++;
      }
      lentmp=lenedg(mesh,met,pt[3]->v[iare[i][0]],pt[3]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[7][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[7][2]++;
      }
    }
#endif
  return(1);
}

/** Split 4 edges in a configuration when 3 lie on the same face */
void split4sf(pMesh mesh,pSol met,int k,int vx[6]) {
  pTetra    pt[6];
  xTetra    xt[6];
  pxTetra   pxt0;
  int       iel,ref0,ref1,ref2,ref4;
  int       newtet[6];
  char      flg,firstxt,isxt[6],imin12,imin23,j,i,tag0,tag1,tag2,tag4;
  unsigned char tau[4],*taued;

  pt[0]  = &mesh->tetra[k];
  flg = pt[0]->flag;
  pt[0]->flag  = 0;
  newtet[0]=k;

#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        lmintmp=lentmp;
      }else if(lentmp>lmaxtmp) {
        lmaxtmp=lentmp;
      }
    }
#endif
  /* Set permutation of vertices : reference configuration : 23 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &permedge[0][0];
  switch(flg){
  case 29:
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &permedge[5][0];
    break;

  case 53:
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    taued = &permedge[9][0];
    break;

  case 60:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &permedge[10][0];
    break;

  case 57:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &permedge[4][0];
    break;

  case 58:
    tau[0] = 2 ; tau[1] = 3 ; tau[2] = 0 ; tau[3] = 1;
    taued = &permedge[8][0];
    break;

  case 27:
    tau[0] = 1 ; tau[1] = 0 ; tau[2] = 3 ; tau[3] = 2;
    taued = &permedge[3][0];
    break;

  case 15:
    tau[0] = 0 ; tau[1] = 2 ; tau[2] = 3 ; tau[3] = 1;
    taued = &permedge[1][0];
    break;

  case 43:
    tau[0] = 2 ; tau[1] = 1 ; tau[2] = 3 ; tau[3] = 0;
    taued = &permedge[7][0];
    break;

  case 39:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &permedge[2][0];
    break;

  case 54:
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    taued = &permedge[11][0];
    break;

  case 46:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &permedge[6][0];
    break;
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[1]],&ref0,&tag0);
  if( ref0 || tag0 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[0]],ref0,tag0);
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[0]],ref0,tag0);
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[2]],&ref1,&tag1);
  if( ref1 || tag1 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[1]],ref1,tag1);
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[1]],ref1,tag1);
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[3]],&ref2,&tag2);
  if( ref2 || tag2 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[2]],ref2,tag2);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[2]],ref2,tag2);
  }

  hPop(&mesh->htab,pt[0]->v[tau[1]],pt[0]->v[tau[3]],&ref4,&tag4);
  if( ref4 || tag4 ){
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[4]],ref4,tag4);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[4]],ref4,tag4);
  }

  imin23 = ((pt[0])->v[tau[2]] < (pt[0])->v[tau[3]]) ? tau[2] : tau[3];
  imin12 = ((pt[0])->v[tau[1]] < (pt[0])->v[tau[2]]) ? tau[1] : tau[2];

  /* create 5 new tetras */
  for(j=1;j<6;j++){
    iel = newElt(mesh);
    assert(iel);
    pt[j] = &mesh->tetra[iel];
    pt[j] = memcpy(pt[j],pt[0],sizeof(Tetra));
    newtet[j]=iel;
  }

  if((pt[0])->xt){
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    for(j=0;j<6;j++){
      memcpy(&xt[j],pxt0, sizeof(xTetra));
    }
  }
  else{
    pxt0 = 0;
    for(j=0;j<6;j++){
      memset(&xt[j],0, sizeof(xTetra));
    }
  }

  /* Generic formulation of split of 4 edges (with 3 on same face) */
  pt[0]->v[tau[1]] = vx[taued[0]] ;   pt[0]->v[tau[2]] = vx[taued[1]] ;   pt[0]->v[tau[3]] = vx[taued[2]];
  xt[0].ref[tau[0]] = 0 ;
  xt[0].ftag[tau[0]] = 0 ;

  pt[1]->v[tau[0]] = vx[taued[2]] ; pt[1]->v[tau[1]] = vx[taued[0]] ;
  pt[1]->v[tau[2]] = vx[taued[1]] ; pt[1]->v[tau[3]] = vx[taued[4]] ;
  xt[1].ref[tau[0]] = 0 ;  xt[1].ref[tau[1]] = 0 ; xt[1].ref[tau[3]] = 0 ;
  xt[1].ftag[tau[0]] = 0 ; xt[1].ftag[tau[1]] = 0 ; xt[1].ftag[tau[3]] = 0 ;

  if(imin12 == tau[1]){
    pt[2]->v[tau[0]] = vx[taued[0]] ; pt[2]->v[tau[2]] = vx[taued[1]] ; pt[2]->v[tau[3]] = vx[taued[4]] ;
    xt[2].ref[tau[0]] = 0 ;   xt[2].ref[tau[1]] = 0 ;
    xt[2].ftag[tau[0]] = 0 ; xt[2].ftag[tau[1]] = 0 ;

    pt[3]->v[tau[0]] = vx[taued[1]] ; pt[3]->v[tau[3]] = vx[taued[4]] ;
    xt[3].ref[tau[1]] = 0 ;   xt[3].ref[tau[2]] = 0 ;
    xt[3].ftag[tau[1]] = 0 ; xt[3].ftag[tau[2]] = 0 ;

  }
  else{
    pt[2]->v[tau[0]] = vx[taued[1]] ; pt[2]->v[tau[1]] = vx[taued[0]] ; pt[2]->v[tau[3]] = vx[taued[4]] ;
    xt[2].ref[tau[0]] = 0 ;   xt[2].ref[tau[1]] = 0 ; xt[2].ref[tau[2]] = 0 ;
    xt[2].ftag[tau[0]] = 0 ; xt[2].ftag[tau[1]] = 0 ; xt[2].ftag[tau[2]] = 0 ;

    pt[3]->v[tau[0]] = vx[taued[0]] ; pt[3]->v[tau[3]] = vx[taued[4]] ;
    xt[3].ref[tau[1]] = 0 ;
    xt[3].ftag[tau[1]] = 0 ;
  }

  if(imin23 == tau[2]){
    pt[4]->v[tau[0]] = vx[taued[2]] ; pt[4]->v[tau[1]] = vx[taued[1]] ; pt[4]->v[tau[3]] = vx[taued[4]] ;
    xt[4].ref[tau[0]] = 0 ;   xt[4].ref[tau[1]] = 0 ; xt[4].ref[tau[2]] = 0;
    xt[4].ftag[tau[0]] = 0 ; xt[4].ftag[tau[1]] = 0 ; xt[4].ftag[tau[2]] = 0;

    pt[5]->v[tau[0]] = vx[taued[2]] ; pt[5]->v[tau[1]] = vx[taued[4]] ;
    xt[5].ref[tau[3]] = 0 ;
    xt[5].ftag[tau[3]] = 0 ;
  }
  else{
    pt[4]->v[tau[0]] = vx[taued[2]] ; pt[4]->v[tau[1]] = vx[taued[4]] ; pt[4]->v[tau[2]] = vx[taued[1]] ;
    xt[4].ref[tau[0]] = 0 ;   xt[4].ref[tau[3]] = 0 ;
    xt[4].ftag[tau[0]] = 0 ; xt[4].ftag[tau[3]] = 0 ;

    pt[5]->v[tau[0]] = vx[taued[1]] ; pt[5]->v[tau[1]] = vx[taued[4]] ;
    xt[5].ref[tau[2]] = 0 ; xt[5].ref[tau[3]] = 0 ;
    xt[5].ftag[tau[2]] = 0 ; xt[5].ftag[tau[3]] = 0 ;
  }

  /* Assignation of the xt fields to the appropriate tets */
  for(j=0;j<6;j++){
    isxt[j] = 0;
  }


  for(i=0;i<4;i++){
    for(j=0;j<6;j++){
      if((xt[j]).ref[i] || xt[j].ftag[i] ) isxt[j] = 1;
    }
  }

  if(pt[0]->xt){
    if(isxt[0]){
      memcpy(pxt0,&xt[0],sizeof(xTetra));
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = 0;

      for(i=1;i<6;i++){
        if(isxt[i]){
          mesh->xt++;
          assert(mesh->xt < mesh->xtmax);
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(xTetra));
        }
      }
    }
    else{
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = 0;

      for(i=1;i<6;i++){
        if(isxt[i]){
          if(firstxt){
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[(pt[i])->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
          else{
            mesh->xt++;
            assert(mesh->xt < mesh->xtmax);
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
        }
      }
      pt[0]->xt = 0;
    }
  }
  for(i=0;i<6;i++){
    pt[i]->qual=orcal(mesh,newtet[i]);
  }

#ifdef DEBUG
  tabtmp[8][3]=lmintmp;  tabtmp[8][4]=lmaxtmp;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[8][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[8][2]++;
      }
      lentmp=lenedg(mesh,met,pt[1]->v[iare[i][0]],pt[1]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[8][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[8][2]++;
      }
      lentmp=lenedg(mesh,met,pt[2]->v[iare[i][0]],pt[2]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[8][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[8][2]++;
      }
      lentmp=lenedg(mesh,met,pt[3]->v[iare[i][0]],pt[3]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[8][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[8][2]++;
      }
      lentmp=lenedg(mesh,met,pt[4]->v[iare[i][0]],pt[4]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[8][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[8][2]++;
      }
      lentmp=lenedg(mesh,met,pt[5]->v[iare[i][0]],pt[5]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[8][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[8][2]++;
      }
    }
#endif
}

/** Split 4 edges in a configuration when no 3 edges lie on the same face */
void split4op(pMesh mesh,pSol met,int k,int vx[6]) {
  pTetra        pt[6];
  xTetra        xt[6];
  pxTetra       pxt0;
  int           iel,ref1,ref2,ref3,ref4;
  int           newtet[6];
  char          flg,firstxt,isxt[6],i,j,imin01,imin23,tag1,tag2,tag3,tag4;
  unsigned char tau[4],*taued;

  pt[0]  = &mesh->tetra[k];
  flg = pt[0]->flag;
  pt[0]->flag  = 0;
  newtet[0]=k;

#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        lmintmp=lentmp;
      }else if(lentmp>lmaxtmp) {
        lmaxtmp=lentmp;
      }
    }
#endif
  /* Set permutation of vertices : reference configuration 30 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &permedge[0][0];

  switch(flg){
  case 45:
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &permedge[5][0];
    break;

  case 51:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &permedge[4][0];
    break;
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[2]],&ref1,&tag1);
  if( ref1 || tag1 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[1]],ref1,tag1);
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[1]],ref1,tag1);
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[3]],&ref2,&tag2);
  if( ref2 || tag2 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[2]],ref2,tag2);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[2]],ref2,tag2);
  }

  hPop(&mesh->htab,pt[0]->v[tau[1]],pt[0]->v[tau[2]],&ref3,&tag3);
  if( ref3 || tag3 ){
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[3]],ref3,tag3);
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[3]],ref3,tag3);
  }

  hPop(&mesh->htab,pt[0]->v[tau[1]],pt[0]->v[tau[3]],&ref4,&tag4);
  if( ref4 || tag4 ){
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[4]],ref4,tag4);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[4]],ref4,tag4);
  }

  imin01 = ((pt[0])->v[tau[0]] < (pt[0])->v[tau[1]]) ? tau[0] : tau[1];
  imin23 = ((pt[0])->v[tau[2]] < (pt[0])->v[tau[3]]) ? tau[2] : tau[3];

  /* create 5 new tetras */
  for(j=1;j<6;j++){
    iel = newElt(mesh);
    assert(iel);
    pt[j] = &mesh->tetra[iel];
    pt[j] = memcpy(pt[j],pt[0],sizeof(Tetra));
    newtet[j]=iel;
  }

  if((pt[0])->xt){
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    for(j=0;j<6;j++){
      memcpy(&xt[j],pxt0, sizeof(xTetra));
    }
  }
  else{
    pxt0 = 0;
    for(j=0;j<6;j++){
      memset(&xt[j],0, sizeof(xTetra));
    }
  }

  /* Generic formulation for split of 4 edges, with no 3 edges lying on the same face */
  if(imin01 == tau[0]){
    pt[0]->v[tau[2]] = vx[taued[3]] ; pt[0]->v[tau[3]] = vx[taued[4]];
    xt[0].ref[tau[1]] = 0;
    xt[0].ftag[tau[1]] = 0;

    pt[1]->v[tau[1]] = vx[taued[4]] ; pt[1]->v[tau[2]] = vx[taued[3]] ; pt[1]->v[tau[3]] = vx[taued[2]];
    xt[1].ref[tau[0]] = 0 ;   xt[1].ref[tau[1]] = 0 ;   xt[1].ref[tau[3]] = 0;
    xt[1].ftag[tau[0]] = 0 ;  xt[1].ftag[tau[1]] = 0 ;  xt[1].ftag[tau[3]] = 0;

    pt[2]->v[tau[1]] = vx[taued[3]] ; pt[2]->v[tau[2]] = vx[taued[1]] ; pt[2]->v[tau[3]] = vx[taued[2]];
    xt[2].ref[tau[0]] = 0 ;   xt[2].ref[tau[2]] = 0 ;
    xt[2].ftag[tau[0]] = 0 ;  xt[2].ftag[tau[2]] = 0 ;
  }
  else{
    pt[0]->v[tau[2]] = vx[taued[1]] ; pt[0]->v[tau[3]] = vx[taued[2]];
    xt[0].ref[tau[0]] = 0;
    xt[0].ftag[tau[0]] = 0;

    pt[1]->v[tau[0]] = vx[taued[1]] ; pt[1]->v[tau[2]] = vx[taued[3]] ; pt[1]->v[tau[3]] = vx[taued[2]];
    xt[1].ref[tau[0]] = 0 ;   xt[1].ref[tau[1]] = 0 ;   xt[1].ref[tau[2]] = 0;
    xt[1].ftag[tau[0]] = 0 ;  xt[1].ftag[tau[1]] = 0 ;  xt[1].ftag[tau[2]] = 0;

    pt[2]->v[tau[0]] = vx[taued[2]] ; pt[2]->v[tau[2]] = vx[taued[3]] ; pt[2]->v[tau[3]] = vx[taued[4]];
    xt[2].ref[tau[1]] = 0 ;   xt[2].ref[tau[3]] = 0 ;
    xt[2].ftag[tau[1]] = 0 ;  xt[2].ftag[tau[3]] = 0 ;
  }

  if(imin23 == tau[2]){
    pt[3]->v[tau[0]] = vx[taued[2]] ; pt[3]->v[tau[1]] = vx[taued[4]];
    xt[3].ref[tau[3]] = 0;
    xt[3].ftag[tau[3]] = 0;

    pt[4]->v[tau[0]] = vx[taued[2]] ; pt[4]->v[tau[1]] = vx[taued[3]] ; pt[4]->v[tau[3]] = vx[taued[4]];
    xt[4].ref[tau[1]] = 0 ;   xt[4].ref[tau[2]] = 0 ;   xt[4].ref[tau[3]] = 0;
    xt[4].ftag[tau[1]] = 0 ;  xt[4].ftag[tau[2]] = 0 ;  xt[4].ftag[tau[3]] = 0;

    pt[5]->v[tau[0]] = vx[taued[1]] ; pt[5]->v[tau[1]] = vx[taued[3]] ; pt[5]->v[tau[3]] = vx[taued[2]];
    xt[5].ref[tau[0]] = 0 ;   xt[5].ref[tau[2]] = 0 ;
    xt[5].ftag[tau[0]] = 0 ;  xt[5].ftag[tau[2]] = 0 ;
  }
  else{
    pt[3]->v[tau[0]] = vx[taued[1]] ; pt[3]->v[tau[1]] = vx[taued[3]];
    xt[3].ref[tau[2]] = 0;
    xt[3].ftag[tau[2]] = 0;

    pt[4]->v[tau[0]] = vx[taued[2]] ; pt[4]->v[tau[1]] = vx[taued[3]] ; pt[4]->v[tau[2]] = vx[taued[1]];
    xt[4].ref[tau[0]] = 0 ;   xt[4].ref[tau[2]] = 0 ;   xt[4].ref[tau[3]] = 0;
    xt[4].ftag[tau[0]] = 0 ;  xt[4].ftag[tau[2]] = 0 ;  xt[4].ftag[tau[3]] = 0;

    pt[5]->v[tau[0]] = vx[taued[2]] ; (pt[5])->v[tau[1]] = vx[taued[4]] ; (pt[5])->v[tau[2]] = vx[taued[3]];
    xt[5].ref[tau[1]] = 0 ;  xt[5].ref[tau[3]] = 0 ;
    xt[5].ftag[tau[1]] = 0 ; xt[5].ftag[tau[3]] = 0 ;
  }

  /* Assignation of the xt fields to the appropriate tets */
  for(j=0;j<6;j++){
    isxt[j] = 0;
  }

  for(i=0;i<4;i++){
    for(j=0;j<6;j++){
      if((xt[j]).ref[i] || xt[j].ftag[i] ) isxt[j] = 1;
    }
  }

  // In this case, at least one of the 4 created tets must have a special field
  if(pt[0]->xt){
    if(isxt[0]){
      memcpy(pxt0,&xt[0],sizeof(xTetra));
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = 0;

      for(i=1;i<6;i++){
        if(isxt[i]){
          mesh->xt++;
          assert(mesh->xt < mesh->xtmax);
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(xTetra));
        }
      }
    }
    else{
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = 0;

      for(i=1;i<6;i++){
        if(isxt[i]){
          if(firstxt){
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[ pt[i]->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
          else{
            mesh->xt++;
            assert(mesh->xt < mesh->xtmax);
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&(xt[i]),sizeof(xTetra));
          }
        }
      }
      pt[0]->xt = 0;

    }
  }
  for(i=0;i<6;i++){
    pt[i]->qual=orcal(mesh,newtet[i]);
  }

#ifdef DEBUG
  tabtmp[9][3]=lmintmp;  tabtmp[9][4]=lmaxtmp;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[9][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[9][2]++;
      }
      lentmp=lenedg(mesh,met,pt[1]->v[iare[i][0]],pt[1]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[9][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[9][2]++;
      }
      lentmp=lenedg(mesh,met,pt[2]->v[iare[i][0]],pt[2]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[9][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[9][2]++;
      }
      lentmp=lenedg(mesh,met,pt[3]->v[iare[i][0]],pt[3]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[9][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[9][2]++;
      }
      lentmp=lenedg(mesh,met,pt[4]->v[iare[i][0]],pt[4]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[9][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[9][2]++;
      }
      lentmp=lenedg(mesh,met,pt[5]->v[iare[i][0]],pt[5]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[9][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[9][2]++;
      }
    }
#endif
}

/** Split 5 edges */
void split5(pMesh mesh,pSol met,int k,int vx[6]) {
  pTetra    pt[7];
  xTetra    xt[7];
  pxTetra   pxt0;
  int       iel,i,j,ref1,ref2,ref3,ref4,ref5;
  int       newtet[7];
  char      flg,firstxt,isxt[7],imin,tag1,tag2,tag3,tag4,tag5;
  unsigned char tau[4],*taued;

  pt[0]  = &mesh->tetra[k];
  flg = pt[0]->flag;
  pt[0]->flag  = 0;
  newtet[0]=k;

#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        lmintmp=lentmp;
      }else if(lentmp>lmaxtmp) {
        lmaxtmp=lentmp;
      }
    }
#endif
  /* create 6 new tetras */
  for(i=1;i<7;i++){
    iel = newElt(mesh);
    assert(iel);
    pt[i] = &mesh->tetra[iel];
    pt[i] = memcpy(pt[i],pt[0],sizeof(Tetra));
    newtet[i]=iel;
  }

  if(pt[0]->xt){
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    for(i=0;i<7;i++){
      memcpy(&xt[i],pxt0, sizeof(xTetra));
    }
  }
  else{
    pxt0 = 0;
    for(i=0;i<7;i++){
      memset(&xt[i],0, sizeof(xTetra));
    }
  }

  /* set permutation of vertices and edges ; reference configuration : 62 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &permedge[0][0];

  switch(flg) {
  case 61:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &permedge[6][0];
    break;

  case 59:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &permedge[2][0];
    break;

  case 55:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &permedge[4][0];
    break;

  case 47:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &permedge[10][0];
    break;

  case 31:
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    taued = &permedge[11][0];
    break;
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[2]],&ref1,&tag1);
  if( ref1 || tag1 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[1]],ref1,tag1);
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[1]],ref1,tag1);
  }

  hPop(&mesh->htab,pt[0]->v[tau[0]],pt[0]->v[tau[3]],&ref2,&tag2);
  if( ref2 || tag2 ){
    hEdge(&mesh->htab,pt[0]->v[tau[0]],vx[taued[2]],ref2,tag2);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[2]],ref2,tag2);
  }

  hPop(&mesh->htab,pt[0]->v[tau[1]],pt[0]->v[tau[2]],&ref3,&tag3);
  if( ref3 || tag3 ){
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[3]],ref3,tag3);
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[3]],ref3,tag3);
  }

  hPop(&mesh->htab,pt[0]->v[tau[1]],pt[0]->v[tau[3]],&ref4,&tag4);
  if( ref4 || tag4 ){
    hEdge(&mesh->htab,pt[0]->v[tau[1]],vx[taued[4]],ref4,tag4);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[4]],ref4,tag4);
  }

  hPop(&mesh->htab,pt[0]->v[tau[2]],pt[0]->v[tau[3]],&ref5,&tag5);
  if( ref5 || tag5 ){
    hEdge(&mesh->htab,pt[0]->v[tau[2]],vx[taued[5]],ref5,tag5);
    hEdge(&mesh->htab,pt[0]->v[tau[3]],vx[taued[5]],ref5,tag5);
  }

  /* Generic formulation of split of 5 edges */
  imin = (pt[0]->v[tau[0]] < pt[0]->v[tau[1]]) ? tau[0] : tau[1];

  pt[0]->v[tau[0]] = vx[taued[2]] ;   pt[0]->v[tau[1]] = vx[taued[4]] ;   pt[0]->v[tau[2]] = vx[taued[5]];
  xt[0].ref[tau[3]] = 0 ;
  xt[0].ftag[tau[3]] = 0 ;

  pt[1]->v[tau[0]] = vx[taued[1]] ; pt[1]->v[tau[1]] = vx[taued[3]] ; pt[1]->v[tau[3]] = vx[taued[5]];
  xt[1].ref[tau[2]] = 0 ;
  xt[1].ftag[tau[2]] = 0 ;

  pt[2]->v[tau[0]] = vx[taued[2]] ; pt[2]->v[tau[1]] = vx[taued[4]];
  pt[2]->v[tau[2]] = vx[taued[3]] ; pt[2]->v[tau[3]] = vx[taued[5]];
  xt[2].ref[tau[1]] = 0 ;   xt[2].ref[tau[2]] = 0 ;   xt[2].ref[tau[3]] = 0 ;
  xt[2].ftag[tau[1]] = 0 ;  xt[2].ftag[tau[2]] = 0 ;  xt[2].ftag[tau[3]] = 0 ;

  pt[3]->v[tau[0]] = vx[taued[2]] ; pt[3]->v[tau[1]] = vx[taued[3]];
  pt[3]->v[tau[2]] = vx[taued[1]] ; pt[3]->v[tau[3]] = vx[taued[5]];
  xt[3].ref[tau[0]] = 0   ; xt[3].ref[tau[2]] = 0   ; xt[3].ref[tau[3]] = 0 ;
  xt[3].ftag[tau[0]] = 0  ; xt[3].ftag[tau[2]] = 0  ; xt[3].ftag[tau[3]] = 0 ;

  if(imin == tau[0]){
    pt[4]->v[tau[2]] = vx[taued[3]] ; pt[4]->v[tau[3]] = vx[taued[4]];
    xt[4].ref[tau[1]] = 0;
    xt[4].ftag[tau[1]] = 0;

    pt[5]->v[tau[1]] = vx[taued[4]] ; pt[5]->v[tau[2]] = vx[taued[3]]; pt[5]->v[tau[3]] = vx[taued[2]];
    xt[5].ref[tau[0]] = 0 ; xt[5].ref[tau[1]] = 0 ; xt[5].ref[tau[3]] = 0 ;
    xt[5].ftag[tau[0]] = 0 ; xt[5].ftag[tau[1]] = 0 ; xt[5].ftag[tau[3]] = 0 ;

    pt[6]->v[tau[1]] = vx[taued[3]] ; pt[6]->v[tau[2]] = vx[taued[1]]; pt[6]->v[tau[3]] = vx[taued[2]];
    xt[6].ref[tau[0]] = 0 ; xt[6].ref[tau[2]] = 0 ;
    xt[6].ftag[tau[0]] = 0 ; xt[6].ftag[tau[2]] = 0 ;

  }
  else{
    pt[4]->v[tau[2]] = vx[taued[1]] ; pt[4]->v[tau[3]] = vx[taued[2]];
    xt[4].ref[tau[0]] = 0;
    xt[4].ftag[tau[0]] = 0;

    pt[5]->v[tau[0]] = vx[taued[2]] ; pt[5]->v[tau[2]] = vx[taued[3]]; pt[5]->v[tau[3]] = vx[taued[4]];
    xt[5].ref[tau[1]] = 0 ; xt[5].ref[tau[3]] = 0 ;
    xt[5].ftag[tau[1]] = 0 ; xt[5].ftag[tau[3]] = 0 ;

    pt[6]->v[tau[0]] = vx[taued[1]] ; pt[6]->v[tau[2]] = vx[taued[3]]; pt[6]->v[tau[3]] = vx[taued[2]];
    xt[6].ref[tau[0]] = 0 ; xt[6].ref[tau[1]] = 0 ;   xt[6].ref[tau[2]] = 0 ;
    xt[6].ftag[tau[0]] = 0 ; xt[6].ftag[tau[1]] = 0 ; xt[6].ftag[tau[2]] = 0 ;
  }

  /* Assignation of the xt fields to the appropriate tets */
  for(j=0;j<7;j++){
    isxt[j] = 0;
  }

  for(i=0;i<4;i++){
    for(j=0;j<7;j++){
      if((xt[j]).ref[i] || xt[j].ftag[i] ) isxt[j] = 1;
    }
  }

  if(pt[0]->xt){
    if(isxt[0]){
      memcpy(pxt0,&xt[0],sizeof(xTetra));
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = pt[6]->xt = 0;

      for(i=1;i<7;i++){
        if(isxt[i]){
          mesh->xt++;
          assert(mesh->xt < mesh->xtmax);
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(xTetra));
        }
      }
    }
    else{
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = pt[6]->xt = 0;

      for(i=1;i<7;i++){
        if(isxt[i]){
          if(firstxt){
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[(pt[i])->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
          else{
            mesh->xt++;
            assert(mesh->xt < mesh->xtmax);
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(xTetra));
          }
        }
      }
      pt[0]->xt = 0;

    }
  }
  for(i=0;i<7;i++){
    pt[i]->qual=orcal(mesh,newtet[i]);
  }
#ifdef DEBUG
  tabtmp[10][3]=lmintmp;  tabtmp[10][4]=lmaxtmp;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[10][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[10][2]++;
      }
      lentmp=lenedg(mesh,met,pt[1]->v[iare[i][0]],pt[1]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[10][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[10][2]++;
      }
      lentmp=lenedg(mesh,met,pt[2]->v[iare[i][0]],pt[2]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[10][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[10][2]++;
      }
      lentmp=lenedg(mesh,met,pt[3]->v[iare[i][0]],pt[3]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[10][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[10][2]++;
      }
      lentmp=lenedg(mesh,met,pt[4]->v[iare[i][0]],pt[4]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[10][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[10][2]++;
      }
      lentmp=lenedg(mesh,met,pt[5]->v[iare[i][0]],pt[5]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[10][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[10][2]++;
      }
      lentmp=lenedg(mesh,met,pt[6]->v[iare[i][0]],pt[6]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[10][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[10][2]++;
      }
    }
#endif
}

/** split all faces (6 edges) */
void split6(pMesh mesh,pSol met,int k,int vx[6]) {
  pTetra    pt[8];
  xTetra    xt0,xt;
  pxTetra   pxt;
  int       i,iel,nxt0,ref0,ref1,ref2,ref3,ref4,ref5;
  int       newtet[8];
  char      isxt0,isxt,tag0,tag1,tag2,tag3,tag4,tag5;

  pt[0]  = &mesh->tetra[k];
  pt[0]->flag  = 0;
  newtet[0]=k;

#ifdef DEBUG
  lmintmp=0.6;lmaxtmp=1.3;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        lmintmp=lentmp;
      }else if(lentmp>lmaxtmp) {
        lmaxtmp=lentmp;
      }
    }
#endif
  nxt0 = pt[0]->xt;
  pxt = &mesh->xtetra[nxt0];
  memcpy(&xt0,pxt,sizeof(xTetra));

  /* create 7 new tetras */
  for(i=1;i<8;i++){
    iel = newElt(mesh);
    assert(iel);
    pt[i] = &mesh->tetra[iel];
    pt[i] = memcpy(pt[i],pt[0],sizeof(Tetra));
    newtet[i]=iel;
  }

  /* Modify edge information */
  hPop(&mesh->htab,pt[0]->v[0],pt[0]->v[1],&ref0,&tag0);
  if( ref0 || tag0 ){
    hEdge(&mesh->htab,pt[0]->v[0],vx[0],ref0,tag0);
    hEdge(&mesh->htab,pt[0]->v[1],vx[0],ref0,tag0);
  }

  hPop(&mesh->htab,pt[0]->v[0],pt[0]->v[2],&ref1,&tag1);
  if( ref1 || tag1 ){
    hEdge(&mesh->htab,pt[0]->v[0],vx[1],ref1,tag1);
    hEdge(&mesh->htab,pt[0]->v[2],vx[1],ref1,tag1);
  }

  hPop(&mesh->htab,pt[0]->v[0],pt[0]->v[3],&ref2,&tag2);
  if( ref2 || tag2 ){
    hEdge(&mesh->htab,pt[0]->v[0],vx[2],ref2,tag2);
    hEdge(&mesh->htab,pt[0]->v[3],vx[2],ref2,tag2);
  }

  hPop(&mesh->htab,pt[0]->v[1],pt[0]->v[2],&ref3,&tag3);
  if( ref3 || tag3 ){
    hEdge(&mesh->htab,pt[0]->v[1],vx[3],ref3,tag3);
    hEdge(&mesh->htab,pt[0]->v[2],vx[3],ref3,tag3);
  }

  hPop(&mesh->htab,pt[0]->v[1],pt[0]->v[3],&ref4,&tag4);
  if( ref4 || tag4 ){
    hEdge(&mesh->htab,pt[0]->v[1],vx[4],ref4,tag4);
    hEdge(&mesh->htab,pt[0]->v[3],vx[4],ref4,tag4);
  }

  hPop(&mesh->htab,pt[0]->v[2],pt[0]->v[3],&ref5,&tag5);
  if( ref5 || tag5 ){
    hEdge(&mesh->htab,pt[0]->v[2],vx[5],ref5,tag5);
    hEdge(&mesh->htab,pt[0]->v[3],vx[5],ref5,tag5);
  }

  /* Modify first tetra */
  pt[0]->v[1] = vx[0] ; pt[0]->v[2] = vx[1]; pt[0]->v[3] = vx[2];
  if(nxt0){
    memcpy(&xt,&xt0,sizeof(xTetra));

    xt.ref[0] = 0;  xt.ftag[0] = 0;
    isxt0 = 0;
    for(i=0;i<4;i++){
      if((xt.ref[i]) || xt.ftag[i] ) isxt0 = 1;
    }

    if(isxt0){
      memcpy(pxt,&xt,sizeof(xTetra));
    }
    else{
      pt[0]->xt = 0;
    }
  }

  /* Modify second tetra */
  pt[1]->v[0] = vx[0] ; pt[1]->v[2] = vx[3]; pt[1]->v[3] = vx[4];

  if(nxt0){
    memcpy(&xt,&xt0,sizeof(xTetra));
    xt.ref[1] = 0;  xt.ftag[1] = 0;

    isxt = 0;

    for(i=0;i<4;i++){
      if((xt.ref[i]) || xt.ftag[i]) isxt = 1;
    }

    pt[1]->xt = 0;
    if(isxt){
      if(!isxt0){
        isxt0 = 1;
        pt[1]->xt = nxt0;
        pxt = &mesh->xtetra[pt[1]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
      else{
        mesh->xt++;
        assert(mesh->xt < mesh->xtmax);
        pt[1]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[1]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
    }
  }

  /* Modify 3rd tetra */
  pt[2]->v[0] = vx[1] ; pt[2]->v[1] = vx[3]; pt[2]->v[3] = vx[5];

  if(nxt0){
    memcpy(&xt,&xt0,sizeof(xTetra));
    xt.ref[2] = 0;  xt.ftag[2] = 0;

    isxt = 0;

    for(i=0;i<4;i++){
      if((xt.ref[i]) || xt.ftag[i]) isxt = 1;
    }

    pt[2]->xt = 0;
    if(isxt){
      if(!isxt0){
        isxt0 = 1;
        pt[2]->xt = nxt0;
        pxt = &mesh->xtetra[pt[2]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
      else{
        mesh->xt++;
        assert(mesh->xt < mesh->xtmax);
        pt[2]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[2]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
    }
  }

  /* Modify 4th tetra */
  pt[3]->v[0] = vx[2] ; pt[3]->v[1] = vx[4]; pt[3]->v[2] = vx[5];

  if(nxt0){
    memcpy(&xt,&xt0,sizeof(xTetra));
    xt.ref[3] = 0;  xt.ftag[3] = 0;

    isxt = 0;

    for(i=0;i<4;i++){
      if((xt.ref[i]) || xt.ftag[i]) isxt = 1;
    }

    pt[3]->xt = 0;
    if(isxt){
      if(!isxt0){
        isxt0 = 1;
        pt[3]->xt = nxt0;
        pxt = &mesh->xtetra[pt[3]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
      else{
        mesh->xt++;
        assert(mesh->xt < mesh->xtmax);
        pt[3]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[3]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
    }
  }

  /* Modify 5th tetra */
  pt[4]->v[0] = vx[0] ; pt[4]->v[1] = vx[3]; pt[4]->v[2] = vx[1] ; pt[4]->v[3] = vx[2];

  if(nxt0){
    memcpy(&xt,&xt0,sizeof(xTetra));
    xt.ref[0] = 0 ; xt.ref[1] = 0 ; xt.ref[2] = 0;
    xt.ftag[0] = 0 ; xt.ftag[1] = 0 ; xt.ftag[2] = 0;

    isxt = 0;

    for(i=0;i<4;i++){
      if((xt.ref[i]) || xt.ftag[i]) isxt = 1;
    }

    pt[4]->xt = 0;
    if(isxt){
      if(!isxt0){
        isxt0 = 1;
        pt[4]->xt = nxt0;
        pxt = &mesh->xtetra[(pt[4])->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
      else{
        mesh->xt++;
        assert(mesh->xt < mesh->xtmax);
        pt[4]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[4]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
    }
  }

  /* Modify 6th tetra */
  pt[5]->v[0] = vx[2] ; pt[5]->v[1] = vx[0]; pt[5]->v[2] = vx[3] ; pt[5]->v[3] = vx[4];

  if(nxt0){
    memcpy(&xt,&xt0,sizeof(xTetra));
    xt.ref[0] = 0 ; xt.ref[1] = 0 ; xt.ref[3] = 0;
    xt.ftag[0] = 0 ; xt.ftag[1] = 0 ; xt.ftag[3] = 0;

    isxt = 0;

    for(i=0;i<4;i++){
      if((xt.ref[i]) || xt.ftag[i]) isxt = 1;
    }

    pt[5]->xt = 0;
    if(isxt){
      if(!isxt0){
        isxt0 = 1;
        pt[5]->xt = nxt0;
        pxt = &mesh->xtetra[pt[5]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
      else{
        mesh->xt++;
        assert(mesh->xt < mesh->xtmax);
        pt[5]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[5]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
    }
  }

  /* Modify 7th tetra */
  pt[6]->v[0] = vx[2] ; pt[6]->v[1] = vx[3]; pt[6]->v[2] = vx[1] ; pt[6]->v[3] = vx[5];

  if(nxt0){
    memcpy(&xt,&xt0,sizeof(xTetra));
    xt.ref[0] = 0 ; xt.ref[2] = 0 ; xt.ref[3] = 0;
    xt.ftag[0] = 0 ; xt.ftag[2] = 0 ; xt.ftag[3] = 0;

    isxt = 0;

    for(i=0;i<4;i++){
      if((xt.ref[i]) || xt.ftag[i]) isxt = 1;
    }

    pt[6]->xt = 0;
    if(isxt){
      if(!isxt0){
        isxt0 = 1;
        pt[6]->xt = nxt0;
        pxt = &mesh->xtetra[pt[6]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
      else{
        mesh->xt++;
        assert(mesh->xt < mesh->xtmax);
        pt[6]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[6]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
    }
  }

  /* Modify last tetra */
  pt[7]->v[0] = vx[2] ; pt[7]->v[1] = vx[3]; pt[7]->v[2] = vx[5] ; pt[7]->v[3] = vx[4];

  if(nxt0){
    memcpy(&xt,&xt0,sizeof(xTetra));

    xt.ref[1] = 0 ; xt.ref[2] = 0 ; xt.ref[3] = 0;
    xt.ftag[1] = 0 ; xt.ftag[2] = 0 ; xt.ftag[3] = 0;

    isxt = 0;

    for(i=0;i<4;i++){
      if((xt.ref[i]) || xt.ftag[i]) isxt = 1;
    }

    pt[7]->xt = 0;
    if(isxt){
      if(!isxt0){
        isxt0 = 1;
        pt[7]->xt = nxt0;
        pxt = &mesh->xtetra[pt[7]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
      else{
        mesh->xt++;
        assert(mesh->xt < mesh->xtmax);
        pt[7]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[7]->xt];
        memcpy(pxt,&xt,sizeof(xTetra));
      }
    }
  }
  for(i=0;i<8;i++){
    pt[i]->qual=orcal(mesh,newtet[i]);
  }
#ifdef DEBUG
  tabtmp[11][3]=lmintmp;  tabtmp[11][4]=lmaxtmp;
  for(i=0;i<6;i++)
    {
      lentmp=lenedg(mesh,met,pt[0]->v[iare[i][0]],pt[0]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[11][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[11][2]++;
      }
      lentmp=lenedg(mesh,met,pt[1]->v[iare[i][0]],pt[1]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[11][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[11][2]++;
      }
      lentmp=lenedg(mesh,met,pt[2]->v[iare[i][0]],pt[2]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[11][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[11][2]++;
      }
      lentmp=lenedg(mesh,met,pt[3]->v[iare[i][0]],pt[3]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[11][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[11][2]++;
      }
      lentmp=lenedg(mesh,met,pt[4]->v[iare[i][0]],pt[4]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[11][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[11][2]++;
      }
      lentmp=lenedg(mesh,met,pt[5]->v[iare[i][0]],pt[5]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[11][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[11][2]++;
      }
      lentmp=lenedg(mesh,met,pt[7]->v[iare[i][0]],pt[7]->v[iare[i][1]]);
      if(lentmp<lmintmp) {
        tabtmp[11][1]++;
      }else if(lentmp>lmaxtmp) {
        tabtmp[11][2]++;
      }
    }
#endif
}
