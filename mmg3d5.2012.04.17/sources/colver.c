#include "mmg3d.h"

extern Info  info;
extern char  ddb;

/* Check whether collapse ip -> iq could be performed, ip internal ;
   'mechanical' tests (positive jacobian) are not performed here */
int chkcol_int(pMesh mesh,int k,char iface,char iedg,int *list) {
	pTetra   pt,pt0;
	double   calold,calnew;
  int      j,iel,ilist,nq;
	char     ip,iq;

  ip = idir[iface][inxt2[iedg]];
  iq = idir[iface][iprv2[iedg]];
	pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
	nq  = pt->v[iq];
  ilist = boulevolp(mesh,k,ip,list);

	calold = calnew = DBL_MAX;
  for (j=0; j<ilist; j++) {
    if ( list[j] < 0 )  continue;
    iel = list[j] / 4;
    ip  = list[j] % 4;  
    pt  = &mesh->tetra[iel];
    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[ip] = nq;

    calold = MG_MIN(calold,orcal(mesh,iel));
    calnew = MG_MIN(calnew,orcal(mesh,0));
  }
  if ( calold < NULKAL && calnew <= calold )  return(0);
  else if ( calnew < 0.3*calold )  return(0);

  return(ilist);
}

/* Check whether collapse ip -> iq could be performed, ip boundary point ;
 'mechanical' tests (positive jacobian) are not performed here ; 
 iface = boundary face on which lie edge iedg - in local face num. 
      (pq, or ia in local tet notation) */
int chkcol_bdy(pMesh mesh,int k,char iface,char iedg,int *listv) {
  pTetra        pt,pt0;
  pPoint        p0;
  Tria          tt;
  double        calold,calnew,nprvold[3],nprvnew[3],ncurold[3],ncurnew[3],ps,devold,devnew;
  int           ipp,ilistv,nump,numq,ilists,lists[LMAX+2],l,iel,ref;
  unsigned char iopp,ia,ip,tag,i,j;
    
  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  ia   = iarf[iface][iedg];
  ip   = idir[iface][inxt2[iedg]];
  nump = pt->v[ip];
  numq = pt->v[idir[iface][iprv2[iedg]]];
  p0   = &mesh->point[nump];
  assert(p0->tag & MG_BDY);
  assert(p0->xp);

  /* collect triangles and tetras around ip */
  boulesurfvolp(mesh,k,ip,iface,listv,&ilistv,lists,&ilists);
  
  /* prevent collapse in case surface ball has 3 triangles */
  if ( ilists <= 2 )  return(0);

  /* Surfacic ball is enumerated with first tet having (pq) as edge n° iprv2[ip] on face iopp */
  startedgsurfball(mesh,nump,numq,lists,ilists);

  /* check tetra quality */
	calold = calnew = DBL_MAX;
  for (l=0; l<ilistv; l++) {
    iel = listv[l] / 4;
    ipp = listv[l] % 4;
    pt  = &mesh->tetra[iel];
    /* do not scan shell */
		for (j=0; j<4; j++)  if ( pt->v[j] == numq )  break;
		if ( j < 4 )  continue;
    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[ipp] = numq;

    calold = MG_MIN(calold,orcal(mesh,iel));
    calnew = MG_MIN(calnew,orcal(mesh,0));
  }
  if ( calold < NULKAL && calnew <= calold )  return(0);
  else if ( calnew < 0.3*calold )  return(0);

  /* analyze surfacic ball of p */
  for (l=1; l<ilists-1; l++) {
    iel  = lists[l] / 4;
    iopp = lists[l] % 4;
    pt   = &mesh->tetra[iel];
	  assert(pt->xt);

    /* retrieve vertex in tetra */
	  for (ip=0; ip<4; ip++)
      if ( pt->v[ip] == nump )  break;
    assert(ip<4);        

    memcpy(pt0,pt,sizeof(Tetra));
    pt0->v[ip] = numq;

	  if ( !norface(mesh,iel,iopp,ncurold) )  return(0);
    if ( !norface(mesh,0,iopp,ncurnew) )    return(0);

    /* check normal flipping */  
    ps = ncurold[0]*ncurnew[0] + ncurold[1]*ncurnew[1] + ncurold[2]*ncurnew[2];
    if ( ps < 0.0 )  return(0);

    /* check normal deviation */
    if ( l > 1 ) {
      ia = idirinv[iopp][ip]; /* index of p in tria iopp */
      ia = iprv2[ia];         /* edge between l-1 and l, in local num of tria */
      ia = iarf[iopp][ia];    /* edge between l-1 and l in local num of tetra */
    
	    hGet(&mesh->htab,pt->v[iare[ia][0]],pt->v[iare[ia][1]],&ref,&tag);
      if ( !(tag & MG_GEO) ) {
        devold = nprvold[0]*ncurold[0] + nprvold[1]*ncurold[1] + nprvold[2]*ncurold[2];
        devnew = nprvnew[0]*ncurnew[0] + nprvnew[1]*ncurnew[1] + nprvnew[2]*ncurnew[2];
        if ( devold < ANGEDG ) {
          if ( devnew < devold )  return(0);
        }
        else if ( devnew < ANGEDG )  return(0);
      }
    }

    /* check Hausdorff distance to geometric support */ 
	  tet2tri(mesh,iel,iopp,&tt);
	  if ( l == 1 ) {
	    for (i=0; i<3; i++) {
	      if ( tt.v[i] == nump )  break;
	    }
	    assert(i<3);
	    /* Index of the third point of the first collapsed triangle */
	    i  = inxt2[i];
  	  ia = inxt2[i];
	    hGet(&mesh->htab,tt.v[i],numq,&ref,&tag);
	    tt.tag[ia] = MG_MAX(tt.tag[ia],tag);
	  }
    else if ( l == ilists-2 ) {
	    for (i=0; i<3; i++) {
	      if ( tt.v[i] == nump )  break;
	    }
	    assert(i<3);
	    /* Index of the third point of the first collapsed triangle */
	    i  = iprv2[i];
	    ia = iprv2[i];
	    hGet(&mesh->htab,tt.v[i],numq,&ref,&tag);
	    tt.tag[ia] = MG_MAX(tt.tag[ia],tag);
 	  }
	 
	  for (i=0; i<3; i++) {
	    if ( tt.v[i] == nump )  break;
    }
	  assert(i<3);
		tt.v[i] = numq;
  	if ( chkedg(mesh,&tt) )  return(0);

    memcpy(nprvold,ncurold,3*sizeof(double));
    memcpy(nprvnew,ncurnew,3*sizeof(double));  
  }

  return(ilistv);
}

/* Collapse vertex p = list[0]%4 of tetra list[0]/4 over vertex indq of tetra list[0]/4. 
   Only physical tests (positive jacobian) are done (i.e. approximation of the surface,
   etc... must be performed outside). */
int colver(pMesh mesh,int *list,int ilist,char indq){
  pTetra          pt,pt1,pt0;
  pxTetra         pxt,pxt1;
  xTetra          xt;
  int             iel,jel,pel,qel,k,np,nq,*adja,op,oq,ref;
  unsigned char   ip,iq,i,j,voy,voyp,voyq;
  char            tag;

  iel = list[0] / 4;
  ip  = list[0] % 4;
  pt  = &mesh->tetra[iel];
  pt0 = &mesh->tetra[0];
  np  = pt->v[ip];
  nq  = pt->v[indq];
                                                              
  /* Mark elements of the shell of edge (pq) */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 4;
    i   = list[k] % 4;
    pt  = &mesh->tetra[iel];  
    for (j=0; j<3; j++) {
      i = inxt3[i];
      if ( pt->v[i] == nq ) {
        list[k] *= -1;
        break;
      }
    }
  }

  /* avoid recreating existing elt */
  for (k=0; k<ilist; k++) {
    if ( list[k] < 0 )  continue;
    iel = list[k] / 4;
    ip  = list[k] % 4;
    pt  = &mesh->tetra[iel];
    adja = &mesh->adja[4*(iel-1)+1];
    jel  = adja[ip] / 4;
    voy  = adja[ip] % 4;
    if ( !jel )  continue;
    pt = &mesh->tetra[jel];
    if ( pt->v[voy] == nq )  return(0);  
  }
   
  /* deal with the shell of edge (pq) and the implied updates */
  for (k=0; k<ilist; k++) {
    if ( list[k] > 0 )  continue;
    iel = (-list[k]) / 4;
    ip  = (-list[k]) % 4;
    pt  = &mesh->tetra[iel];
    iq  = ip;
    for (j=0; j<3; j++) {
      iq = inxt3[iq];
      if ( pt->v[iq] == nq )  break;
    }
    assert(j<3);

    adja = &mesh->adja[4*(iel-1)+1];
    /* pel = neighbour of iel that belongs to ball of p \setminus shell, same for qel */
    pel  = adja[iq] / 4;
    voyp = adja[iq] % 4;
    qel  = adja[ip] / 4;
    voyq = adja[ip] % 4;
		op = 0;
    if ( pel ) {
      pt1 = &mesh->tetra[pel];
      op  = pt1->v[voyp];
    }
    if ( qel ) {
      pt1 = &mesh->tetra[qel];
      oq  = pt1->v[voyq];
      //assert(op != oq);
    }

    /* Update adjacency relations */
    if ( pel ) {
      adja = &mesh->adja[4*(pel-1)+1];
      adja[voyp] = 4*qel+voyq; 
    }
    if ( qel ) {
      adja = &mesh->adja[4*(qel-1)+1];
      adja[voyq] = 4*pel+voyp; 
    }
    
    /* Update references for edges (pa)->(qa) when pqa is a face of the mesh */    
    for (j=0; j<4; j++) {
      if ( j == ip || j == iq )  continue;
      hPop(&mesh->htab,np,pt->v[j],&ref,&tag);
      if( tag || ref )
        hEdge(&mesh->htab,nq,pt->v[j],ref,tag);
    }
    
    /* Update references for faces (one in pel) ; 
       possibly, creation of a new field pxt for pel must be carried out */
    if ( pel ) {
      pt1 = &mesh->tetra[pel];
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        if ( pt1->xt > 0 ) {    
          pxt1 = &mesh->xtetra[pt1->xt];
          pxt1->ref[voyp] = pxt->ref[ip];
          pxt1->ftag[voyp] = pxt->ftag[ip];
        }
        else {   
          pxt1 = &xt;
          memset(pxt1,0,sizeof(xTetra));
          pxt1->ref[voyp] = pxt->ref[ip];
          pxt1->ftag[voyp] = pxt->ftag[ip];
          /* Recover the already used place by pxt */
          pt1->xt = pt->xt;
          memcpy(pxt,pxt1,sizeof(xTetra));
        }
      }
      else {
        /* Only the values corresponding to pt become 0 */
        if ( pt1->xt > 0 ) {
          pxt1 = &mesh->xtetra[pt1->xt];
          pxt1->ref[voyp]  = 0;
          pxt1->ftag[voyp] = 0;
        }
      }

      if ( qel ) {
        pt1 = &mesh->tetra[qel];
        if ( pt->xt ) {
          pxt = &mesh->xtetra[pt->xt];
          if ( pt1->xt > 0 ) {    
            pxt1 = &mesh->xtetra[pt1->xt];
            pxt1->ref[voyq] = MG_MAX(pxt1->ref[voyq],pxt->ref[iq]);
            pxt1->ftag[voyq] = (pxt1->ftag[voyq] | pxt->ftag[iq]);
          }
          else {   
            pxt1 = &xt;
            memset(pxt1,0,sizeof(xTetra));
            pxt1->ref[voyp] = pxt->ref[ip];
            pxt1->ftag[voyp] = pxt->ftag[ip];
            /* Create new field xt */
            mesh->xt++;
            pt1->xt = mesh->xt;
            memcpy(pxt,pxt1,sizeof(xTetra));
          }
        }
        else {
          /* Only the values corresponding to pt become 0 */
          if ( pt1->xt > 0 ) {
            pxt1 = &mesh->xtetra[pt1->xt];
            pxt1->ref[voyq]  = 0;
            pxt1->ftag[voyq] = 0;
          }
        }
      }
    }
    else {
      assert(pt->xt);
      pxt = &mesh->xtetra[pt->xt];
      if ( qel ) {
        pt1 = &mesh->tetra[qel];
        if ( pt1->xt > 0 ) {    
          pxt1 = &mesh->xtetra[pt1->xt];
          pxt1->ref[voyq] = pxt->ref[iq];
          pxt1->ftag[voyq] = pxt->ftag[iq];
        }
        else {   
          pxt1 = &xt;
          memset(pxt1,0,sizeof(xTetra));
          pxt1->ref[voyq] = pxt->ref[iq];
          pxt1->ftag[voyq] = pxt->ftag[iq]; 
          /* Recover the already used place by pxt */
          pt1->xt = pt->xt;
          memcpy(pxt,pxt1,sizeof(xTetra));
        }
      }
    }
    delElt(mesh,iel);
  }

  /* Update vertices coordinates for elements that do not belong to the shell of (pq) */
  for (k=0; k<ilist;  k++) {
    if ( list[k] < 0 )  continue;
    iel = list[k] / 4;
    ip  = list[k] % 4;
    pt  = &mesh->tetra[iel];
    for (j=0; j<4; j++) {
      if ( j==ip )  continue;
      hPop(&mesh->htab,np,pt->v[j],&ref,&tag);
      if ( tag || ref )
        hEdge(&mesh->htab,nq,pt->v[j],ref,tag);
    }
    pt->v[ip] = nq;
  }

  if ( mesh->point[np].tag & MG_BDY )
    hPop(&mesh->htab,np,nq,&ref,&tag);

  delPt(mesh,np);  
  return(1);
}