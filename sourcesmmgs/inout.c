#include "mmgs.h"
#include <math.h>

extern Info   info;


int loadMesh(pMesh mesh) {
  pTria      pt1,pt2;
  pPoint     ppt;
	double    *norm,*n,dd;
  float      fp1,fp2,fp3;
  int        k,ia,nq,nc,nrv,nri,nr,inm,num,ip,idn,ng;
  char      *ptr,*name,i,data[256];

  name = mesh->namein;
  strcpy(data,name);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
    }
  }
  else if( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  mesh->npi = GmfStatKwd(inm,GmfVertices);
  mesh->nti = GmfStatKwd(inm,GmfTriangles);
  nq = GmfStatKwd(inm,GmfQuadrilaterals);
  if ( !mesh->npi || !mesh->nti ) {
    fprintf(stdout,"  ** MISSING DATA\n");
    return(0);
  }
  mesh->np = mesh->npi;
  mesh->nt = mesh->nti + 2*nq;

  /* mem alloc */
  if ( !zaldy(mesh) )  return(0);

  /* read vertices */
  GmfGotoKwd(inm,GmfVertices);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( mesh->ver == GmfFloat ) {
      GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
      ppt->c[0] = fp1;
      ppt->c[1] = fp2;
      ppt->c[2] = fp3;
    }
    else
      GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
    ppt->tag = MS_NUL;
  }

  /* read triangles and set seed */
  GmfGotoKwd(inm,GmfTriangles);
  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
    for (i=0; i<3; i++) {    
      ppt = &mesh->point[pt1->v[i]];
      ppt->tag &= ~MS_NUL;
    }
  }

  /* read quads */
  if ( nq > 0 ) {
    GmfGotoKwd(inm,GmfQuadrilaterals);
    for (k=1; k<=nq; k++) {
      mesh->nti++;
      pt1 = &mesh->tria[mesh->nti];
      mesh->nti++;
      pt2 = &mesh->tria[mesh->nti];
      GmfGetLin(inm,GmfQuadrilaterals,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt2->v[2],&pt1->ref);
      pt2->v[0] = pt1->v[0];
      pt2->v[1] = pt1->v[2];
      pt2->ref  = pt1->ref;
      for (i=0; i<3; i++) {    
        ppt = &mesh->point[pt1->v[i]];
        ppt->tag &= ~MS_NUL;
        ppt = &mesh->point[pt2->v[i]];
        ppt->tag &= ~MS_NUL;
      }
    }
    mesh->nt = mesh->nti; 
  }

  /* read corners */
  nc = GmfStatKwd(inm,GmfCorners);
  if ( nc ) {
    GmfGotoKwd(inm,GmfCorners);
    for (k=1; k<=nc; k++) {
      GmfGetLin(inm,GmfCorners,&num);
      mesh->point[num].tag |= MS_CRN;
    }
  }
  /* read required vertices */
  nrv = GmfStatKwd(inm,GmfRequiredVertices);
  if ( nrv ) {
    GmfGotoKwd(inm,GmfRequiredVertices);
    for (k=1; k<=nrv; k++) {
      GmfGetLin(inm,GmfRequiredVertices,&num);
      mesh->point[num].tag |= MS_REQ;
    }
  }

  mesh->na = GmfStatKwd(inm,GmfEdges);
  nr = nri = 0;
  if ( mesh->na ) {
    mesh->edge = (pEdge)calloc(mesh->na+1,sizeof(Edge));
    assert(mesh->edge);
    GmfGotoKwd(inm,GmfEdges);
    for (k=1; k<=mesh->na; k++) {
      GmfGetLin(inm,GmfEdges,&mesh->edge[k].a,&mesh->edge[k].b,&mesh->edge[k].ref);
      mesh->edge[k].tag |= MS_REF;
      mesh->point[mesh->edge[k].a].tag |= MS_REF;
      mesh->point[mesh->edge[k].b].tag |= MS_REF;
    }

    nri = GmfStatKwd(inm,GmfRidges);
    if ( nri ) {
      GmfGotoKwd(inm,GmfRidges);
      for (k=1; k<=nri; k++) {
        GmfGetLin(inm,GmfRidges,&ia);
        if ( ia > 0 && ia <= mesh->na )  mesh->edge[ia].tag |= MS_GEO;  
      }
    }
    nr = GmfStatKwd(inm,GmfRequiredEdges);
    if ( nr ) {
      GmfGotoKwd(inm,GmfRequiredEdges);
      for (k=1; k<=nr; k++) {
        GmfGetLin(inm,GmfRequiredEdges,&ia);
        if ( ia > 0 && ia <= mesh->na )   mesh->edge[ia].tag |= MS_REQ;  
      }
    }
    /* set tria edges tags */
    assignEdge(mesh);
  }

  /* read geometric entities */
  ng = GmfStatKwd(inm,GmfNormals);
  if ( ng > 0 ) {
		norm = (double*)calloc(3*ng+1,sizeof(double));
		assert(norm);

    GmfGotoKwd(inm,GmfNormals);
    for (k=1; k<=ng; k++) {
      n = &norm[3*(k-1)+1];
      if ( mesh->ver == GmfFloat ) {
        GmfGetLin(inm,GmfNormals,&fp1,&fp2,&fp3);
				n[0] = fp1;  n[1] = fp2;  n[2] = fp3;
			}
      else {
				GmfGetLin(inm,GmfNormals,&n[0],&n[1],&n[2]);
			}
      dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
      if ( dd > EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        n[0] *= dd;
        n[1] *= dd;
        n[2] *= dd;
      }
    }

    mesh->nc1 = GmfStatKwd(inm,GmfNormalAtVertices);
    GmfGotoKwd(inm,GmfNormalAtVertices);
    for (k=1; k<=mesh->nc1; k++) {
	    GmfGetLin(inm,GmfNormalAtVertices,&ip,&idn);
      if ( idn > 0 && ip < mesh->np+1 )
				memcpy(&mesh->point[ip].n,&norm[3*(idn-1)+1],3*sizeof(double));
	  }
		free(norm);
  }

  if ( abs(info.imprim) > 4 ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d / %8d   CORNERS/REQ. %d / %d\n",mesh->npi,mesh->npmax,nc,nrv);
    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d  RIDGES %6d\n",mesh->na,nri);
    fprintf(stdout,"     NUMBER OF TRIANGLES  %8d / %8d\n",mesh->nti,mesh->ntmax);
  }
  GmfCloseMesh(inm);
  return(1);
}

int saveMesh(pMesh mesh) {
  pPoint       ppt;
  pTria        pt;
  pEdge        edge;
  pGeom        go;
  int         *adja,k,jel,outm,np,nt,na,nc,ng,nn,nr,nre;
  char         data[128],i,i1,i2;

  edge = 0;
  mesh->ver = GmfDouble;
  strcpy(data,mesh->nameout);
  if ( !(outm = GmfOpenMesh(data,GmfWrite,mesh->ver,mesh->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* vertices */
  np = nc = ng = nn = nre = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->tmp = 0;
    if ( MS_VOK(ppt) ) {
      np++;
      ppt->tmp = np;
      if ( ppt->tag & MS_CRN )  nc++;
      if ( ppt->tag & MS_REQ )  nre++;
      if ( MS_EDG(ppt->tag) )   ng++;
    }
  }
    
  GmfSetKwd(outm,GmfVertices,np);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MS_VOK(ppt) ) {
      GmfSetLin(outm,GmfVertices,ppt->c[0],ppt->c[1],ppt->c[2],ppt->ref);
      if ( !(ppt->tag & MS_GEO) )  nn++;
    }
  }

  nt = na = nr = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( MS_EOK(pt) ) { 
      nt++;
      for (i=0; i<3; i++)
        if ( MS_EDG(pt->tag[i]) )  na++;
    }
  }

  /* memory alloc */
  if ( na ) {
    edge = (pEdge)calloc(na+1,sizeof(Edge));
    assert(edge);
  }

  /* write triangles */
  GmfSetKwd(outm,GmfTriangles,nt);
  na = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( MS_EOK(pt) ) {
      GmfSetLin(outm,GmfTriangles,mesh->point[pt->v[0]].tmp,mesh->point[pt->v[1]].tmp,
                                  mesh->point[pt->v[2]].tmp,pt->ref);
      
      for (i=0; i<3; i++) {
        if ( !MS_EDG(pt->tag[i]) )  continue;
        
        adja = &mesh->adja[3*(k-1)+1];
        jel  = adja[i] / 3;
        if ( !jel || jel > k ) {
          i1 = inxt[i];
          i2 = inxt[i1];
          na++;
          edge[na].a    = mesh->point[pt->v[i1]].tmp;
          edge[na].b    = mesh->point[pt->v[i2]].tmp;
          edge[na].ref  = pt->edg[i];
          edge[na].tag |= pt->tag[i];
          if ( pt->tag[i] & MS_GEO )  nr++;
        }
      }
    }
  }

  /* write corners */
  if ( nc ) {
    GmfSetKwd(outm,GmfCorners,nc);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MS_VOK(ppt) && ppt->tag & MS_CRN )
        GmfSetLin(outm,GmfCorners,ppt->tmp);
    }
  }
  if ( nre ) {
    GmfSetKwd(outm,GmfRequiredVertices,nre);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MS_VOK(ppt) && ppt->tag & MS_REQ )
        GmfSetLin(outm,GmfRequiredVertices,ppt->tmp);
    }
  }

  /* write edges, ridges */
  if ( na ) {
    GmfSetKwd(outm,GmfEdges,na);
    nre = 0;
    for (k=1; k<=na; k++) {
      GmfSetLin(outm,GmfEdges,edge[k].a,edge[k].b,edge[k].ref);
      if ( edge[k].tag & MS_REQ )  nre++;
    }
    /* ridges */
    if ( nr ) {
      GmfSetKwd(outm,GmfRidges,nr);
      for (k=1; k<=na; k++) {
        if ( edge[k].tag & MS_GEO )  GmfSetLin(outm,GmfRidges,k);
      }
    }
    if ( nre ) {
      GmfSetKwd(outm,GmfRequiredEdges,nre);
      for (k=1; k<=na; k++)
        if ( edge[k].tag & MS_REQ )  GmfSetLin(outm,GmfRequiredEdges,k);
    }
    free(edge);
  }

  /* write normals */
  if ( nn ) {
    GmfSetKwd(outm,GmfNormals,nn);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MS_VOK(ppt) )  continue;
      else if ( !(ppt->tag & MS_GEO) ) {
        if ( ppt->tag & MS_REF ) {
          go = &mesh->geom[ppt->ig];
          GmfSetLin(outm,GmfNormals,go->n1[0],go->n1[1],go->n1[2]);
        }
        else
          GmfSetLin(outm,GmfNormals,ppt->n[0],ppt->n[1],ppt->n[2]);
      }
    }
    GmfSetKwd(outm,GmfNormalAtVertices,nn);
    nn = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MS_VOK(ppt) && !(ppt->tag & MS_GEO) )
        GmfSetLin(outm,GmfNormalAtVertices,ppt->tmp,++nn);
    }
/*
		nn = 0;
		for (k=1; k<=mesh->nt; k++) {
			pt = &mesh->tria[k];
			if ( !MS_EOK(pt) )  continue;
			for (i=0; i<3; i++) {
				ppt = &mesh->point[pt->v[i]];
				if ( ppt->tag & MS_GEO )  nn++;
			}
		}
    GmfSetKwd(outm,GmfNormalAtTriangleVertices,nn);
		for (k=1; k<=mesh->nt; k++) {
			pt = &mesh->tria[k];
			if ( !MS_EOK(pt) )  continue;
			for (i=0; i<3; i++) {
				ppt = &mesh->point[pt->v[i]];
				if ( ppt->tag & MS_GEO )  nn++;
			}
*/
  }

  /* write tangents */
  if ( ng ) {
    GmfSetKwd(outm,GmfTangents,ng);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MS_VOK(ppt) && MS_EDG(ppt->tag) )
        GmfSetLin(outm,GmfTangents,ppt->n[0],ppt->n[1],ppt->n[2]);
    }
    GmfSetKwd(outm,GmfTangentAtVertices,ng);
    ng = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MS_VOK(ppt) && MS_EDG(ppt->tag) )
        GmfSetLin(outm,GmfTangentAtVertices,ppt->tmp,++ng);
    }
  }

  if ( abs(info.imprim) > 4 ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d  CORNERS    %6d\n",np,nc);
    if ( na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d  RIDGES     %6d\n",na,nr);
    fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",nt);
    if ( nn+ng )
      fprintf(stdout,"     NUMBER OF NORMALS    %8d  TANGENTS   %6d\n",nn,ng);
  }
  GmfCloseMesh(outm);
  return(1);
}

/* load metric field */
int loadMet(pSol met) {
  double       tmpd,dbuf[GmfMaxTyp];
  float        tmpf,fbuf[GmfMaxTyp];
  int          i,k,inm,type,typtab[GmfMaxTyp];
  char        *ptr,data[128];

	if ( !met->namein )  return(0);
  strcpy(data,met->namein);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".solb");
  if (!(inm = GmfOpenMesh(data,GmfRead,&met->ver,&met->dim)) ) {
    ptr  = strstr(data,".sol");
    *ptr = '\0';
    strcat(data,".sol");
    if (!(inm = GmfOpenMesh(data,GmfRead,&met->ver,&met->dim)) ) {
      if ( info.imprim < 0 ) 
        fprintf(stderr,"  ** %s  NOT FOUND. USE DEFAULT METRIC.\n",data);
      return(-1);
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* read solution or metric */
  met->np = GmfStatKwd(inm,GmfSolAtVertices,&type,&met->size,typtab);
  if ( !met->np ) {
    fprintf(stdout,"  ** MISSING DATA.\n");
    return(1);
  }
  if ( (type != 1) || (typtab[0] != 1 && typtab[0] != 3) ) {
    fprintf(stdout,"  ** DATA IGNORED %d  %d\n",type,typtab[0]);
    met->np = met->npmax = 0;
    return(-1);
  }

  /* mem alloc */
  met->m = (double*)calloc(met->size*(met->npmax+1)+1,sizeof(double));
  assert(met->m);

  GmfGotoKwd(inm,GmfSolAtVertices);
  /* isotropic metric */
  if ( met->size == 1 ) {
		if ( met->ver == GmfFloat ) {
	    for (k=1; k<=met->np; k++) {
	      GmfGetLin(inm,GmfSolAtVertices,fbuf);
				met->m[k] = fbuf[0];
	    }
	  }
	  else {      
	    for (k=1; k<=met->np; k++) {
	      GmfGetLin(inm,GmfSolAtVertices,dbuf);
				met->m[k] = fbuf[0];
			}
	  }
	}
  /* anisotropic metric */
  else {
    if ( met->ver == GmfFloat ) {
	    for (k=1; k<=met->np; k++) {
	      GmfGetLin(inm,GmfSolAtVertices,fbuf);
			  tmpf    = fbuf[2];
			  fbuf[2] = fbuf[3];
			  fbuf[3] = tmpf;
        for (i=0; i<6; i++)  met->m[6*k+1+i] = fbuf[i];
      }
    }
    else {
	    for (k=1; k<=met->np; k++) {
        GmfGetLin(inm,GmfSolAtVertices,dbuf);
			  tmpd    = dbuf[2];
			  dbuf[2] = dbuf[3];
			  dbuf[3] = tmpd;
        for (i=0; i<met->size; i++)  met->m[6*k+1+i] = dbuf[i];
      }
	  }
  }
  GmfCloseMesh(inm);
  return(1);  
}

/* write iso or aniso metric */
int saveMet(pMesh mesh,pSol met) {
  pPoint     ppt;
  double     dbuf[GmfMaxTyp],tmp;
  int        k,i,np,outm,nbm,typtab[GmfMaxTyp];
  char      *ptr,data[128];

  met->ver = GmfDouble;
  strcpy(data,met->nameout);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  ptr = strstr(data,".sol");
  if ( !ptr )  strcat(data,".sol");
  if ( !(outm = GmfOpenMesh(data,GmfWrite,met->ver,met->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

	np = 0;
	for (k=1; k<=mesh->np; k++) {
		ppt = &mesh->point[k];
		if ( MS_VOK(ppt) )  np++;
	}

  /* write isotropic metric */
  if ( met->size == 1 ) {
		typtab[0] = 1;
		nbm = 1;
    GmfSetKwd(outm,GmfSolAtVertices,np,nbm,typtab);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MS_VOK(ppt) ) {
				dbuf[0] = met->m[k];
				GmfSetLin(outm,GmfSolAtVertices,dbuf);
	    }
	  }
	}
  /* write anisotropic metric */
	else {
    typtab[0] = 3;
		nbm = 1;
    GmfSetKwd(outm,GmfSolAtVertices,np,nbm,typtab);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MS_VOK(ppt) ) {
        for (i=0; i<met->size; i++)  dbuf[i] = met->m[met->size*(k)+1+i];
				tmp = dbuf[2];
				dbuf[2] = dbuf[3];
				dbuf[3] = tmp;
        GmfSetLin(outm,GmfSolAtVertices,dbuf);
      }
    }
  }
  GmfCloseMesh(outm);
  return(1);
}


