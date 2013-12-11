#include "mmg3d.h"

extern Info  info;

/** read mesh data */
int loadMesh(pMesh mesh) {
  pTetra       pt;
  pTria        pt1;
  pEdge        pa;
  pPoint       ppt;
  float        fp1,fp2,fp3;
  int          i,k,inm,ia,np,nr,nre,nc,aux,nt,ref,v[3],na,*ina,nereq;
  char        *ptr,*name,data[128];

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
  mesh->nt  = GmfStatKwd(inm,GmfTriangles);
  mesh->nei = GmfStatKwd(inm,GmfTetrahedra);
  mesh->nai = GmfStatKwd(inm,GmfEdges);
  if ( !mesh->npi || !mesh->nei ) {
    fprintf(stdout,"  ** MISSING DATA. Exit program.\n");
    return(0);
  }
  /* memory allocation */
  mesh->np = mesh->npi;
  mesh->ne = mesh->nei;
  mesh->na = mesh->nai;
  if ( !zaldy(mesh) )  return(0);

  if ( (mesh->np > mesh->npmax) || (mesh->nt > mesh->ntmax) ||
       (mesh->ne > mesh->nemax) ) {
    fprintf(stdout,"  ## Error: lack of memory\n");
    fprintf(stdout,"  ## Increase the allocated ");
    fprintf(stdout,"memory with the -m option.\n");
    return(0);
  }
  /* read mesh vertices */
  GmfGotoKwd(inm,GmfVertices);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if (mesh->ver == GmfFloat) {
      GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
      ppt->c[0] = fp1;
      ppt->c[1] = fp2;
      ppt->c[2] = fp3;
    } else {
      GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
    }
    ppt->tag  = MG_NUL;
  }
  /* get required vertices */
  np = GmfStatKwd(inm,GmfRequiredVertices);
  if ( np ) {
    GmfGotoKwd(inm,GmfRequiredVertices);
    for (k=1; k<=np; k++) {
      GmfGetLin(inm,GmfRequiredVertices,&i);
      assert(i <= mesh->np);
      ppt = &mesh->point[i];
      ppt->tag |= MG_REQ;
    }
  }
  /* get corners */
  nc = GmfStatKwd(inm,GmfCorners);
  if ( nc ) {
    GmfGotoKwd(inm,GmfCorners);
    for (k=1; k<=nc; k++) {
      GmfGetLin(inm,GmfCorners,&i);
      assert(i <= mesh->np);
      ppt = &mesh->point[i];
      ppt->tag |= MG_CRN;
    }
  }

  /* read mesh triangles */
  nt = 0;
  if ( mesh->nt ) {
    /* Skip triangles with negative refs */
    if( info.iso ) {
      GmfGotoKwd(inm,GmfTriangles);
      nt = mesh->nt;
      mesh->nt = 0;
      ina = (int*)calloc(nt+1,sizeof(int));
      if ( !ina ) {
        perror("  ## Memory problem: calloc");
        exit(EXIT_FAILURE);
      }
      for (k=1; k<=nt; k++) {
        GmfGetLin(inm,GmfTriangles,&v[0],&v[1],&v[2],&ref);
        if( abs(ref) != MG_ISO ) {
          pt1 = &mesh->tria[++mesh->nt];
          pt1->v[0] = v[0];
          pt1->v[1] = v[1];
          pt1->v[2] = v[2];
          pt1->ref = ref;
          ina[k]=mesh->nt;
        }
      }
      if( !mesh->nt ){
        free(mesh->tria);
        mesh->tria=NULL;
      }
      else if ( mesh->nt < nt )
        mesh->tria = (pTria)realloc(mesh->tria,(mesh->nt+1)*sizeof(Tria));
    }
    else {
      GmfGotoKwd(inm,GmfTriangles);
      for (k=1; k<=mesh->nt; k++) {
        pt1 = &mesh->tria[k];
        GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
      }
    }
    /* get required triangles */
    nt = GmfStatKwd(inm,GmfRequiredTriangles);
    if ( nt ) {
      GmfGotoKwd(inm,GmfRequiredTriangles);
      for (k=1; k<=nt; k++) {
        GmfGetLin(inm,GmfRequiredTriangles,&i);
        assert(i <= mesh->nt);
        if( info.iso ){
          if( ina[i] == 0 ) continue;
          else {
            pt1 = &mesh->tria[ina[i]];
            pt1->tag[0] |= MG_REQ;
            pt1->tag[1] |= MG_REQ;
            pt1->tag[2] |= MG_REQ;
          }
        }
        else{
          pt1 = &mesh->tria[i];
          pt1->tag[0] |= MG_REQ;
          pt1->tag[1] |= MG_REQ;
          pt1->tag[2] |= MG_REQ;
        }

      }
    }
    if ( info.iso ) {
      free(ina);
      ina=NULL;
    }
  }

  /* read mesh edges */
  nr = nre = 0;
  if ( mesh->na ) {
    na = mesh->na;
    if (info.iso ) {
      mesh->na = 0;
      ina = (int*)calloc(na+1,sizeof(int));
      if ( !ina ) {
        perror("  ## Memory problem: calloc");
        exit(EXIT_FAILURE);
      }
    }

    GmfGotoKwd(inm,GmfEdges);

    for (k=1; k<=na; k++) {
      pa = &mesh->edge[k];
      GmfGetLin(inm,GmfEdges,&pa->a,&pa->b,&pa->ref);
      pa->tag |= MG_REF;
      if ( info.iso ) {
        if( abs(pa->ref) != MG_ISO ) {
          ++mesh->na;
          pa->ref = abs(pa->ref);
          memmove(&mesh->edge[mesh->na],&mesh->edge[k],sizeof(Edge));
          ina[k] = mesh->na;
        }
      }
    }

    /* get ridges */
    nr = GmfStatKwd(inm,GmfRidges);
    if ( nr ) {
      GmfGotoKwd(inm,GmfRidges);
      for (k=1; k<=nr; k++) {
        GmfGetLin(inm,GmfRidges,&ia);
        assert(ia <= na);
        if( info.iso ){
          if( ina[ia] == 0 )
						continue;
          else {
            pa = &mesh->edge[ina[ia]];
            pa->tag |= MG_GEO;
          }
        }
        else{
          pa = &mesh->edge[ia];
          pa->tag |= MG_GEO;
        }
      }
    }
    /* get required edges */
    nre = GmfStatKwd(inm,GmfRequiredEdges);
    if ( nre ) {
      GmfGotoKwd(inm,GmfRequiredEdges);
      for (k=1; k<=nre; k++) {
        GmfGetLin(inm,GmfRequiredEdges,&ia);
        assert(ia <= na);
        if( info.iso ){
          if( ina[ia] == 0 ) continue;
          else {
            pa = &mesh->edge[ina[ia]];
            pa->tag |= MG_REQ;
          }
        }
        else{
          pa = &mesh->edge[ia];
          pa->tag |= MG_REQ;
        }

      }
    }
    if (info.iso ) {
      free(ina);
      ina = NULL;
    }
  }

  /* read mesh tetrahedra */
  GmfGotoKwd(inm,GmfTetrahedra);
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&pt->ref);
    pt->qual = orcal(mesh,k);
    pt->mark = 0;
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pt->v[i]];
      ppt->tag &= ~MG_NUL;
    }

    if ( info.iso )  pt->ref = 0;

    /* Possibly switch 2 vertices number so that each tet is positively oriented */
    if ( orvol(mesh->point,pt->v) < 0.0 ) {
      aux = pt->v[2];
      pt->v[2] = pt->v[3];
      pt->v[3] = aux;
    }
  }
  /* get required tetrahedra */
  nereq = GmfStatKwd(inm,GmfRequiredTetrahedra);
  if ( nereq ) {
    GmfGotoKwd(inm,GmfRequiredTetrahedra);
    for (k=1; k<=nereq; k++) {
      GmfGetLin(inm,GmfRequiredTetrahedra,&i);
      assert(i <= mesh->ne);
      pt = &mesh->tetra[i];
      pt->tag |= MG_REQ;
    }
  }


  /* stats */
  if ( abs(info.imprim) > 4 ) {
    fprintf(stdout,"     NUMBER OF VERTICES     %8d / %8d\n",mesh->np,mesh->npmax);
    if ( mesh->na ) {
      fprintf(stdout,"     NUMBER OF EDGES        %8d\n",mesh->na);
      if ( nr )
        fprintf(stdout,"     NUMBER OF RIDGES        %8d\n",nr);
        }
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES    %8d\n",mesh->nt);
    fprintf(stdout,"     NUMBER OF ELEMENTS     %8d / %8d\n",mesh->ne,mesh->nemax);

    if ( np || nre || nt || nereq ) {
      fprintf(stdout,"     NUMBER OF REQUIRED ENTITIES: \n");
      if ( np )
        fprintf(stdout,"                  VERTICES    %8d \n",np);
      if ( nre )
        fprintf(stdout,"                  EDGES       %8d \n",nre);
      if ( nt )
        fprintf(stdout,"                  TRIANGLES   %8d \n",nt);
      if ( nereq )
        fprintf(stdout,"                  TETRAHEDRAS %8d \n",nereq);
    }
    if(nc) fprintf(stdout,"     NUMBER OF CORNERS        %8d \n",nc);
  }
  GmfCloseMesh(inm);
  return(1);
}

/** Save mesh data */
int saveMesh(pMesh mesh) {
  pPoint       ppt;
  pTetra       pt;
  pTria        ptt;
  xPoint      *pxp;
  hgeom       *ph;
  int          k,i,na,nc,np,ne,nn,nr,nre,nedreq,ntreq,nt,outm,nereq;
  char         data[128];

  mesh->ver = GmfDouble;
  strcpy(data,mesh->nameout);
  if ( !(outm = GmfOpenMesh(data,GmfWrite,mesh->ver,mesh->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* vertices */
  np = nc = na = nr = nre = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      ppt->tmp = ++np;
      if ( ppt->tag & MG_CRN )  nc++;
      if ( ppt->tag & MG_REQ )  nre++;
    }
  }
  GmfSetKwd(outm,GmfVertices,np);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) )
      GmfSetLin(outm,GmfVertices,ppt->c[0],ppt->c[1],ppt->c[2],ppt->ref);
  }

  /* corners+required */
  if ( nc ) {
    GmfSetKwd(outm,GmfCorners,nc);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && ppt->tag & MG_CRN )
        GmfSetLin(outm,GmfCorners,ppt->tmp);
    }
  }
  if ( nre ) {
    GmfSetKwd(outm,GmfRequiredVertices,nre);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && ppt->tag & MG_REQ )
        GmfSetLin(outm,GmfRequiredVertices,ppt->tmp);
    }
  }

  nn = nt = 0;
  if ( mesh->xp ) {
    /* Count tangents and normals */
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) || MG_SIN(ppt->tag) )  continue;
      else if ( (ppt->tag & MG_BDY)
                && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM)) )
        nn++;
      if ( MG_EDG(ppt->tag) || (ppt->tag & MG_NOM) ) nt++;
    }

    /* write normals */
    GmfSetKwd(outm,GmfNormals,nn);

    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) || MG_SIN(ppt->tag) )  continue;
      else if ( (ppt->tag & MG_BDY)
                && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM)) ) {
        pxp = &mesh->xpoint[ppt->xp];
        GmfSetLin(outm,GmfNormals,pxp->n1[0],pxp->n1[1],pxp->n1[2]);
      }
    }

    GmfSetKwd(outm,GmfNormalAtVertices,nn);
    nn = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) || MG_SIN(ppt->tag) )  continue;
      else if ( (ppt->tag & MG_BDY)
                && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM)) )
        GmfSetLin(outm,GmfNormalAtVertices,ppt->tmp,++nn);
    }

    if ( nt ) {
      /* Write tangents */
      GmfSetKwd(outm,GmfTangents,nt);
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || MG_SIN(ppt->tag) )  continue;
        else if ( MG_EDG(ppt->tag) || (ppt->tag & MG_NOM) ) {
          pxp = &mesh->xpoint[ppt->xp];
          GmfSetLin(outm,GmfTangents,pxp->t[0],pxp->t[1],pxp->t[2]);
        }
      }
      GmfSetKwd(outm,GmfTangentAtVertices,nt);
      nt = 0;
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || MG_SIN(ppt->tag) )  continue;
        else if ( MG_EDG(ppt->tag) || (ppt->tag & MG_NOM) )
          GmfSetLin(outm,GmfTangentAtVertices,ppt->tmp,++nt);
      }
    }
  }
  free(mesh->xpoint);
  mesh->xpoint = NULL;
  mesh->xp = 0;

  /* boundary mesh */
  /* tria + required tria */
  mesh->nt = ntreq = 0;
  if ( mesh->tria){
    free(mesh->tria);
    mesh->tria=NULL;
  }
  chkNumberOfTri(mesh);
  if ( bdryTria(mesh) ) {
    GmfSetKwd(outm,GmfTriangles,mesh->nt);
    for (k=1; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];
      if ( ptt->tag[0] & MG_REQ && ptt->tag[1] & MG_REQ && ptt->tag[2] & MG_REQ )  ntreq++;
      GmfSetLin(outm,GmfTriangles,mesh->point[ptt->v[0]].tmp,mesh->point[ptt->v[1]].tmp, \
                mesh->point[ptt->v[2]].tmp,ptt->ref);
    }
    if ( ntreq ) {
      GmfSetKwd(outm,GmfRequiredTriangles,ntreq);
      for (k=0; k<=mesh->nt; k++) {
        ptt = &mesh->tria[k];
        if ( ptt->tag[0] & MG_REQ && ptt->tag[1] & MG_REQ && ptt->tag[2] & MG_REQ )
          GmfSetLin(outm,GmfRequiredTriangles,k);
      }
    }
    free(mesh->adjt);
    mesh->adjt=NULL;
    free(mesh->adja);
    mesh->adja = NULL;

    /* build hash table for edges */
    if ( mesh->htab.geom ) {
      free(mesh->htab.geom);
      mesh->htab.geom=NULL;
    }
    /* in the wost case (all edges are marked), we will have around 1 edge per *
     * triangle (we count edges only one time) */
    na = nr = nedreq = 0;
    if ( hNew(&mesh->htab,mesh->nt,3*(mesh->nt),0) ) {
      for (k=1; k<=mesh->ne; k++) {
        pt   = &mesh->tetra[k];
        if ( MG_EOK(pt) &&  pt->xt ) {
          for (i=0; i<6; i++) {
            if ( mesh->xtetra[pt->xt].edg[i] ||
                 ( MG_EDG(mesh->xtetra[pt->xt].tag[i] ) ||
                   (mesh->xtetra[pt->xt].tag[i] & MG_REQ) ) )
              hEdge(&mesh->htab,pt->v[iare[i][0]],pt->v[iare[i][1]],
                    mesh->xtetra[pt->xt].edg[i],mesh->xtetra[pt->xt].tag[i]);
          }
        }
      }
      /* edges + ridges + required edges */
      for (k=0; k<=mesh->htab.max; k++) {
        ph = &mesh->htab.geom[k];
        if ( !ph->a )  continue;
        na++;
        if ( ph->tag & MG_GEO )  nr++;
        if ( ph->tag & MG_REQ )  nedreq++;
      }
      if ( na ) {
        GmfSetKwd(outm,GmfEdges,na);
        for (k=0; k<=mesh->htab.max; k++) {
          ph = &mesh->htab.geom[k];
          if ( !ph->a )  continue;
          GmfSetLin(outm,GmfEdges,mesh->point[ph->a].tmp,mesh->point[ph->b].tmp,ph->ref);
        }
        if ( nr ) {
          GmfSetKwd(outm,GmfRidges,nr);
          na = 0;
          for (k=0; k<=mesh->htab.max; k++) {
            ph = &mesh->htab.geom[k];
            if ( !ph->a )  continue;
            na++;
            if ( ph->tag & MG_GEO )  GmfSetLin(outm,GmfRidges,na);
          }
        }
        if ( nedreq ) {
          GmfSetKwd(outm,GmfRequiredEdges,nedreq);
          na = 0;
          for (k=0; k<=mesh->htab.max; k++) {
            ph = &mesh->htab.geom[k];
            if ( !ph->a )  continue;
            na++;
            if ( ph->tag & MG_REQ )  GmfSetLin(outm,GmfRequiredEdges,na);
          }
        }
      }
      //freeXTets(mesh);
      free(mesh->htab.geom);
      mesh->htab.geom = NULL;
    }
  }

  /* tetrahedra */
  ne = nereq = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    ne++;
    if ( pt->tag & MG_REQ ){
      nereq++;
    }
  }

  GmfSetKwd(outm,GmfTetrahedra,ne);
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( MG_EOK(pt) ) GmfSetLin(outm,GmfTetrahedra,mesh->point[pt->v[0]].tmp,mesh->point[pt->v[1]].tmp, \
                                mesh->point[pt->v[2]].tmp,mesh->point[pt->v[3]].tmp,pt->ref);
  }

  if ( nereq ) {
    GmfSetKwd(outm,GmfRequiredTetrahedra,nereq);
    ne = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;
      ne++;
      if ( pt->tag & MG_REQ ) GmfSetLin(outm,GmfRequiredTetrahedra,ne);
    }
  }

  if ( info.imprim ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",np,nc+nre);
    if ( na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d   RIDGES  %8d\n",na,nr);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",mesh->nt);
    fprintf(stdout,"     NUMBER OF ELEMENTS   %8d\n",mesh->ne);
  }

  GmfCloseMesh(outm);
  return(1);
}

/** load metric field */
int loadMet(pSol met) {
  double       dbuf[GmfMaxTyp];
  float        fbuf[GmfMaxTyp];
  int          k,inm,typtab[GmfMaxTyp];
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
  met->np = GmfStatKwd(inm,GmfSolAtVertices,&met->type,&met->size,typtab);
  if ( !met->np ) {
    fprintf(stdout,"  ** MISSING DATA.\n");
    return(1);
  }
  if ( (met->type != 1) || (typtab[0] != 1 && typtab[0] != 3) ) {
    fprintf(stdout,"  ** DATA IGNORED %d  %d\n",met->type,typtab[0]);
    met->np = met->npmax = 0;
    return(-1);
  }

  /* mem alloc */
  met->m = (double*)calloc(met->size*met->npmax+1,sizeof(double));
  if ( !met->m ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }

  /* read mesh solutions */
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
        met->m[k] = dbuf[0];
      }
    }
  }
  /* anisotropic metric */
  /*else {
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
    }*/

  GmfCloseMesh(inm);
  return(1);
}

/** write iso or aniso metric */
int saveMet(pMesh mesh,pSol met) {
  pPoint     ppt;
  double     dbuf[GmfMaxTyp]/*,tmp*/;
  int        k/*,i*/,np,outm,nbm,typtab[GmfMaxTyp];
  char      *ptr,data[128];

  if ( !met->m )  return(-1);
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
    if ( MG_VOK(ppt) )  np++;
  }

  /* write isotropic metric */
  if ( met->size == 1 ) {
    typtab[0] = 1;
    nbm = 1;
    GmfSetKwd(outm,GmfSolAtVertices,np,nbm,typtab);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) ) {
        dbuf[0] = met->m[k];
        GmfSetLin(outm,GmfSolAtVertices,dbuf);
      }
    }
  }
  /* write anisotropic metric */
  /*else {
    typtab[0] = 3;
    nbm = 1;
    GmfSetKwd(outm,GmfSolAtVertices,np,nbm,typtab);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) ) {
        for (i=0; i<met->size; i++)  dbuf[i] = met->m[met->size*(k)+1+i];
        tmp = dbuf[2];
        dbuf[2] = dbuf[3];
        dbuf[3] = tmp;
        GmfSetLin(outm,GmfSolAtVertices,dbuf);
      }
    }
    }*/
  GmfCloseMesh(outm);
  return(1);
}

#ifdef SINGUL
/** Read singul data. Here we suppose that the file contains the singularities *
 *  (corner, required, ridges....) */
int loadSingul(pSingul singul) {
  Mesh         mesh;
  pEdge        pa,pas;
  pPoint       ppt;
  psPoint      ppts;
  float        fp1,fp2,fp3;
  int          i,k,inm,nr,nre,nc,npr,na,ns;
  char         *ptr,data[128],*filein;

  filein = singul->namein;
  strcpy(data,filein);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if( !(inm = GmfOpenMesh(data,GmfRead,&mesh.ver,&mesh.dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if( !(inm = GmfOpenMesh(data,GmfRead,&mesh.ver,&mesh.dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
    }
  }
  else if( !(inm = GmfOpenMesh(data,GmfRead,&mesh.ver,&mesh.dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  mesh.np = GmfStatKwd(inm,GmfVertices);
  mesh.nt = GmfStatKwd(inm,GmfTriangles);
  mesh.na = GmfStatKwd(inm,GmfEdges);
  if ( !mesh.np ) {
    fprintf(stdout,"  ** MISSING DATA. Exit program.\n");
    return(0);
  }

  /* memory allocation */
  mesh.point = (pPoint)calloc(mesh.np+1,sizeof(Point));
  if ( !mesh.point ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }

  if ( mesh.nt ) {
    mesh.tria = (pTria)calloc(mesh.nt+1,sizeof(Tria));
    if ( !mesh.tria ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
  }
  if ( mesh.na ) {
    mesh.edge = (pEdge)calloc(mesh.na+1,sizeof(Edge));
    if ( !mesh.edge ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
  }

  /* find bounding box */
  for (i=0; i<mesh.dim; i++) {
    singul->min[i] =  1.e30;
    singul->max[i] = -1.e30;
  }

  /* read mesh vertices */
  GmfGotoKwd(inm,GmfVertices);
  for (k=1; k<=mesh.np; k++) {
    ppt = &mesh.point[k];
    if (mesh.ver == GmfFloat) {
      GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
      ppt->c[0] = fp1;
      ppt->c[1] = fp2;
      ppt->c[2] = fp3;
    } else {
      GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
    }
    for (i=0; i<mesh.dim; i++) {
      if ( ppt->c[i] > singul->max[i] )  singul->max[i] = ppt->c[i];
      if ( ppt->c[i] < singul->min[i] )  singul->min[i] = ppt->c[i];
    }
    ppt->tag  = MG_NOTAG;
  }

  /* fill singul */
  /* get required vertices */
  npr = GmfStatKwd(inm,GmfRequiredVertices);
  if ( npr ) {
    GmfGotoKwd(inm,GmfRequiredVertices);
    for (k=1; k<=npr; k++) {
      GmfGetLin(inm,GmfRequiredVertices,&i);
      assert(i <= mesh.np);
      ppt = &mesh.point[i];
      ppt->tag |= MG_REQ;
    }
  }
  /* get corners */
  nc = GmfStatKwd(inm,GmfCorners);
  if ( nc ) {
    GmfGotoKwd(inm,GmfCorners);
    for (k=1; k<=nc; k++) {
      GmfGetLin(inm,GmfCorners,&i);
      assert(i <= mesh.np);
      ppt = &mesh.point[i];
      if ( !MG_SIN(ppt->tag) ){
        npr++;
      }
      ppt->tag |= MG_CRN;
    }
  }

  /* read mesh edges */
  if ( mesh.na ) {
    GmfGotoKwd(inm,GmfEdges);

    for (k=1; k<=mesh.na; k++) {
      pa = &mesh.edge[k];
      GmfGetLin(inm,GmfEdges,&pa->a,&pa->b,&pa->ref);
      pa->tag = MG_NOTAG;
    }
  }

  /* get ridges */
  nr = GmfStatKwd(inm,GmfRidges);
  if ( nr ) {
    GmfGotoKwd(inm,GmfRidges);
    for (k=1; k<=nr; k++) {
      GmfGetLin(inm,GmfRidges,&i);
      assert(i <= mesh.na);
      pa = &mesh.edge[i];
      pa->tag |= MG_GEO;
      ppt = &mesh.point[pa->a];
      if ( !(ppt->tag & MG_GEO) ){
        ppt->tag |= MG_GEO;
        if ( !MG_SIN(ppt->tag) )  npr++;
      }
      ppt = &mesh.point[pa->b];
      if ( !(ppt->tag & MG_GEO) ){
        ppt->tag |= MG_GEO;
        if ( !MG_SIN(ppt->tag) )  npr++;
      }
    }
  }
  /* get required edges */
  nre = GmfStatKwd(inm,GmfRequiredEdges);
  na  = 0;
  if ( nre ) {
    GmfGotoKwd(inm,GmfRequiredEdges);
    for (k=1; k<=nre; k++) {
      GmfGetLin(inm,GmfRequiredEdges,&i);
      assert(i <= mesh.na);
      pa = &mesh.edge[i];
      if ( !(pa->tag & MG_GEO) ) na++;
      pa->tag |= MG_REQ;
      ppt = &mesh.point[pa->a];
      if ( !(ppt->tag & MG_REQ) ){
        ppt->tag |= MG_REQ;
        if ( !MG_SIN(ppt->tag) )  npr++;
      }
      ppt = &mesh.point[pa->b];
      if ( !(ppt->tag & MG_REQ) ){
        ppt->tag |= MG_REQ;
        if ( !MG_SIN(ppt->tag) )  npr++;
      }
    }
  }

  singul->ns = npr;
  ns = 1;
  if ( singul->ns ) {
    singul->point = (psPoint)calloc(singul->ns+1,sizeof(sPoint));
    if ( !singul->point ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    for ( k=1; k<=mesh.np; k++ ) {
      ppt = &mesh.point[k];
      if ( (ppt->tag & MG_REQ) || (ppt->tag & MG_GEO) ) {
        ppts = &singul->point[ns];
        ppts->c[0] = ppt->c[0];
        ppts->c[1] = ppt->c[1];
        ppts->c[2] = ppt->c[2];
        ppts->tag  = ppt->tag | MG_SGL;
        ppt->tmp   = ns;
        ns++;
      }
    }
  }


  singul->na = nr+na;
  na = 1;
  if ( singul->na ) {
    singul->edge = (pEdge)calloc(singul->na+1,sizeof(Edge));
    if ( !singul->edge ) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    for ( k=1; k<=mesh.na; k++ ) {
      pa = &mesh.edge[k];
      if ( (pa->tag & MG_REQ) || (pa->tag & MG_GEO) ) {
        pas = &singul->edge[na];
        pas->a = mesh.point[pa->a].tmp;
        pas->b = mesh.point[pa->b].tmp;
        pas->tag  = pa->tag | MG_SGL;
        na++;
      }
    }
  }

  /* stats */
  if ( singul->ns )
    fprintf(stdout,"     NUMBER OF REQUIRED VERTICES : %8d \n",singul->ns);
  if ( singul->na )
    fprintf(stdout,"     NUMBER OF REQUIRED EDGES    : %8d \n",singul->na);

  GmfCloseMesh(inm);

  /* memory free */
  free (mesh.point);
  mesh.point=NULL;
  if ( mesh.na ) {
    free(mesh.edge);
    mesh.edge=NULL;
  }
  if ( mesh.nt ) {
    free(mesh.tria);
    mesh.tria=NULL;
  }
  return(1);
}
#endif
