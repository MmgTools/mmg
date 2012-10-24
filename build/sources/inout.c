#include "mmg3d.h"

extern Info  info;

/* read mesh data */
int loadMesh(pMesh mesh) {
  pTetra       pt;
  pTria        pt1;
  pEdge        pa;
  pPoint       ppt;
  float        fp1,fp2,fp3;
  int          i,k,inm,ia,nr,aux,nt,ref,v[3],na,*ina;
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
  /* read mesh triangles */
  if ( mesh->nt ) {
    /* Skip triangles with negative refs */
    if( info.iso ){
      GmfGotoKwd(inm,GmfTriangles);
      nt = mesh->nt;
      mesh->nt = 0;
      for (k=1; k<=nt; k++) {
        GmfGetLin(inm,GmfTriangles,&v[0],&v[1],&v[2],&ref);
        if( ref >= 0 ) {
          pt1 = &mesh->tria[++mesh->nt];
          pt1->v[0] = v[0];
          pt1->v[1] = v[1];
          pt1->v[2] = v[2];
          pt1->ref = ref;
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
  }

  /* read mesh edges */
  nr = 0;
  if ( mesh->na ) {
    na = mesh->na;
    if (info.iso ) {
      mesh->na = 0;
      ina = (int*)calloc(na+1,sizeof(int));
    }

    GmfGotoKwd(inm,GmfEdges);

    for (k=1; k<=na; k++) {
      pa = &mesh->edge[k];
      GmfGetLin(inm,GmfEdges,&pa->a,&pa->b,&pa->ref);
      pa->tag |= MG_REF;
      if ( info.iso ) {
        if( pa->ref != MG_ISO ) {
          ++mesh->na;
          pa->ref = abs(pa->ref);
          memcpy(&mesh->edge[mesh->na],&mesh->edge[k],sizeof(Edge));
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
          if( ina[ia] == 0) continue;
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
    nr = GmfStatKwd(inm,GmfRequiredEdges);
    if ( nr ) {
      GmfGotoKwd(inm,GmfRequiredEdges);
      for (k=1; k<=nr; k++) {
        GmfGetLin(inm,GmfRequiredEdges,&ia);
        assert(ia <= na);
        if( info.iso ){
          if( ina[ia] == 0) continue;
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
      ina=NULL;
    }
  }

  /* read mesh tetrahedra */
  GmfGotoKwd(inm,GmfTetrahedra);
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&pt->ref);
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pt->v[i]];
      ppt->tag &= ~MG_NUL;
    }

    if (info.iso ) pt->ref = 0;

    /* Possibly switch 2 vertices number so that each tet is positively oriented */
    if ( orvol(mesh->point,pt->v) < 0.0 ) {
      aux = pt->v[2];
      pt->v[2] = pt->v[3];
      pt->v[3] = aux;
    }
  }

  /* stats */
  if ( abs(info.imprim) > 4 ) {
    fprintf(stdout,"     NUMBER OF VERTICES     %8d / %8d\n",mesh->np,mesh->npmax);
    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES/RIDGES %8d / %8d\n",mesh->na,nr);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES    %8d\n",mesh->nt);
    fprintf(stdout,"     NUMBER OF ELEMENTS     %8d / %8d\n",mesh->ne,mesh->nemax);
  }
  GmfCloseMesh(inm);
  return(1);
}

/* Save mesh data */
int saveMesh(pMesh mesh) {
  pPoint       ppt;
  pTetra       pt;
  pTria        ptt;
  xPoint      *pxp;
  hgeom       *ph;
  int          k,na,nc,np,ne,nn,nr,nre,nt,outm;
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

  /* boundary mesh */
  mesh->nt = 0;
  if ( mesh->tria){
    free(mesh->tria);
    mesh->tria=NULL;
  }
  if ( bdryTria(mesh) ) {
    GmfSetKwd(outm,GmfTriangles,mesh->nt);
    for (k=1; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];
      GmfSetLin(outm,GmfTriangles,mesh->point[ptt->v[0]].tmp,mesh->point[ptt->v[1]].tmp,\
                mesh->point[ptt->v[2]].tmp,ptt->ref);
    }
    free(mesh->adjt);
    free(mesh->tria);
    mesh->adjt=NULL;
    mesh->tria=NULL;

    /* edges + ridges */
    na = nr = 0;
    for (k=0; k<=mesh->htab.max; k++) {
      ph = &mesh->htab.geom[k];
      if ( !ph->a )  continue;
      na++;
      if ( ph->tag & MG_GEO )  nr++;
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
    }
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

  /* tetrahedra */
  ne = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( MG_EOK(pt) ) ne++;
  }
  GmfSetKwd(outm,GmfTetrahedra,ne);
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( MG_EOK(pt) ) GmfSetLin(outm,GmfTetrahedra,mesh->point[pt->v[0]].tmp,mesh->point[pt->v[1]].tmp, \
                                mesh->point[pt->v[2]].tmp,mesh->point[pt->v[3]].tmp,pt->ref);
  }

  /* write normals */
  nn = nt = 0;
  GmfSetKwd(outm,GmfNormals,mesh->xp);

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_SIN(ppt->tag) )  continue;
    else if ( MG_VOK(ppt) && (ppt->tag & MG_BDY) && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM))) {
      pxp = &mesh->xpoint[ppt->xp];
      GmfSetLin(outm,GmfNormals,pxp->n1[0],pxp->n1[1],pxp->n1[2]);
      nn++;
    }
    if ( MG_EDG(ppt->tag) ) nt++;
  }
  GmfSetKwd(outm,GmfNormalAtVertices,nn);
  nn = 0;

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_SIN(ppt->tag) )  continue;
    else if ( MG_VOK(ppt) && (ppt->tag & MG_BDY) && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM)) )
      GmfSetLin(outm,GmfNormalAtVertices,ppt->tmp,++nn);
  }

  if ( nt ) {
    GmfSetKwd(outm,GmfTangents,nt);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_SIN(ppt->tag) )  continue;
      else if ( MG_VOK(ppt) && (MG_EDG(ppt->tag) || (ppt->tag & MG_NOM) )) {
        pxp = &mesh->xpoint[ppt->xp];
        GmfSetLin(outm,GmfTangents,pxp->t[0],pxp->t[1],pxp->t[2]);
      }
    }
    GmfSetKwd(outm,GmfTangentAtVertices,nt);
    nt = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_SIN(ppt->tag) )  continue;
      else if ( MG_VOK(ppt) && (MG_EDG(ppt->tag) || (ppt->tag & MG_NOM) ) )
        GmfSetLin(outm,GmfTangentAtVertices,ppt->tmp,++nt);
    }
  }

  if ( info.imprim ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",np,nc+nre);
    if ( na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d   RIDGES  %8d\n",na,nr);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",mesh->nt);
    fprintf(stdout,"     NUMBER OF ELEMENTS   %8d\n",ne);
  }

  GmfCloseMesh(outm);
  return(1);
}

/* load metric field */
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
  assert(met->m);

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

/* write iso or aniso metric */
int saveMet(pMesh mesh,pSol met) {
  pPoint     ppt;
  double     dbuf[GmfMaxTyp],tmp;
  int        k,i,np,outm,nbm,typtab[GmfMaxTyp];
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
  else {
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
  }
  GmfCloseMesh(outm);
  return(1);
}
