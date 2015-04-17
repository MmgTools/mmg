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
 * \file mmgs/inout.c
 * \brief Input / Output Functions.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"
#include <math.h>

#define sw 4
#define sd 8

int swapbin(int sbin)
{
  int inv;
  char *p_in = (char *) &sbin;
  char *p = (char *)&inv;


  p[0] = p_in[3];
  p[1] = p_in[2];
  p[2] = p_in[1];
  p[3] = p_in[0];

  return(inv);
  /*unsigned char c1, c2, c3, c4;

    c1 = sbin & 255;
    c2 = (sbin >> 8) & 255;
    c3 = (sbin >> 16) & 255;
    c4 = (sbin >> 24) & 255;

    return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;   */

}
float swapf(float sbin)
{
  float out;
  char *p_in = (char *) &sbin;
  char *p_out = (char *) &out;
  p_out[0] = p_in[3];
  p_out[1] = p_in[2];
  p_out[2] = p_in[1];
  p_out[3] = p_in[0];

  return(out);
}
double swapd(double sbin)
{
  float out;
  char *p_in = (char *) &sbin;
  char *p_out = (char *) &out;
  int i;

  for(i=0;i<8;i++)
  {
    p_out[i] = p_in[7-i];
  }
  //printf("CONVERTION DOUBLE\n");
  return(out);
}
int loadMesh(MMG5_pMesh mesh) {
  FILE        *inm;
  MMG5_pTria  pt1,pt2;
  MMG5_pPoint ppt;
  double      *norm,*n,dd;
  float       fc;
  int         i,k,ia,nq,nri,nr,ip,idn,ng;
  char        *ptr,*name,data[256],chaine[128];
  int         posnp,posnt,posne,posncor,posnq,posned,posnr;
  int         posnpreq,npreq,ntreq,posnormal,posnc1,posntreq;
  int         ncor,bin,iswp,nedreq,posnedreq,bdim,binch,bpos;

  posnp = posnt = posne = posncor = posnq = posntreq = 0;
  posned = posnr = posnpreq = posnc1 = npreq = 0;
  posnedreq = posnormal = 0;
  ncor = nr = ng = nedreq = nq = ntreq = 0;
  bin = 0;
  iswp = 0;
  mesh->np = mesh->nt = mesh->nti = mesh->npi = 0;

  name = mesh->namein;
  strcpy(data,name);
  ptr = strstr(data,".mesh");
    if ( !ptr ) {
    /* data contains the filename without extension */
    strcat(data,".meshb");
    if( !(inm = fopen(data,"rb")) ) {
      /* our file is not a .meshb file, try with .mesh ext */
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if( !(inm = fopen(data,"r")) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
      else {
        if ( !strstr(mesh->nameout,".mesh") ) {
          _MMG5_ADD_MEM(mesh,5*sizeof(char),"output file name",
                        printf("  Exit program.\n");
                        exit(EXIT_FAILURE));
          _MMG5_SAFE_REALLOC(mesh->nameout,strlen(mesh->nameout)+6,char,"output mesh name");
          strcat(mesh->nameout,".mesh");
        }
      }
    }
    else {
      bin = 1;
      if ( !strstr(mesh->nameout,".mesh") ) {
        _MMG5_ADD_MEM(mesh,6*sizeof(char),"input file name",
                      printf("  Exit program.\n");
                      exit(EXIT_FAILURE));
        _MMG5_SAFE_REALLOC(mesh->nameout,strlen(mesh->nameout)+7,char,"input file name");
        strcat(mesh->nameout,".meshb");
      }
    }
  }
  else {
    ptr = strstr(data,".meshb");
    if( !ptr ) {
      if( !(inm = fopen(data,"r")) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
    } else {
      bin = 1;
      if( !(inm = fopen(data,"rb")) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if (!bin) {
    strcpy(chaine,"D");
    while(fscanf(inm,"%s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
      if(!strncmp(chaine,"MeshVersionFormatted",strlen("MeshVersionFormatted"))) {
        fscanf(inm,"%d",&mesh->ver);
        continue;
      } else if(!strncmp(chaine,"Dimension",strlen("Dimension"))) {
        fscanf(inm,"%d",&mesh->dim);
        if(mesh->dim!=3) {
          fprintf(stdout,"BAD DIMENSION : %d\n",mesh->dim);
          return(0);
        }
        continue;
      } else if(!strncmp(chaine,"Vertices",strlen("Vertices"))) {
        fscanf(inm,"%d",&mesh->npi);
        posnp = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredVertices",strlen("RequiredVertices"))) {
        fscanf(inm,"%d",&npreq);
        posnpreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Triangles",strlen("Triangles"))) {
        fscanf(inm,"%d",&mesh->nti);
        posnt = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredTriangles",strlen("RequiredTriangles"))) {
        fscanf(inm,"%d",&ntreq);
        posntreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Quadrilaterals",strlen("Quadrilaterals"))) {
        fscanf(inm,"%d",&nq);
        posnq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Corners",strlen("Corners"))) {
        fscanf(inm,"%d",&ncor);
        posncor = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Edges",strlen("Edges"))) {
        fscanf(inm,"%d",&mesh->na);
        posned = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredEdges",strlen("RequiredEdges"))) {
        fscanf(inm,"%d",&nedreq);
        posnedreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Ridges",strlen("Ridges"))) {
        fscanf(inm,"%d",&nri);
        posnr = ftell(inm);
        continue;
      } else if(!ng && !strncmp(chaine,"Normals",strlen("Normals"))) {
        fscanf(inm,"%d",&ng);
        posnormal = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"NormalsAtVertices",strlen("NormalsAtVertices"))) {
        fscanf(inm,"%d",&mesh->nc1);
        posnc1 = ftell(inm);
        continue;
      }
    }
  } else { //binary file
    bdim = 0;
    fread(&mesh->ver,sw,1,inm);
    iswp=0;
    if(mesh->ver==16777216)
      iswp=1;
    else if(mesh->ver!=1) {
      fprintf(stdout,"BAD FILE ENCODING\n");
    }
    fread(&mesh->ver,sw,1,inm);
    if(iswp) mesh->ver = swapbin(mesh->ver);
    while(fread(&binch,sw,1,inm)!=0 && binch!=54 ) {
      if(iswp) binch=swapbin(binch);
      if(binch==54) break;
      if(!bdim && binch==3) {  //Dimension
        fread(&bdim,sw,1,inm);  //NulPos=>20
        if(iswp) bdim=swapbin(bdim);
        fread(&bdim,sw,1,inm);
        if(iswp) bdim=swapbin(bdim);
        mesh->dim = bdim;
        if(bdim!=3) {
          fprintf(stdout,"BAD SOL DIMENSION : %d\n",mesh->dim);
          exit(0);
          return(1);
        }
        continue;
      } else if(!mesh->npi && binch==4) {  //Vertices
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&mesh->npi,sw,1,inm);
        if(iswp) mesh->npi=swapbin(mesh->npi);
        posnp = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==15) {  //RequiredVertices
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&npreq,sw,1,inm);
        if(iswp) npreq=swapbin(npreq);
        posnpreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!mesh->nti && binch==6) {//Triangles
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&mesh->nti,sw,1,inm);
        if(iswp) mesh->nti=swapbin(mesh->nti);
        posnt = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==17) {  //RequiredTriangles
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&ntreq,sw,1,inm);
        if(iswp) ntreq=swapbin(ntreq);
        posntreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==7) {//Quadrilaterals
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&nq,sw,1,inm);
        if(iswp) nq=swapbin(nq);
        posnq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!ncor && binch==13) { //Corners
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&ncor,sw,1,inm);
        if(iswp) ncor=swapbin(ncor);
        posncor = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!mesh->na && binch==5) { //Edges
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&mesh->na,sw,1,inm);
        if(iswp) mesh->na=swapbin(mesh->na);
        posned = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==16) {  //RequiredEdges
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&nedreq,sw,1,inm);
        if(iswp) nedreq=swapbin(nedreq);
        posnedreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==14) {  //Ridges
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&nri,sw,1,inm);
        if(iswp) nri=swapbin(nri);
        posnr = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!ng && binch==60) {  //Normals
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&ng,sw,1,inm);
        if(iswp) ng=swapbin(ng);
        posnormal = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==20) {  //NormalsAtVertices
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        fread(&mesh->nc1,sw,1,inm);
        if(iswp) mesh->nc1=swapbin(mesh->nc1);
        posnc1 = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else {
        //printf("on traite ? %d\n",binch);
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=swapbin(bpos);
        //printf("on avance... Nulpos %d\n",bpos);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
      }
    }
  }

  if ( !mesh->npi || !mesh->nti ) {
    fprintf(stdout,"  ** MISSING DATA\n");
    return(0);
  }
  mesh->np = mesh->npi;
  mesh->nt = mesh->nti + 2*nq;

  /* mem alloc */
  if ( !zaldy(mesh) )  return(0);

  /* read vertices */

  rewind(inm);
  fseek(inm,posnp,SEEK_SET);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if (mesh->ver < 2) { /*float*/
      if (!bin) {
        for (i=0 ; i<3 ; i++) {
          fscanf(inm,"%f",&fc);
          ppt->c[i] = (double) fc;
        }
        fscanf(inm,"%d",&ppt->ref);
      } else {
        for (i=0 ; i<3 ; i++) {
          fread(&fc,sw,1,inm);
          if(iswp) fc=swapf(fc);
          ppt->c[i] = (double) fc;
        }
        fread(&ppt->ref,sw,1,inm);
        if(iswp) ppt->ref=swapbin(ppt->ref);
      }
    } else {
      if (!bin)
        fscanf(inm,"%lf %lf %lf %d",&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
      else {
        for (i=0 ; i<3 ; i++) {
          fread(&ppt->c[i],sd,1,inm);
          if(iswp) ppt->c[i]=swapd(ppt->c[i]);
        }
        fread(&ppt->ref,sw,1,inm);
        if(iswp) ppt->ref=swapbin(ppt->ref);
      }
    }
    ppt->tag = MG_NUL;
  }

  /* read triangles and set seed */
  rewind(inm);
  fseek(inm,posnt,SEEK_SET);
  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    if (!bin)
      fscanf(inm,"%d %d %d %d",&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
    else {
      for (i=0 ; i<3 ; i++) {
        fread(&pt1->v[i],sw,1,inm);
        if(iswp) pt1->v[i]=swapbin(pt1->v[i]);
      }
      fread(&pt1->ref,sw,1,inm);
      if(iswp) pt1->ref=swapbin(pt1->ref);
    }
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt1->v[i]];
      ppt->tag &= ~MG_NUL;
    }
  }
  /* read quads */
  if ( nq > 0 ) {
    rewind(inm);
    fseek(inm,posnq,SEEK_SET);

    for (k=1; k<=nq; k++) {
      mesh->nti++;
      pt1 = &mesh->tria[mesh->nti];
      mesh->nti++;
      pt2 = &mesh->tria[mesh->nti];

      if (!bin)
        fscanf(inm,"%d %d %d %d %d",&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt2->v[2],&pt1->ref);
      else {
        for (i=0 ; i<3 ; i++) {
          fread(&pt1->v[i],sw,1,inm);
          if(iswp) pt1->v[i]=swapbin(pt1->v[i]);
        }
        fread(&pt2->v[2],sw,1,inm);
        if(iswp) pt2->v[2]=swapbin(pt2->v[2]);
        fread(&pt1->ref,sw,1,inm);
        if(iswp) pt1->ref=swapbin(pt1->ref);
      }
      pt2->v[0] = pt1->v[0];
      pt2->v[1] = pt1->v[2];
      pt2->ref  = pt1->ref;
      for (i=0; i<3; i++) {
        ppt = &mesh->point[pt1->v[i]];
        ppt->tag &= ~MG_NUL;
        ppt = &mesh->point[pt2->v[i]];
        ppt->tag &= ~MG_NUL;
      }
    }
    mesh->nt = mesh->nti;
  }

  if(ncor) {
    rewind(inm);
    fseek(inm,posncor,SEEK_SET);
    for (k=1; k<=ncor; k++) {
      if(!bin)
        fscanf(inm,"%d",&i);
      else {
        fread(&i,sw,1,inm);
        if(iswp) i=swapbin(i);
      }
      if(i>mesh->np) {
        fprintf(stdout,"   Warning Corner number %8d IGNORED\n",i);
      } else {
        ppt = &mesh->point[i];
        ppt->tag |= MG_CRN;
      }
    }
  }

  /* read required vertices */
  if(npreq) {
    rewind(inm);
    fseek(inm,posnpreq,SEEK_SET);
    for (k=1; k<=npreq; k++) {
      if(!bin)
        fscanf(inm,"%d",&i);
      else {
        fread(&i,sw,1,inm);
        if(iswp) i=swapbin(i);
      }
      if(i>mesh->np) {
        fprintf(stdout,"   Warning Required Vertices number %8d IGNORED\n",i);
      } else {
        ppt = &mesh->point[i];
        ppt->tag |= MG_REQ;
      }
    }
  }

  nri = 0;
  if ( mesh->na ) {
    _MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"initial edges",return(0));
    _MMG5_SAFE_CALLOC(mesh->edge,mesh->na+1,MMG5_Edge);
    rewind(inm);
    fseek(inm,posned,SEEK_SET);
    for (k=1; k<=mesh->na; k++) {
      if (!bin)
        fscanf(inm,"%d %d %d",&mesh->edge[k].a,&mesh->edge[k].b,&mesh->edge[k].ref);
      else {
        fread(&mesh->edge[k].a,sw,1,inm);
        if(iswp) mesh->edge[k].a=swapbin(mesh->edge[k].a);
        fread(&mesh->edge[k].b,sw,1,inm);
        if(iswp) mesh->edge[k].b=swapbin(mesh->edge[k].b);
        fread(&mesh->edge[k].ref,sw,1,inm);
        if(iswp) mesh->edge[k].ref=swapbin(mesh->edge[k].ref);
      }
      mesh->edge[k].tag |= MG_REF;
      mesh->point[mesh->edge[k].a].tag |= MG_REF;
      mesh->point[mesh->edge[k].b].tag |= MG_REF;
    }

    if ( nri ) {
      rewind(inm);
      fseek(inm,posnr,SEEK_SET);
      for (k=1; k<=nri; k++) {
        if(!bin)
          fscanf(inm,"%d",&ia);
        else {
          fread(&ia,sw,1,inm);
          if(iswp) ia=swapbin(ia);
        }
        if ( ia > 0 && ia <= mesh->na )  mesh->edge[ia].tag |= MG_GEO;
      }
    }
    if ( nedreq ) {
      rewind(inm);
      fseek(inm,posnedreq,SEEK_SET);
      for (k=1; k<=nedreq; k++) {
        if(!bin)
          fscanf(inm,"%d",&ia);
        else {
          fread(&ia,sw,1,inm);
          if(iswp) ia=swapbin(ia);
        }
        if ( ia > 0 && ia <= mesh->na )   mesh->edge[ia].tag |= MG_REQ;
      }
    }
    /* set tria edges tags */
    assignEdge(mesh);
  }

  /* read geometric entities */
  if ( ng > 0 ) {
    _MMG5_SAFE_CALLOC(norm,3*ng+1,double);

    rewind(inm);
    fseek(inm,posnormal,SEEK_SET);
    for (k=1; k<=ng; k++) {
      n = &norm[3*(k-1)+1];
      if ( mesh->ver == 1 ) {
        if (!bin) {
          for (i=0 ; i<3 ; i++) {
            fscanf(inm,"%f",&fc);
            n[i] = (double) fc;
          }
        } else {
          for (i=0 ; i<3 ; i++) {
            fread(&fc,sw,1,inm);
            if(iswp) fc=swapf(fc);
            n[i] = (double) fc;
          }
        }
      }
      else {
        if (!bin)
          fscanf(inm,"%lf %lf %lf",&n[0],&n[1],&n[2]);
        else {
          for (i=0 ; i<3 ; i++) {
            fread(&n[i],sd,1,inm);
            if(iswp) n[i]=swapd(n[i]);
          }
        }
      }
      dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
      if ( dd > _MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        n[0] *= dd;
        n[1] *= dd;
        n[2] *= dd;
      }
    }

    rewind(inm);
    fseek(inm,posnc1,SEEK_SET);
    for (k=1; k<=mesh->nc1; k++) {
      if (!bin)
        fscanf(inm,"%d %d",&ip,&idn);
      else {
        fread(&ip,sw,1,inm);
        if(iswp) ip=swapbin(ip);
        fread(&idn,sw,1,inm);
        if(iswp) idn=swapbin(idn);
      }
      if ( idn > 0 && ip < mesh->np+1 )
        memcpy(&mesh->point[ip].n,&norm[3*(idn-1)+1],3*sizeof(double));
    }
    _MMG5_SAFE_FREE(norm);
  }

  if ( abs(mesh->info.imprim) > 4 ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d / %8d   CORNERS/REQ. %d / %d\n",mesh->npi,mesh->npmax,ncor,npreq);
    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d  RIDGES %6d\n",mesh->na,nri);
    fprintf(stdout,"     NUMBER OF TRIANGLES  %8d / %8d\n",mesh->nti,mesh->ntmax);
  }
  fclose(inm);
  return(1);
}

int MMG5_saveMesh(MMG5_pMesh mesh) {
  FILE         *inm;
  MMG5_pPoint  ppt;
  MMG5_pTria   pt;
  MMG5_pEdge   edge;
  MMG5_pxPoint go;
  int          *adja,k,jel,np,nt,na,nc,ng,nn,nr,nre;
  int          bin,binch,bpos;
  // int          outm;
  char         data[128],*ptr,chaine[128],i,i1,i2;

  edge = 0;
  mesh->ver = 2;
  strcpy(data,mesh->nameout);
  bin=0;
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if( !(inm = fopen(data,"wb")) ) {
      ptr  = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if( !(inm = fopen(data,"w")) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        return(0);
      }
    } else {
      bin = 1;
    }
  }
  else {
    ptr = strstr(data,".meshb");
    if( ptr ) bin = 1;
    if( !(inm = fopen(data,"w")) ) {
      fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
      return(0);
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /*entete fichier*/
  if(!bin) {
    strcpy(&chaine[0],"MeshVersionFormatted 2\n");
    fprintf(inm,"%s",chaine);
    strcpy(&chaine[0],"\n\nDimension 3\n");
    fprintf(inm,"%s ",chaine);
  } else {
    binch = 1; //MeshVersionFormatted
    fwrite(&binch,sw,1,inm);
    binch = 2; //version
    fwrite(&binch,sw,1,inm);
    binch = 3; //Dimension
    fwrite(&binch,sw,1,inm);
    bpos = 20; //Pos
    fwrite(&bpos,sw,1,inm);
    binch = 3;
    fwrite(&binch,sw,1,inm);

  }
  /* vertices */
  np = nc = ng = nn = nre = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->tmp = 0;
    if ( MG_VOK(ppt) ) {
      np++;
      ppt->tmp = np;
      if ( ppt->tag & MG_CRN )  nc++;
      if ( ppt->tag & MG_REQ )  nre++;
      if ( MG_EDG(ppt->tag) )   ng++;
    }
  }

  if(!bin) {
    strcpy(&chaine[0],"\n\nVertices\n");
    fprintf(inm,"%s",chaine);
    fprintf(inm,"%d\n",np);
  } else {
    binch = 4; //Vertices
    fwrite(&binch,sw,1,inm);
    bpos += 12+(1+3*mesh->ver)*4*np; //NullPos
    fwrite(&bpos,sw,1,inm);
    fwrite(&np,sw,1,inm);
  }
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      if(!bin) {
        fprintf(inm,"%.15lg %.15lg %.15lg %d\n",ppt->c[0],ppt->c[1],ppt->c[2],abs(ppt->ref));
      } else {
        fwrite((unsigned char*)&ppt->c[0],sd,1,inm);
        fwrite((unsigned char*)&ppt->c[1],sd,1,inm);
        fwrite((unsigned char*)&ppt->c[2],sd,1,inm);
        ppt->ref = abs(ppt->ref);
        fwrite((unsigned char*)&ppt->ref,sw,1,inm);
      }
      if ( !(ppt->tag & MG_GEO) )  nn++;
    }
  }

  nt = na = nr = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( MG_EOK(pt) ) {
      nt++;
      for (i=0; i<3; i++)
        if ( MG_EDG(pt->tag[i]) )  na++;
    }
  }

  /* memory alloc */
  if ( na ) {
    _MMG5_SAFE_CALLOC(edge,na+1,MMG5_Edge);
  }

  /* write triangles */
  if(!bin) {
    strcpy(&chaine[0],"\n\nTriangles\n");
    fprintf(inm,"%s",chaine);
    fprintf(inm,"%d \n",nt);
  } else {
    binch = 6; //Triangles
    fwrite(&binch,sw,1,inm);
    bpos += 12+16*nt; //Pos
    fwrite(&bpos,sw,1,inm);
    fwrite(&nt,sw,1,inm);
  }
  na = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( MG_EOK(pt) ) {
      if(!bin) {
        fprintf(inm,"%d %d %d %d\n",mesh->point[pt->v[0]].tmp,mesh->point[pt->v[1]].tmp
                ,mesh->point[pt->v[2]].tmp,abs(pt->ref));
      } else {
        fwrite(&mesh->point[pt->v[0]].tmp,sw,1,inm);
        fwrite(&mesh->point[pt->v[1]].tmp,sw,1,inm);
        fwrite(&mesh->point[pt->v[2]].tmp,sw,1,inm);
	pt->ref = abs(pt->ref);
        fwrite(&pt->ref,sw,1,inm);
      }
      for (i=0; i<3; i++) {
        if ( !MG_EDG(pt->tag[i]) )  continue;

        adja = &mesh->adja[3*(k-1)+1];
        jel  = adja[i] / 3;
        if ( !jel || jel > k ) {
          i1 = _MMG5_inxt2[i];
          i2 = _MMG5_inxt2[i1];
          na++;
          edge[na].a    = mesh->point[pt->v[i1]].tmp;
          edge[na].b    = mesh->point[pt->v[i2]].tmp;
          edge[na].ref  = pt->edg[i];
          edge[na].tag |= pt->tag[i];
          if ( pt->tag[i] & MG_GEO )  nr++;
        }
      }
    }
  }

  /* write corners */
  if ( nc ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nCorners\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",nc);
    } else {
      binch = 13; //
      fwrite(&binch,sw,1,inm);
      bpos += 12+4*nc; //NullPos
      fwrite(&bpos,sw,1,inm);
      fwrite(&nc,sw,1,inm);
    }
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && ppt->tag & MG_CRN ) {
        if(!bin) {
          fprintf(inm,"%d \n",ppt->tmp);
        } else {
          fwrite(&ppt->tmp,sw,1,inm);
        }
      }
    }
  }
  if ( nre ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nRequiredVertices\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",nre);
    } else {
      binch = 15; //
      fwrite(&binch,sw,1,inm);
      bpos += 12+4*nre; //NullPos
      fwrite(&bpos,sw,1,inm);
      fwrite(&nre,sw,1,inm);
    }
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && ppt->tag & MG_REQ ) {
        if(!bin) {
          fprintf(inm,"%d \n",ppt->tmp);
        } else {
          fwrite(&ppt->tmp,sw,1,inm);
        }
      }
    }
  }

  /* write edges, ridges */
  if ( na ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nEdges\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",na);
    } else {
      binch = 5; //Edges
      fwrite(&binch,sw,1,inm);
      bpos += 12 + 3*4*na;//Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&na,sw,1,inm);
    }
    nre = 0;
    for (k=1; k<=na; k++) {
      if(!bin) {
        fprintf(inm,"%d %d %d \n",edge[k].a,edge[k].b,edge[k].ref);
      } else {
        fwrite(&edge[k].a,sw,1,inm);
        fwrite(&edge[k].b,sw,1,inm);
        fwrite(&edge[k].ref,sw,1,inm);
      }
      if ( edge[k].tag & MG_REQ )  nre++;
    }
    /* ridges */
    if ( nr ) {
      if(!bin) {
        strcpy(&chaine[0],"\n\nRidges\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%d\n",nr);
      } else {
        binch = 14; //Ridges
        fwrite(&binch,sw,1,inm);
        bpos += 12 + 4*nr;//Pos
        fwrite(&bpos,sw,1,inm);
        fwrite(&nr,sw,1,inm);
      }
      for (k=1; k<=na; k++) {
        if ( edge[k].tag & MG_GEO ) {
          if(!bin) {
            fprintf(inm,"%d \n",k);
          } else {
            fwrite(&k,sw,1,inm);
          }
        }
      }
    }
    if ( nre ) {
      if(!bin) {
        strcpy(&chaine[0],"\n\nRequiredEdges\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%d\n",nre);
      } else {
        binch = 16; //RequiredEdges
        fwrite(&binch,sw,1,inm);
        bpos += 12 + 4*nre;//Pos
        fwrite(&bpos,sw,1,inm);
        fwrite(&nre,sw,1,inm);
      }
      for (k=1; k<=na; k++)
        if ( edge[k].tag & MG_REQ )  {
          if(!bin) {
            fprintf(inm,"%d \n",k);
          } else {
            fwrite(&k,sw,1,inm);
          }
        }
    }
    _MMG5_SAFE_FREE(edge);
  }

  /* write normals */
  if ( nn ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nNormals\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",nn);
    } else {
      binch = 60; //normals
      fwrite(&binch,sw,1,inm);
      bpos += 12+(3*mesh->ver)*4*nn; //Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&nn,sw,1,inm);
    }
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      else if ( !(ppt->tag & MG_GEO) ) {
        if ( ppt->tag & MG_REF ) {
          go = &mesh->xpoint[ppt->xp];
          if(!bin) {
            fprintf(inm,"%.15lg %.15lg %.15lg \n",go->n1[0],go->n1[1],go->n1[2]);
          } else {
            fwrite((unsigned char*)&go->n1[0],sd,1,inm);
            fwrite((unsigned char*)&go->n1[1],sd,1,inm);
            fwrite((unsigned char*)&go->n1[2],sd,1,inm);
          }
        }
        else {
          if(!bin) {
            fprintf(inm,"%.15lg %.15lg %.15lg \n",ppt->n[0],ppt->n[1],ppt->n[2]);
          } else {
            fwrite((unsigned char*)&ppt->n[0],sd,1,inm);
            fwrite((unsigned char*)&ppt->n[1],sd,1,inm);
            fwrite((unsigned char*)&ppt->n[2],sd,1,inm);
          }
        }
      }
    }
    if(!bin) {
      strcpy(&chaine[0],"\n\nNormalAtVertices\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",nn);
    } else {
      binch = 20; //normalatvertices
      fwrite(&binch,sw,1,inm);
      bpos += 12 + 2*4*nn;//Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&nn,sw,1,inm);
    }
    nn = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && !(ppt->tag & MG_GEO) ) {
        if(!bin) {
          fprintf(inm,"%d %d\n",ppt->tmp,++nn);
        } else {
          fwrite(&ppt->tmp,sw,1,inm);
          ++nn;
          fwrite(&nn,sw,1,inm);
        }
      }
    }
/*
  nn = 0;
  for (k=1; k<=mesh->nt; k++) {
  pt = &mesh->tria[k];
  if ( !MG_EOK(pt) )  continue;
  for (i=0; i<3; i++) {
  ppt = &mesh->point[pt->v[i]];
  if ( ppt->tag & MG_GEO )  nn++;
  }
  }
  GmfSetKwd(outm,GmfNormalAtTriangleVertices,nn);
  for (k=1; k<=mesh->nt; k++) {
  pt = &mesh->tria[k];
  if ( !MG_EOK(pt) )  continue;
  for (i=0; i<3; i++) {
  ppt = &mesh->point[pt->v[i]];
  if ( ppt->tag & MG_GEO )  nn++;
  }
*/
  }

  /* write tangents */
  if ( ng ) {
    /* Write tangents */
    if(!bin) {
      strcpy(&chaine[0],"\n\nTangents\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",ng);
    } else {
      binch = 59; //tangent
      fwrite(&binch,sw,1,inm);
      bpos += 12+(3*mesh->ver)*4*ng; //Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&ng,sw,1,inm);
    }
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && MG_EDG(ppt->tag) ) {
        if(!bin) {
          fprintf(inm,"%.15lg %.15lg %.15lg \n",ppt->n[0],ppt->n[1],ppt->n[2]);
        } else {
          fwrite((unsigned char*)&ppt->n[0],sd,1,inm);
          fwrite((unsigned char*)&ppt->n[1],sd,1,inm);
          fwrite((unsigned char*)&ppt->n[2],sd,1,inm);
        }
      }
    }
    if(!bin) {
      strcpy(&chaine[0],"\n\nTangentAtVertices\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",ng);
    } else {
      binch = 61; //tangentatvertices
      fwrite(&binch,sw,1,inm);
      bpos += 12 + 2*4*ng;//Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&ng,sw,1,inm);
    }
    ng = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && MG_EDG(ppt->tag) ) {
        if(!bin) {
          fprintf(inm,"%d %d\n",ppt->tmp,++ng);
        } else {
          fwrite(&ppt->tmp,sw,1,inm);
          ++ng;
          fwrite(&(ng),sw,1,inm);
        }
      }
    }
  }

  if ( abs(mesh->info.imprim) > 4 ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d  CORNERS    %6d\n",np,nc);
    if ( na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d  RIDGES     %6d\n",na,nr);
    fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",nt);
    if ( nn+ng )
      fprintf(stdout,"     NUMBER OF NORMALS    %8d  TANGENTS   %6d\n",nn,ng);
  }
  /*fin fichier*/
  if(!bin) {
    strcpy(&chaine[0],"\n\nEnd\n");
    fprintf(inm,"%s",chaine);
  } else {
    binch = 54; //End
    fwrite(&binch,sw,1,inm);
  }
  fclose(inm);
  return(1);
}

/* load metric field */
int MMG5_loadMet(MMG5_pMesh mesh,MMG5_pSol met) {
  FILE       *inm;
  float       fbuf[6],tmpf;
  double      dbuf[6],tmpd;
  int         binch,bdim,iswp;
  int         k,i,type,bin,bpos;
  long        posnp;
  char        *ptr,data[128],chaine[128];
  int         compute_hmin, compute_hmax;
  double      lambda[3],eigenv[3][3];

  if ( !met->namein )  return(0);
  posnp = 0;
  bin   = 0;
  iswp  = 0;

  strcpy(data,met->namein);
  ptr = strstr(data,".sol");
  if ( !ptr ) {
    /* data contains the filename without extension */
    strcat(data,".solb");
    if (!(inm = fopen(data,"rb"))  ) {
      /* our file is not a .solb file, try with .sol ext */
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      if (!(inm = fopen(data,"r"))  ) {
        fprintf(stderr,"  ** %s  NOT FOUND. USE DEFAULT METRIC.\n",data);
        return(-1);
      }
    } else {
      bin = 1;
    }
  }
  else {
    if (!(inm = fopen(data,"r")) ) {
      fprintf(stderr,"  ** %s  NOT FOUND. USE DEFAULT METRIC.\n",data);
      return(-1);
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* read solution or metric */
  if(!bin) {
    strcpy(chaine,"DDD");
    while(fscanf(inm,"%s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
      if(!strncmp(chaine,"MeshVersionFormatted",strlen("MeshVersionFormatted"))) {
        fscanf(inm,"%d",&met->ver);
        continue;
      } else if(!strncmp(chaine,"Dimension",strlen("Dimension"))) {
        fscanf(inm,"%d",&met->dim);
        if(met->dim!=3) {
          fprintf(stdout,"BAD SOL DIMENSION : %d\n",met->dim);
          return(1);
        }
        continue;
      } else if(!strncmp(chaine,"SolAtVertices",strlen("SolAtVertices"))) {
        fscanf(inm,"%d",&met->np);
        fscanf(inm,"%d",&type);
        if(type!=1) {
          fprintf(stdout,"SEVERAL SOLUTION => IGNORED : %d\n",type);
          return(1);
        }
        fscanf(inm,"%d",&met->size);
        posnp = ftell(inm);
        break;
      }
    }
  } else {
    fread(&binch,sw,1,inm);
    iswp=0;
    if(binch==16777216) iswp=1;
    else if(binch!=1) {
      fprintf(stdout,"BAD FILE ENCODING\n");
    }
    fread(&met->ver,sw,1,inm);
    if(iswp) met->ver = swapbin(met->ver);
    while(fread(&binch,sw,1,inm)!=EOF && binch!=54 ) {
      if(iswp) binch=swapbin(binch);
      if(binch==54) break;
      if(binch==3) {  //Dimension
        fread(&bdim,sw,1,inm);  //NulPos=>20
        if(iswp) bdim=swapbin(bdim);
        fread(&met->dim,sw,1,inm);
        if(iswp) met->dim=swapbin(met->dim);
        if(met->dim!=3) {
          fprintf(stdout,"BAD SOL DIMENSION : %d\n",met->dim);
          exit(0);
          return(1);
        }
        continue;
      } else if(binch==62) {  //SolAtVertices
        fread(&binch,sw,1,inm); //NulPos
        if(iswp) binch=swapbin(binch);
        fread(&met->np,sw,1,inm);
        if(iswp) met->np=swapbin(met->np);
        fread(&type,sw,1,inm); //nb sol
        if(iswp) type=swapbin(type);
        if(type!=1) {
          fprintf(stdout,"SEVERAL SOLUTION => IGNORED : %d\n",type);
          return(1);
        }
        fread(&met->size,sw,1,inm); //typsol
        if(iswp) met->size=swapbin(met->size);
        posnp = ftell(inm);
        break;
      } else {
        fread(&bpos,sw,1,inm); //Pos
        if(iswp) bpos=swapbin(bpos);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
      }
    }

  }

  if ( !met->np ) {
    fprintf(stdout,"  ** MISSING DATA.\n");
    return(1);
  }
  if ( (type != 1) || (met->size != 1 && met->size != 3) ) {
    fprintf(stdout,"  ** DATA IGNORED %d  %d\n",type,met->size);
    met->np = met->npmax = 0;
    return(-1);
  }

  if(met->size == 3) met->size = 6;

  /* mem alloc */
  if ( met->m )
    _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));

  _MMG5_ADD_MEM(mesh,(met->size*(met->npmax+1))*sizeof(double),
                "initial solution",return(0));
  _MMG5_SAFE_CALLOC(met->m,met->size*(met->npmax+1),double);

  rewind(inm);
  fseek(inm,posnp,SEEK_SET);

  /* If they are not provided by the user, enforce default values for hmin and
   * hmax:
   *    - for hmin we take 0.1 \times the minimum of the metric sizes.
   *    - for hmax we take 10 \times the max of the metric sizes. */
  compute_hmin = compute_hmax = 0;

  if ( mesh->info.hmin < 0. ) {
    compute_hmin=1;
    mesh->info.hmin = FLT_MAX;
  }
  if ( mesh->info.hmax < 0. ){
    compute_hmax=1;
    mesh->info.hmax = 0.;
  }

  /* isotropic metric */
  if ( met->size == 1 ) {
    if ( met->ver == 1 ) {
      for (k=1; k<=met->np; k++) {
        if(!bin){
          fscanf(inm,"%f",&fbuf[0]);
        } else {
          fread(&fbuf[0],sw,1,inm);
          if(iswp) fbuf[0]=swapf(fbuf[0]);
        }
        met->m[k] = fbuf[0];
      }
    }
    else {
      for (k=1; k<=met->np; k++) {
        if(!bin){
          fscanf(inm,"%lf",&dbuf[0]);
        } else {
          fread(&dbuf[0],sd,1,inm);
          if(iswp) dbuf[0]=swapd(dbuf[0]);
        }
        met->m[k] = dbuf[0];
      }
    }
    /* Find the minimum and maximum of the metric sizes. */
    if ( compute_hmin ) {
      for (k=1; k<=met->np; k++) {
        mesh->info.hmin = MG_MIN(mesh->info.hmin,met->m[k]);
      }
    }
    if ( compute_hmax ) {
      for (k=1; k<=met->np; k++) {
        mesh->info.hmax = MG_MAX(mesh->info.hmax,met->m[k]);
      }
    }
  }

  /* anisotropic metric */
  else {
    if ( met->ver == 1 ) {
      for (k=1; k<=met->np; k++) {
        if(!bin){
          for(i=0 ; i<met->size ; i++)
            fscanf(inm,"%f",&fbuf[i]);
        } else {
          for(i=0 ; i<met->size ; i++) {
            fread(&fbuf[i],sw,1,inm);
            if(iswp) fbuf[i]=swapf(fbuf[i]);
          }
        }
        tmpf    = fbuf[2];
        fbuf[2] = fbuf[3];
        fbuf[3] = tmpf;
        for (i=0; i<6; i++)  met->m[6*k+i] = fbuf[i];
      }
    }
    else {
      for (k=1; k<=met->np; k++) {
        if(!bin){
          for(i=0 ; i<met->size ; i++)
            fscanf(inm,"%lf",&dbuf[i]);
        } else {
          for(i=0 ; i<met->size ; i++) {
            fread(&dbuf[i],sw,1,inm);
            if(iswp) dbuf[i]=swapf(dbuf[i]);
          }
        }
        tmpd    = dbuf[2];
        dbuf[2] = dbuf[3];
        dbuf[3] = tmpd;
        for (i=0; i<met->size; i++)  met->m[6*k+i] = dbuf[i];
      }
    }

    for (k=1; k<=met->np; k++) {
      if ( ! _MMG5_eigenv(1,&met->m[6*k],lambda,eigenv ) ) {
        printf("Error: metric diagonalisation fail,"
               " unable to compute the sizes associated to the vertex %d.\n",
               k);
        return(0);
      }
      if ( compute_hmin ) {
        mesh->info.hmin = MG_MIN(mesh->info.hmin,1./sqrt(lambda[0]));
        mesh->info.hmin = MG_MIN(mesh->info.hmin,1./sqrt(lambda[1]));
        mesh->info.hmin = MG_MIN(mesh->info.hmin,1./sqrt(lambda[2]));
      }
      if ( compute_hmax ) {
        mesh->info.hmax = MG_MAX(mesh->info.hmax,1./sqrt(lambda[0]));
        mesh->info.hmax = MG_MAX(mesh->info.hmax,1./sqrt(lambda[1]));
        mesh->info.hmax = MG_MAX(mesh->info.hmax,1./sqrt(lambda[2]));
      }
    }
  }

  if ( compute_hmin ) {
    mesh->info.hmin *=.1;
    /* Check that user has not given a hmax value lower that the founded
     * hmin. */
    if ( mesh->info.hmin > mesh->info.hmax ) {
      mesh->info.hmin = 0.1*mesh->info.hmax;
    }
  }
  if ( compute_hmax ) {
    mesh->info.hmax *=10.;
    /* Check that user has not given a hmin value bigger that the founded
     * hmax. */
    if ( mesh->info.hmax < mesh->info.hmin ) {
      mesh->info.hmax = 10.*mesh->info.hmin;
    }
  }

  fclose(inm);
  return(1);
}

/* write iso or aniso metric */
int saveMet(MMG5_pMesh mesh,MMG5_pSol met) {
  FILE*        inm;
  MMG5_pPoint  ppt;
  double       dbuf[6],tmp;
  char        *ptr,data[128],chaine[128];
  int          binch,bpos,bin,np,k,typ,i;

  met->ver = 2;
  bin = 0;
  strcpy(data,met->nameout);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  ptr = strstr(data,".sol");
  if ( !ptr )  strcat(data,".sol");
  if (  !(inm = fopen(data,"w")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /*entete fichier*/
  if(!bin) {
    strcpy(&chaine[0],"MeshVersionFormatted 2\n");
    fprintf(inm,"%s",chaine);
    strcpy(&chaine[0],"\n\nDimension 3\n");
    fprintf(inm,"%s ",chaine);
  } else {
    binch = 1; //MeshVersionFormatted
    fwrite(&binch,sw,1,inm);
    binch = 2; //version
    fwrite(&binch,sw,1,inm);
    binch = 3; //Dimension
    fwrite(&binch,sw,1,inm);
    bpos = 20; //Pos
    fwrite(&bpos,sw,1,inm);
    binch = 3;
    fwrite(&binch,sw,1,inm);

  }


  np = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) )  np++;
  }

  if(met->size==1) {
    typ = 1;
  } else {
    typ = 3;
  }

  if(!bin) {
    strcpy(&chaine[0],"\n\nSolAtVertices\n");
    fprintf(inm,"%s",chaine);
    fprintf(inm,"%d\n",np);
    fprintf(inm,"%d %d\n",1,typ);
  } else {
    binch = 62; //Vertices
    fwrite(&binch,sw,1,inm);
    bpos += 20+(met->size*met->ver)*4*np; //Pos
    fwrite(&bpos,sw,1,inm);
    fwrite(&np,sw,1,inm);
    binch = 1; //nb sol
    fwrite(&binch,sw,1,inm);
    binch = typ; //typ sol
    fwrite(&binch,sw,1,inm);
  }

  /* write isotropic metric */
  if ( met->size == 1 ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) ) {
        dbuf[0] = met->m[k];
        if(!bin) {
          fprintf(inm,"%.15lg \n",dbuf[0]);
        } else {
          fwrite((unsigned char*)&dbuf[0],sd,1,inm);
        }
      }
    }
  }
  /* write anisotropic metric */
  else {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) ) {
        for (i=0; i<met->size; i++)  dbuf[i] = met->m[met->size*(k)+i];
        tmp = dbuf[2];
        dbuf[2] = dbuf[3];
        dbuf[3] = tmp;
        if(!bin) {
          for(i=0; i<met->size; i++)
            fprintf(inm,"%.15lg  ",dbuf[i]);
          fprintf(inm,"\n");
        } else {
          for(i=0; i<met->size; i++)
            fwrite((unsigned char*)&dbuf[i],sd,1,inm);
        }
      }
    }
  }
  /*fin fichier*/
  if(!bin) {
    strcpy(&chaine[0],"\n\nEnd\n");
    fprintf(inm,"%s",chaine);
  } else {
    binch = 54; //End
    fwrite(&binch,sw,1,inm);
  }
  fclose(inm);
  return(1);
}


