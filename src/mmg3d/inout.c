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
 * \file mmg3d/inout.c
 * \brief Input / Output Functions.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmg3d.h"


#define sw 4
#define sd 8

static int _MMG5_swapbin(int sbin)
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
static float _MMG5_swapf(float sbin)
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
static double _MMG5_swapd(double sbin)
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

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh data.
 *
 */
int MMG5_loadMesh(MMG5_pMesh mesh) {
  FILE*       inm;
  MMG5_pTetra pt;
  MMG5_pTria  pt1;
  MMG5_pEdge  pa;
  MMG5_pPoint ppt;
  int         posnp,posnt,posne,posned,posncor,posnpreq,posntreq,posnereq,posnedreq;
  int         posnr;
  int         npreq,ntreq,nereq,nedreq,ncor,ned,bin,iswp;
  int         binch,bdim,bpos,i,k;
  int         *ina,v[3],ref,nt,na,nr,ia,aux;
  float            fc;
  char        *ptr,*name,data[128],chaine[128];

  posnp = posnt = posne = posncor = 0;
  posnpreq = posntreq = posnereq = posned = posnedreq = posnr = 0;
  ncor = ned = npreq = ntreq = nereq = nedreq = nr = 0;
  bin = 0;
  iswp = 0;
  ina = NULL;
  mesh->np = mesh->nt = mesh->ne = 0;

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
          return(-1);
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
      } else if(!strncmp(chaine,"Tetrahedra",strlen("Tetrahedra"))) {
        fscanf(inm,"%d",&mesh->nei);
        posne = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredTetrahedra",strlen("RequiredTetrahedra"))) {
        fscanf(inm,"%d",&nereq);
        posnereq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Corners",strlen("Corners"))) {
        fscanf(inm,"%d",&ncor);
        posncor = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Edges",strlen("Edges"))) {
        fscanf(inm,"%d",&mesh->nai);
        posned = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredEdges",strlen("RequiredEdges"))) {
        fscanf(inm,"%d",&nedreq);
        posnedreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Ridges",strlen("Ridges"))) {
        fscanf(inm,"%d",&nr);
        posnr = ftell(inm);
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
    if(iswp) mesh->ver = _MMG5_swapbin(mesh->ver);
    while(fread(&binch,sw,1,inm)!=0 && binch!=54 ) {
      if(iswp) binch=_MMG5_swapbin(binch);
      if(binch==54) break;
      if(!bdim && binch==3) {  //Dimension
        fread(&bdim,sw,1,inm);  //NulPos=>20
        if(iswp) bdim=_MMG5_swapbin(bdim);
        fread(&bdim,sw,1,inm);
        if(iswp) bdim=_MMG5_swapbin(bdim);
        mesh->dim = bdim;
        if(bdim!=3) {
          fprintf(stdout,"BAD SOL DIMENSION : %d\n",mesh->dim);
          fprintf(stdout," Exit program.\n");
          return(-1);
        }
        continue;
      } else if(!mesh->npi && binch==4) {  //Vertices
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        fread(&mesh->npi,sw,1,inm);
        if(iswp) mesh->npi=_MMG5_swapbin(mesh->npi);
        posnp = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==15) {  //RequiredVertices
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        fread(&npreq,sw,1,inm);
        if(iswp) npreq=_MMG5_swapbin(npreq);
        posnpreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!mesh->nti && binch==6) {//Triangles
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        fread(&mesh->nti,sw,1,inm);
        if(iswp) mesh->nti=_MMG5_swapbin(mesh->nti);
        posnt = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==17) {  //RequiredTriangles
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        fread(&ntreq,sw,1,inm);
        if(iswp) ntreq=_MMG5_swapbin(ntreq);
        posntreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!mesh->nei && binch==8) {//Tetra
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        fread(&mesh->nei,sw,1,inm);
        if(iswp) mesh->nei=_MMG5_swapbin(mesh->nei);
        posne = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==12) {  //RequiredTetra
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        fread(&nereq,sw,1,inm);
        if(iswp) nereq=_MMG5_swapbin(nereq);
        posnereq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!ncor && binch==13) { //Corners
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        fread(&ncor,sw,1,inm);
        if(iswp) ncor=_MMG5_swapbin(ncor);
        posncor = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!mesh->nai && binch==5) { //Edges
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        fread(&mesh->nai,sw,1,inm);
        if(iswp) mesh->nai=_MMG5_swapbin(mesh->nai);
        posned = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==16) {  //RequiredEdges
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        fread(&nedreq,sw,1,inm);
        if(iswp) nedreq=_MMG5_swapbin(nedreq);
        posnedreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      }  else if(binch==14) {  //Ridges
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        fread(&nr,sw,1,inm);
        if(iswp) nr=_MMG5_swapbin(nr);
        posnr = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else {
        //printf("on traite ? %d\n",binch);
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        //printf("on avance... Nulpos %d\n",bpos);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
      }
    }
  }

  if ( !mesh->npi || !mesh->nei ) {
    fprintf(stdout,"  ** MISSING DATA.\n");
    fprintf(stdout," Check that your mesh contains points and tetrahedra.\n");
    fprintf(stdout," Exit program.\n");
    return(-1);
  }
  /* memory allocation */
  mesh->np = mesh->npi;
  mesh->nt = mesh->nti;
  mesh->ne = mesh->nei;
  mesh->na = mesh->nai;
  if ( !_MMG5_zaldy(mesh) )  return(0);
  if (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne) {
    return(-1);
  }

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
          if(iswp) fc=_MMG5_swapf(fc);
          ppt->c[i] = (double) fc;
        }
        fread(&ppt->ref,sw,1,inm);
        if(iswp) ppt->ref=_MMG5_swapbin(ppt->ref);
      }
    } else {
      if (!bin)
        fscanf(inm,"%lf %lf %lf %d",&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
      else {
        for (i=0 ; i<3 ; i++) {
          fread(&ppt->c[i],sd,1,inm);
          if(iswp) ppt->c[i]=_MMG5_swapd(ppt->c[i]);
        }
        fread(&ppt->ref,sw,1,inm);
        if(iswp) ppt->ref=_MMG5_swapbin(ppt->ref);
      }
    }
    ppt->tag  = MG_NUL;
    ppt->tmp  = 0;
  }
  /* get required vertices */
  if(npreq) {
    rewind(inm);
    fseek(inm,posnpreq,SEEK_SET);
    for (k=1; k<=npreq; k++) {
      if(!bin)
        fscanf(inm,"%d",&i);
      else {
        fread(&i,sw,1,inm);
        if(iswp) i=_MMG5_swapbin(i);
      }
      if(i>mesh->np) {
        fprintf(stdout,"   Warning: Required Vertices number %8d IGNORED\n",i);
      } else {
        ppt = &mesh->point[i];
        ppt->tag |= MG_REQ;
      }
    }
  }

  /* get corners */
  if(ncor) {
    rewind(inm);
    fseek(inm,posncor,SEEK_SET);
    for (k=1; k<=ncor; k++) {
      if(!bin)
        fscanf(inm,"%d",&i);
      else {
        fread(&i,sw,1,inm);
        if(iswp) i=_MMG5_swapbin(i);
      }
      if(i>mesh->np) {
        fprintf(stdout,"   Warning: Corner number %8d IGNORED\n",i);
      } else {
        ppt = &mesh->point[i];
        ppt->tag |= MG_CRN;
      }
    }
  }

  /* read mesh triangles */
  nt = 0;
  if ( mesh->nt ) {
    rewind(inm);
    fseek(inm,posnt,SEEK_SET);
    /* Skip triangles with MG_ISO refs */
    if( mesh->info.iso ) {
      nt = mesh->nt;
      mesh->nt = 0;
      _MMG5_SAFE_CALLOC(ina,nt+1,int);

      for (k=1; k<=nt; k++) {
        if (!bin)
          fscanf(inm,"%d %d %d %d",&v[0],&v[1],&v[2],&ref);
        else {
          for (i=0 ; i<3 ; i++) {
            fread(&v[i],sw,1,inm);
            if(iswp) v[i]=_MMG5_swapbin(v[i]);
          }
          fread(&ref,sw,1,inm);
          if(iswp) ref=_MMG5_swapbin(ref);
        }
        if( abs(ref) != MG_ISO ) {
          pt1 = &mesh->tria[++mesh->nt];
          pt1->v[0] = v[0];
          pt1->v[1] = v[1];
          pt1->v[2] = v[2];
          pt1->ref = ref;
          ina[k]=mesh->nt;
        }
      }
      if( !mesh->nt )
        _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(MMG5_Tria));

      else if ( mesh->nt < nt ) {
        _MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(MMG5_Tria),"triangles",
                      printf("  Exit program.\n");
                      exit(EXIT_FAILURE));
        _MMG5_SAFE_RECALLOC(mesh->tria,nt+1,(mesh->nt+1),MMG5_Tria,"triangles");
      }
    }
    else {
      for (k=1; k<=mesh->nt; k++) {
        pt1 = &mesh->tria[k];
        if (!bin)
          fscanf(inm,"%d %d %d %d",&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
        else {
          for (i=0 ; i<3 ; i++) {
            fread(&pt1->v[i],sw,1,inm);
            if(iswp) pt1->v[i]=_MMG5_swapbin(pt1->v[i]);
          }
          fread(&pt1->ref,sw,1,inm);
          if(iswp) pt1->ref=_MMG5_swapbin(pt1->ref);
        }
      }
    }
    /* get required triangles */
    if(ntreq) {
      rewind(inm);
      fseek(inm,posntreq,SEEK_SET);
      for (k=1; k<=ntreq; k++) {
        if(!bin)
          fscanf(inm,"%d",&i);
        else {
          fread(&i,sw,1,inm);
          if(iswp) i=_MMG5_swapbin(i);
        }
        if ( i>mesh->nt ) {
          fprintf(stdout,"   Warning: Required MMG5_Triangles number %8d IGNORED\n",i);
        } else {
          if( mesh->info.iso ){
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
    }
    if ( mesh->info.iso )
      _MMG5_SAFE_FREE(ina);
  } //end if mesh->Nt

    /* read mesh edges */
  if ( mesh->na ) {
    na = mesh->na;
    if (mesh->info.iso ) {
      mesh->na = 0;
      _MMG5_SAFE_CALLOC(ina,na+1,int);
    }

    rewind(inm);
    fseek(inm,posned,SEEK_SET);

    for (k=1; k<=na; k++) {
      pa = &mesh->edge[k];
      if (!bin)
        fscanf(inm,"%d %d %d",&pa->a,&pa->b,&pa->ref);
      else {
        fread(&pa->a,sw,1,inm);
        if(iswp) pa->a=_MMG5_swapbin(pa->a);
        fread(&pa->b,sw,1,inm);
        if(iswp) pa->b=_MMG5_swapbin(pa->b);
        fread(&pa->ref,sw,1,inm);
        if(iswp) pa->ref=_MMG5_swapbin(pa->ref);
      }
      pa->tag |= MG_REF;
      if ( mesh->info.iso ) {
        if( abs(pa->ref) != MG_ISO ) {
          ++mesh->na;
          pa->ref = abs(pa->ref);
          memmove(&mesh->edge[mesh->na],&mesh->edge[k],sizeof(MMG5_Edge));
          ina[k] = mesh->na;
        }
      }
    }

    /* get ridges */
    if ( nr ) {
      rewind(inm);
      fseek(inm,posnr,SEEK_SET);
      for (k=1; k<=nr; k++) {
        if(!bin)
          fscanf(inm,"%d",&ia);
        else {
          fread(&ia,sw,1,inm);
          if(iswp) ia=_MMG5_swapbin(ia);
        }
        if(ia>na) {
          fprintf(stdout,"   Warning Ridge number %8d IGNORED\n",ia);
          continue;
        }
        if( mesh->info.iso ){
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
    if ( nedreq ) {
      rewind(inm);
      fseek(inm,posnedreq,SEEK_SET);
      for (k=1; k<=nedreq; k++) {
        if(!bin)
          fscanf(inm,"%d",&ia);
        else {
          fread(&ia,sw,1,inm);
          if(iswp) ia=_MMG5_swapbin(ia);
        }
        if(ia>na) {
          fprintf(stdout,"   Warning Required Edges number %8d/%8d IGNORED\n",ia,na);
          continue;
        }
        if( mesh->info.iso ){
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
    if (mesh->info.iso )
      _MMG5_SAFE_FREE(ina);
  }

  /* read mesh tetrahedra */
  rewind(inm);
  fseek(inm,posne,SEEK_SET);
  mesh->xt = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if (!bin)
      fscanf(inm,"%d %d %d %d %d",&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&ref);
    else {
      for (i=0 ; i<4 ; i++) {
        fread(&pt->v[i],sw,1,inm);
        if(iswp) pt->v[i]=_MMG5_swapbin(pt->v[i]);
      }
      fread(&ref,sw,1,inm);
      if(iswp) ref=_MMG5_swapbin(ref);
    }
    pt->ref  = ref;//0;//ref ;
    pt->qual = _MMG5_orcal(mesh,k);
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pt->v[i]];
      ppt->tag &= ~MG_NUL;
    }

    if ( mesh->info.iso )  pt->ref = 0;

    /* Possibly switch 2 vertices number so that each tet is positively oriented */
    if ( _MMG5_orvol(mesh->point,pt->v) < 0.0 ) {
      /* mesh->xt temporary used to count reoriented tetra */
      mesh->xt++;
      aux = pt->v[2];
      pt->v[2] = pt->v[3];
      pt->v[3] = aux;
    }
  }
  if(mesh->xt) {
    fprintf(stdout,"\n     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n");
    fprintf(stdout,"         BAD ORIENTATION : vol < 0 -- %8d tetra reoriented\n",mesh->xt);
    fprintf(stdout,"     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n\n");
  }
  mesh->xt = 0;
  /* get required tetrahedra */
  if(nereq) {
    rewind(inm);
    fseek(inm,posnereq,SEEK_SET);
    for (k=1; k<=nereq; k++) {
      if(!bin)
        fscanf(inm,"%d",&i);
      else {
        fread(&i,sw,1,inm);
        if(iswp) i=_MMG5_swapbin(i);
      }
      if(i>mesh->ne) {
        fprintf(stdout,"   Warning Required Tetra number %8d IGNORED\n",i);
        continue;
      }
      pt = &mesh->tetra[i];
      pt->tag |= MG_REQ;
    }
  }


  /* stats */
  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"     NUMBER OF VERTICES     %8d\n",mesh->np);
    if ( mesh->na ) {
      fprintf(stdout,"     NUMBER OF EDGES        %8d\n",mesh->na);
      if ( nr )
        fprintf(stdout,"     NUMBER OF RIDGES       %8d\n",nr);
    }
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES    %8d\n",mesh->nt);
    fprintf(stdout,"     NUMBER OF ELEMENTS     %8d\n",mesh->ne);

    if ( npreq || nedreq || ntreq || nereq ) {
      fprintf(stdout,"     NUMBER OF REQUIRED ENTITIES: \n");
      if ( npreq )
        fprintf(stdout,"                  VERTICES    %8d \n",npreq);
      if ( nedreq )
        fprintf(stdout,"                  EDGES       %8d \n",nedreq);
      if ( ntreq )
        fprintf(stdout,"                  TRIANGLES   %8d \n",ntreq);
      if ( nereq )
        fprintf(stdout,"                  TETRAHEDRAS %8d \n",nereq);
    }
    if(ncor) fprintf(stdout,"     NUMBER OF CORNERS        %8d \n",ncor);
  }
  fclose(inm);
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 otherwise.
 * \remark Function call through the function pointer \ref MMG5_saveMesh.
 *
 * Save mesh data.
 *
 */
int _MMG5_saveAllMesh(MMG5_pMesh mesh) {
  FILE*        inm;
  MMG5_pPoint  ppt;
  MMG5_pTetra  pt;
  MMG5_pTria   ptt;
  MMG5_xPoint *pxp;
  MMG5_hgeom  *ph;
  int          k,i,na,nc,np,ne,nn,nr,nre,nedreq,ntreq,nt,nereq;
  int          bin,binch,bpos;
  char         data[128],chaine[128],*ptr;

  mesh->ver = 2;
  bin = 0;
  strcpy(data,mesh->nameout);
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
  binch=0; bpos=10;
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
  np = nc = na = nr = nre = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      ppt->tmp = ++np;
      if ( ppt->tag & MG_CRN )  nc++;
      if ( ppt->tag & MG_REQ )  nre++;
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
    }
  }

  /* corners+required */
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
  if ( nre && !mesh->info.nosurf ) {
    /* Don't save the required vertices when no surface remeshing (because all
     * the surface vertices are required). */
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
      if ( !MG_VOK(ppt) || MG_SIN(ppt->tag) )  continue;
      else if ( (ppt->tag & MG_BDY)
                && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM)) ) {
        pxp = &mesh->xpoint[ppt->xp];
        if(!bin) {
          fprintf(inm,"%.15lg %.15lg %.15lg \n",pxp->n1[0],pxp->n1[1],pxp->n1[2]);
        } else {
          fwrite((unsigned char*)&pxp->n1[0],sd,1,inm);
          fwrite((unsigned char*)&pxp->n1[1],sd,1,inm);
          fwrite((unsigned char*)&pxp->n1[2],sd,1,inm);
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
      if ( !MG_VOK(ppt) || MG_SIN(ppt->tag) )  continue;
      else if ( (ppt->tag & MG_BDY)
                && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM)) ) {
        if(!bin) {
          fprintf(inm,"%d %d\n",ppt->tmp,++nn);
        } else {
          fwrite(&ppt->tmp,sw,1,inm);
          ++nn;
          fwrite(&nn,sw,1,inm);
        }
      }
    }

    if ( nt ) {
      /* Write tangents */
      if(!bin) {
        strcpy(&chaine[0],"\n\nTangents\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%d\n",nt);
      } else {
        binch = 59; //tangent
        fwrite(&binch,sw,1,inm);
        bpos += 12+(3*mesh->ver)*4*nt; //Pos
        fwrite(&bpos,sw,1,inm);
        fwrite(&nt,sw,1,inm);
      }

      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || MG_SIN(ppt->tag) )  continue;
        else if ( MG_EDG(ppt->tag) || (ppt->tag & MG_NOM) ) {
          pxp = &mesh->xpoint[ppt->xp];
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
        fprintf(inm,"%d\n",nt);
      } else {
        binch = 61; //tangentatvertices
        fwrite(&binch,sw,1,inm);
        bpos += 12 + 2*4*nt;//Pos
        fwrite(&bpos,sw,1,inm);
        fwrite(&nt,sw,1,inm);
      }
      nt = 0;
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || MG_SIN(ppt->tag) )  continue;
        else if ( MG_EDG(ppt->tag) || (ppt->tag & MG_NOM) ) {
          if(!bin) {
            fprintf(inm,"%d %d\n",ppt->tmp,++nn);
          } else {
            fwrite(&ppt->tmp,sw,1,inm);
            ++nn;
            fwrite(&(nn),sw,1,inm);
          }
        }
      }
    }
  }
  _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
  mesh->xp = 0;

  /* boundary mesh */
  /* tria + required tria */
  mesh->nt = ntreq = 0;
  if ( mesh->tria )
    _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(MMG5_Tria));

  _MMG5_chkNumberOfTri(mesh);
  if ( _MMG5_bdryTria(mesh) ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nTriangles\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d \n",mesh->nt);
    } else {
      binch = 6; //Triangles
      fwrite(&binch,sw,1,inm);
      bpos += 12+16*mesh->nt; //Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&mesh->nt,sw,1,inm);
    }
    for (k=1; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];
      if ( ptt->tag[0] & MG_REQ && ptt->tag[1] & MG_REQ && ptt->tag[2] & MG_REQ )  ntreq++;
      if(!bin) {
        fprintf(inm,"%d %d %d %d\n",mesh->point[ptt->v[0]].tmp,mesh->point[ptt->v[1]].tmp
                ,mesh->point[ptt->v[2]].tmp,ptt->ref);
      } else {
        fwrite(&mesh->point[ptt->v[0]].tmp,sw,1,inm);
        fwrite(&mesh->point[ptt->v[1]].tmp,sw,1,inm);
        fwrite(&mesh->point[ptt->v[2]].tmp,sw,1,inm);
        fwrite(&ptt->ref,sw,1,inm);
      }
    }
    if ( ntreq  && !mesh->info.nosurf ) {
      /* Don't save the required triangles when no surface remeshing (because
       * all the surface triangles are required). */
      if(!bin) {
        strcpy(&chaine[0],"\n\nRequiredTriangles\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%d \n",ntreq);
      } else {
        binch = 17; //ReqTriangles
        fwrite(&binch,sw,1,inm);
        bpos += 12+4*ntreq; //Pos
        fwrite(&bpos,sw,1,inm);
        fwrite(&ntreq,sw,1,inm);
      }
      for (k=0; k<=mesh->nt; k++) {
        ptt = &mesh->tria[k];
        if ( ptt->tag[0] & MG_REQ && ptt->tag[1] & MG_REQ && ptt->tag[2] & MG_REQ ) {
          if(!bin) {
            fprintf(inm,"%d \n",k);
          } else {
            fwrite(&k,sw,1,inm);
          }
        }
      }
    }

    /* Release memory before htab allocation */
    if ( mesh->adjt )
      _MMG5_DEL_MEM(mesh,mesh->adjt,(3*mesh->nt+4)*sizeof(int));
    if ( mesh->adja )
      _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

    /* build hash table for edges */
    if ( mesh->htab.geom )
      _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));

    /* in the wost case (all edges are marked), we will have around 1 edge per *
     * triangle (we count edges only one time) */
    na = nr = nedreq = 0;
    mesh->memCur += (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));
    if ( (mesh->memCur) > (mesh->memMax) ) {
      mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));
      fprintf(stdout,"  ## Warning:");
      fprintf(stdout," unable to allocate htab.\n");
    }
    else if ( _MMG5_hNew(&mesh->htab,mesh->nt,3*(mesh->nt),0) ) {
      for (k=1; k<=mesh->ne; k++) {
        pt   = &mesh->tetra[k];
        if ( MG_EOK(pt) &&  pt->xt ) {
          for (i=0; i<6; i++) {
            if ( mesh->xtetra[pt->xt].edg[i] ||
                 ( MG_EDG(mesh->xtetra[pt->xt].tag[i] ) ||
                   (mesh->xtetra[pt->xt].tag[i] & MG_REQ) ) )
              _MMG5_hEdge(mesh,pt->v[_MMG5_iare[i][0]],pt->v[_MMG5_iare[i][1]],
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
        for (k=0; k<=mesh->htab.max; k++) {
          ph = &mesh->htab.geom[k];
          if ( !ph->a )  continue;
          if(!bin) {
            fprintf(inm,"%d %d %d \n",mesh->point[ph->a].tmp,mesh->point[ph->b].tmp,ph->ref);
          } else {
            fwrite(&mesh->point[ph->a].tmp,sw,1,inm);
            fwrite(&mesh->point[ph->b].tmp,sw,1,inm);
            fwrite(&ph->ref,sw,1,inm);
          }
        }
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
          na = 0;
          for (k=0; k<=mesh->htab.max; k++) {
            ph = &mesh->htab.geom[k];
            if ( !ph->a )  continue;
            na++;
            if ( ph->tag & MG_GEO )  {
              if(!bin) {
                fprintf(inm,"%d \n",na);
              } else {
                fwrite(&na,sw,1,inm);
              }
            }
          }
        }
        if ( nedreq  && !mesh->info.nosurf ) {
          /* Don't save the required edges when no surface remeshing (because
           * all the surface edges are required). */
          if(!bin) {
            strcpy(&chaine[0],"\n\nRequiredEdges\n");
            fprintf(inm,"%s",chaine);
            fprintf(inm,"%d\n",nedreq);
          } else {
            binch = 16; //RequiredEdges
            fwrite(&binch,sw,1,inm);
            bpos += 12 + 4*nedreq;//Pos
            fwrite(&bpos,sw,1,inm);
            fwrite(&nedreq,sw,1,inm);
          }
          na = 0;
          for (k=0; k<=mesh->htab.max; k++) {
            ph = &mesh->htab.geom[k];
            if ( !ph->a )  continue;
            na++;
            if ( ph->tag & MG_REQ )  {
              if(!bin) {
                fprintf(inm,"%d \n",na);
              } else {
                fwrite(&na,sw,1,inm);
              }
            }
          }
        }
      }
      //_MMG5_freeXTets(mesh);
      _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));
    }
    else
      mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));
  } //fin if bdrytria....

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

  if(!bin) {
    strcpy(&chaine[0],"\n\nTetrahedra\n");
    fprintf(inm,"%s",chaine);
    fprintf(inm,"%d\n",ne);
  } else {
    binch = 8; //Tetra
    fwrite(&binch,sw,1,inm);
    bpos += 12 + 20*ne;//Pos
    fwrite(&bpos,sw,1,inm);
    fwrite((unsigned char*)&ne,sw,1,inm);
  }
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( MG_EOK(pt) ) {
      if(!bin) {
        fprintf(inm,"%d %d %d %d %d\n",mesh->point[pt->v[0]].tmp,mesh->point[pt->v[1]].tmp
                ,mesh->point[pt->v[2]].tmp,mesh->point[pt->v[3]].tmp,pt->ref);
      } else {
        fwrite(&mesh->point[pt->v[0]].tmp,sw,1,inm);
        fwrite(&mesh->point[pt->v[1]].tmp,sw,1,inm);
        fwrite(&mesh->point[pt->v[2]].tmp,sw,1,inm);
        fwrite(&mesh->point[pt->v[3]].tmp,sw,1,inm);
        fwrite(&pt->ref,sw,1,inm);
      }
    }
  }

  if ( nereq ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nRequiredTetrahedra\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",nereq);
    } else {
      binch = 12; //RequiredTetra
      fwrite(&binch,sw,1,inm);
      bpos += 12 + 4*nereq;//Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&nereq,sw,1,inm);
    }
    ne = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;
      ne++;
      if ( pt->tag & MG_REQ ) {
        if(!bin) {
          fprintf(inm,"%d \n",ne);
        } else {
          fwrite(&ne,sw,1,inm);
        }
      }
    }
  }

  if ( mesh->info.imprim ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",np,nc+nre);
    if ( na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d   RIDGES  %8d\n",na,nr);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",mesh->nt);
    fprintf(stdout,"     NUMBER OF ELEMENTS   %8d\n",mesh->ne);
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

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 otherwise.
 * \remark Function call through the function pointer \ref MMG5_saveMesh.
 *
 * Save mesh data without adjacency and xtetra tables (for library version).
 *
 */
int _MMG5_saveLibraryMesh(MMG5_pMesh mesh) {
  FILE        *inm;
  MMG5_pPoint  ppt;
  MMG5_pTetra  pt;
  MMG5_pTria   ptt;
  int          k,nc,np,ne,nr,nre,nedreq,ntreq,nereq;
  int          bin,binch,bpos,na;
  char         data[128],chaine[128],*ptr;

  mesh->ver = 2;
  bin = 0;
  strcpy(data,mesh->nameout);
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
  np = nc = nre = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      ppt->tmp = ++np;
      if ( ppt->tag & MG_CRN )  nc++;
      if ( ppt->tag & MG_REQ )  nre++;
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
    }
  }

  /* corners+required */
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
  if ( nre  && !mesh->info.nosurf ) {
    /* Don't save the required vertices when no surface remeshing (because all
     * the surface vertices are required). */
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
  /* boundary mesh */
  /* tria + required tria */
  ntreq = 0;
  if ( mesh->nt ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nTriangles\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d \n",mesh->nt);
    } else {
      binch = 6; //Triangles
      fwrite(&binch,sw,1,inm);
      bpos += 12+16*mesh->nt; //Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&mesh->nt,sw,1,inm);
    }
    for (k=1; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];
      if ( ptt->tag[0] & MG_REQ && ptt->tag[1] & MG_REQ && ptt->tag[2] & MG_REQ )  ntreq++;
      if(!bin) {
        fprintf(inm,"%d %d %d %d\n",mesh->point[ptt->v[0]].tmp,mesh->point[ptt->v[1]].tmp
                ,mesh->point[ptt->v[2]].tmp,ptt->ref);
      } else {
        fwrite(&mesh->point[ptt->v[0]].tmp,sw,1,inm);
        fwrite(&mesh->point[ptt->v[1]].tmp,sw,1,inm);
        fwrite(&mesh->point[ptt->v[2]].tmp,sw,1,inm);
        fwrite(&ptt->ref,sw,1,inm);
      }
    }
    if ( ntreq && !mesh->info.nosurf ) {
      /* Don't save the required triangles when no surface remeshing (because
       * all the surface triangles are required). */
      if(!bin) {
        strcpy(&chaine[0],"\n\nRequiredTriangles\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%d \n",ntreq);
      } else {
        binch = 17; //ReqTriangles
        fwrite(&binch,sw,1,inm);
        bpos += 12+4*ntreq; //Pos
        fwrite(&bpos,sw,1,inm);
        fwrite(&ntreq,sw,1,inm);
      }
      for (k=0; k<=mesh->nt; k++) {
        ptt = &mesh->tria[k];
        if ( ptt->tag[0] & MG_REQ && ptt->tag[1] & MG_REQ && ptt->tag[2] & MG_REQ ) {
          if(!bin) {
            fprintf(inm,"%d \n",k);
          } else {
            fwrite(&k,sw,1,inm);
          }
        }
      }
    }
  }

  nr = nedreq = 0;
  if ( mesh->na ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nEdges\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",mesh->na);
    } else {
      binch = 5; //Edges
      fwrite(&binch,sw,1,inm);
      bpos += 12 + 3*4*mesh->na;//Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&mesh->na,sw,1,inm);
    }
    for (k=1; k<=mesh->na; k++) {
      if(!bin) {
        fprintf(inm,"%d %d %d \n",mesh->point[mesh->edge[k].a].tmp,
                mesh->point[mesh->edge[k].b].tmp,mesh->edge[k].ref);
      } else {
        fwrite(&mesh->point[mesh->edge[k].a].tmp,sw,1,inm);
        fwrite(&mesh->point[mesh->edge[k].b].tmp,sw,1,inm);
        fwrite(&mesh->edge[k].ref,sw,1,inm);
      }
      if ( mesh->edge[k].tag & MG_GEO ) nr++;
      if ( mesh->edge[k].tag & MG_REQ ) nedreq++;
    }

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
      na = 0;
      for (k=1; k<=mesh->na; k++) {
        na++;
        if ( mesh->edge[k].tag & MG_GEO ) {
          if(!bin) {
            fprintf(inm,"%d \n",na);
          } else {
            fwrite(&na,sw,1,inm);
          }
        }
      }
    }
    if ( nedreq && !mesh->info.nosurf ) {
      /* Don't save the required edges when no surface remeshing (because
       * all the surface edges are required). */
      if(!bin) {
        strcpy(&chaine[0],"\n\nRequiredEdges\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%d\n",nedreq);
      } else {
        binch = 16; //RequiredEdges
        fwrite(&binch,sw,1,inm);
        bpos += 12 + 4*nedreq;//Pos
        fwrite(&bpos,sw,1,inm);
        fwrite(&nedreq,sw,1,inm);
      }
      na = 0;
      for (k=1; k<=mesh->na; k++) {
        na++;
        if (  mesh->edge[k].tag & MG_REQ ) {
          if(!bin) {
            fprintf(inm,"%d \n",na);
          } else {
            fwrite(&na,sw,1,inm);
          }
        }
      }
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

  if(!bin) {
    strcpy(&chaine[0],"\n\nTetrahedra\n");
    fprintf(inm,"%s",chaine);
    fprintf(inm,"%d\n",ne);
  } else {
    binch = 8; //Tetra
    fwrite(&binch,sw,1,inm);
    bpos += 12 + 20*ne;//Pos
    fwrite(&bpos,sw,1,inm);
    fwrite((unsigned char*)&ne,sw,1,inm);
  }
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( MG_EOK(pt) ) {
      if(!bin) {
        fprintf(inm,"%d %d %d %d %d\n",mesh->point[pt->v[0]].tmp,mesh->point[pt->v[1]].tmp
                ,mesh->point[pt->v[2]].tmp,mesh->point[pt->v[3]].tmp,pt->ref);
      } else {
        fwrite(&mesh->point[pt->v[0]].tmp,sw,1,inm);
        fwrite(&mesh->point[pt->v[1]].tmp,sw,1,inm);
        fwrite(&mesh->point[pt->v[2]].tmp,sw,1,inm);
        fwrite(&mesh->point[pt->v[3]].tmp,sw,1,inm);
        fwrite(&pt->ref,sw,1,inm);
      }
    }
  }

  if ( nereq ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nRequiredTetrahedra\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",nereq);
    } else {
      binch = 12; //RequiredTetra
      fwrite(&binch,sw,1,inm);
      bpos += 12 + 4*nereq;//Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&nereq,sw,1,inm);
    }
    ne = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;
      ne++;
      if ( pt->tag & MG_REQ ) {
        if(!bin) {
          fprintf(inm,"%d \n",ne);
        } else {
          fwrite(&ne,sw,1,inm);
        }
      }
    }
  }

  if ( mesh->info.imprim ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",np,nc+nre);
    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d   RIDGES  %8d\n",mesh->na,nr);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",mesh->nt);
    fprintf(stdout,"     NUMBER OF ELEMENTS   %8d\n",mesh->ne);
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

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load metric field.
 *
 */
int MMG5_loadMet(MMG5_pMesh mesh,MMG5_pSol met) {
  FILE       *inm;
  float       fbuf[6],tmpf;
  double      dbuf[6],tmpd,lambda[3],eigenv[3][3];
  int         binch,bdim,iswp;
  int         i,k,bin,bpos;
  int         compute_hmin, compute_hmax;
  long        posnp;
  char        *ptr,data[128],chaine[128];

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
        return(0);
      }
    } else {
      bin = 1;
    }
  }
  else {
    if (!(inm = fopen(data,"r")) ) {
      fprintf(stderr,"  ** %s  NOT FOUND. USE DEFAULT METRIC.\n",data);
      return(0);
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);


  /* read solution or metric */
  if(!bin) {
    strcpy(chaine,"DDD");
    while(fscanf(inm,"%s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
      if(!strncmp(chaine,"Dimension",strlen("Dimension"))) {
        fscanf(inm,"%d",&met->dim);
        if(met->dim!=3) {
          fprintf(stdout,"BAD SOL DIMENSION : %d\n",met->dim);
          return(-1);
        }
        continue;
      } else if(!strncmp(chaine,"SolAtVertices",strlen("SolAtVertices"))) {
        fscanf(inm,"%d",&met->np);
        fscanf(inm,"%d",&met->type);
        if(met->type!=1) {
          fprintf(stdout,"SEVERAL SOLUTION => IGNORED : %d\n",met->type);
          return(-1);
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
    if(iswp) met->ver = _MMG5_swapbin(met->ver);
    while(fread(&binch,sw,1,inm)!=EOF && binch!=54 ) {
      if(iswp) binch=_MMG5_swapbin(binch);
      if(binch==54) break;
      if(binch==3) {  //Dimension
        fread(&bdim,sw,1,inm);  //NulPos=>20
        if(iswp) bdim=_MMG5_swapbin(bdim);
        fread(&met->dim,sw,1,inm);
        if(iswp) met->dim=_MMG5_swapbin(met->dim);
        if(met->dim!=3) {
          fprintf(stdout,"BAD SOL DIMENSION : %d\n",met->dim);
          printf("  Exit program.\n");
          return(-1);
        }
        continue;
      } else if(binch==62) {  //SolAtVertices
        fread(&binch,sw,1,inm); //NulPos
        if(iswp) binch=_MMG5_swapbin(binch);
        fread(&met->np,sw,1,inm);
        if(iswp) met->np=_MMG5_swapbin(met->np);
        fread(&met->type,sw,1,inm); //nb sol
        if(iswp) met->type=_MMG5_swapbin(met->type);
        if(met->type!=1) {
          fprintf(stdout,"SEVERAL SOLUTION => IGNORED : %d\n",met->type);
          return(-1);
        }
        fread(&met->size,sw,1,inm); //typsol
        if(iswp) met->size=_MMG5_swapbin(met->size);
        posnp = ftell(inm);
        break;
      } else {
        fread(&bpos,sw,1,inm); //Pos
        if(iswp) bpos=_MMG5_swapbin(bpos);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
      }
    }

  }
  if ( mesh->np != met->np ) {
    return(-1);
  }
  if ( mesh->info.lag == -1 ) {
    if(met->size!=1 && met->size!=3) {
      fprintf(stdout,"  ** DATA TYPE IGNORED %d \n",met->size);
      return(-1);
    }
    if(met->size > 1) met->size = 6;
  }
  else if ( met->size != 2 ) {
    return(-1);
  }

  met->npi = met->np;

  /* mem alloc */
  if ( met->m )  _MMG5_DEL_MEM(mesh,met->m,(met->size*met->npmax+1)*sizeof(double));
  met->npmax = mesh->npmax;

  _MMG5_ADD_MEM(mesh,(met->size*met->npmax+1)*sizeof(double),"initial solution",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(met->m,met->size*met->npmax+1,double);

  /* read mesh solutions */
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
          if(iswp) fbuf[0]=_MMG5_swapf(fbuf[0]);
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
          if(iswp) dbuf[0]=_MMG5_swapd(dbuf[0]);
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
  /* vector displacement only */
  else if(met->size==2) {
    met->size = 3;
    if ( met->ver == 1 ) {
      for (k=1; k<=met->np; k++) {
        for (i=1; i<=3; i++) {
          if(!bin){
            fscanf(inm,"%f",&fbuf[0]);
          } else {
            fread(&fbuf[0],sw,1,inm);
            if(iswp) fbuf[0]=_MMG5_swapf(fbuf[0]);
          }
          met->m[3*(k-1)+i] = fbuf[0];
        }
      }
    }
    else {
      for (k=1; k<=met->np; k++) {
        for (i=1; i<=3; i++) {
          if(!bin){
            fscanf(inm,"%lf",&dbuf[0]);
          } else {
            fread(&dbuf[0],sd,1,inm);
            if(iswp) dbuf[0]=_MMG5_swapd(dbuf[0]);
          }
          met->m[3*(k-1)+i] = dbuf[0];
        }
      }
    }
  }
  /* anisotropic metric */
  else { //met->size==6
    if ( met->ver == 1 ) {
      for (k=1; k<=met->np; k++) {
        if(!bin){
          for(i=0 ; i<met->size ; i++)
            fscanf(inm,"%f",&fbuf[i]);
        } else {
          for(i=0 ; i<met->size ; i++) {
            fread(&fbuf[i],sw,1,inm);
            if(iswp) fbuf[i]=_MMG5_swapf(fbuf[i]);
          }
        }
        tmpf    = fbuf[2];
        fbuf[2] = fbuf[3];
        fbuf[3] = tmpf;
        for (i=0; i<6; i++)  met->m[6*k+1+i] = fbuf[i];
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
            if(iswp) dbuf[i]=_MMG5_swapf(dbuf[i]);
          }
        }
        tmpd    = dbuf[2];
        dbuf[2] = dbuf[3];
        dbuf[3] = tmpd;
        for (i=0; i<met->size; i++)  met->m[6*k+1+i] = dbuf[i];
      }
    }
    for (k=1; k<=met->np; k++) {
      if ( ! _MMG5_eigenv(1,&met->m[6*k+1],lambda,eigenv ) ) {
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

  met->npi = met->np;
  fclose(inm);

  /* stats */
  if ( abs(mesh->info.imprim) > 3 ) {
    if ( met->size == 1 ) 
      fprintf(stdout,"     NUMBER OF SCALAR VALUES %8d\n",met->np);
    else if ( met->size == 3 )
      fprintf(stdout,"     NUMBER OF VECTOR VALUES %8d\n",met->np);
    else
      fprintf(stdout,"     NUMBER OF TENSOR VALUES %8d\n",met->np);
  }
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric.
 *
 */
int MMG5_saveMet(MMG5_pMesh mesh,MMG5_pSol met) {
  FILE*        inm;
  MMG5_pPoint  ppt;
  char        *ptr,data[128],chaine[128];
  int          binch,bpos,bin,np,k,typ;

  if ( !met->m || !met->nameout )  return(-1);
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
  binch=bpos=0;
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
        if(!bin) {
          fprintf(inm,"%.15lg \n ",met->m[k]);
        } else {
          fwrite((unsigned char*)&met->m[k],sd,1,inm);
        }
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
