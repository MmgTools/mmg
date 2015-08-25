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
 * \file mmg3d/mmgtols.c
 * \brief Tools for interfacing mmg with LS.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmg3d.h"

/** Conversion of a MMG5_pPoint (Pt ipm) into a LS_pPoint (Pt ipl) */
int _MMG5_ptMMGtoLS(MMG5_pMesh mesh,int ipm,LS_pMesh lsmesh,int ipl) {
  MMG5_pPoint    p0m;
  LS_pPoint      p0l;
  char           i;
  
  p0m = &mesh->point[ipm];
  p0l = &lsmesh->point[ipl];
  
  p0l->ref = p0m->ref;
  p0l->old = ipm;
  
  for(i=0; i<3; i++)
    p0l->c[i] = p0m->c[i];
  
  return(1);
}

/** Conversion of a MMG5_pTetra into a LS_pTetra */
int _MMG5_tetMMGtoLS(MMG5_pMesh mesh,int ielm,LS_pMesh lsmesh,int iell,int *perm) {
  MMG5_pTetra   ptm;
  LS_pTetra     ptl;
  char          i;
  
  ptm = &mesh->tetra[ielm];
  ptl = &lsmesh->tetra[iell];
  
  ptl->ref = 0; //ptm->ref;
  
  for(i=0; i<4; i++)
    ptl->v[i] = perm[ptm->v[i]];
  
  return(1);
}

/** Creation of a LS_pTria (number nt in lsmesh->tria) from face i in tetra k in mesh; 
    the numbers of vertices are turned into those of lsmesh thanks to table perm */
int _MMG5_creaLStria(MMG5_pMesh mesh,int iel,char i,LS_pMesh lsmesh,int nt,int ref,int *perm) {
  MMG5_pTetra   ptm;
  LS_pTria      pttl;
  char          j,i0;
  
  ptm  = &mesh->tetra[iel];
  pttl = &lsmesh->tria[nt];
  
  pttl->ref = ref;
  for(j=0; j<3; j++) {
    i0 = _MMG5_idir[i][j];
    pttl->v[j] = perm[ptm->v[i0]];
  }

  return(1);
}

/** Create lsmesh, suitable for the ELASTIC library, from the datum of mesh 
    ne = number of elements in the (new) mesh, 
    list = list of the tetras to retain, 
    perm = perm table for the points to retain, 
    np = number of points */
int _MMG5_iniLSmesh(MMG5_pMesh mesh,LS_pMesh lsmesh,int ne,int *list,int *perm,int np) {
  MMG5_pTetra    ptm;
  MMG5_pxTetra   pxt;
  MMG5_pPoint    p0m;
  LS_pTetra      ptl;
  LS_pPoint      p0l;
  int            k,ip,iel,jel,nt,*adja;
  char           i;
  
  nt = 0;
  
  lsmesh->dim = 3;
  lsmesh->np = np;
  lsmesh->npi= np;
  lsmesh->np2 = 0;
  lsmesh->ne = ne;
  lsmesh->na = 0;
  
  lsmesh->name = (char*)calloc(128,sizeof(char));
  lsmesh->point = (LS_Point*)calloc(lsmesh->np+1,sizeof(LS_Point));
  lsmesh->tetra = (LS_Tetra*)calloc(lsmesh->ne+1,sizeof(LS_Tetra));
  
  /* Setting of the name of lsmesh */
  strcat(lsmesh->name,"LS");
  strcat(lsmesh->name,mesh->nameout);
  
  /* Creation of points */
  for(k=1; k<=mesh->np; k++) {
    ip = perm[k];
    if ( !ip ) continue;
    if ( !_MMG5_ptMMGtoLS(mesh,k,lsmesh,ip) ) return(0);
  }
  
  /* Creation of tetras */
  for(k=1; k<=ne; k++) {
    iel = list[k];
    if ( !_MMG5_tetMMGtoLS(mesh,iel,lsmesh,k,perm) ) return(0);
  }
  
  /* Creation of the surface triangles for boundary conditions */
  /* Void loop for counting purposes */
  for(k=1; k<=ne; k++) {
    iel = list[k];
    ptm = &mesh->tetra[iel];
    adja = &mesh->adja[4*(iel-1)+1];
    if (ptm->xt) pxt = &mesh->xtetra[ptm->xt];
    
    for(i=0; i<4; i++) {
      if ( ptm->xt && (pxt->ftag[i] & MG_BDY) && (pxt->ref[i] == _MMG5_DISPREF) ) nt++;
      else {
        jel = adja[i] / 4;
        if ( !jel || (!mesh->tetra[jel].mark) ) nt++;
      }
    }
  }
  
  lsmesh->nt = nt;
  lsmesh->tria = (LS_Tria*)calloc(lsmesh->nt,sizeof(LS_Tria));

  /* Effective creation of triangles */
  nt = 0;
  
  for(k=1; k<=ne; k++) {
    iel = list[k];
    ptm = &mesh->tetra[iel];
    adja = &mesh->adja[4*(iel-1)+1];
    if (ptm->xt) pxt = &mesh->xtetra[ptm->xt];
    
    for(i=0; i<4; i++) {
      if ( ptm->xt && (pxt->ftag[i] & MG_BDY) && (pxt->ref[i] == _MMG5_DISPREF) ) {
        nt++;
        _MMG5_creaLStria(mesh,iel,i,lsmesh,nt,_MMG5_DISPREF,perm);
      }
      else {
        jel = adja[i] / 4;
        if ( !jel || (!mesh->tetra[jel].mark) ) {
          nt++;
          _MMG5_creaLStria(mesh,iel,i,lsmesh,nt,_LS_REFDIR,perm);
        }
      }
    }
  }
  
  return(1);
}

/** Create a LS_pSol file from the MMG5_pSol disp */
int _MMG5_creaLSdisp(MMG5_pSol disp,LS_pSol lsdisp,int np,int *perm) {
  LS_pCl   pcl;
  LS_pMat  pmat;
  int      k,ip;
  char     j;
  
  lsdisp->np    = np;
  lsdisp->na    = 0;
  lsdisp->dim   = 3;
  lsdisp->ver   = disp->ver;
  lsdisp->nbcl  = 2;
  lsdisp->nmat  = 1;
  lsdisp->cltyp = _LS_Tri;
  lsdisp->err = 1.e-6;
  lsdisp->nit = 10000;
  
  lsdisp->u = (double*)calloc(3*np,sizeof(double));
  
  /* TEMPORARY: same numbering procedure for LS and MG... maybe not to keep !! */
  for(k=1; k<=disp->np; k++) {
    ip = perm[k];
    if ( !ip ) continue;
    for (j=0; j<3; j++) lsdisp->u[3*(ip-1)+j] = disp->m[3*k+j];
  }
  
  lsdisp->cl  = (LS_pCl)calloc(2,sizeof(LS_Cl));
  lsdisp->mat = (LS_pMat)calloc(1,sizeof(LS_Mat));
  
  pcl  = &lsdisp->cl[0];
  
  pcl->typ = Dirichlet;
  pcl->ref = _MMG5_DISPREF;
  pcl->att = 'f';
  pcl->elt = _LS_Tri;
  
  pcl  = &lsdisp->cl[1];
  pcl->typ = Dirichlet;
  pcl->ref = _LS_REFDIR;
  pcl->att = 'v';
  pcl->u[0] = pcl->u[1] = pcl->u[2] = 0.0;
  pcl->elt = _LS_Tri;
  
  pmat = &lsdisp->mat[0];

  pmat->ref = 0;
  pmat->lambda = _LS_LAMBDA;
  pmat->mu  = _LS_MU;
  
  return(1);
}

/** Extracts a layer around the surface triangles tagged dispref, 
    and converts it into a LS_pMesh - truncate the MMG5_sol file into a LS_sol file in the meantime */
int _MMG5_packLS(MMG5_pMesh mesh,LS_pMesh lsmesh,MMG5_pSol disp,LS_pSol lsdisp) {
  MMG5_pTetra   pt,pt1;
  MMG5_pxTetra  pxt;
  int           *list,ilist,ilisto,ilistck,k,kk,nlay,n,iel,jel,*adja,*perm,np,npf;
  char          i,j;
  
  nlay = 5;
  npf = 0;
  list = (int*)calloc(mesh->ne+1,sizeof(int));
  perm = (int*)calloc(mesh->np+1,sizeof(int));
  ilist = ilisto = ilistck = 0;
  
  for (k=1; k<=mesh->ne; k++)
    mesh->tetra[k].mark = 0;      // A faire fusionner avec celui de mmg3d3.c
  
  /* Step 1: pile all the tetras containing a triangle with ref DISPREF */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];
    
    for(i=0; i<4; i++) {
      if ( (pxt->ftag[i] & MG_BDY) && (pxt->ref[i] == _MMG5_DISPREF) ) {
        ilist++;
        list[ilist] = k;
        pt->mark = 1;
        
        for(j=0; j<4; j++) {
          np = pt->v[j];
          if ( !perm[np] ) {
            npf++;
            perm[np] = npf;
          }
        }
        
        break;
      }
    }
  }
  
  /* Step 2: create a layer around these tetras */
  for(n=0; n<nlay; n++) {
    ilistck = ilisto;
    ilisto = ilist;
    
    for(k=ilistck+1; k<=ilisto; k++) {
      iel = list[k];
      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];
     
      for(i=0; i<4; i++) {
        jel = adja[i] / 4;
        if ( !jel ) continue;
        pt1 = &mesh->tetra[jel];
        if ( MG_EOK(pt1) && (!pt1->mark ) ) {
          ilist++;
          assert( ilist <= mesh->ne );
          pt1->mark = 1;
          list[ilist] = jel;
          
          for(j=0; j<4; j++) {
            np = pt1->v[j];
            if ( !perm[np] ) {
              npf++;
              perm[np] = npf;
            }
          }
        }
      }
    }
  }
  
  /* Step 3: conversion into a LS_pMesh -- creation of the surface triangles for boundary conditions */
  if ( !_MMG5_iniLSmesh(mesh,lsmesh,ilist,list,perm,npf) ) {
    fprintf(stdout,"  ## Problem in fn MMG5_iniLSmesh. Exiting.\n");
    return(0);
  }
  
  /* Step 4: truncate .sol file */
  if ( !_MMG5_creaLSdisp(disp,lsdisp,np,perm) ) {
    fprintf(stdout,"  ## Problem in fn MMG5_creaLSdisp. Exiting.\n");
    return(0);
  }
  
  /* Save mesh */
  /*if ( !_MMG5_saveLSmesh(lsmesh) ) {
    fprintf(stdout,"  ## Problem in fn MMG5_saveLSmesh. Exiting.\n");
    return(0);
  }*/
  
  free(list);
  free(perm);
  
  return(1);
}

/** Redistribute the (extended) displacement on the LS mesh to the whole MMG mesh */
int _MMG5_unpackLS(MMG5_pMesh mesh,LS_pMesh lsmesh,MMG5_pSol disp,LS_pSol lsdisp) {
  LS_pPoint   p0;
  int         k,ip;
  char        j;
  
  for(k=1; k<=mesh->np; k++) {
    for(j=0; j<3; j++)
      disp->m[3*k+j] = 0.0;
  }
  
  for(k=1; k<=lsmesh->np; k++) {
    p0 = &lsmesh->point[k];
    ip = p0->old;

    for(j=0; j<3; j++)
      disp->m[3*ip+j] = lsdisp->u[3*(k-1)+j];
  }

  return(1);
}

/** For debugging purposes: save lsmesh */
int _MMG5_saveLSmesh(LS_pMesh lsmesh) {
  FILE        *out;
  LS_pTetra   pt;
  LS_pTria    ptt;
  LS_pPoint   p0;
  int         k;
  char        i;
  
  out = fopen(lsmesh->name,"w");
  
  fprintf(out,"MeshVersionFormatted 1\n\nDimension\n%d\n\n",lsmesh->dim);
  fprintf(out,"Vertices\n%d\n",lsmesh->np);
  
  /* Print points */
  for(k=1; k<= lsmesh->np; k++) {
    p0 = &lsmesh->point[k];
    fprintf(out,"%f %f %f %d\n",p0->c[0],p0->c[1],p0->c[2],p0->ref);
  }
  
  /* Print Tetrahedra */
  fprintf(out,"\nTetrahedra\n%d\n",lsmesh->ne);
  
  for(k=1; k<= lsmesh->ne; k++) {
    pt = &lsmesh->tetra[k];
    fprintf(out,"%d %d %d %d %d\n",pt->v[0],pt->v[1],pt->v[2],pt->v[3],pt->ref);
  }
  
  /* Print Triangles */
  fprintf(out,"\nTriangles\n%d\n",lsmesh->nt);
  
  for(k=1; k<= lsmesh->nt; k++) {
    ptt = &lsmesh->tria[k];
    fprintf(out,"%d %d %d %d\n",ptt->v[0],ptt->v[1],ptt->v[2],ptt->ref);
  }
  
  fprintf(out,"\nEnd");
  fclose(out);
  
  return(1);
}

/** For debugging purposes: save disp */
int _MMG5_saveDisp(MMG5_pMesh mesh,MMG5_pSol disp) {
  FILE        *out;
  int         k;
  char        j,data[256],*ptr;

  strcpy(data,disp->namein);
  ptr = strstr(data,".sol");
  *ptr = '\0';
  strcat(data,".o.disp.sol");
  
  out = fopen(data,"w");
  
  fprintf(out,"MeshVersionFormatted 1\n\nDimension\n%d\n\n",disp->dim);
  fprintf(out,"SolAtVertices\n%d\n 1 2\n",disp->np);
  
  /* Print solutions */
  for(k=1; k<= disp->np; k++) {
    fprintf(out,"%f %f %f\n",disp->m[3*k+0],disp->m[3*k+1],disp->m[3*k+2]);
  }
  
  fprintf(out,"\nEnd");
  fclose(out);
  
  return(1);
}