/* =============================================================================
**  This file is part of the MMG3D 5 software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
**
**  MMG3D 5 is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  MMG3D 5 is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with MMG3D 5 (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the MMG3D 5 distribution only if you accept them.
** =============================================================================
*/

/**
 * \file mmg3d/mmg3d1_delone.c
 * \brief Perform volume and surface mesh adaptation in delaunay mode.
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 *
 * Perform volume and surface mesh adaptation in delaunay mode (\a
 * PATTERN preprocessor flag set to OFF).
 *
 */
#include "mmg3d.h"

char  ddb;

#define LOPTLDEL     1.41//1.41
#define LOPTSDEL     0.6
int MMG_npuiss,MMG_nvol,MMG_npres;

/** Internal edge flipping */
/*static*/ int swptetdel(pMesh mesh,pSol met,double crit,pBucket bucket) {
  pTetra   pt;
  pxTetra  pxt;
  int      list[LMAX+2],ilist,k,it,nconf,maxit,ns,nns,ier;
  char     i;
  double critloc;

  maxit = 2;
  it = nns = 0;
  critloc = crit;
  do {
    ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
      if ( pt->qual > 0.0288675 /*0.6/ALPHAD*/ )  continue;
      //if ( pt->qual < 0.009) critloc = 1.01;
      //else critloc = crit;

      for (i=0; i<6; i++) {
        /* Prevent swap of a ref or tagged edge */
        if ( pt->xt ) {
          pxt = &mesh->xtetra[pt->xt];
          if ( pxt->edg[i] || pxt->tag[i] ) continue;
        }

        nconf = chkswpgen(mesh,k,i,&ilist,list,critloc);
        if ( nconf ) {
#ifdef PATTERN
          printf("ERROR: function not available in pattern mode. Exiting\n");
          exit(EXIT_FAILURE);
#else
          ier = swpgen(mesh,met,nconf,ilist,list,bucket);
#endif
          if ( ier > 0 )  ns++;
          else if ( ier < 0 ) return(-1);
          break;
        }
      }
    }
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8d edge swapped\n",nns);

  return(nns);
}

/** Analyze tetrahedra and move points so as to make mesh more uniform */
/*static*/ int movtetdel(pMesh mesh,pSol met,int maxitin) {
  pTetra        pt;
  pPoint        ppt;
  pxTetra       pxt;
  double        *n;
  int           i,k,ier,nm,nnm,ns,lists[LMAX+2],listv[LMAX+2],ilists,ilistv,it;
  int           improve;
  unsigned char j,i0,base;
  int internal,maxit;

  if(maxitin<0) {
    internal = 0;
    maxit = abs(maxitin);
  } else {
    maxit = maxitin;
    internal=1;
  }
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** OPTIMIZING MESH\n");

  base = 1;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = base;

  it = nnm = 0;
  do {
    base++;
    nm = ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;

      /* point j on face i */
      for (i=0; i<4; i++) {
        for (j=0; j<3; j++) {
          if ( pt->xt ) {
            pxt = &mesh->xtetra[pt->xt];
            if ( pxt->tag[iarf[i][j]] & MG_REQ )  continue;
          }
          else  pxt = 0;
          i0  = idir[i][j];
          ppt = &mesh->point[pt->v[i0]];
          if ( ppt->flag == base )  continue;
          else if ( MG_SIN(ppt->tag) )  continue;
#ifdef SINGUL
          else if ( ppt->tag & MG_SGL )  continue;
          else if ( mesh->info.sing && pt->xt && (pxt->tag[iarf[i][j]] & MG_SGL) )
            continue;
#endif
          if ( maxit != 1 ) {
            ppt->flag = base;
            improve   = 1;
          }
          else {
            improve = 0;
          }
          ier = 0;
          if ( ppt->tag & MG_BDY ) {
            /* Catch a boundary point by a boundary face */
            if ( !pt->xt || !(MG_BDY & pxt->ftag[i]) )  continue;
            else if( ppt->tag & MG_NOM ){
              if( mesh->adja[4*(k-1)+1+i] ) continue;
              if( !(ier=bouleext(mesh,k,i0,i,listv,&ilistv,lists,&ilists)) )
                continue;
              else if ( ier>0 )
                ier = movbdynompt(mesh,listv,ilistv,lists,ilists);
              else  return(-1);
            }
            else if ( ppt->tag & MG_GEO ) {
              if ( !(ier=boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists)) )
                continue;
              else if ( ier>0)
                ier = movbdyridpt(mesh,listv,ilistv,lists,ilists);
              else
                return(-1);
            }
            else if ( ppt->tag & MG_REF ) {
              if ( !(ier=boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists)) )
                continue;
              else if ( ier>0 )
                ier = movbdyrefpt(mesh,listv,ilistv,lists,ilists);
              else
                return(-1);
            }
            else {
              if ( !(ier=boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists)) )
                continue;
              else if ( ier<0 )
                return(-1);

              n = &(mesh->xpoint[ppt->xp].n1[0]);
              if ( !directsurfball(mesh, pt->v[i0],lists,ilists,n) )  continue;
              ier = movbdyregpt(mesh,listv,ilistv,lists,ilists);
              if ( ier )  ns++;
            }
          }
          else if(internal){
            ilistv = boulevolp(mesh,k,i0,listv);
            if ( !ilistv )  continue;
            ier = movintpt(mesh,listv,ilistv,improve);
          }
          if ( ier ) {
            nm++;
            if(maxit==1){
              ppt->flag = base;
            }
          }
        }
      }
    }
    nnm += nm;
    if ( mesh->info.ddebug )  fprintf(stdout,"     %8d moved, %d geometry\n",nm,ns);

    outqua(mesh,met);
  }
  while( ++it < maxit && nm > 0 );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nnm )
    fprintf(stdout,"     %8d vertices moved, %d iter.\n",nnm,it);

  return(nnm);
}

static inline int boucle_for(pMesh mesh, pSol met,pBucket bucket,int ne,int* ifilt,int* ns,int* nc,int* warn,int it) {
  pTetra     pt;
  pxTetra    pxt;
  Tria       ptt;
  pPoint     p0,p1,ppt;
  pxPoint    pxp;
  double     dd,len,lmax,o[3],to[3],ro[3],no1[3],no2[3],v[3];
  int        k,ip,ip1,ip2,list[LMAX+2],ilist,ref;
  char       imax,tag,j,i,i1,i2,ifa0,ifa1;
  int        lon,ret,ier;
  double     lmin;
  int        imin,iq,nnc,nns,nnf,nnm;
  int        ii,MMG_npd;
  double     maxgap,lmaxtet,lmintet;
  int nconf,imaxtet,imintet;
  double critloc;

  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt)  || (pt->tag & MG_REQ) )   continue;
    //  if(it>1)
    //if(pt->qual > 0.038/*0.0288675*/ /*0.6/ALPHAD*/) continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /* 1) find longest and shortest edge  and try to manage it*/
    imax = -1; lmax = 0.0;
    imin = -1; lmin = DBL_MAX;
    for (ii=0; ii<6; ii++) {
      if ( pt->xt && (pxt->tag[ii] & MG_REQ) )  continue;
      ip1  = iare[ii][0];
      ip2  = iare[ii][1];
      len = lenedg(mesh,met,pt->v[ip1],pt->v[ip2]);
      if ( len > lmax ) {
	lmax = len;
	imax = ii;
      }
      if ( len < lmin ) {
	lmin = len;
	imin = ii;
      }
    }
    if ( imax==-1 )
      fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are boundary and required\n",
	      __FILE__,__LINE__,k);
    if ( imin==-1 )
      fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are boundary and required\n",
	      __FILE__,__LINE__,k);
    /* imax = ii; */
    /* lmax = len; */
    /* imin = ii; */
    /* lmin = len; */
    if ( lmax >= LOPTLDEL )  {
      /* proceed edges according to lengths */
      ifa0 = ifar[imax][0];
      ifa1 = ifar[imax][1];
      i  = (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
      j  = iarfinv[i][imax];
      i1 = idir[i][inxt2[j]];
      i2 = idir[i][iprv2[j]];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      p0  = &mesh->point[ip1];
      p1  = &mesh->point[ip2];

      /* Case of a boundary face */
      if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
	if ( !(MG_GET(pxt->ori,i)) ) continue;
	ref = pxt->edg[iarf[i][j]];
	tag = pxt->tag[iarf[i][j]];
	if ( tag & MG_REQ )  continue;
	tag |= MG_BDY;
	ilist = coquil(mesh,k,imax,list);
	if ( !ilist )  continue;
	else if ( ilist<0 ) return(-1);
	if ( tag & MG_NOM ){
	  if( !BezierNom(mesh,ip1,ip2,0.5,o,no1,to) )
	    continue;
	  else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
	    tet2tri(mesh,k,i,&ptt);
	    nortri(mesh,&ptt,no1);
	    if ( !MG_GET(pxt->ori,i) ) {
	      no1[0] *= -1.0;
	      no1[1] *= -1.0;
	      no1[2] *= -1.0;
	    }
	  }
	}
	else if ( tag & MG_GEO ) {
	  if ( !BezierRidge(mesh,ip1,ip2,0.5,o,no1,no2,to) )
	    continue;
	  if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
	    tet2tri(mesh,k,i,&ptt);
	    nortri(mesh,&ptt,no1);
	    no2[0] = to[1]*no1[2] - to[2]*no1[1];
	    no2[1] = to[2]*no1[0] - to[0]*no1[2];
	    no2[2] = to[0]*no1[1] - to[1]*no1[0];
	    dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
	    if ( dd > EPSD2 ) {
	      dd = 1.0 / sqrt(dd);
	      no2[0] *= dd;
	      no2[1] *= dd;
	      no2[2] *= dd;
	    }
	  }
	}
	else if ( tag & MG_REF ) {
	  if ( !BezierRef(mesh,ip1,ip2,0.5,o,no1,to) )
	    goto collapse;//continue;
	}
	else {
	  //CECILE : je comprend pas pourquoi la normale est mauvaise a la fin
	  //goto collapse;
	  if ( !norface(mesh,k,i,v) )  goto collapse;//continue;
	  if ( !BezierReg(mesh,ip1,ip2,0.5,v,o,no1) ) goto collapse;

	}
	ier = simbulgept(mesh,list,ilist,o);
	if ( !ier ) {
	  ier = dichoto1b(mesh,list,ilist,o,ro);
	  memcpy(o,ro,3*sizeof(double));
	}
	ip = newPt(mesh,o,tag);

	if ( !ip ){
	  /* reallocation of point table */
#ifndef PATTERN

	  POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
				   *warn=1;
				   goto collapse//break
				   ,o,tag);

#else
	  printf("ERROR: function not available in delaunay mode. Exiting\n");
	  exit(EXIT_FAILURE);
#endif
	}
	//CECILE
	if ( met->m )
	  met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
	//CECILE
	ier = split1b(mesh,met,list,ilist,ip,1);
	/* if we realloc memory in split1b pt and pxt pointers are not valid */
	pt = &mesh->tetra[k];
	pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

	if ( ier < 0 ) {
	  fprintf(stdout,"  ## Error: unable to split.\n");
	  return(-1);
	}
	else if ( !ier ) {
	  delPt(mesh,ip);
	  goto collapse;//continue;
	} else {
	  (*ns)++;
	  //addBucket(mesh,bucket,ip);

	  ppt = &mesh->point[ip];
	  if ( MG_EDG(tag) || (tag & MG_NOM) )
	    ppt->ref = ref;
	  else
	    ppt->ref = pxt->ref[i];
	  ppt->tag = tag;
	  if ( met->m )
	    met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

	  pxp = &mesh->xpoint[ppt->xp];
	  if ( tag & MG_NOM ){
	    memcpy(pxp->n1,no1,3*sizeof(double));
	    memcpy(pxp->t,to,3*sizeof(double));
	  }
	  else if ( tag & MG_GEO ) {
	    memcpy(pxp->n1,no1,3*sizeof(double));
	    memcpy(pxp->n2,no2,3*sizeof(double));
	    memcpy(pxp->t,to,3*sizeof(double));
	  }
	  else if ( tag & MG_REF ) {
	    memcpy(pxp->n1,no1,3*sizeof(double));
	    memcpy(pxp->t,to,3*sizeof(double));
	  }
	  else
	    memcpy(pxp->n1,no1,3*sizeof(double));
	}
	continue;//break;//imax continue;
      }
      else if(pt->xt){
	if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) {
	  continue;
	}
	ilist = coquil(mesh,k,imax,list);
	if ( !ilist )    continue;
	else if ( ilist<0 ) return(-1);
	o[0] = 0.5*(p0->c[0] + p1->c[0]);
	o[1] = 0.5*(p0->c[1] + p1->c[1]);
	o[2] = 0.5*(p0->c[2] + p1->c[2]);
	ip = newPt(mesh,o,MG_NOTAG);

	if ( !ip )  {
	  /* reallocation of point table */
#ifndef PATTERN
	  POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
				   *warn=1;
				   goto collapse//break
				   ,o,MG_NOTAG);
#else
	  printf("ERROR: function not available in delaunay mode. Exiting\n");
	  exit(EXIT_FAILURE);
#endif
	}
	//CECILE
	if ( met->m )
	  met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
	//CECILE
	ier = split1b(mesh,met,list,ilist,ip,1);
	if ( ier < 0 ) {
	  fprintf(stdout,"  ## Error: unable to split.\n");
	  return(-1);
	}
	else if ( !ier ) { //Et on teste pas du tout les qualités ici ?
	  delPt(mesh,ip);
	  goto collapse;//continue;
	}
	else {
	  ppt = &mesh->point[ip];
	  met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
	  addBucket(mesh,bucket,ip);
	  (*ns)++;
	  continue;//break;//imax continue;
	}
	printf("on doit pas passer la\n");
	/* Case of an internal face */
      } else {
	/*TEST POUR LES ARETES AYANT DEUX POINTS DE PEAU...*/
	/*            if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) { */
	/*        //tente le swap si arete pas frontiere */
	/*        ilist = 0; */
	/*        critloc = 1.053; */
	/*        nconf = chkswpgen(mesh,k,imax,&ilist,list,critloc); */
	/*        if ( nconf ) { */
	/* #ifdef PATTERN */
	/*    printf("ERROR: function not available in pattern mode. Exiting\n"); */
	/*    exit(EXIT_FAILURE); */
	/* #else */
	/*    ier = swpgen(mesh,met,nconf,ilist,list,bucket); */
	/* #endif */
	/*    if ( ier > 0 )  { */
	/*      //printf("on a swappe"); */
	/*      continue; */
	/*    } */
	/*    else if ( ier < 0 ) { */
	/*      //printf("pas swap"); */
	/*      goto collapse; */
	/*    } */
	/*        } else { */
	/*    goto collapse; */
	/*        } */
	/*      } */
	/*FIN DU TEST*/
	ilist = coquil(mesh,k,imax,list);
	if ( !ilist )    continue;
	else if ( ilist<0 ) return(-1);
	else if(ilist%2) goto collapse; //bdry edge
	o[0] = 0.5*(p0->c[0] + p1->c[0]);
	o[1] = 0.5*(p0->c[1] + p1->c[1]);
	o[2] = 0.5*(p0->c[2] + p1->c[2]);
	ip = newPt(mesh,o,MG_NOTAG);

	if ( !ip )  {
	  /* reallocation of point table */
#ifndef PATTERN
	  POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
				   *warn=1;
				   goto collapse//break
				   ,o,MG_NOTAG);

#else
	  printf("ERROR: function not available in delaunay mode. Exiting\n");
	  exit(EXIT_FAILURE);
#endif
	}
	//CECILE
	if ( met->m )
	  met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
	//CECILE
	//LA DELONE

	if ( /*it &&*/ !buckin_iso(mesh,met,bucket,ip) ) {
	  delPt(mesh,ip);
	  (*ifilt)++;
	  goto collapse;////continue;
	} else {
	  lon = cavity(mesh,met,k,ip,list,ilist/2);
	  if ( lon < 1 ) {
	    MMG_npd++;
	    delPt(mesh,ip);
	    goto collapse;//continue;
	  } else {
	    ret = delone(mesh,met,ip,list,lon);
	    if ( ret > 0 ) {
	      ppt = &mesh->point[ip];
	      met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
	      //chkmsh(mesh,0,0);
	      addBucket(mesh,bucket,ip);
	      (*ns)++;
	      continue;//break;//imax continue;
	    }
	    else if ( ret == 0 ) {
	      MMG_npd++;
	      delPt(mesh,ip);
	      goto collapse;//continue;
	    }
	    else { /*allocation problem ==> saveMesh*/
	      return(0);
	      /* MMG_npd++; */
	      /* delPt(mesh,ip); */
	      /* goto collapse;//continue; */
	    }
	    printf("on passe pas la1\n");
	  }
	  printf("on passe pas la2\n");
	}
	printf("on passe pas la3\n");
      }
      printf("on passe pas la3\n");
    }
  collapse:
    if(lmin <= LOPTSDEL) {// continue;
      ifa0 = ifar[imin][0];
      ifa1 = ifar[imin][1];
      i  =  (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
      j  = iarfinv[i][imin];
      i1 = idir[i][inxt2[j]];
      i2 = idir[i][iprv2[j]];
      ip = pt->v[i1];
      iq = pt->v[i2];
      p0 = &mesh->point[ip];
      p1 = &mesh->point[iq];

      if ( (p0->tag > p1->tag) || (p0->tag & MG_REQ) )  continue;

      /* Case of a boundary face */
      ilist = 0;
      if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
	tag = pxt->tag[iarf[i][j]];
	if ( tag & MG_REQ )  continue;
	tag |= MG_BDY;
	if ( p0->tag > tag )   continue;
	if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) continue;
	ilist = chkcol_bdy(mesh,k,i,j,list);
	if ( ilist > 0 ) {
	  ier = colver(mesh,list,ilist,i2);
	  //nc += ier;
	  if ( ier < 0 ) return(-1);
	  else if(ier) {
	    //delBucket(mesh,bucket,ier);
	    delPt(mesh,ier);
	    (*nc)++;
	    continue;//break;//imax continue;
	  }
	}
	else if (ilist < 0 )  return(-1);
      }
      /* Case of an internal face */
      else {
	if ( p0->tag & MG_BDY )  continue;
	ilist = chkcol_int(mesh,met,k,i,j,list,2);
	if ( ilist > 0 ) {
	  ier = colver(mesh,list,ilist,i2);
	  if ( ilist < 0 ) continue;
	  //nc += ier;
	  if ( ier < 0 ) return(-1);
	  else if(ier) {
	    delBucket(mesh,bucket,ier);
	    delPt(mesh,ier);
	    (*nc)++;
	    continue;//break;//imax continue;
	  }
	}
	else if (ilist < 0 )  return(-1);
      }

      // }//end for ii
    } //end if lmin < LOPTSDEL
    /*2) longest and shortest edges are stucked => try another edges*/
    imaxtet = imax;
    imintet = imin;
    lmaxtet = lmax;
    lmintet = lmin;
    for (ii=0; ii<6; ii++) {
      if ( pt->xt && (pxt->tag[ii] & MG_REQ) )  continue;
      if ( (ii==imintet) && (lmintet < LOPTSDEL)) continue;
      if ( (ii==imaxtet) && (lmaxtet > LOPTLDEL) ) continue;

      ip1  = iare[ii][0];
      ip2  = iare[ii][1];
      len = lenedg(mesh,met,pt->v[ip1],pt->v[ip2]);
	  
      imax = ii;
      lmax = len;
      imin = ii;
      lmin = len;
      if ( lmax >= LOPTLDEL )  {
	/* proceed edges according to lengths */
	ifa0 = ifar[imax][0];
	ifa1 = ifar[imax][1];
	i  = (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
	j  = iarfinv[i][imax];
	i1 = idir[i][inxt2[j]];
	i2 = idir[i][iprv2[j]];
	ip1 = pt->v[i1];
	ip2 = pt->v[i2];
	p0  = &mesh->point[ip1];
	p1  = &mesh->point[ip2];

	/* Case of a boundary face */
	if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
	  if ( !(MG_GET(pxt->ori,i)) ) continue;
	  ref = pxt->edg[iarf[i][j]];
	  tag = pxt->tag[iarf[i][j]];
	  if ( tag & MG_REQ )  continue;
	  tag |= MG_BDY;
	  ilist = coquil(mesh,k,imax,list);
	  if ( !ilist )  continue;
	  else if ( ilist<0 ) return(-1);
	  if ( tag & MG_NOM ){
	    if( !BezierNom(mesh,ip1,ip2,0.5,o,no1,to) )
	      continue;
	    else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
	      tet2tri(mesh,k,i,&ptt);
	      nortri(mesh,&ptt,no1);
	      if ( !MG_GET(pxt->ori,i) ) {
		no1[0] *= -1.0;
		no1[1] *= -1.0;
		no1[2] *= -1.0;
	      }
	    }
	  }
	  else if ( tag & MG_GEO ) {
	    if ( !BezierRidge(mesh,ip1,ip2,0.5,o,no1,no2,to) )
	      continue;
	    if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
	      tet2tri(mesh,k,i,&ptt);
	      nortri(mesh,&ptt,no1);
	      no2[0] = to[1]*no1[2] - to[2]*no1[1];
	      no2[1] = to[2]*no1[0] - to[0]*no1[2];
	      no2[2] = to[0]*no1[1] - to[1]*no1[0];
	      dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
	      if ( dd > EPSD2 ) {
		dd = 1.0 / sqrt(dd);
		no2[0] *= dd;
		no2[1] *= dd;
		no2[2] *= dd;
	      }
	    }
	  }
	  else if ( tag & MG_REF ) {
	    if ( !BezierRef(mesh,ip1,ip2,0.5,o,no1,to) )
	      goto collapse2;//continue;
	  }
	  else {
	    //CECILE : je comprend pas pourquoi la normale est mauvaise a la fin
	    //goto collapse;
	    if ( !norface(mesh,k,i,v) )  goto collapse2;//continue;
	    if ( !BezierReg(mesh,ip1,ip2,0.5,v,o,no1) ) goto collapse2;

	  }
	  ier = simbulgept(mesh,list,ilist,o);
	  if ( !ier ) {
	    ier = dichoto1b(mesh,list,ilist,o,ro);
	    memcpy(o,ro,3*sizeof(double));
	  }
	  ip = newPt(mesh,o,tag);

	  if ( !ip ){
	    /* reallocation of point table */
#ifndef PATTERN
	    POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
				     *warn=1;
				     goto collapse2//break
				     ,o,tag);
#else
	    printf("ERROR: function not available in delaunay mode. Exiting\n");
	    exit(EXIT_FAILURE);
#endif
	  }
	  //CECILE
	  if ( met->m )
	    met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
	  //CECILE
	  ier = split1b(mesh,met,list,ilist,ip,1);
	  /* if we realloc memory in split1b pt and pxt pointers are not valid */
	  pt = &mesh->tetra[k];
	  pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

	  if ( ier < 0 ) {
	    fprintf(stdout,"  ## Error: unable to split.\n");
	    return(-1);
	  }
	  else if ( !ier ) {
	    delPt(mesh,ip);
	    goto collapse2;//continue;
	  } else {
	    (*ns)++;
	    //addBucket(mesh,bucket,ip);

	    ppt = &mesh->point[ip];
	    if ( MG_EDG(tag) || (tag & MG_NOM) )
	      ppt->ref = ref;
	    else
	      ppt->ref = pxt->ref[i];
	    ppt->tag = tag;
	    if ( met->m )
	      met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);

	    pxp = &mesh->xpoint[ppt->xp];
	    if ( tag & MG_NOM ){
	      memcpy(pxp->n1,no1,3*sizeof(double));
	      memcpy(pxp->t,to,3*sizeof(double));
	    }
	    else if ( tag & MG_GEO ) {
	      memcpy(pxp->n1,no1,3*sizeof(double));
	      memcpy(pxp->n2,no2,3*sizeof(double));
	      memcpy(pxp->t,to,3*sizeof(double));
	    }
	    else if ( tag & MG_REF ) {
	      memcpy(pxp->n1,no1,3*sizeof(double));
	      memcpy(pxp->t,to,3*sizeof(double));
	    }
	    else
	      memcpy(pxp->n1,no1,3*sizeof(double));
	  }
	  break;//imax continue;
	}
	else if(pt->xt){
	  if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) {
	    continue;
	  }
	  ilist = coquil(mesh,k,imax,list);
	  if ( !ilist )    continue;
	  else if ( ilist<0 ) return(-1);
	  o[0] = 0.5*(p0->c[0] + p1->c[0]);
	  o[1] = 0.5*(p0->c[1] + p1->c[1]);
	  o[2] = 0.5*(p0->c[2] + p1->c[2]);
	  ip = newPt(mesh,o,MG_NOTAG);

	  if ( !ip )  {
	    /* reallocation of point table */
#ifndef PATTERN
	    POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
				     *warn=1;
				     goto collapse2//break
				     ,o,MG_NOTAG);
#else
	    printf("ERROR: function not available in delaunay mode. Exiting\n");
	    exit(EXIT_FAILURE);
#endif
	  }
	  //CECILE
	  if ( met->m )
	    met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
	  //CECILE
	  ier = split1b(mesh,met,list,ilist,ip,1);
	  if ( ier < 0 ) {
	    fprintf(stdout,"  ## Error: unable to split.\n");
	    return(-1);
	  }
	  else if ( !ier ) { //Et on teste pas du tout les qualités ici ?
	    delPt(mesh,ip);
	    goto collapse2;//continue;
	  }
	  else {
	    ppt = &mesh->point[ip];
	    met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
	    addBucket(mesh,bucket,ip);
	    (*ns)++;
	    break;//imax continue;
	  }
	  printf("on doit pas passer la\n");
	  /* Case of an internal face */
	} else {
	  /*TEST POUR LES ARETES AYANT DEUX POINTS DE PEAU...*/
	  /*            if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) { */
	  /*        //tente le swap si arete pas frontiere */
	  /*        ilist = 0; */
	  /*        critloc = 1.053; */
	  /*        nconf = chkswpgen(mesh,k,imax,&ilist,list,critloc); */
	  /*        if ( nconf ) { */
	  /* #ifdef PATTERN */
	  /*    printf("ERROR: function not available in pattern mode. Exiting\n"); */
	  /*    exit(EXIT_FAILURE); */
	  /* #else */
	  /*    ier = swpgen(mesh,met,nconf,ilist,list,bucket); */
	  /* #endif */
	  /*    if ( ier > 0 )  { */
	  /*      //printf("on a swappe"); */
	  /*      continue; */
	  /*    } */
	  /*    else if ( ier < 0 ) { */
	  /*      //printf("pas swap"); */
	  /*      goto collapse; */
	  /*    } */
	  /*        } else { */
	  /*    goto collapse; */
	  /*        } */
	  /*      } */
	  /*FIN DU TEST*/
	  ilist = coquil(mesh,k,imax,list);
	  if ( !ilist )    continue;
	  else if ( ilist<0 ) return(-1);
	  else if(ilist%2) goto collapse2; //bdry edge
	  o[0] = 0.5*(p0->c[0] + p1->c[0]);
	  o[1] = 0.5*(p0->c[1] + p1->c[1]);
	  o[2] = 0.5*(p0->c[2] + p1->c[2]);
	  ip = newPt(mesh,o,MG_NOTAG);

	  if ( !ip )  {
	    /* reallocation of point table */
#ifndef PATTERN

	    POINT_AND_BUCKET_REALLOC(mesh,met,ip,mesh->gap,
				     *warn=1;
				     goto collapse2//break
				     ,o,MG_NOTAG);
#else
	    printf("ERROR: function not available in delaunay mode. Exiting\n");
	    exit(EXIT_FAILURE);
#endif
	  }
	  //CECILE
	  if ( met->m )
	    met->m[ip] = 0.5 * (met->m[ip1]+met->m[ip2]);
	  //CECILE
	  //LA DELONE
	  if ( /*lmax>4 &&*/ /*it &&*/  !buckin_iso(mesh,met,bucket,ip) ) {
	    delPt(mesh,ip);
	    (*ifilt)++;
	    goto collapse2;////continue;
	  } else {
	    lon = cavity(mesh,met,k,ip,list,ilist/2);
	    if ( lon < 1 ) {
	      MMG_npd++;
	      delPt(mesh,ip);
	      goto collapse2;//continue;
	    } else {
	      ret = delone(mesh,met,ip,list,lon);
	      if ( ret > 0 ) {
		ppt = &mesh->point[ip];
		met->m[ip] = 0.5 * (met->m[ip1] + met->m[ip2]);
		//chkmsh(mesh,0,0);
		addBucket(mesh,bucket,ip);
		(*ns)++;
		break;//imax continue;
	      }
	      else if ( ret == 0 ) {
		MMG_npd++;
		delPt(mesh,ip);
		goto collapse2;//continue;
	      }
	      else { /*allocation problem ==> savemesh*/
		return(0);
		/* MMG_npd++; */
		/* delPt(mesh,ip); */
		/* goto collapse2;//continue; */
	      }
	      printf("on passe pas la1\n");
	    }
	    printf("on passe pas la2\n");
	  }
	  printf("on passe pas la3\n");
	}
	printf("on passe pas la3\n");
      }
    collapse2:
      if(lmin > LOPTSDEL) continue;
      ifa0 = ifar[imin][0];
      ifa1 = ifar[imin][1];
      i  =  (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
      j  = iarfinv[i][imin];
      i1 = idir[i][inxt2[j]];
      i2 = idir[i][iprv2[j]];
      ip = pt->v[i1];
      iq = pt->v[i2];
      p0 = &mesh->point[ip];
      p1 = &mesh->point[iq];

      if ( (p0->tag > p1->tag) || (p0->tag & MG_REQ) )  continue;

      /* Case of a boundary face */
      ilist = 0;
      if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
	tag = pxt->tag[iarf[i][j]];
	if ( tag & MG_REQ )  continue;
	tag |= MG_BDY;
	if ( p0->tag > tag )   continue;
	if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) continue;
	ilist = chkcol_bdy(mesh,k,i,j,list);
	if ( ilist > 0 ) {
	  ier = colver(mesh,list,ilist,i2);
	  //nc += ier;
	  if ( ier < 0 ) return(-1);
	  else if(ier) {
	    //delBucket(mesh,bucket,ier);
	    delPt(mesh,ier);
	    (*nc)++;
	    break;//imax continue;
	  }
	}
	else if (ilist < 0 )  return(-1);
      }
      /* Case of an internal face */
      else {
	if ( p0->tag & MG_BDY )  continue;
	ilist = chkcol_int(mesh,met,k,i,j,list,2);
	if ( ilist > 0 ) {
	  ier = colver(mesh,list,ilist,i2);
	  if ( ilist < 0 ) continue;
	  //nc += ier;
	  if ( ier < 0 ) return(-1);
	  else if(ier) {
	    delBucket(mesh,bucket,ier);
	    delPt(mesh,ier);
	    (*nc)++;
	    break;//imax continue;
	  }
	}
	else if (ilist < 0 )  return(-1);
      }

    }//end for ii

    // } //end else de if(it>=3)
  }

  return(1);
}


/** Split edges of length bigger than LOPTL */
int adpsplcol(pMesh mesh,pSol met,pBucket bucket, int* warn) {
  pTetra     pt;
  pxTetra    pxt;
  Tria       ptt;
  pPoint     p0,p1,ppt;
  pxPoint    pxp;
  double     dd,len,lmax,o[3],to[3],ro[3],no1[3],no2[3],v[3];
  int        k,ip,ip1,ip2,list[LMAX+2],ilist,ns,ref;
  char       imax,tag,j,i,i1,i2,ifa0,ifa1;
  int        ifilt,lon,ret,ne,ier;
  double     lmin;
  int        imin,iq,nc,it,nnc,nns,nnf,nnm,maxit,nf,nm;
  int        ii,MMG_npd;
  double     maxgap;
  int nconf,imaxtet,imintet;
  double critloc;

  /* Iterative mesh modifications */
  it = nnc = nns = nnf = nnm = 0;
  maxit = 10;
  mesh->gap = maxgap = 0.5;
  MMG_npuiss=MMG_nvol=MMG_npres =MMG_npd=0 ;
  do {
    if ( !mesh->info.noinsert ) {
      *warn=0;
      ns = nc = 0;
      nf = nm = 0;
      ifilt = 0;
      ne = mesh->ne;
      ier = boucle_for(mesh,met,bucket,ne,&ifilt,&ns,&nc,warn,it);
      if(ier<0) exit(EXIT_FAILURE);
      else if(!ier) return(-1);
    } /* End conditional loop on mesh->info.noinsert */
    else  ns = nc = ifilt = 0;

    /* prilen(mesh,met); */
    /*     fprintf(stdout,"    REJECTED : %5d\n",MMG_npd); */
    /* fprintf(stdout,"          VOL      : %6.2f %%    %5d \n", */
    /*         100*(MMG_nvol/(float) */
    /*              MMG_npd),MMG_nvol); */
    /* fprintf(stdout,"          PUISS    : %6.2f %%    %5d \n", */
    /*         100*(MMG_npuiss/(float) MMG_npd),MMG_npuiss); */
    /* fprintf(stdout,"         PROCHE    : %6.2f %%    %5d \n", */
    /*         100*(MMG_npres/(float) MMG_npuiss),MMG_npres); */
    /* MMG_npd=0; */
    /* MMG_npuiss=0; */
    /* MMG_nvol=0; */
    /* MMG_npres=0; */
    if ( !mesh->info.noswap ) {
      nf = swpmsh(mesh,met
#ifndef PATTERN
                  ,bucket
#endif
                  );
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;
      if(it==2 || it==6/*&& it==1 || it==3 || it==5 || it > 8*/) {
        nf += swptetdel(mesh,met,1.053,bucket);
      } else {
        nf += 0;
      }
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;
    nnf+=nf;

    //fprintf(stdout,"$$$$$$$$$$$$$$$$$$ ITER SWAP %7d\n",nnf);

    if ( !mesh->info.nomove ) {
      nm = movtetdel(mesh,met,-1);
      if ( nm < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh.\n");
        return(0);
      }
    }
    else  nm = 0;
    nnm += nm;
    nnc += nc;
    nns += ns;

    /* decrease size of gap for reallocation */

    if ( mesh->gap > maxgap/(double)maxit )
      mesh->gap -= maxgap/(double)maxit;
    else
      mesh->gap -= mesh->gap/(double)maxit;


    if ( 1 || ((abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && ns+nc > 0) )
      fprintf(stdout,"     %8d filtered %8d splitted, %8d collapsed, %8d swapped, %8d moved\n",ifilt,ns,nc,nf,nm);
    if ( ns < 10 && abs(nc-ns) < 3 )  break;
    else if ( it > 3 && abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
/* #ifdef USE_SCOTCH */
/*     /\*check enough vertex to renum*\/ */
/*     if ( mesh->info.renum && (mesh->np/2. > BOXSIZE) ) { */
/*       /\* renumbering begin *\/ */
/*       if ( mesh->info.imprim > 5 ) */
/* 	fprintf(stdout,"renumbering"); */
/*       renumbering(BOXSIZE,mesh, met); */

/*       if ( mesh->info.imprim > 5) { */
/* 	fprintf(stdout,"  -- PHASE RENUMBERING COMPLETED. \n"); */
/*       } */
/*       if ( mesh->info.ddebug )  chkmsh(mesh,1,0); */
/*       /\* renumbering end *\/ */
/*     } */
/* #else */
/*     // free(mesh->adja); */
/*     //mesh->adja=NULL; */
/*     //hashTetra(mesh,1); */
/* #endif */
/*     /\*free bucket*\/ */
/*     DEL_MEM(mesh,bucket->head,(bucket->size*bucket->size*bucket->size+1)*sizeof(int)); */
/*     DEL_MEM(mesh,bucket->link,(mesh->npmax+1)*sizeof(int)); */
/*     DEL_MEM(mesh,bucket,sizeof(Bucket)); */
    
/*     bucket = newBucket(mesh,mesh->info.bucket);//M_MAX(mesh->mesh->info.bucksiz,BUCKSIZ)); */
/*     if ( !bucket )  return(0); */


  }
  while( ++it < maxit && nc+ns > 0 );

  return(1);
}

int optet(pMesh mesh, pSol met,pBucket bucket) {
  int it,nnm,nnf,maxit,nm,nf;
  double declic;

  /*shape optim*/
  it = nnm = nnf = 0;
  maxit = 4;
  declic = 1.053;
  do {
    /* badly shaped process */
    if ( !mesh->info.noswap ) {
      nf = swpmsh(mesh,met
#ifndef PATTERN
                  ,bucket
#endif
                  );
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = swptetdel(mesh,met,declic,bucket);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      nm = movtetdel(mesh,met,0);
      if ( nm < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh.\n");
        return(0);
      }
    }
    else  nm = 0;
    nnm += nm;

    /* ier = badelt(mesh,met);
    if ( ier < 0 ) {
      fprintf(stdout,"  ## Unable to remove bad elements.\n");
      return(0);
      }*/





    if ( (abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && nf+nm > 0 ){
      fprintf(stdout,"                                            ");
      fprintf(stdout,"%8d swapped, %8d moved\n",nf,nm);
    }
  }
  while( ++it < maxit && nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    nm = movtetdel(mesh,met,3);
    if ( nm < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh.\n");
      return(0);
    }
  }
  else  nm = 0;
  nnm += nm;
  if ( (abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && nm > 0 )
    fprintf(stdout,"                                            ");
  fprintf(stdout,"                  %8d moved\n",nm);

  return(1);
}

/** Analyze tetrahedra and split long / collapse short, according to prescribed metric */
/*static*/ int adptet1(pMesh mesh,pSol met,pBucket bucket) {
  int      nnf,maxit,ns,nf;
  int      warn;

  /*initial swap*/
  if ( !mesh->info.noswap ) {
    nf = swpmsh(mesh,met
#ifndef PATTERN
                ,bucket
#endif
                );
    if ( nf < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
    nnf = nf;
    nf = swptetdel(mesh,met,1.053,bucket);
    if ( nf < 0 ) {
      fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
    nnf+=nf;
  } else  nnf = nf = 0;
  fprintf(stdout,"$$$$$$$$$$$$$$$$$$ INITIAL SWAP %7d\n",nnf);
  outqua(mesh,met);

  /* Iterative mesh modifications */
  warn = 0;
#ifdef DEBUG
  tabtmp[0][0]=0;
  tabtmp[0][1]=0;     tabtmp[0][2]=0;
#endif
  ns = adpsplcol(mesh,met,bucket,&warn);
#ifdef DEBUG
  if ( ns ) { printf("APS ADPSPLCOL == %d\n",ns);
    prilen(mesh,met);
    printf(" histo %5.1f  %5.1f %5.1f\n", tabtmp[0][0],tabtmp[0][1],tabtmp[0][2]);}
#endif
  if ( ns < 0 ) {
    fprintf(stdout,"  ## Unable to complete mesh. Exit program.\n");
    return(0);
  }

  if ( warn ) {
    fprintf(stdout,"  ## Error:");
    fprintf(stdout," unable to allocate a new point in last call of adpspl.\n");
    fprintf(stdout,"  ## Check the mesh size or ");
    fprintf(stdout,"increase the maximal authorized memory with the -m option.\n");
    fprintf(stdout,"  ## Uncomplete mesh. Exiting\n" );
    return(0);
  }

#ifdef USE_SCOTCH
  /*check enough vertex to renum*/
  if ( mesh->info.renum && (mesh->np/2. > BOXSIZE) ) {
    /* renumbering begin */
    if ( mesh->info.imprim > 5 )
      fprintf(stdout,"renumbering");
    renumbering(BOXSIZE,mesh, met);

    if ( mesh->info.imprim > 5) {
      fprintf(stdout,"  -- PHASE RENUMBERING COMPLETED. \n");
    }
    if ( mesh->info.ddebug )  chkmsh(mesh,1,0);
    /* renumbering end */
  }
#endif
  outqua(mesh,met);

  if(!optet(mesh,met,bucket)) return(0);

  return(1);
}




/** analyze tetrahedra and split if needed */
/*static*/ int anatetdel(pMesh mesh,pSol met,char typchk) {
  int     ier,nc,ns,nf,nnc,nns,nnf,it,maxit;

  /* analyze tetras : initial splitting */
  nns = nnc = nnf = it = 0;
  maxit = 5;
  mesh->gap = 0.5;
  do {
    /* memory free */
    DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));
#ifdef DEBUG
    puts("AVT ANATET4");
    prilen(mesh,met);
#endif
    if ( !mesh->info.noinsert ) {
      /* split tetra with more than 2 bdry faces */
      ier = anatet4(mesh,met);
#ifdef DEBUG
      if ( ier ) { printf("APS ANATET4 == %d\n",ier);
        prilen(mesh,met);}
#endif
      if ( ier < 0 )  return(0);
      ns = ier;

      /* analyze surface tetras */
      ier = anatets(mesh,met,typchk);
#ifdef DEBUG
      if ( ier ) { printf("APS ANATETS == %d\n",ier);
        prilen(mesh,met);}
#endif
      if ( ier < 0 ) {
        fprintf(stdout,"  ## Unable to complete surface mesh. Exit program.\n");
        return(0);
      }
      ns += ier;

      /* analyze internal tetras */
      ier = 0;
#ifdef DEBUG
      if(ier){ printf("APS ANATETV == %d\n",ier);
        prilen(mesh,met);}
#endif
      if ( ier < 0 ) {
        fprintf(stdout,"  ## Unable to complete volume mesh. Exit program.\n");
        return(0);
      }
      ns += ier;
    }
    else  ns = 0;

    if ( !hashTetra(mesh,1) ) {
      fprintf(stdout,"  ## Hashing problem. Exit program.\n");
      return(0);
    }
    if ( typchk == 2 && it == maxit-1 )  mesh->info.fem = 1;

    /* collapse short edges */
    if ( !mesh->info.noinsert ) {
      nc = coltet(mesh,met,typchk);
      if ( nc < 0 ) {
        fprintf(stdout,"  ## Unable to collapse mesh. Exiting.\n");
        return(0);
      }
    }
    else  nc = 0;

    /* attempt to swap */
    if ( !mesh->info.noswap ) {
      nf = swpmsh(mesh,met
#ifndef PATTERN
                  ,NULL
#endif
                  );
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = swptetdel(mesh,met,1.1,NULL);
      if ( nf < 0 ) {
        fprintf(stdout,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    nnc += nc;
    nns += ns;
    nnf += nf;
    if ( (abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && ns+nc+nf > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped\n",ns,nc,nf);
    if ( it > 3 && abs(nc-ns) < 0.1 * MG_MAX(nc,ns) )  break;
  }
  while ( ++it < maxit && ns+nc+nf > 0 );

  if ( (abs(mesh->info.imprim) < 4 || mesh->info.ddebug ) && nns+nnc > 0 )
    fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %d iter.\n",nns,nnc,nnf,it);
#ifdef DEBUG
  puts("FIN ANATET");
  prilen(mesh,met);
#endif
  return(1);
}

/** main adaptation routine */
int mmg3d1_delone(pMesh mesh,pSol met) {
  pBucket bucket;

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  if ( mesh->info.iso && !chkmani(mesh) ) {
    fprintf(stdout,"  ## Non orientable implicit surface. Exit program.\n");
    return(0);
  }

  /**--- stage 1: geometric mesh */
  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !anatetdel(mesh,met,1) ) {
    fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }
  outqua(mesh,met);

#ifdef DEBUG
  outqua(mesh,met);
#endif

  /**--- stage 2: computational mesh */
  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");

  /* define metric map */
  if ( !defsiz(mesh,met) ) {
    fprintf(stdout,"  ## Metric undefined. Exit program.\n");
    return(0);
  }

  if ( mesh->info.hgrad > 0. && !gradsiz(mesh,met) ) {
    fprintf(stdout,"  ## Gradation problem. Exit program.\n");
    return(0);
  }
  if ( !anatetdel(mesh,met,2) ) {
    fprintf(stdout,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }
  outqua(mesh,met);

#ifdef DEBUG
  puts("---------------------------Fin anatet---------------------");
  outqua(mesh,met);
#endif
#ifdef USE_SCOTCH
    /*check enough vertex to renum*/
    if ( mesh->info.renum  && (mesh->np/2. > BOXSIZE) ) {
      /* renumbering begin */
      if ( mesh->info.imprim > 5 )
        fprintf(stdout,"renumbering");
      renumbering(BOXSIZE,mesh, met);

      if ( mesh->info.imprim > 5) {
        fprintf(stdout,"  -- PHASE RENUMBERING COMPLETED. \n");
      }

      if ( mesh->info.ddebug )  chkmsh(mesh,1,0);
      /* renumbering end */
    }
#endif

  /* CEC : create filter */
  bucket = newBucket(mesh,mesh->info.bucket);//M_MAX(mesh->mesh->info.bucksiz,BUCKSIZ));
  if ( !bucket )  return(0);

  if ( !adptet1(mesh,met,bucket) ) {
    fprintf(stdout,"  ## Unable to adapt. Exit program.\n");
    return(0);
  }
  outqua(mesh,met);

#ifdef DEBUG
  puts("---------------------Fin adptet-----------------");
  outqua(mesh,met);
#endif
  /* in test phase: check if no element with 2 bdry faces */
  if ( !chkfemtopo(mesh) ) {
    fprintf(stdout,"  ## Topology of mesh unsuited for fem computations. Exit program.\n");
    return(0);
  }

  if ( mesh->info.iso && !chkmani(mesh) ) {
    fprintf(stdout,"  ## Non orientable implicit surface. Exit program.\n");
    return(0);
  }

  /*free bucket*/
  DEL_MEM(mesh,bucket->head,(bucket->size*bucket->size*bucket->size+1)*sizeof(int));
  DEL_MEM(mesh,bucket->link,(mesh->npmax+1)*sizeof(int));
  DEL_MEM(mesh,bucket,sizeof(Bucket));

  return(1);
}
