#include "mmg3d.h"

extern char  ddb;

/** solve 3*3 non symmetric system Ar = b */
static inline int invsl(double A[3][3],double b[3],double r[3]) {
  double detA;

  detA = A[0][0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2]) \
    - A[0][1]*(A[1][0]*A[2][2] - A[2][0]*A[1][2]) \
    + A[0][2]*(A[1][0]*A[2][1] - A[2][0]*A[1][1]);
  if ( detA < EPSD )  return(0);
  detA = 1.0 / detA;

  r[0] =  b[0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2]) \
    - A[0][1]*(b[1]*A[2][2] - b[2]*A[1][2]) \
    + A[0][2]*(b[1]*A[2][1] - b[2]*A[1][1]);

  r[1] = A[0][0]*(b[1]*A[2][2] - b[2]*A[1][2]) \
    - b[0]*(A[1][0]*A[2][2] - A[2][0]*A[1][2]) \
    + A[0][2]*(A[1][0]*b[2] - A[2][0]*b[1]);

  r[2] = A[0][0]*(A[1][1]*b[2] - A[2][1]*b[1]) \
    - A[0][1]*(A[1][0]*b[2] - A[2][0]*b[1]) \
    + b[0]*(A[1][0]*A[2][1] - A[2][0]*A[1][1]);

  r[0] *= detA;
  r[1] *= detA;
  r[2] *= detA;

  return(1);
}

/** Check ball of point np, return 0 if a nonmanifold situation is created */
static int ismaniball(pMesh mesh,pSol sol,int k,int indp) {
  pTetra   pt,pt1;
  double   v,v0,v1,v2;
  int      *adja,list[LMAX+1],bdy[LMAX+1],ibdy,np,ilist,base,cur,iel,jel,res,l;
  char     i,i0,i1,i2,j0,j1,j2,j,ip,nzeros,nopp,nsame;

  pt = &mesh->tetra[k];
  np = pt->v[indp];
  if ( fabs(sol->m[np]) > EPSD2 )  return(1);

  memset(bdy,0,(LMAX+1)*sizeof(int));
#ifdef SINGUL
 restart:
#endif
  memset(list,0,(LMAX+1)*sizeof(int));

  /* Sign of a starting point in ball of np */
  for (j=0; j<3; j++) {
    ip = idir[indp][j];
    if ( sol->m[pt->v[ip]] != 0.0 )  break;
  }
  if ( j == 3) {
    fprintf(stdout,"  *** Problem in function ismaniball : tetra %d : 4 null values",k);
    exit(EXIT_FAILURE);
  }

  v = sol->m[pt->v[ip]];
  base = ++mesh->base;
  pt->flag = base;
  ilist = 0;
  list[ilist] = 4*k+indp;
  ilist++;

  /* travel list and pile up, by adjacency, faces of ball of np while they have at least
     a vertex with same sign as v */
  res = cur = 0;
  while ( cur < ilist ) {
    iel = list[cur] / 4;
    i   = list[cur] % 4;
    pt  = &mesh->tetra[iel];
    adja = &mesh->adja[4*(iel-1)+1];

    /* Store a face for starting back enumeration with the opposite sign */
    if ( !res ) {
      for (j=0; j<3; j++) {
        i1 = idir[i][j];
        v1 = sol->m[pt->v[i1]];
        if ( ( v1 != 0.0 ) && !MG_SMSGN(v,v1) ) {
          res = 4*iel + i;
          break;
        }
      }
    }

    /* Pile up faces sharing a vertex with same sign as v */
    for (j=0; j<3; j++) {
      i1 = idir[i][inxt2[j]];
      i2 = idir[i][iprv2[j]];
      v1 = sol->m[pt->v[i1]];
      v2 = sol->m[pt->v[i2]];

      if ( ( ( v1 != 0.0 ) && MG_SMSGN(v,v1) ) ||
           ( ( v2 != 0.0 ) && MG_SMSGN(v,v2) ) ) {
        jel = adja[idir[i][j]] / 4;
        if( !jel ) continue;
        pt1 = &mesh->tetra[jel];

        if ( pt1->flag == base )  continue;
        for (ip=0; ip<4; ip++) {
          if ( pt1->v[ip] == np )  break;
        }
        assert( ip < 4 );
        pt1->flag   = base;
        list[ilist] = 4*jel + ip;
        ilist++;
        assert(ilist < LMAX);
      }
    }
    cur++;
  }
  /* 0 value has been snapped accidentally */
  if ( !res ) {
#ifdef SINGUL
    fprintf(stdout,"%s:%d: Warning:\n",__FILE__,__LINE__);
    fprintf(stdout,"Point with value 0 arounded by points of");
    fprintf(stdout," same sign: elt %d (%d), pt %d (%d)\n",
	    k,indElt(mesh,k),np,indPt(mesh,np));
    /* Try to find a tetrahedra to swap */
    for(l=0; l<ilist; l++) {
      iel  = list[l] / 4;
      i    = list[l] % 4;
      pt   = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];
      jel  = adja[i]/4;
      if ( !jel )  continue;
      j    = adja[i]%4;
      pt1  = &mesh->tetra[jel];
      v1   = sol->m[pt1->v[j]];
      if ( ( v1 != 0.0 ) && !MG_SMSGN(v,v1) ) {
        if ( swap23(mesh,iel,i) ) {
          if ( l==0 ) {
            for (indp=0; indp<4; indp++ ) {
              if ( pt->v[indp] == np ) break;
            }
            assert(indp<4);
          }
          goto restart;
        }
      }
    }
#endif
    return(0);
  }


  /* Fill in list bdy, corresponding to the support tetras of the boundary to be created */
  ibdy = 0;
  for(l=0; l<ilist; l++) {
    iel = list[l] / 4;
    i   = list[l] % 4;
    pt  = &mesh->tetra[iel];

    nzeros = nsame = nopp = 0;

    i0 = idir[i][0];
    i1 = idir[i][1];
    i2 = idir[i][2];

    v0 = sol->m[pt->v[i0]];
    v1 = sol->m[pt->v[i1]];
    v2 = sol->m[pt->v[i2]];

    if ( v0 == 0.0 )
      nzeros++;
    else if ( MG_SMSGN(v,v0) )
      nsame++;
    else
      nopp++;

    if ( v1 == 0.0 )
      nzeros++;
    else if ( MG_SMSGN(v,v1) )
      nsame++;
    else
      nopp++;

    if ( v2 == 0.0 )
      nzeros++;
    else if ( MG_SMSGN(v,v2) )
      nsame++;
    else
      nopp++;

    if ( ( nzeros == 2 && nsame == 1 ) || ( nsame >= 1 && nopp >= 1 ) )  {
      bdy[ibdy] = list[l];
      ibdy++;
    }
  }

  /* Reset the current part of the ball, and start back the process with the other sign */
  iel = res / 4;
  pt = &mesh->tetra[iel];
  base = ++mesh->base;
  pt->flag = base;

  memset(list,0,(LMAX+1)*sizeof(int));
  ilist = cur = 0;
  list[ilist] = res;
  ilist++;
  while ( cur < ilist ) {
    iel = list[cur] / 4;
    i   = list[cur] % 4;
    pt  = &mesh->tetra[iel];
    adja = &mesh->adja[4*(iel-1)+1];

    /* Pile up faces sharing a vertex with opposite sign to v */
    for (j=0; j<3; j++) {
      i1 = idir[i][inxt2[j]];
      i2 = idir[i][iprv2[j]];
      v1 = sol->m[pt->v[i1]];
      v2 = sol->m[pt->v[i2]];

      if ( v1 == 0.0 && v2 == 0.0 ) {
        jel = adja[idir[i][j]] / 4;
        if( !jel ) continue;
        pt1 = &mesh->tetra[jel];
        pt1->flag = base;
      }

      else if ( ( ( v1 != 0.0 ) && (!MG_SMSGN(v,v1)) ) || ( ( v2 != 0.0 ) && (!MG_SMSGN(v,v2)) ) ) {
        jel = adja[idir[i][j]] / 4;
        if( !jel ) continue;
        pt1 = &mesh->tetra[jel];

        j0 = idir[i][0];
        j1 = idir[i][1];
        j2 = idir[i][2];

        v0 = sol->m[pt1->v[j0]];
        v1 = sol->m[pt1->v[j1]];
        v2 = sol->m[pt1->v[j2]];

        nzeros = nsame = nopp = 0;

        if ( v0 == 0.0 )
          nzeros++;
        else if ( MG_SMSGN(v,v0) )
          nsame++;
        else
          nopp++;

        if ( v1 == 0.0 )
          nzeros++;
        else if ( MG_SMSGN(v,v1) )
          nsame++;
        else
          nopp++;

        if ( v2 == 0.0 )
          nzeros++;
        else if ( MG_SMSGN(v,v2) )
          nsame++;
        else
          nopp++;

        if ( ( nzeros == 2 && nsame == 1 ) || ( nsame >= 1 && nopp >= 1 ) )  {
          if ( pt1->flag < base - 1 ) return(0);
        }

        if ( pt1->flag == base ) continue;
        for (ip=0; ip<4; ip++) {
          if ( pt1->v[ip] == np )  break;
        }
        assert( ip < 4 );
        pt1->flag   = base;
        list[ilist] = 4*jel + ip;
        ilist++;
        assert(ilist < LMAX);
      }
    }
    cur++;
  }

  /* Now, all elements of bdy should have been marked by a flag base + 1 */
  for (l=0; l<ibdy; l++) {
    iel = bdy[l] / 4;
    pt = &mesh->tetra[iel];
    if ( pt->flag != base ) return(0);
  }

  return(1);
}

/** Snap values of the level set function very close to 0 to exactly 0,
    and prevent nonmanifold patterns from being generated */
static int snpval_ls(pMesh mesh,pSol sol,double *tmp) {
  pTetra   pt;
  pPoint   p0;
  int      k,nc,ns,ip;
  char     i;

  /* create tetra adjacency */
  if ( !hashTetra(mesh,1) ) {
    fprintf(stdout,"  ## Hashing problem (1). Exit program.\n");
    return(0);
  }

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Snap values of sol that are close to 0 to 0 exactly */
  ns = nc = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) ) continue;
    if ( fabs(sol->m[k]) < EPS ) {
      if ( mesh->info.ddebug )  fprintf(stdout,"  Snapping value %d ; previous value : %E\n",k,fabs(sol->m[k]));
      tmp[k] = ( fabs(sol->m[k]) < EPSD ) ? (-100.0*EPS) : sol->m[k];
      p0->flag = 1;
      sol->m[k] = 0.0;
      ns++;
    }
  }

  /* Check snapping did not lead to a nonmanifold situation */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    for (i=0; i<4; i++) {
      ip = pt->v[i];
      p0 = &mesh->point[ip];
      if ( p0->flag ) {
        if ( !ismaniball(mesh,sol,k,i) ) {
          sol->m[ip] = tmp[ip];
          nc++;
        }
        p0->flag = 0;
        tmp[ip]  = 0.0;
      }
    }
  }
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && ns+nc > 0 )
    fprintf(stdout,"     %8d points snapped, %d corrected\n",ns,nc);

  /* memory free */
  DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));
  return(1);
}

/** Proceed to discretization of the implicit function carried by sol into mesh, once values
    of sol have been snapped/checked */
static int cuttet_ls(pMesh mesh, pSol sol/*,double *tmp*/){
  pTetra   pt;
  pPoint   p0,p1;
  Hash     hash;
  double   c[3],v0,v1,s;
  int      vx[6],nb,k,ip0,ip1,np,ns,ne;
  char     ia;
  /* Commented because unused */
  /*pPoint  p[4];*/
  /*double   *grad,A[3][3],b[3],*g0,*g1,area,a,d,dd,s1,s2;*/
  /*int       ip[4],ng*/
  /*char    i,ier;*/

  /* reset point flags and h */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* compute the number nb of intersection points on edges */
  nb = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[iare[ia][0]];
      ip1 = pt->v[iare[ia][1]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];
      if ( p0->flag && p1->flag )  continue;
      v0  = sol->m[ip0];
      v1  = sol->m[ip1];
      if ( fabs(v0) > EPSD2 && fabs(v1) > EPSD2 && v0*v1 < 0.0 ) {
        if ( !p0->flag ) {
          p0->flag = nb;
          nb++;
        }
        if ( !p1->flag ) {
          p1->flag = nb;
          nb++;
        }
      }
    }
  }
  if ( ! nb )  return(1);

  /* Store gradients of level set function at those points */
  /* Commented because unused */
  /* grad = (double*)calloc(3*nb+1,sizeof(double)); */
  /* assert(grad); */

  /* for (k=1; k<=mesh->ne; k++) { */
  /*   pt = &mesh->tetra[k]; */
  /*   ia = 0; */
  /*   for (i=0; i<4; i++) { */
  /*     ip[i] = pt->v[i]; */
  /*     p[i]  = &mesh->point[ip[i]]; */
  /*     if ( p[i]->flag == 0 )  ia++; */
  /*   } */
  /*   if ( ia == 4 )  continue; */

  /*   A[0][0] = p[1]->c[0] - p[0]->c[0];  A[0][1] = p[1]->c[1] - p[0]->c[1];  A[0][2] = p[1]->c[2] - p[0]->c[2]; */
  /*   A[1][0] = p[2]->c[0] - p[0]->c[0];  A[1][1] = p[2]->c[1] - p[0]->c[1];  A[1][2] = p[2]->c[2] - p[0]->c[2]; */
  /*   A[2][0] = p[3]->c[0] - p[0]->c[0];  A[2][1] = p[3]->c[1] - p[0]->c[1];  A[2][2] = p[3]->c[2] - p[0]->c[2]; */

  /*   b[0] = sol->m[ip[1]] - sol->m[ip[0]]; */
  /*   b[1] = sol->m[ip[2]] - sol->m[ip[0]]; */
  /*   b[2] = sol->m[ip[3]] - sol->m[ip[0]]; */

  /*   area = det4pt(p[0]->c,p[1]->c,p[2]->c,p[3]->c); */
  /*   ier  = invsl(A,b,c); */
  /*   if ( !ier )  continue; */

  /*   for (i=0; i<4; i++) { */
  /*     if ( p[i]->flag ) { */
  /*       ng = p[i]->flag; */
  /*       tmp[ip[i]] += fabs(area); */
  /*       grad[3*(ng-1)+1] += (area*c[0]); */
  /*       grad[3*(ng-1)+2] += (area*c[1]); */
  /*       grad[3*(ng-1)+3] += (area*c[2]); */
  /*     } */
  /*   } */
  /* } */
  /* for (k=1; k<=mesh->np; k++) { */
  /*   p0 = &mesh->point[k]; */
  /*   if ( p0->flag ) { */
  /*     area = MG_MAX(EPSD2,tmp[k]); */
  /*     area = 1.0 / area; */
  /*     ng   = p0->flag; */
  /*     grad[3*(ng-1)+1] *= area; */
  /*     grad[3*(ng-1)+2] *= area; */
  /*     grad[3*(ng-1)+3] *= area; */
  /*   } */
  /* } */

  /* Create intersection points at 0 isovalue and set flags to tetras */
  if ( !hashNew(mesh,&hash,nb,7*nb) ) return(0);
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[iare[ia][0]];
      ip1 = pt->v[iare[ia][1]];
      np  = hashGet(&hash,ip0,ip1);
      if ( np )  continue;

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      v0 = sol->m[ip0];
      v1 = sol->m[ip1];
      if ( fabs(v0) < EPSD2 || fabs(v1) < EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;

      /* g0 = &grad[3*(p0->flag -1)+1]; */
      /* g1 = &grad[3*(p1->flag -1)+1]; */
      /* a = 0.5 * ((g1[0]-g0[0])*(p1->c[0]-p0->c[0]) + (g1[1]-g0[1])*(p1->c[1]-p0->c[1]) \ */
      /*            + (g1[2]-g0[2])*(p1->c[2]-p0->c[2])); */
      /* d  = v1 - v0 - a; */
      /* dd = d*d - 4.0*a*v0; */
      /* dd = MG_MAX(EPSD2,dd); */
      /* dd = sqrt(dd); */
      /* if ( fabs(a) < EPSD2 ) */
      /*   s = v0 / (v0-v1); */
      /* else { */
      /*   s1 = 0.5*( dd -d) / a; */
      /*   s2 = 0.5*(-dd -d) / a; */
      /*   if ( s1 > 0.0 && s1 < 1.0 ) */
      /*     s = s1; */
      /*   else if (s2 > 0.0 && s2 < 1.0) */
      /*     s = s2; */
      /*   else */
      /*     s = MG_MIN(fabs(s1),fabs(s1-1.0)) < MG_MIN(fabs(s2),fabs(s2-1.0)) ? s1 : s2 ; */
      /* } */
      // IMPORTANT A REGARDER
      s = v0 / (v0-v1);

      s = MG_MAX(MG_MIN(s,1.0-EPS),EPS);
      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
      c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

      np = newPt(mesh,c,0);
      if ( !np ) {
        POINT_REALLOC(mesh,sol,np,0.2,
                      printf("  ## Error: unable to allocate a new point\n");
                      printf("  ## Check the mesh size or increase");
                      printf(" the allocated memory with the -m option.\n");
                      return(0)
                      ,c,0);
      }
      sol->m[np] = 0.0;
      hashEdge(mesh,&hash,ip0,ip1,np);
    }
  }

  /* Proceed to splitting, according to flags to tets */
  ne = mesh->ne;
  ns = 0;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    pt->flag = 0;
    memset(vx,0,6*sizeof(int));
    for (ia=0; ia<6; ia++) {
      vx[ia] = hashGet(&hash,pt->v[iare[ia][0]],pt->v[iare[ia][1]]);
      if ( vx[ia] )  MG_SET(pt->flag,ia);
    }
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
      split1(mesh,sol,k,vx);
      ns++;
      break;

    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      split2sf(mesh,sol,k,vx);
      ns++;
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      split3cone(mesh,sol,k,vx);
      ns++;
      break;

    case 30: case 45: case 51:
      split4op(mesh,sol,k,vx);
      ns++;
      break;

    default :
      assert(pt->flag == 0);
      break;
    }
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",ns);

  DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
  return(ns);
}

#ifdef SINGUL
/** Check if all singular edges will appear after split on the level-set function *
 *  (the shell of singular edges must have vertices with opposite sign). */
static inline
int chkedg_ls(pMesh mesh, pSol sol){
  pTetra   pt,pt1;
  pxTetra  pxt;
  double   v0,v1,v;
  int      k,l,iel,i,ia,piv,ipa,ipb,na,nb,ilist;
  int      *adja,adj,list[LMAX+2],nf,nc;

  nf = nc = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( (!MG_EOK(pt)) || (!pt->xt) ) continue;
    pxt = &mesh->xtetra[pt->xt];

    for (ia=0; ia<6; ia++) {
    nextEdg:
      if ( pxt->tag[ia] & MG_SGL ) {

        /* Sign of a starting point in coquil of i */
        v0 = sol->m[pt->v[ifar[ia][0]]];
        v1 = sol->m[pt->v[ifar[ia][1]]];
        if ( v0 != 0.0 ) {
          if ( v1 != 0.0 && !MG_SMSGN(v0,v1) ) continue;
          v  = v0;
        }
        else if ( v1 != 0.0 ) {
          if ( v1 != 0.0 && !MG_SMSGN(v0,v1) ) continue;
          v = v1;
        } else {
          fprintf(stdout,"  *** Problem in function chkedg_ls :");
          fprintf(stdout," tetra %d : 4 null values",k);
          return(0);
        }
        na   = pt->v[ iare[ia][0] ];
        nb   = pt->v[ iare[ia][1] ];
        memset(list,0,(LMAX+2)*sizeof(int));
        ilist = 0;
        list[ilist] = 6*k+ia;
        ilist++;

        /* travel list and pile up, by adjacency, faces of coquil of edge ia
           while they have at least a vertex with same sign as v */
        adja = &mesh->adja[4*(k-1)+1];
        adj  = adja[ifar[ia][0]] / 4; // start travelling by face (ia,0)
        piv  = pt->v[ifar[ia][1]];

        while ( adj && (adj != k) ) {
          pt1 = &mesh->tetra[adj];
          /* identification of edge number in tetra adj */
          for (i=0; i<6; i++) {
            ipa = iare[i][0];
            ipb = iare[i][1];
            if ( (pt1->v[ipa] == na && pt1->v[ipb] == nb) ||
                 (pt1->v[ipa] == nb && pt1->v[ipb] == na))  break;
          }
          assert(i<6);
          v0 = sol->m[pt1->v[ifar[ia][0]]];
          v1 = sol->m[pt1->v[ifar[ia][1]]];
          if ( ( ( v0 != 0.0 ) && !MG_SMSGN(v,v0) ) ||
               ( ( v1 != 0.0 ) && !MG_SMSGN(v,v1) ) ) {
            ++ia;
            goto nextEdg;
          }

          list[ilist] = 6*adj +i;
          ilist++;
          /* overflow */
          assert( ilist <= LMAX-3 );

          /* set new triangle for travel */
          adja = &mesh->adja[4*(adj-1)+1];
          if ( pt1->v[ ifar[i][0] ] == piv ) {
            adj = adja[ ifar[i][0] ] / 4;
            piv = pt1->v[ ifar[i][1] ];
          }
          else {
            assert(pt1->v[ ifar[i][1] ] == piv );
            adj = adja[ ifar[i][1] ] /4;
            piv = pt1->v[ ifar[i][0] ];
          }
        }

        /* At this point, the first travel, in one direction, *
         * of the shell is complete. Now, analyze why the travel ended. */
        if ( adj != k ) {
          /* Theoretically, singular edges are far enough to domain boundaries so *
           * we hope to have only closed shell and this conditional loop *
           * is useless. */

          assert( !adj );
          adj = list[ilist-1] / 6;
          i   = list[ilist-1] % 6;
          ilist = 0;

          /* Start back everything from this tetra adj */
          list[ilist] = 6*adj + i;
          ilist++;
          /* overflow */
          assert( ilist <= LMAX-3 );

          adja = &mesh->adja[4*(adj-1)+1];
          if ( pt1->v[ ifar[i][0] ] == piv ) {
            adj = adja[ ifar[i][0] ] / 4;
            piv = pt1->v[ ifar[i][1] ];
          }
          else {
            adj = adja[ ifar[i][1] ] /4;
            piv = pt1->v[ ifar[i][0] ];
          }
          while ( adj ) {
            pt1 = &mesh->tetra[adj];
            if ( pt1->tag & MG_REQ )  return(0);
            /* identification of edge number in tetra adj */
            for (i=0; i<6; i++) {
              ipa = iare[i][0];
              ipb = iare[i][1];
              if ( (pt1->v[ipa] == na && pt1->v[ipb] == nb) ||
                   (pt1->v[ipa] == nb && pt1->v[ipb] == na))  break;
            }
            assert(i<6);
            v1 = sol->m[pt1->v[ifar[ia][1]]];
            if ( ( ( v0 != 0.0 ) && !MG_SMSGN(v,v0) ) ||
                 ( ( v1 != 0.0 ) && !MG_SMSGN(v,v1) ) ) {
              ++ia;
              goto nextEdg;
            }
            list[ilist] = 6*adj +i;
            ilist++;
            /* overflow */
            assert( ilist <= LMAX-2 );

            /* set new triangle for travel */
            adja = &mesh->adja[4*(adj-1)+1];
            if ( pt1->v[ ifar[i][0] ] == piv ) {
              adj = adja[ ifar[i][0] ] / 4;
              piv = pt1->v[ ifar[i][1] ];
            }
            else {
              adj = adja[ ifar[i][1] ] /4;
              piv = pt1->v[ ifar[i][0] ];
            }
          }
          assert(!adj);
        }

        /* shell is complete and we have'nt found a vertex with different sign
         * for the level-set function. */
        nf++;
        fprintf(stdout,"%s:%d: Warning:\n",__FILE__,__LINE__);
        fprintf(stdout,"Singular edge arounded by elements of");
        fprintf(stdout," same sign: elt %d (%d), edge %d--%d (%d--%d)\n",
                k, indElt(mesh,k),na,nb,indPt(mesh,na),indPt(mesh,nb));
        for(l=0; l<ilist; l++) {
          iel  = list[l] / 6;
          i    = list[l] % 6;
          adja = &mesh->adja[4*(iel-1)+1];
          adj  = adja[ iare[i][0] ]/4;
          if ( adj ) {
            piv  = adja[ iare[i][0] ]%4;
            pt1  = &mesh->tetra[adj];
            v1   = sol->m[pt1->v[piv]];
            if ( ( v1 != 0.0 ) && !MG_SMSGN(v,v1) ) {
              if ( swap23(mesh,iel,iare[i][0]) ) {
                nc++;
                break;
              }
            }
          }
          adj  = adja[ iare[i][1] ]/4;
          if ( adj ) {
            piv  = adja[ iare[i][1] ]%4;
            pt1  = &mesh->tetra[adj];
            v1   = sol->m[pt1->v[piv]];
            if ( ( v1 != 0.0 ) && !MG_SMSGN(v,v1) ) {
              if ( swap23(mesh,iel,iare[i][1]) ) {
                nc++;
                break;
              }
            }
          }
        }
      }
    }
  }
  if ( (abs(mesh->info.imprim) > 3 || mesh->info.ddebug) && nf > 0 )
    fprintf(stdout,"    %8d problematic edges, %d corrected\n",nf,nc);
  return(1);
}
#endif

/** Set references to tets according to the sign of the level set function */
static int setref_ls(pMesh mesh, pSol sol) {
  pTetra   pt;
  int      k,ip;
  char     nmns,npls,nz,i;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    nmns = npls = nz = 0;
    for (i=0; i<4; i++) {
      ip = pt->v[i];
      if ( sol->m[ip] > 0.0 )
        npls++;
      else if ( sol->m[ip] < 0.0 )
        nmns++;
      else
        nz ++;
    }
    assert(nz < 4);
    if ( npls ) {
      assert(!nmns);
      pt->ref = MG_PLUS;
    }
    else {
      assert(nmns);
      pt->ref = MG_MINUS;
    }
  }
  return(1);
}

/** Check whether implicit surface is orientable in ball of point ip in tet iel ;
    Beware : may return 0 when implicit boundary is tangent to outer boundary */
int chkmaniball(pMesh mesh, int start, char ip){
  pTetra    pt,pt1;
  int       ref,base,ilist,nump,k,cur,k1,nref;
  int       *adja,list[LMAX+2];
  char      i,l,j;

  base = ++mesh->base;
  ilist = 0;

  pt = &mesh->tetra[start];
  nump = pt->v[ip];
  ref = pt->ref;

  /* Store initial tetrahedron */
  pt->flag = base;
  list[ilist] = 4*start+ip;
  ilist++;

  /* explore list, and find all tets in ball of p belonging to the component ref */
  cur = 0;
  while( cur < ilist ) {
    k = list[cur] / 4;
    i = list[cur] % 4;
    pt = &mesh->tetra[k];

    adja = &mesh->adja[4*(k-1)+1];
    for(l=0; l<3; l++){
      i = inxt3[i];

      /* Travel only through non boundary faces. */
      k1 = adja[i] / 4;
      if(!k1) continue;
      pt1 = &mesh->tetra[k1];

      if( pt1 ->ref != ref ) continue;

      if( pt1->flag == base ) continue;
      pt1->flag = base;

      for(j=0; j<4 ; j++){
        if(pt1->v[j] == nump)
          break;
      }
      assert(j<4);

      /* overflow */
      assert ( ilist <= LMAX-3 );
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }

  /* Number of caught tets with ref ptstart->ref*/
  nref = ilist;

  /* Complete ball of point */
  cur = 0;
  while(cur < ilist){
    k = list[cur] / 4;
    i = list[cur] % 4;
    pt = &mesh->tetra[k];

    adja = &mesh->adja[4*(k-1)+1];
    for(l=0; l<3; l++){
      i = inxt3[i];

      k1 = adja[i]/4;
      if(!k1) continue;
      pt1 = &mesh->tetra[k1];
      if(pt1->flag == base) continue;
      pt1->flag = base;

      for(j=0; j<4 ; j++){
        if(pt1->v[j] == nump)
          break;
      }
      assert(j<4);

      /* overflow */
      assert ( ilist <= LMAX-3 );
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }

  /* Elements from nref to ilist-1 must not have ref ptstart->ref */
  for(cur=nref; cur<ilist; cur++) {
    k = list[cur] / 4;
    pt = &mesh->tetra[k];
    if( pt->ref == ref ) {
      fprintf(stdout,"   *** Topological problem:");
      fprintf(stdout," non manifold surface at point %d \n",nump);
      return(0);
    }
  }

  return(1);
}

/** Check whether implicit surface enclosed in volume is orientable */
int chkmani(pMesh mesh){
  pTetra    pt,pt1;
  int       k,iel,ref;
  int       *adja;
  char      i,j,ip,cnt;

  for(k=1; k<=mesh->np; k++){
    mesh->point[k].flag = 0;
  }

  /** First test : check whether a tetra has 4 boundary faces */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )   continue;
    adja = &mesh->adja[4*(k-1)+1];

    ref = pt->ref;
    cnt = 0;
    for(i=0; i<4; i++) {
      if( !adja[i] ) {
        cnt++;
      }
      else {
        pt1 = &mesh->tetra[adja[i]/4];
        if ( pt1->ref != ref ) cnt++;
      }
    }
    if ( cnt == 4 ) {
      fprintf(stdout,"Tetra %d : 4 boundary faces \n",k);
      //return(0);
    }
  }

  /** Second test : Check whether configuration is manifold in each ball */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ))   continue;
    adja = &mesh->adja[4*(k-1)+1];

    for(i=0; i<4; i++){
      if(!adja[i]) continue;
      iel = adja[i] / 4;
      pt1 = &mesh->tetra[iel];
      if(pt1->ref == pt->ref) continue;

      for(j=0; j<3; j++){
        ip = idir[i][j];

        if(!chkmaniball(mesh,k,ip))
          return(0);
      }
    }
  }

  fprintf(stdout,"  *** Manifold implicit surface.\n");
  return(1);
}

/** Check whether implicit surface enclosed in volume is orientable (perform an additionnal
    test w.r.t. chkmani) */
int chkmani2(pMesh mesh,pSol sol) {
  pTetra    pt,pt1;
  int       k,iel;
  int       *adja;
  char      i,j,ip,cnt;

  for(k=1; k<=mesh->np; k++){
    mesh->point[k].flag = 0;
  }

  /** First test : assure no tetra has its 4 vertices on implicit boundary */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ))   continue;

    cnt = 0;
    for(j=0; j<4; j++) {
      if( sol->m[pt->v[j]] == 0.0 ) cnt++;
    }
    if(cnt == 4) {
      fprintf(stdout,"Problem in tetra %d : 4 vertices on implicit boundary",k);
      exit(EXIT_FAILURE);
    }
  }

  /** Second test : check whether configuration is manifold in each ball */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ))   continue;
    adja = &mesh->adja[4*(k-1)+1];

    for(i=0; i<4; i++){
      if(!adja[i]) continue;
      iel = adja[i] / 4;
      pt1 = &mesh->tetra[iel];
      if(pt1->ref == pt->ref) continue;

      for(j=0; j<3; j++){
        ip = idir[i][j];

        if(!chkmaniball(mesh,k,ip)){
          fprintf(stdout,"Non orientable implicit surface : ball of point %d\n",pt->v[ip]);
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if ( mesh->info.ddebug )  fprintf(stdout,"  *** Manifold implicit surface.\n");
  return(1);
}

/** Check whether collapse of point np to nq does not create a non manifold situation at nq
    ndepmin, ndepplus = tetra of ref minus, plus in ball of np, not in shell of (np,nq). */
int chkmanicoll(pMesh mesh,int k,int iface,int iedg,int ndepmin,int ndepplus,char isminp,char isplp) {
  pTetra    pt,pt1;
  int       nump,numq,ilist,ref,cur,stor,iel,jel,base,ndepmq,ndeppq;
  int       list[LMAX+2],*adja,*adja1;
  char      i,j,ip,jp,iq,jq,voy,indp,indq,isminq,isplq,ismin,ispl;

  ilist = 0;
  ndepmq = ndeppq = 0;
  isplq = isminq = 0;

  pt    = &mesh->tetra[k];
  ip    = idir[iface][inxt2[iedg]];
  iq    = idir[iface][iprv2[iedg]];
  nump  = pt->v[ip];
  numq  = pt->v[iq];

  /* Case when nump does not have any interior (resp. ext.) tetra which will not
     disappear : search for start in ball of q */
  if ( !ndepmin || !ndepplus ) {
    base  = ++mesh->base;

    pt = &mesh->tetra[k];
    for(j=0; j<4; j++) {
      if ( pt->v[j] == numq ) break;
    }
    assert( j < 4 );
    list[ilist] = 4*k+j;
    ilist++;
    assert( ilist < LMAX+2 );
    pt->flag = base;

    if ( pt->ref == MG_MINUS ) isminq = 1;
    else if ( pt->ref == MG_PLUS ) isplq = 1;

    cur = 0;
    while( cur < ilist ) {
      iel = list[cur] / 4;
      i = list[cur] % 4;
      adja = &mesh->adja[4*(iel-1)+1];

      for (j=0; j<3; j++) {
        i = inxt3[i];
        jel = adja[i] / 4;
        if ( !jel ) continue;

        pt1 = &mesh->tetra[jel];

        if ( pt1->ref == MG_MINUS ) isminq = 1;
        else if ( pt1->ref == MG_PLUS ) isplq = 1;

        if ( pt1->flag == base ) continue;
        pt1->flag = base;
        for(iq=0; iq<4; iq++)
          if ( pt1->v[iq] == numq ) break;
        assert( iq < 4 );
        /* overflow */
        assert ( ilist < LMAX+2 );

        list[ilist] = 4*jel+iq;
        ilist++;

        /* check if jel is an available starting tetra for further enumeration */
        if ( !ndeppq && pt1->ref == MG_PLUS ) {
          for(ip=0; ip<4; ip++)
            if ( pt1->v[ip] == nump ) break;
          if( ip == 4 ) ndeppq = jel;
        }
        if( !ndepmq && pt1->ref == MG_MINUS ) {
          for(ip=0; ip<4; ip++)
            if ( pt1->v[ip] == nump ) break;
          if( ip == 4 ) ndepmq = jel;
        }
      }
      cur++;
    }

    memset(list,0,(LMAX+2)*sizeof(int));
    ilist = 0;
  }

  ispl = ( isplp || isplq ) ? 1 : 0;
  ismin = ( isminp || isminq ) ? 1 : 0;

  /** First step : pile up tetras of future ball of nq, crossing
      through the shell of (np,nq), as long as they have same ref as ndepmin
      list[l] <= 0 if element of ball of np, >= 0, if element of ball of nq */
  base  = ++mesh->base;

  if( ndepmin ) {
    pt = &mesh->tetra[ndepmin];
    ref = pt->ref;

    for(j=0; j<4; j++) {
      if ( pt->v[j] == nump ) break;
    }
    assert( j < 4 );

    pt->flag = base;
    list[ilist] = - (4*ndepmin+j);
    ilist++;
  }
  else if ( ndepmq ) {
    pt = &mesh->tetra[ndepmq];
    ref = pt->ref;

    for(j=0; j<4; j++) {
      if ( pt->v[j] == numq ) break;
    }
    assert( j < 4 );

    pt->flag = base;
    list[ilist] = 4*ndepmq+j;
    ilist++;
  }
  else {
    if ( ismin && ispl )
      return(0);
    else
      return(1);
  }


  cur = 0;
  while ( cur < ilist ) {
    stor = list[cur];
    /* Element belongs to the ball of np */
    if ( stor <= 0 ) {
      stor *= -1;
      iel = stor / 4;
      ip  = stor % 4;

      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jp = ip;
      for (i=0; i<3; i++) {
        jp = inxt3[jp];
        jel = adja[jp] / 4;
        voy = adja[jp] % 4;
        if ( !jel ) continue;

        pt1 = &mesh->tetra[jel];
        if ( pt1->ref != ref ) continue;

        /* Current tetra is neighbour of a tetra of the shell of (np,nq) */
        if( pt1->v[voy] == numq ) {
          adja1 = &mesh->adja[4*(jel-1)+1];
          for(j=0; j<4; j++)
            if (pt1->v[j] == nump ) break;
          assert( j< 4);

          jel = adja1[j] / 4;
          if (!jel ) continue;

          pt1 = &mesh->tetra[jel];

          if ( pt1->ref != ref) continue;   // ICI, il ne faut pas autoriser à passer si on a à nouveau un tet de la coquille (avant de marquer)
          if ( pt1->flag == base ) continue;

          /* New tetra to be added must not be itself an element of the shell */
          for(j=0; j<4; j++) {
            if ( pt1->v[j] == nump ) break;
          }
          if ( j<4 ) continue;

          pt1->flag = base;

          for(j=0; j<4; j++)
            if ( pt1->v[j] == numq ) break;
          assert( j< 4);

          list[ilist] = 4*jel+j;
          ilist++;
          assert( ilist < LMAX+1 );
        }
        else {
          if ( pt1->flag == base ) continue;
          pt1->flag = base;
          for(j=0; j<4; j++)
            if ( pt1->v[j] == nump ) break;
          assert( j< 4 );

          list[ilist] = - (4*jel+j);
          ilist++;
          assert( ilist < LMAX+1 );
        }
      }
    }
    /* Element belongs to the ball of nq */
    else {
      iel = stor / 4;
      iq  = stor % 4;

      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jq = iq;
      for (i=0; i<3; i++) {
        jq = inxt3[jq];
        jel = adja[jq] / 4;
        voy = adja[jq] % 4;
        if ( !jel ) continue;

        pt1 = &mesh->tetra[jel];
        if ( pt1->ref != ref ) continue;

        /* Current tetra is neighbour of a tetra of the shell of (np,nq) */
        if( pt1->v[voy] == nump ) {
          adja1 = &mesh->adja[4*(jel-1)+1];
          for(j=0; j<4; j++)
            if (pt1->v[j] == numq ) break;
          assert( j< 4);

          jel = adja1[j] / 4;
          if (!jel ) continue;

          pt1 = &mesh->tetra[jel];
          if ( pt1->ref != ref) continue;
          if ( pt1->flag == base ) continue;

          /* New tetra to be added must not be itself an element of the shell */
          for(j=0; j<4; j++) {
            if ( pt1->v[j] == numq ) break;
          }
          if ( j<4 ) continue;

          pt1->flag = base;
          for(j=0; j<4; j++)
            if (pt1->v[j] == nump ) break;
          assert( j< 4);

          list[ilist] = -(4*jel+j);
          ilist++;
          assert( ilist < LMAX+1 );
        }
        else {
          if ( pt1->flag == base ) continue;
          pt1->flag = base;
          for(j=0; j<4; j++)
            if (pt1->v[j] == numq ) break;
          assert( j< 4);

          list[ilist] = 4*jel+j;
          ilist++;
          assert( ilist < LMAX+1 );
        }
      }
    }
    cur++;
  }

  assert( cur == ilist );

  /** Second step : same process, starting with a tetra of different reference, in the ball of np */
  if( ndepplus ) {
    pt = &mesh->tetra[ndepplus];
    for(j=0; j<4; j++) {
      if ( pt->v[j] == nump ) break;
    }
    assert( j < 4 );

    pt->flag = base;
    list[ilist] = - (4*ndepplus+j);
    ilist++;
    ref = pt->ref;
  }
  else if ( ndeppq ) {
    pt = &mesh->tetra[ndeppq];
    for(j=0; j<4; j++) {
      if ( pt->v[j] == numq ) break;
    }
    assert( j < 4 );

    pt->flag = base;
    list[ilist] = 4*ndeppq+j;
    ilist++;
    ref = pt->ref;
  }
  else {
    if ( ismin && ispl )
      return(0);
    else
      return(1);
  }

  while ( cur < ilist ) {
    stor = list[cur];
    /* Element belongs to the ball of np */
    if ( stor <= 0 ) {
      stor *= -1;
      iel = stor / 4;
      ip  = stor % 4;

      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jp = ip;

      for (i=0; i<3; i++) {
        jp = inxt3[jp];
        jel = adja[jp] / 4;
        voy = adja[jp] % 4;
        if ( !jel ) continue;

        pt1 = &mesh->tetra[jel];
        if ( pt1->ref != ref ) continue;

        /* Current tetra is neighbour of a tetra of the shell of (np,nq) */
        if( pt1->v[voy] == numq ) {
          adja1 = &mesh->adja[4*(jel-1)+1];
          for(j=0; j<4; j++)
            if (pt1->v[j] == nump ) break;
          assert( j< 4);

          jel = adja1[j] / 4;
          if (!jel ) continue;

          pt1 = &mesh->tetra[jel];
          if ( pt1->ref != ref) continue;
          if ( pt1->flag == base ) continue;

          /* New tetra to be added must not be itself an element of the shell */
          for(j=0; j<4; j++) {
            if ( pt1->v[j] == nump ) break;
          }
          if ( j<4 ) continue;

          pt1->flag = base;
          for(j=0; j<4; j++)
            if ( pt1->v[j] == numq ) break;
          assert( j< 4 );

          list[ilist] = 4*jel+j;
          ilist++;
          assert( ilist < LMAX+1 );
        }
        else {
          if ( pt1->flag == base ) continue;
          pt1->flag = base;
          for(j=0; j<4; j++)
            if (pt1->v[j] == nump ) break;
          assert( j< 4);

          list[ilist] = - (4*jel+j);
          ilist++;
          assert( ilist < LMAX+1 );
        }
      }
    }
    /* Element belongs to the ball of nq */
    else {
      iel = stor / 4;
      iq  = stor % 4;

      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jq = iq;

      for (i=0; i<3; i++) {
        jq = inxt3[jq];
        jel = adja[jq] / 4;
        voy = adja[jq] % 4;

        if ( !jel ) continue;

        pt1 = &mesh->tetra[jel];

        if ( pt1->ref != ref ) continue;

        /* Current tetra is neighbour of a tetra of the shell of (np,nq) */
        if( pt1->v[voy] == nump ) {
          adja1 = &mesh->adja[4*(jel-1)+1];
          for(j=0; j<4; j++)
            if (pt1->v[j] == numq ) break;
          assert( j< 4);

          jel = adja1[j] / 4;
          if (!jel ) continue;

          pt1 = &mesh->tetra[jel];
          if ( pt1->ref != ref) continue;
          if ( pt1->flag == base ) continue;

          /* New tetra to be added must not be itself an element of the shell */
          for(j=0; j<4; j++) {
            if ( pt1->v[j] == numq ) break;
          }
          if ( j<4 ) continue;

          pt1->flag = base;
          for(j=0; j<4; j++)
            if (pt1->v[j] == nump ) break;
          assert( j< 4);

          list[ilist] = -(4*jel+j);
          ilist++;
          assert( ilist < LMAX+1 );
        }
        else {
          if ( pt1->flag == base ) continue;
          pt1->flag = base;
          for(j=0; j<4; j++)
            if (pt1->v[j] == numq ) break;
          assert( j< 4 );

          list[ilist] = 4*jel+j;
          ilist++;
          assert( ilist < LMAX+1 );
        }
      }
    }
    cur++;
  }
  assert( cur == ilist );

  /* At this point, all elements of ball np \cup ball nq \setminus shell have been tagged
     unless the future ball of nq, ending up from collapse is non manifold */
  cur = 0;
  while ( cur < ilist ) {
    stor = list[cur];

    if ( stor <= 0 ) {
      stor *= -1;
      iel = stor / 4;
      ip  = stor % 4;
      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jp = ip;
      for(i=0; i<3; i++) {
        jp = inxt3[jp];
        jel = adja[jp] / 4;

        if ( !jel ) continue;
        pt1 = &mesh->tetra[jel];
        if (pt1->flag == base ) continue;
        pt1->flag = base;

        indp = -1;
        indq = -1;
        for(j=0; j<4; j++) {
          if ( pt1->v[j] == nump )
            indp = j;
          else if ( pt1->v[j] == numq )
            indq = j;
        }
        assert( indp >= 0 && indp < 4 );

        /* Only tets of the shell of (np,nq) can be added, unless future ball is non manifold */
        if ( indq == -1 ) {
          fprintf(stdout,"On devrait passer ici TRES rarement : ");
          fprintf(stdout,"tetra numero %d =  %d %d %d %d , sa ref = %d\n",
                  jel,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3],pt1->ref);
          return(0);
        }

        list[ilist] = -(4*jel+indp);
        ilist++;
        assert( ilist < LMAX +1 );
      }
    }
    else {
      iel = stor / 4;
      iq  = stor % 4;
      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jq = iq;
      for(i=0; i<3; i++) {
        jq = inxt3[jq];
        jel = adja[jq] / 4;

        if ( !jel ) continue;
        pt1 = &mesh->tetra[jel];
        if (pt1->flag == base ) continue;
        pt1->flag = base;

        indp = -1;
        indq = -1;

        for(j=0; j<4; j++) {
          if ( pt1->v[j] == nump )
            indp = j;
          else if ( pt1->v[j] == numq )
            indq = j;
        }
        assert( indq >= 0 && indq < 4 );

        /* Only tets of the shell of (np,nq) can be added, unless future ball is non manifold */
        if ( indp == -1 ) {
          fprintf(stdout,"On devrait passer ici TRES rarement : ");
          fprintf(stdout,"tetra numero %d =  %d %d %d %d , sa ref = %d\n",
                  jel,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3],pt1->ref);
          return(0);
        }

        list[ilist] = 4*jel+indq;
        ilist++;
        assert( ilist < LMAX +1 );
      }
    }
    cur++;
  }

  return(1);
}

/** Create implicit surface in mesh */
int mmg3d2(pMesh mesh,pSol sol) {
  double   *tmp;

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION\n");

  ADD_MEM(mesh,(mesh->npmax+1)*sizeof(double),"temporary table",
          printf("  Exit program.\n");
          exit(EXIT_FAILURE));
  SAFE_CALLOC(tmp,mesh->npmax+1,double);

  /* Snap values of level set function if need be, then discretize it */
  if ( !snpval_ls(mesh,sol,tmp) ) {
    fprintf(stdout,"  ## Problem with implicit function. Exit program.\n");
    return(0);
  }
  DEL_MEM(mesh,tmp,(mesh->npmax+1)*sizeof(double));

  if ( !hashTetra(mesh,1) ) {
    fprintf(stdout,"  ## Hashing problem. Exit program.\n");
    return(0);
  }
  if ( !chkNumberOfTri(mesh) ) {
    if ( !bdryTria(mesh) ) {
      fprintf(stdout,"  ## Boundary problem. Exit program.\n");
      return(0);
    }
    freeXTets(mesh);
  }
  else if ( !bdryPerm(mesh) ) {
    fprintf(stdout,"  ## Boundary orientation problem. Exit program.\n");
    return(0);
  }

  /* build hash table for initial edges */
  if ( !hGeom(mesh) ) {
    fprintf(stdout,"  ## Hashing problem (0). Exit program.\n");
    return(0);
  }

  if ( !bdrySet(mesh) ) {
    fprintf(stdout,"  ## Problem in setting boundary. Exit program.\n");
    return(0);
  }

#ifdef SINGUL
  if ( !chkedg_ls(mesh,sol) ) {
    fprintf(stdout,"  ## Warning: some singular edges will be lost.\n");
  }
#endif

  if ( !cuttet_ls(mesh,sol/*,tmp*/) ) {
    fprintf(stdout,"  ## Problem in discretizing implicit function. Exit program.\n");
    return(0);
  }

  DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));
  DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(Tria));
  mesh->nt = 0;

  if ( !setref_ls(mesh,sol) ) {
    fprintf(stdout,"  ## Problem in setting references. Exit program.\n");
    return(0);
  }

  /* Clean memory (but not pointer) */
  DEL_MEM(mesh,sol->m,(sol->size*sol->npmax+1)*sizeof(double));
  sol->np = 0;

  return(1);
}
