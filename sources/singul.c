#ifdef SINGUL
#include "mmg3d.h"

extern Info  info;
int basetet = 0;
int basept  = 1;

/** Create a key to know where is the point of barycentric
 *  coordinates cb in the tetra:
 *  key=0 for an inside point,
 *  key=1,2,3,4 if the point is on the face 0,1,2,3 (resp.),
 *  key=5,6,7,8,9,10 if the point is on edge 0,1,2,3,4,5 (resp.),
 *  key=11,12,13,14 if the point is on vertex 0,1,2,3 (resp.).*/
static inline
int locate(double cb[4]) {
  int i, j, key;

  /* the point must be in the tet */
  assert ( (-EPS<=cb[0]) && (cb[0]<=1+EPS) );
  assert ( (-EPS<=cb[1]) && (cb[1]<=1+EPS) );
  assert ( (-EPS<=cb[2]) && (cb[2]<=1+EPS) );
  assert ( (-EPS<=cb[3]) && (cb[3]<=1+EPS) );
  assert ( fabs((cb[0] + cb[1] + cb[2] + cb[3])-1.0)<EPS );

  for ( i=0; i<4; i++ ) {
    if ( fabs(cb[i]-1.0) < EPS ) {
      /* the point is the i^th vertex of the tet */
      key = 11+i;
      return(key);
    }
  }

  for ( i=0; i<4; i++ ) {
    if ( fabs(cb[i]) < EPS ) {
      for ( j=0; j<3; j++ ) {
        if ( fabs(cb[idir[i][j]]) < EPS ) {
          /* the point is on edge key */
          if ( !i )
            key = 11-idir[i][j];
          else
            key = 10-i-idir[i][j];
          return(key);
        }
        continue;
      }
      /* the point is on the i^th tri */
      key = i+1;
      return(key);
    }
    continue;
  }
  /* the point is inside the tet */
  return(0);
}

/** update cb, key, ptt and p values for the next visited tet
    when we travel through a vertex is using the tri ifac */
static inline
int updatevertex(pMesh mesh, pTetra *ptt, int nsfin, pPoint p[4],
                 int is, int ifac, double cb[4], int *key) {
  int i ,ip, vois;

  if ( !nsfin ) {
    fprintf(stdout,"%s:%d: Error: no tetra\n",__FILE__,__LINE__);
    return(0);
  }

  ip = (*ptt)->v[is];

  (*ptt)  = &mesh->tetra[nsfin];
  p[0]    = &mesh->point[(*ptt)->v[0]];
  p[1]    = &mesh->point[(*ptt)->v[1]];
  p[2]    = &mesh->point[(*ptt)->v[2]];
  p[3]    = &mesh->point[(*ptt)->v[3]];

  for ( i=0; i<3; i++ ) {
    if ( (*ptt)->v[idir[ifac][i]] == ip )  break;
  }
  assert(i<3);

  vois    = idir[ifac][i];
  (*key)  = 11 + vois;

  cb[0]    = 0.0;
  cb[1]    = 0.0;
  cb[2]    = 0.0;
  cb[3]    = 0.0;
  cb[vois] = 1.0;

  return(1);
}

/** update cb, key, ptt and p values for the next visited tet
    when we travel through an edge ia using the tri ifac */
static inline
int updateedge(pMesh mesh, pTetra *ptt, int nsfin, pPoint p[4], int ia,
               int ifac, double cb[4], int *key) {
  int    i, ip0, ip1, ja;
  double cbtmp[4];

  if ( !nsfin ) {
    fprintf(stdout,"%s:%d: Error: no tetra\n",__FILE__,__LINE__);
    return(0);
  }

  ip0 = (*ptt)->v[iare[ia][0]];
  ip1 = (*ptt)->v[iare[ia][1]];

  (*ptt)   = &mesh->tetra[nsfin];
  p[0]     = &mesh->point[(*ptt)->v[0]];
  p[1]     = &mesh->point[(*ptt)->v[1]];
  p[2]     = &mesh->point[(*ptt)->v[2]];
  p[3]     = &mesh->point[(*ptt)->v[3]];
  cbtmp[0] = cb[0];
  cbtmp[1] = cb[1];
  cbtmp[2] = cb[2];
  cbtmp[3] = cb[3];
  cb[0]    = 0.0;
  cb[1]    = 0.0;
  cb[2]    = 0.0;
  cb[3]    = 0.0;

  for ( i=0; i<3; i++ ) {
    ja = iarf[ifac][i];
    if ( ((*ptt)->v[iare[ja][0]] == ip0 && (*ptt)->v[iare[ja][1]] == ip1) ) {
      cb[iare[ja][0]] = cbtmp[iare[ia][0]];
      cb[iare[ja][1]] = cbtmp[iare[ia][1]];
      (*key) = 5+ja;
      break;
    } else if ( (*ptt)->v[iare[ja][0]] == ip1 && (*ptt)->v[iare[ja][1]] == ip0 ) {
      cb[iare[ja][1]] = cbtmp[iare[ia][0]];
      cb[iare[ja][0]] = cbtmp[iare[ia][1]];
      (*key) = 5+ja;
      break;
    }
  }
  assert( i<3 );
  return(1);
}

/** update cb, key, ptt and p values for the next visited tet
    when we travel through a tri ifac whose vertices are i0, i1 and i2 */
static inline
int updatefac(pMesh mesh, pTetra *ptt, int nsfin, pPoint p[4],
              int i0, int i1, int i2, int ifac, double cb[4], int *key) {
  int    i,ip0,ip1,ip2;
  double cbtmp[4];

  if ( !nsfin ) {
    fprintf(stdout,"%s:%d: Error: no tetra\n",__FILE__,__LINE__);
    return(0);
  }

  ip0  = (*ptt)->v[i0];
  ip1  = (*ptt)->v[i1];
  ip2  = (*ptt)->v[i2];

  (*ptt) = &mesh->tetra[nsfin];
  p[0]   = &mesh->point[(*ptt)->v[0]];
  p[1]   = &mesh->point[(*ptt)->v[1]];
  p[2]   = &mesh->point[(*ptt)->v[2]];
  p[3]   = &mesh->point[(*ptt)->v[3]];

  cbtmp[0] = cb[0];
  cbtmp[1] = cb[1];
  cbtmp[2] = cb[2];
  cbtmp[3] = cb[3];

  for ( i=0; i<3; i++ ) {
    if ( (*ptt)->v[idir[ifac][i]] == ip0 )  break;
  }
  assert( i<3 );
  cb[idir[ifac][i]] = cbtmp[i0];

  if ( (*ptt)->v[idir[ifac][iprv2[i]]] == ip1 ) {
    i = iprv2[i];
    cb[idir[ifac][i]] = cbtmp[i1];

    i = iprv2[i];
    assert( (*ptt)->v[idir[ifac][i]] == ip2 );
    cb[idir[ifac][i]] = cbtmp[i2];
  }
  else {
    i = inxt2[i];
    assert( (*ptt)->v[idir[ifac][i]] == ip1 );
    cb[idir[ifac][i]] = cbtmp[i1];

    i = inxt2[i];
    assert( (*ptt)->v[idir[ifac][i]] == ip2 );
    cb[idir[ifac][i]] = cbtmp[i2];
  }

  cb[ifac] = 0.0;
  (*key) = ifac+1;

  return(1);
}

/** update cb, key, ptt and p values for the next visited tet
    when we don't know where we travel (vertex, edge, tri...)
    except that we come from the previous tet by tri ifac */
static inline
int updateunknown(pMesh mesh, pTetra *ptt, int nsfin, pPoint p[4],
                  int ifac, double cb[4], int *key) {
  int loc;

  if ( (*key)/11 ) {
    loc = (*key)%11;
    if ( !updatevertex(mesh, ptt, nsfin, p, loc, ifac, cb, key) )
      return(0);
    return(1);
  }

  else if ( (*key)/5 ) {
    loc = (*key)-5;
    if ( !updateedge(mesh, ptt, nsfin, p, loc, ifac, cb, key) )
      return(0);
    return(1);
  }

  else if ( *key ) {
    loc  = (*key)-1;
    if ( !updatefac(mesh, ptt, nsfin, p, idir[loc][0], idir[loc][1],
                    idir[loc][2], ifac, cb, key) )  return(0);
    return(1);
  }

  else {
    fprintf(stdout,"%s:%d: Error: unexpected case: %d \n",
            __FILE__,__LINE__,*key);
    return(0);
  }
}

/** fill mat and ap */
static inline
int barycentric(int is, int i0, int i1, int i2, pPoint p[4], double e2[3],
                double mat[3][3], double ap[3]) {

  /* barycentric */
  mat[0][0]  = p[i0]->c[0] - p[is]->c[0];
  mat[1][0]  = p[i0]->c[1] - p[is]->c[1];
  mat[2][0]  = p[i0]->c[2] - p[is]->c[2];
  mat[0][1]  = p[i1]->c[0] - p[is]->c[0];
  mat[1][1]  = p[i1]->c[1] - p[is]->c[1];
  mat[2][1]  = p[i1]->c[2] - p[is]->c[2];
  mat[0][2]  = p[i2]->c[0] - p[is]->c[0];
  mat[1][2]  = p[i2]->c[1] - p[is]->c[1];
  mat[2][2]  = p[i2]->c[2] - p[is]->c[2];

  ap[0] = e2[0] - p[is]->c[0];
  ap[1] = e2[1] - p[is]->c[1];
  ap[2] = e2[2] - p[is]->c[2];

  return(1);
}

/** compute and check the volume of the tet */
static inline
int calcvol(double mat[3][3], double *dd) {

  /* test volume */
  (*dd) =  mat[0][0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])
    +      mat[1][0]*(mat[2][1]*mat[0][2] - mat[0][1]*mat[2][2])
    +      mat[2][0]*(mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2]);

  if ( (*dd) > 0.0 ) return(1);
  if ( info.ddebug ) {
    fprintf(stdout,"%s:%d: Error: wrong orientation of element\n",
            __FILE__,__LINE__);
    exit(1);
  }
  return(0);
}



/** Compute the intersection point between tri i0 and seg e1-e2 */
static inline
int intersegtri_i0(int is, int i0, int i1, int i2, double ee[3], double e1[3],
                   double ep[3], double dd, double pis[3], double mat[3][3],
                   double u[3], double c[3], double cb[4],
                   int *key, Travel *travcell, int* nsfin, int adji0) {

  double t, cp[3], v[3], vol[4];

  assert(ee[0]*u[0]+ee[1]*u[1]+ee[2]*u[2]);
  t = (-ep[0]*u[0]-ep[1]*u[1]-ep[2]*u[2])/(ee[0]*u[0]+ee[1]*u[1]+ee[2]*u[2]);

  /* coordinates of the intersection point of (e1,e2) with tri */
  c[0]  = ee[0]*t+e1[0];
  c[1]  = ee[1]*t+e1[1];
  c[2]  = ee[2]*t+e1[2];
  cp[0] = c[0]-pis[0];
  cp[1] = c[1]-pis[1];
  cp[2] = c[2]-pis[2];
  v[0]  = mat[1][0]*cp[2] - mat[2][0]*cp[1];
  v[1]  = mat[2][0]*cp[0] - mat[0][0]*cp[2];
  v[2]  = mat[0][0]*cp[1] - mat[1][0]*cp[0];
  /* is the point inside the triangle?
   * (if yes, vol[is], vol[i1] and vol[i2]>0) */
  vol[i1] =  (mat[0][2]*v[0]+mat[1][2]*v[1]+mat[2][2]*v[2]);

  if ( -EPS2 <= vol[i1] ) {
    vol[i2] = (-mat[0][1]*v[0]-mat[1][1]*v[1]-mat[2][1]*v[2]);

    if ( -EPS2 <= vol[i2] ) {
      vol[is] = dd-vol[i1]-vol[i2];

      if ( -EPS2 <= vol[is] )  {
        /* We have found our intersection */
        dd = 1.0 / dd;

        cb[i0] = 0.0;
        cb[is] = vol[is] * dd;
        cb[i1] = vol[i1] * dd;
        cb[i2] = vol[i2] * dd;
        (*key) = locate(cb);
        travcell->cb[0] = cb[0];
        travcell->cb[1] = cb[1];
        travcell->cb[2] = cb[2];
        travcell->cb[3] = cb[3];
        travcell->c[0]  = c[0];
        travcell->c[1]  = c[1];
        travcell->c[2]  = c[2];
        travcell->key   = *key;

        (*nsfin) = adji0;
        return(1);
      }
    }
  }
  return(0);
}

/** Compute the intersection point between tri i1 and seg e1-e2 */
static inline
int intersegtri_i1(int is, int i0, int i1, int i2, double ee[3], double e1[3],
                   double ep[3], double dd, double pis[3], double mat[3][3],
                   double c[3], double cb[4], int *key,
                   Travel *travcell, int* nsfin, int adji1) {

  double t, cp[3], v[3], vol[4];

  v[0]  = mat[1][0]*mat[2][2] - mat[2][0]*mat[1][2];
  v[1]  = mat[2][0]*mat[0][2] - mat[0][0]*mat[2][2];
  v[2]  = mat[0][0]*mat[1][2] - mat[1][0]*mat[0][2];

  assert(ee[0]*v[0]+ee[1]*v[1]+ee[2]*v[2]);

  t = (-ep[0]*v[0]-ep[1]*v[1]-ep[2]*v[2])/(ee[0]*v[0]+ee[1]*v[1]+ee[2]*v[2]);
  /* coordinates of the intersection point of (e1,e2) with tri */
  c[0]  = ee[0]*t+e1[0];
  c[1]  = ee[1]*t+e1[1];
  c[2]  = ee[2]*t+e1[2];
  cp[0] = c[0]-pis[0];
  cp[1] = c[1]-pis[1];
  cp[2] = c[2]-pis[2];

  /* is the point inside the triangle?
   * (if yes, vol[is], vol[i0] and vol[i2]>0) */
  v[0]  = mat[1][1]*cp[2] - mat[2][1]*cp[1];
  v[1]  = mat[2][1]*cp[0] - mat[0][1]*cp[2];
  v[2]  = mat[0][1]*cp[1] - mat[1][1]*cp[0];

  vol[i0] = (-mat[0][2]*v[0]-mat[1][2]*v[1]-mat[2][2]*v[2]);

  if ( -EPS2 <= vol[i0] ) {
    vol[i2] = (mat[0][0]*v[0]+mat[1][0]*v[1]+mat[2][0]*v[2]);

    if ( -EPS2 <= vol[i2] ) {
      vol[is] = dd-vol[i0]-vol[i2];

      if ( -EPS2 <= vol[is] )  {

        /* We have found our intersection */
        dd = 1.0 / dd;

        cb[i1] = 0.0;
        cb[i0] = vol[i0] * dd;
        cb[is] = vol[is] * dd;
        cb[i2] = vol[i2] * dd;
        (*key) = locate(cb);
        travcell->cb[0] = cb[0];
        travcell->cb[1] = cb[1];
        travcell->cb[2] = cb[2];
        travcell->cb[3] = cb[3];
        travcell->c[0]  = c[0];
        travcell->c[1]  = c[1];
        travcell->c[2]  = c[2];
        travcell->key   = (*key);

        (*nsfin) = adji1;
        return(1);
      }
    }
  }
  return(0);
}

/** Compute the intersection point between tri i2 and seg e1-e2 */
static inline
int intersegtri_i2(int is, int i0, int i1, int i2, double ee[3], double e1[3],
                   double ep[3], double dd, double pis[3], double mat[3][3],
                   double c[3], double cb[4], int *key,
                   Travel *travcell, int* nsfin, int adji2) {

  double t, cp[3], v[3], vol[4];

  v[0]  = mat[1][1]*mat[2][0] - mat[2][1]*mat[1][0];
  v[1]  = mat[2][1]*mat[0][0] - mat[0][1]*mat[2][0];
  v[2]  = mat[0][1]*mat[1][0] - mat[1][1]*mat[0][0];

  assert(ee[0]*v[0]+ee[1]*v[1]+ee[2]*v[2]);

  t = (-ep[0]*v[0]-ep[1]*v[1]-ep[2]*v[2])/(ee[0]*v[0]+ee[1]*v[1]+ee[2]*v[2]);
  /* coordinates of the intersection point of (e1,e2) with tri */
  c[0]  = ee[0]*t+e1[0];
  c[1]  = ee[1]*t+e1[1];
  c[2]  = ee[2]*t+e1[2];
  cp[0] = c[0]-pis[0];
  cp[1] = c[1]-pis[1];
  cp[2] = c[2]-pis[2];

  /* is the point inside the triangle?
   * (if yes, vol[is], vol[i0] and vol[i1]>0) */
  v[0]  = mat[1][2]*cp[2] - mat[2][2]*cp[1];
  v[1]  = mat[2][2]*cp[0] - mat[0][2]*cp[2];
  v[2]  = mat[0][2]*cp[1] - mat[1][2]*cp[0];

  vol[i0] = (mat[0][1]*v[0]+mat[1][1]*v[1]+mat[2][1]*v[2]);

  if ( -EPS2 <= vol[i0] ) {
    vol[i1] = (-mat[0][0]*v[0]-mat[1][0]*v[1]-mat[2][0]*v[2]);

    if ( -EPS2 <= vol[i1] ) {
      vol[is] = dd-vol[i0]-vol[i1];

      if ( -EPS2 <= vol[is] )  {

        /* We have found our intersection */
        dd = 1.0 / dd;

        cb[i2] = 0.0;
        cb[i0] = vol[i0] * dd;
        cb[i1] = vol[i1] * dd;
        cb[is] = vol[is] * dd;
        (*key) = locate(cb);
        travcell->cb[0] = cb[0];
        travcell->cb[1] = cb[1];
        travcell->cb[2] = cb[2];
        travcell->cb[3] = cb[3];
        travcell->c[0]  = c[0];
        travcell->c[1]  = c[1];
        travcell->c[2]  = c[2];
        travcell->key   = (*key);

        (*nsfin) = adji2;
        return(1);
      }
    }
  }
  return(0);
}


/** Create a point singular in tetra iel and if needed split the tetra and its neighbour */
int creaPoint(pMesh mesh, pSol met, int iel,double c[3], double cb[4], char tag){
  pTetra pt;
  double hnew;
  int    ia,ip,ilist,key,i;
  int    list[LMAX+2];

  key=locate(cb);
  pt = &mesh->tetra[iel];

  if ( key/11 ) {
    /* the point is a vertex of the tet */
    mesh->point[pt->v[key%11]].tag |= tag;
    return(1);
  }

  if ( key/5 ) {
    /* the point is on the edge ia */
    hnew = 0.0;
    if ( met->m ) {
      for (i=0; i<4; i++) {
        ip    = mesh->tetra[iel].v[i];
        hnew += met->m[ip]*cb[i];
      }
    }
    ip = newPt(mesh,c,tag );
    if ( !ip ) {
      fprintf(stdout,"  ## Allocation problem (point), not enough memory.\n");
      fprintf(stdout,"  ## Increase the allocated memory with the -m option.\n");
      fprintf(stdout,"  Exit program.\n");
      return(0);
    }
    if ( met->m )  met->m[ip] = hnew;

    ia = key-5;
    ilist = coquil(mesh,iel,ia,list);
    if ( (pt->tag & MG_REQ) || !ilist ) {
      fprintf(stdout,"  ## Unable to insert singularity: element required.\n");
      fprintf(stdout,"  ## Delete required elements.\n");
      fprintf(stdout,"  Exit program.\n");
      return(0);
    }
    if ( !split1b(mesh, met, list, ilist, ip, 0) ) return(0);

    /* Update of barycentric coordinate in tet list[0]/6 */
    switch(list[0]%6){
    case 0:
      cb[0] = cb[2] = cb[3] = 0;
      cb[1] = 1.;
      break;
    case 1:
      cb[1] = cb[2] = cb[3] = 0;
      cb[0] = 1.;
      break;
    case 2:
      cb[0] = cb[1] = cb[2] = 0;
      cb[3] = 1.;
      break;
    case 3:
      cb[0] = cb[1] = cb[3] = 0;
      cb[2] = 1.;
      break;
    case 4:
      cb[0] = cb[2] = cb[3] = 0;
      cb[1] = 1.;
      break;
    case 5:
      cb[0] = cb[1] = cb[3] = 0;
      cb[2] = 1.;
      break;
    }
    return(1);
  }

  /* the point is on a tri */
  if ( key ) {
    if ( (pt->tag & MG_REQ) ||
         (pt->xt && (mesh->xtetra[pt->xt].ftag[key-1] & MG_REQ) ) ) {
      fprintf(stdout,"  ## Unable to insert singularity: element or");
      fprintf(stdout," face required.\n");
      fprintf(stdout,"  ## Delete required elements.\n");
      fprintf(stdout,"  Exit program.\n");
      return(0);
    }
    if ( split3cb(mesh, met, iel, key-1, c, cb, &ip)<0 )  return(0);
    mesh->point[ip].tag |= tag;
    return(1);
  }

  /* The point is inside the tet */
  if ( !split4cb(mesh, met, iel, c, cb, &ip) )  return(0);
  mesh->point[ip].tag |= tag;
  return(1);
}


/** Create a required or ridge edge in tetras of list travel */
int creaEdge(pMesh mesh, pSol met, Travel *trav, char tag){
  pTetra  pt;
  double  c[3],cb[4],hnew;
  int     ia,ip,ilist,key,kel,i;
  int     list[LMAX+2];
  char    majcb[6] = {1,0,3,2,1,2};

  /* first tet */
  kel    = trav->kel;
  pt     = &mesh->tetra[kel];
  key    = trav->key;
  c[0]   = trav->c[0];
  c[1]   = trav->c[1];
  c[2]   = trav->c[2];
  cb[0]  = trav->cb[0];
  cb[1]  = trav->cb[1];
  cb[2]  = trav->cb[2];
  cb[3]  = trav->cb[3];

  if ( key/11 ) {
    /* split on a vertex */
    ip = pt->v[key%11];
    mesh->point[ip].tag |= tag;

    /* Update htab */
    hEdge(&mesh->htab,trav->np,pt->v[key%11],0,trav->tag);

    return(1);
  }

  else if ( key/5 ) {
    /* split on an edge */
    hnew = 0.0;
    if ( met->m ) {
      for (i=0; i<4; i++) {
        ip    = mesh->tetra[kel].v[i];
        hnew += met->m[ip]*cb[i];
      }
    }
    ip = newPt(mesh,c,tag);
    if ( !ip ) {
      fprintf(stdout,"  ## Allocation problem (point), not enough memory.\n");
      fprintf(stdout,"  ## Increase the allocated memory with the -m option.\n");
      fprintf(stdout,"  Exit program.\n");
      return(0);
    }
    if ( met->m )  met->m[ip] = hnew;

    ia = key-5;
    ilist = coquil(mesh,kel,ia,list);
    if ( (pt->tag & MG_REQ) || !ilist ) {
      fprintf(stdout,"  ## Unable to insert singularity: element required.\n");
      fprintf(stdout,"  ## Delete required elements.\n");
      fprintf(stdout,"  Exit program.\n");
      return(0);
    }
    if ( !split1b(mesh, met, list, ilist, ip, 0) )  return(0);

    /* Update of barycentric coordinate in tet list[0]/6 */
    i = majcb[list[0]%6];
    trav->cb[i] = 1.;
    trav->cb[iprv3[i]] = 0.;
    trav->cb[inxt3[i]] = 0.;
    trav->cb[inxt3[inxt3[i]]] = 0.;

    /* Update htab */
    hEdge(&mesh->htab,trav->np,ip,0,trav->tag);
    return(1);
  }

  else if ( key ) {
    /* split on a tri */
    if ( (pt->tag & MG_REQ) ||
         (pt->xt && (mesh->xtetra[pt->xt].ftag[key-1] & MG_REQ) ) ) {
      fprintf(stdout,"  ## Unable to insert singularity: element or");
      fprintf(stdout," face required.\n");
      fprintf(stdout,"  ## Delete required elements.\n");
      fprintf(stdout,"  Exit program.\n");
      return(0);
    }
    i = split3cb(mesh,met,kel,key-1,c,trav->cb,&ip);
    if ( i < 0 ) return(0);

    /* update barycentric coor of new point in kel */
    trav->cb[i] = 1.;
    trav->cb[iprv3[i]] = 0.;
    trav->cb[inxt3[i]] = 0.;
    trav->cb[inxt3[inxt3[i]]] = 0.;
    mesh->point[ip].tag |= tag;

    /* Update htab */
    hEdge(&mesh->htab,trav->np,ip,0,trav->tag);
    return(1);
  }

  else {
    /* split inside the tet */
    if ( !split4cb(mesh, met, kel, c,trav->cb,&ip) )  return(0);

    /* update barycentric coor of new point in kel */
    trav->cb[1] = trav->cb[2] = trav->cb[3] = 0.;
    trav->cb[0] = 1.;
    mesh->point[ip].tag |= tag;

    /* Update htab */
    hEdge(&mesh->htab,trav->np,ip,0,trav->tag);
    return(1);
  }
  return(0);
}

/** Return the index of element to which belong the point of coordinates point[3] and compute
 *  the barycentric coordinates of the point  */
int seekPoint(pMesh mesh, psPoint ppt, double cb[4]) {
  int      nsfin,nstart;
  pTetra   pt;
  pPoint   p0,p1,p2,p3;
  double   point[3];
  double   bx,by,bz,cx,cy,cz,dx,dy,dz,ux,uy,uz,apx,apy,apz;
  double   vol0,vol1,vol2,vol3,dd;
  int     *adj,iadr;

  point[0] = ppt->c[0];
  point[1] = ppt->c[1];
  point[2] = ppt->c[2];

  if ( ppt->tet ) {
    nsfin  = ppt->tet;
    nstart = nsfin;
  }
  else {
    nsfin  = 1;
    nstart = 1;
  }
  basetet++;

  do {
    if ( !nsfin ) {
      fprintf(stdout,"%s:%d: Error: no tetra\n",__FILE__,__LINE__);
      return(0);
    }
    pt = &mesh->tetra[nsfin];
    if ( !pt->v[0] ) {
      fprintf(stdout,"%s:%d: Error: wrong element\n",__FILE__,__LINE__);
      return(0);
    }
    if ( pt->flag == basetet ) {
      nstart  = nstart%(mesh->ne)+1;
      if ( nstart == ppt->tet ) {
        fprintf(stdout,"%s:%d: Error: all elements already checked\n",
                __FILE__,__LINE__);
        return(0);
      }

      nsfin = nstart;
    }
    pt->flag = basetet;

    iadr = 4*(nsfin-1)+1;
    adj  = &mesh->adja[iadr];
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    p3 = &mesh->point[pt->v[3]];

    /* barycentric */
    bx  = p1->c[0] - p0->c[0];
    by  = p1->c[1] - p0->c[1];
    bz  = p1->c[2] - p0->c[2];
    cx  = p2->c[0] - p0->c[0];
    cy  = p2->c[1] - p0->c[1];
    cz  = p2->c[2] - p0->c[2];
    dx  = p3->c[0] - p0->c[0];
    dy  = p3->c[1] - p0->c[1];
    dz  = p3->c[2] - p0->c[2];

    /* test volume */
    ux  = cy*dz - cz*dy;
    uy  = cz*dx - cx*dz;
    uz  = cx*dy - cy*dx;

    dd =  bx*ux + by*uy + bz*uz;

    if ( info.ddebug && dd < 0.0 )
      fprintf(stdout,"%s:%d: Warning: wrong orientation of tetra %d\n",
              __FILE__,__LINE__,nsfin);

    apx = point[0] - p0->c[0];
    apy = point[1] - p0->c[1];
    apz = point[2] - p0->c[2];

    /* p in 1 */
    vol1  = apx*ux + apy*uy + apz*uz;
    if ( -EPS2 > vol1 ) {
      nsfin = adj[1] / 4;
      if ( nsfin )  continue;
      else goto end;
    }

    /* p in 2 */
    ux  = by*apz - bz*apy;
    uy  = bz*apx - bx*apz;
    uz  = bx*apy - by*apx;
    vol2 = dx*ux + dy*uy + dz*uz;
    if ( -EPS2 > vol2 ) {
      nsfin = adj[2] / 4;
      if ( nsfin )  continue;
      else goto end;
    }

    /* p in 3 */
    vol3 = -cx*ux - cy*uy - cz*uz;
    if ( -EPS2 > vol3 ) {
      nsfin = adj[3] / 4;
      if ( nsfin )  continue;
      else goto end;
    }

    /* p in 0 */
    vol0 = dd - vol1 - vol2 - vol3;
    if ( -EPS2 > vol0 ) {
      nsfin = adj[0] / 4;
      if ( nsfin )  continue;
    }

  end:
    if ( vol0 < -EPS2 || vol1 < -EPS2 || vol2 < -EPS2 || vol3 < -EPS2 ) {
      nsfin=adj[0]/4;
      if ( nsfin && !mesh->tetra[nsfin].flag ) continue;
      nsfin=adj[1]/4;
      if ( nsfin && !mesh->tetra[nsfin].flag ) continue;
      nsfin=adj[2]/4;
      if ( nsfin && !mesh->tetra[nsfin].flag ) continue;
      nsfin=adj[3]/4;
      if ( nsfin && !mesh->tetra[nsfin].flag ) continue;

      /* All the neighbours of our tet are marked so we restart from another tet */
      nsfin  = nstart%(mesh->ne)+1;
      nstart = nsfin;
      continue;
    }

    dd = vol0+vol1+vol2+vol3;
    if ( dd != 0.0 )  dd = 1.0 / dd;
    else {
      fprintf(stdout,"%s:%d: Error: tetrahedra volume is null\n",
             __FILE__,__LINE__);
      return(0);
    }

    cb[0]     = vol0 * dd;
    cb[1]     = vol1 * dd;
    cb[2]     = vol2 * dd;
    cb[3]     = vol3 * dd;
    ppt->tet  = nsfin;
    ppt->flag = basept;

    return(1);
  }while(1);

  return(0);
}

/** Return the index of element to which belong the point of coordinates point[3] and compute
 *  cb, the barycentric coordinates of the point  */
int seekEdge(pMesh mesh, pSol met, psPoint ppt0, psPoint ppt1,
             Travel *trav, int *lastet){
  pTetra pt;
  pPoint p[4];
  double mat[3][3],u[3],v[3],e0[3],e1[3];
  double ap[3],ep[3],ee[3],c[3],vol[4];
  double dd,cb[4];
  int    it,is,i0,i1,i2,ind,nsfin,*adj,iadr,key;

  it       = 0;

  if ( (*lastet) ) {

    e0[0] = ppt0->c[0];
    e0[1] = ppt0->c[1];
    e0[2] = ppt0->c[2];
    e1[0] = ppt1->c[0];
    e1[1] = ppt1->c[1];
    e1[2] = ppt1->c[2];

    if ( !seekPoint(mesh,ppt0,cb) ) {
      fprintf(stdout,"%s:%d: Error: point not found %e %e %e\n",
             __FILE__,__LINE__,e0[0],e0[1],e0[2]);
      return(0);
    }
    if ( !creaPoint(mesh, met, ppt0->tet, ppt0->c, cb, ppt0->tag) ) {
      fprintf(stdout,"%s:%d: Error: unable to create point %e %e %e\n",
             __FILE__,__LINE__,e0[0],e0[1],e0[2]);
      return(0);
    }
    nsfin = ppt0->tet;
    key   = locate(cb);

  }
  else {
    e0[0] = trav->c[0];
    e0[1] = trav->c[1];
    e0[2] = trav->c[2];

    cb[0] = trav->cb[0];
    cb[1] = trav->cb[1];
    cb[2] = trav->cb[2];
    cb[3] = trav->cb[3];

    e1[0] = ppt1->c[0];
    e1[1] = ppt1->c[1];
    e1[2] = ppt1->c[2];

    nsfin = trav->kel;
    key   = locate(cb);
  }
  c[0] = e0[0];
  c[1] = e0[1];
  c[2] = e0[2];

  basetet++;

  ee[0] = e1[0]-e0[0];
  ee[1] = e1[1]-e0[1];
  ee[2] = e1[2]-e0[2];

  pt   = &mesh->tetra[nsfin];
  p[0] = &mesh->point[pt->v[0]];
  p[1] = &mesh->point[pt->v[1]];
  p[2] = &mesh->point[pt->v[2]];
  p[3] = &mesh->point[pt->v[3]];

  do {

    if ( !pt->v[0] ) {
      fprintf(stdout,"%s:%d: Error: wrong element\n",__FILE__,__LINE__);
      return(0);
    }

    if ( pt->flag == basetet ) {
      fprintf(stdout,"%s:%d: Error: we pass through an element already seen, tet %d \n",
             __FILE__,__LINE__,nsfin);
      return(0);
    }
    pt->flag = basetet;

    iadr = 4*(nsfin-1)+1;
    adj  = &mesh->adja[iadr];

    /* if, for j!=key%11, lambda_j<0, then we are not on the wanted tet */
    ind = key%11;
    is  = ind;
    i0  = idir[is][0];
    i1  = idir[is][1];
    i2  = idir[is][2];

    barycentric(is,i0,i1,i2,p,e1,mat,ap);
    calcvol(mat, &dd);

    if ( adj[i0] && (mesh->tetra[adj[i0]/4].flag != basetet) ) {
      /* lambda_i0 computation */
      u[0] = mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2];
      u[1] = mat[2][1]*mat[0][2] - mat[0][1]*mat[2][2];
      u[2] = mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2];

      vol[i0] = (u[0]*ap[0]+u[1]*ap[1]+u[2]*ap[2]);

      if ( -EPS2 > vol[i0] ) {
        nsfin = adj[i0]/4;

        /* update of key, cb, pt and p */
        if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i0]%4, cb, &key) )
          return(0);
        continue;
      }

      u[0]  = mat[1][0]*ap[2] - mat[2][0]*ap[1];
      u[1]  = mat[2][0]*ap[0] - mat[0][0]*ap[2];
      u[2]  = mat[0][0]*ap[1] - mat[1][0]*ap[0];
      if ( adj[i1] && (mesh->tetra[adj[i1]/4].flag != basetet) ) {
        /* lambda_i1 computation */
        vol[i1] = (mat[0][2]*u[0] + mat[1][2]*u[1] + mat[2][2]*u[2]);

        if ( -EPS2 > vol[i1] ) {
          nsfin = adj[i1]/4;

          /* update of key, cb, pt and p */
          if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i1]%4, cb, &key) )
            return(0);
          continue;
        }

        /* lambda_i2 computation */
        vol[i2] = (-mat[0][1]*u[0] - mat[1][1]*u[1] - mat[2][1]*u[2]);

        if ( -EPS2 > vol[i2] ) {
          if ( adj[i2] && (mesh->tetra[adj[i2]/4].flag != basetet) ) {
            nsfin = adj[i2]/4;

            /* update of key, cb, pt and p */
            if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i2]%4, cb, &key) )
              return(0);
            continue;
          }
          fprintf(stdout,"%s:%d: Warning: tet %d:",__FILE__,__LINE__,nsfin);
          fprintf(stdout," vol[%d]>0, vol[%d]>0, but vol[%d]<0 \n",i0,i1,i2);
          fprintf(stdout," adjacent by tri %d marked => ",i2);
          fprintf(stdout,"we travel through the wrong tri %d\n",i0);
          nsfin = adj[i0]/4;

          /* update of key, cb, pt and p */
          if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i0]%4, cb, &key) )
            return(0);
          continue;
        }
      }
      else {
        /* lambda_i2 computation */
        vol[i2] = (-mat[0][1]*u[0] - mat[1][1]*u[1] - mat[2][1]*u[2]);
        if ( -EPS2 > vol[i2] ) {
          if ( adj[i2] && (mesh->tetra[adj[i2]/4].flag != basetet) ) {
            nsfin = adj[i2]/4;

            /* update of key, cb, pt and p */
            if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i2]%4, cb, &key) )
              return(0);
            continue;
          }
          fprintf(stdout,"%s:%d: Warning: tet %d:",__FILE__,__LINE__,nsfin);
          fprintf(stdout," vol[%d]>0 but vol[%d]<0 \n",i0,i2);
          fprintf(stdout," adjacent by tri %d marked => ",i2);
          fprintf(stdout,"we travel through the wrong tri %d\n",i0);
          nsfin = adj[i0]/4;

          /* update of key, cb, pt and p */
          if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i0]%4, cb, &key) )
            return(0);
          continue;
        }

        /* lambda_i1 computation */
        vol[i1] = (mat[0][2]*u[0] + mat[1][2]*u[1] + mat[2][2]*u[2]);
        if ( -EPS2 > vol[i1] ) {
          fprintf(stdout,"%s:%d: Warning:",__FILE__,__LINE__);
          fprintf(stdout," vol[%d]>0,vol[%d]>0, but vol[%d]<0 \n",i0,i1,i2);
          fprintf(stdout," adjacent by tri %d marked => ",i1);
          fprintf(stdout,"we travel through the wrong tri %d\n",i0);
          nsfin = adj[i0]/4;

          /* update of key, cb, pt and p */
          if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i0]%4, cb, &key) )
            return(0);
          continue;
        }
      }
    }
    else if ( adj[i1] && (mesh->tetra[adj[i1]/4].flag != basetet) ) {
      /* lambda_i1 computation */
      u[0] = mat[2][0]*mat[1][2] - mat[1][0]*mat[2][2];
      u[1] = mat[0][0]*mat[2][2] - mat[2][0]*mat[0][2];
      u[2] = mat[1][0]*mat[0][2] - mat[0][0]*mat[1][2];

      vol[i1]  = (u[0]*ap[0]+u[1]*ap[1]+u[2]*ap[2]);

      if ( -EPS2 > vol[i1] ) {
        nsfin = adj[i1]/4;

        /* update of key, cb, pt and p */
        if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i1]%4, cb, &key) )
          return(0);
        continue;
      }

      /* lambda_i2 computation */
      u[0]  = mat[1][1]*ap[2] - mat[2][1]*ap[1];
      u[1]  = mat[2][1]*ap[0] - mat[0][1]*ap[2];
      u[2]  = mat[0][1]*ap[1] - mat[1][1]*ap[0];

      vol[i2] = (mat[0][0]*u[0] + mat[1][0]*u[1] + mat[2][0]*u[2]);

      if ( -EPS2 > vol[i2] ) {
        if ( adj[i2] && (mesh->tetra[adj[i2]/4].flag != basetet) ) {
          nsfin = adj[i2]/4;

          /* update of key, cb, pt and p */
          if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i2]%4, cb, &key) )
            return(0);
          continue;
        }
        fprintf(stdout,"%s:%d: Warning: vol[%d]>0, but vol[%d]<0 \n",
               __FILE__,__LINE__,i1,i2);
        fprintf(stdout," adjacent by tri %d marked => ",i2);
        fprintf(stdout,"we travel through the wrong tri %d\n",i1);
        nsfin = adj[i1]/4;

        /* update of key, cb, pt and p */
        if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i1]%4, cb, &key) )
          return(0);
        continue;
      }

      /* lambda_i0 computation */
      vol[i0] = (-mat[0][2]*u[0] - mat[1][2]*u[1] - mat[2][2]*u[2]);

      if ( -EPS2 > vol[i0] ) {
        fprintf(stdout,"%s:%d: Warning:",__FILE__,__LINE__);
        fprintf(stdout," vol[%d]>0, vol[%d]>0 and vol[%d]<0 \n",i1,i2,i0);
        fprintf(stdout," adjacent by tri %d marked => ",i0);
        fprintf(stdout,"we travel through the wrong tri %d\n",i1);
        nsfin = adj[i1]/4;

        /* update of key, cb, pt and p */
        if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i1]%4, cb, &key) )
          return(0);
        continue;
      }
    }
    else if ( adj[i2] && (mesh->tetra[adj[i2]/4].flag != basetet) ) {
      /* lambda_i2 computation */
      u[0] = mat[1][0]*mat[2][1] - mat[2][0]*mat[1][1];
      u[1] = mat[2][0]*mat[0][1] - mat[0][0]*mat[2][1];
      u[2] = mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1];

      vol[i2]  = (u[0]*ap[0]+u[1]*ap[1]+u[2]*ap[2]);

      if ( -EPS2 > vol[i2] ) {
        nsfin = adj[i2]/4;

        /* update of key, cb, pt and p */
        if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i2]%4, cb, &key) )
          return(0);
        continue;
      }
      /* lambda_i0 computation */
      u[0]  = mat[1][2]*ap[2] - mat[2][2]*ap[1];
      u[1]  = mat[2][2]*ap[0] - mat[0][2]*ap[2];
      u[2]  = mat[0][2]*ap[1] - mat[1][2]*ap[0];

      vol[i0] = (mat[0][1]*u[0] + mat[1][1]*u[1] + mat[2][1]*u[2]);

      if ( -EPS2 > vol[i0] ) {
        fprintf(stdout,"%s:%d: Warning: vol[%d]>0 but vol[%d]<0 \n",
               __FILE__,__LINE__,i2,i0);
        fprintf(stdout," adjacent by tri %d marked => ",i0);
        fprintf(stdout,"we travel through the wrong tri %d\n",i2);
        nsfin = adj[i2]/4;

        /* update of key, cb, pt and p */
        if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i2]%4, cb, &key) )
          return(0);
        continue;
      }

      /* lambda_i1 computation */
      vol[i1] = (-mat[0][0]*u[0] - mat[1][0]*u[1] - mat[2][0]*u[2]);

      if ( -EPS2 > vol[i1] ) {
        fprintf(stdout,"%s:%d: Warning:",__FILE__,__LINE__);
        fprintf(stdout," vol[%d]>0, vol[%d]>0 but vol[%d]<0 \n",i2,i0,i1);
        fprintf(stdout," adjacent by tri %d marked => ",i1);
        fprintf(stdout,"we travel through the wrong tri %d\n",i2);
        nsfin = adj[i2]/4;

        /* update of key, cb, pt and p */
        if ( !updatevertex(mesh, &pt, nsfin, p, is, adj[i2]%4, cb, &key) )
          return(0);
        continue;
      }
    }
    else{
      /* lambda_i0 computation */
      u[0] = mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2];
      u[1] = mat[2][1]*mat[0][2] - mat[0][1]*mat[2][2];
      u[2] = mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2];

      vol[i0] = (u[0]*ap[0]+u[1]*ap[1]+u[2]*ap[2]);

      if ( -EPS2 > vol[i0] ) {
        fprintf(stdout,"%s:%d: Error: Element %d, ",__FILE__,__LINE__,nsfin);
        fprintf(stdout,"sommet %d --> vol %d: no-win situation\n",ind,i0);
        return(0);
      }

      /* lambda_i1 computation */
      v[0]  = mat[1][0]*ap[2] - mat[2][0]*ap[1];
      v[1]  = mat[2][0]*ap[0] - mat[0][0]*ap[2];
      v[2]  = mat[0][0]*ap[1] - mat[1][0]*ap[0];

      vol[i1] = (mat[0][2]*v[0] + mat[1][2]*v[1] + mat[2][2]*v[2]);

      if ( -EPS2 > vol[i1] ) {
        fprintf(stdout,"%s:%d: Error: Element %d, ",__FILE__,__LINE__,nsfin);
        fprintf(stdout,"sommet %d --> vol %d: no-win situation\n",ind,i1);
        return(0);
      }

      /* lambda_i2 computation */
      vol[i2] = (-mat[0][1]*v[0] - mat[1][1]*v[1] - mat[2][1]*v[2]);

      if ( -EPS2 > vol[i2] ) {
        fprintf(stdout,"%s:%d: Error: Element %d, ",__FILE__,__LINE__,nsfin);
        fprintf(stdout,"sommet %d --> vol %d: no-win situation\n",ind,i2);
        return(0);
      }
    }

    /* we intersect the tri key%11 */
    trav->kel    = nsfin;
    vol[is]  = dd-vol[i0]-vol[i1]-vol[i2];

    dd = vol[i0]+vol[i1]+vol[i2]+vol[is];
    if ( dd == 0.0 ) {
      fprintf(stdout,"%s:%d: Error: tetrahedra volume is null\n",
             __FILE__,__LINE__);
      return(0);
    }

    if ( -EPS2 > vol[is] ) {
      (*lastet) = 0;

      /* barycentric */
      is = inxt3[is];
      i0 = idir[is][0];
      i1 = idir[is][1];
      i2 = idir[is][2];

      barycentric(is,i0,i1,i2,p,e0,mat,ep);

      ep[0] = e0[0]-p[is]->c[0];
      ep[1] = e0[1]-p[is]->c[1];
      ep[2] = e0[2]-p[is]->c[2];

      switch (ind) {
      case(0):
        u[0]  = mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2];
        u[1]  = mat[2][1]*mat[0][2] - mat[0][1]*mat[2][2];
        u[2]  = mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2];
        if ( !intersegtri_i0(is, i0, i1, i2, ee, e0, ep, dd, p[is]->c,
                             mat, u, c, cb, &key, trav,
                             &nsfin, adj[ind]/4) ) {
          fprintf(stdout,"%s:%d: Error: We didn't intersect the tri %d\n",
                 __FILE__,__LINE__,ind);
          return(0);
        }
        break;
      case(1):
      case(2):
        if ( !intersegtri_i1(is, i0, i1, i2, ee, e0, ep, dd, p[is]->c,
                             mat, c, cb, &key, trav,
                             &nsfin, adj[ind]/4) ) {
          fprintf(stdout,"%s:%d: Error: We didn't intersect the tri %d\n",
                 __FILE__,__LINE__,ind);
          return(0);
        }
        break;

      case(3):
        if ( !intersegtri_i2(is, i0, i1, i2, ee, e0, ep, dd, p[is]->c,
                             mat, c, cb, &key, trav,
                             &nsfin, adj[ind]/4) ) {
          fprintf(stdout,"%s:%d: Error: We didn't intersect the tri %d\n",
                 __FILE__,__LINE__,ind);
          return(0);
        }
        break;
      default:
        fprintf(stdout,"%s:%d: Error: unexpected case: %d \n",
               __FILE__,__LINE__,ind);
        return(0);
      }

      trav->np = pt->v[ind];
      if ( !creaEdge(mesh,met,trav,trav->tag) ) {
        fprintf(stdout,"%s:%d: Error: not able to create edge\n",__FILE__,__LINE__);
        return(0);
      }

      if ( nsfin && (mesh->tetra[nsfin].flag != basetet) )  return(1);

      fprintf(stdout,"%s:%d: Warning: we can't travel through the wanted tri %d\n",
             __FILE__,__LINE__,ind);

      if ( key/11 ) {
        ind = key%11;
        i0  = inxt3[ind];
        i1  = inxt3[i0];
        i2  = inxt3[i1];

        if ( adj[i0] && (mesh->tetra[adj[i0]/4].flag != basetet) ) {
          nsfin = adj[i0]/4;

          return(1);
        }
        if ( adj[i1] && (mesh->tetra[adj[i1]/4].flag != basetet) ) {
          nsfin = adj[i1]/4;

          return(1);
        }
        if ( adj[i2] && (mesh->tetra[adj[i2]/4].flag != basetet) ) {
          nsfin = adj[i2]/4;

          return(1);
        }
      }
      else if ( key/5 ) {
        ind = key-5;
        i0  = ifar[ind][0];
        i1  = ifar[ind][1];
        if ( adj[i0] && (mesh->tetra[adj[i0]/4].flag != basetet) ) {
          nsfin = adj[i0]/4;

          return(1);
        }
        if ( adj[i1] && (mesh->tetra[adj[i1]/4].flag != basetet) ) {
          nsfin = adj[i1]/4;

          return(1);
        }
      }
      fprintf(stdout,"%s:%d: Error: no-win situation\n",__FILE__,__LINE__);
      return(0);
    }

    /* the point is inside the tet */
    dd = 1.0 / dd;

    cb[is]      = vol[is] * dd;
    cb[i0]      = vol[i0] * dd;
    cb[i1]      = vol[i1] * dd;
    cb[i2]      = vol[i2] * dd;
    key         = locate(cb);
    trav->cb[0] = cb[0];
    trav->cb[1] = cb[1];
    trav->cb[2] = cb[2];
    trav->cb[3] = cb[3];
    trav->c[0]  = e1[0];
    trav->c[1]  = e1[1];
    trav->c[2]  = e1[2];
    trav->np    = pt->v[ind];
    trav->key   = key;
    ppt1->tet   = nsfin;
    ppt1->flag  = basept;
    (*lastet)   = 1;

    if ( !creaEdge(mesh,met,trav,trav->tag) ) {
      fprintf(stdout,"%s:%d: Error: not able to create edge\n",
              __FILE__,__LINE__);
      return(0);
    }

    return(1);

  } while( ++it<=mesh->ne );
  fprintf(stdout,"%s:%d: Error: we reach the end",__FILE__,__LINE__);
  fprintf(stdout," of the function, strange... %d %d\n",it,mesh->ne);
  return(0);
}

/** put singularities stored in singul in the mesh */
int inserSingul(pMesh mesh, pSol met, pSingul singul){
  int     lastet, k, k0, k1;
  double  cb[4];
  Travel  trav;

  if ( (!singul->na) && (!singul->ns) ) return(-1);

  if ( abs(info.imprim) > 4 || info.ddebug )
    fprintf(stdout,"\n  ** INSERTION OF SINGULARITIES\n");

  if ( !hashTetra(mesh) ) {
    fprintf(stdout,"  ## Hashing problem (1). Exit program.\n");
    return(0);
  }

  /* Alloc of htab for singular edges */
  /* if edges exist in mesh, hash special edges from existing field */
  if ( mesh->na || singul->na ) {
    mesh->namax = MG_MAX(1.5*(mesh->na+singul->na),NAMAX);
    hNew(&mesh->htab,mesh->na+singul->na,3*mesh->namax);
  }

  /* put edge singularities in mesh */
  for ( k=1; k<=singul->na; k++) {
    lastet = 1;
    if ( singul->point[singul->edge[k].a].tet ||
         !(singul->point[singul->edge[k].b].tet) ) {
      k0 = singul->edge[k].a;
      k1 = singul->edge[k].b;
    }
    else {
      k0 = singul->edge[k].a;
      k1 = singul->edge[k].b;
    }
    assert( k0 && k1 );

    trav.tag = singul->edge[k].tag;
    do {
      if ( !seekEdge(mesh,met,&singul->point[k0],&singul->point[k1],
                     &trav,&lastet) ) {
        fprintf(stdout,"%s:%d: Error: edge %d not found in the mesh\n",
               __FILE__,__LINE__,k);
        return(0);
      }
    } while( !lastet );
  }

  /* put point singularities in mesh */
  for ( k=1; k<=singul->ns; k++) {
    if ( singul->point[k].flag != basept ) {
      if ( !seekPoint(mesh,&singul->point[k],cb) ) {
        fprintf(stdout,"%s:%d: Error: point %d not found in the mesh\n",
                __FILE__,__LINE__,k);
        return(0);
      }
      if ( !creaPoint(mesh,met,singul->point[k].tet,singul->point[k].c,
                      cb,singul->point[k].tag) ) {
        fprintf(stdout,"%s:%d: Error: unable to create point %d\n",
                __FILE__,__LINE__,k);
        return(0);
      }
    }
  }

  /* clean markers */
  for ( k=1; k<=mesh->np; k++) {
    mesh->point[k].flag = 0;
  }
  for ( k=1; k<=mesh->ne; k++) {
    mesh->tetra[k].flag = 0;
  }

  return(1);
}

/** collapse of singularities */
int colSing(pMesh mesh,pSol met) {
  pTetra  pt;
  pxTetra pxt;
  pPoint  p0,p1;
  int     k,nc,nnc,list[LMAX+2],ilist;
  int     it,maxit,i,ifac,jseg,ier;

  if ( abs(info.imprim) > 4 )
    fprintf(stdout,"  ** SINGULARITIES PRE-REMESHING\n");
  nnc = it = 0;
  maxit = 5;
  do {
    nc = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( (!MG_EOK(pt)) || (!pt->xt) )  continue;
      pxt = &mesh->xtetra[pt->xt];

      /* Attempt to collapse */
      for (i=0; i<6; i++) {
        if ( (!(pxt->tag[i] & MG_GEO)) || (pxt->tag[i] & MG_BDY) )  continue;
        p0 = &mesh->point[pt->v[iare[i][0]]];
        p1 = &mesh->point[pt->v[iare[i][1]]];

        if ( !MG_SIN(p0->tag) ) {
          ifac  = isar[i][0];
          jseg  = iarfinv[ifac][i];
          ilist = chkcol_int(mesh,met,k,ifac,jseg,list,2);
          if ( ilist ) {
            ier = colver(mesh,list,ilist,iare[i][1]);
            if ( ier < 0 ) return(-1);
            else if ( ier ) {
              delPt(mesh,ier);
              nc++;
              break;
            }
          }
        }
        if ( !MG_SIN(p1->tag) ) {
          ifac  = isar[i][1];
          jseg  = iarfinv[ifac][i];
          ilist = chkcol_int(mesh,met,k,ifac,jseg,list,2);
          if ( ilist ) {
            ier = colver(mesh,list,ilist,iare[i][0]);
            if ( ier < 0 ) return(-1);
            else if ( ier ) {
              delPt(mesh,ier);
              nc++;
              break;
            }
          }
        }
      }
    }
    nnc += nc;
    if ( nc > 0 && (abs(info.imprim) > 5 || info.ddebug) )
      fprintf(stdout,"     %8d vertices removed\n",nc);
  }
  while ( ++it < maxit && nc > 0 );
  return(1);
}

/* Check that we don't have a tetra with 4 singular points otherwise *
 * try to swap or split the tetra. *
 * ( warning: here we don't perform quality tests) */
int solveUnsignedTet(pMesh mesh,pSol met) {
  pTetra  pt;
  int     k,nf,ns;
  int     *adja;
  int     i,ip;

  if ( abs(info.imprim) > 4 )
    fprintf(stdout,"  ** SINGULARITIES POST-REMESHING\n");

  ns = nf = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( (!MG_EOK(pt)) || (!pt->xt) )  continue;

    if ( (mesh->point[pt->v[0]].tag & MG_SGL) &&
         (mesh->point[pt->v[1]].tag & MG_SGL) &&
         (mesh->point[pt->v[2]].tag & MG_SGL) &&
         (mesh->point[pt->v[3]].tag & MG_SGL) ) {

      /* First: attempt to swap (swap 2->3) */
      adja =  &mesh->adja[4*(k-1)+1];
      for (i=0; i<4; i++) {
        ip     = mesh->tetra[adja[i]/4].v[adja[i]%4];
        if ( !(mesh->point[ip].tag & MG_SGL ) && swap23(mesh,k,i) ) {
          nf++;
          break;
        }
      }

      /* Second: split on barycenter of tetra */
      if ( i == 4 ) {
        /* Swapping is useless so we add a new degree of freedom */
        if ( !split4bar(mesh,met,k) ) {
          fprintf(stdout,"%s:%d: Error: we can't split element %d",
                  __FILE__,__LINE__,k);
          fprintf(stdout," whose all vertices are on inserted singularities\n");
          return(0);
        }
        ns++;
        continue;
      }
    }
  }
  if ( abs(info.imprim) < 5 && (ns+nf) > 0 ) {
    fprintf(stdout,"     %8d corrected ",ns+nf);
    fprintf(stdout,"(%8d splitted, %8d swapped)\n",ns,nf);
  }
  return(1);
}

#endif
