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
#include "mmg2d.h"

#define BUCKSIZ    64
#define M_LONG     1.4//1.85//1.4//1.421
#define M_SHORT    0.65//0.8//0.65//0.707
int ddebug = 0;

int MMG2_invmat(double *m,double *minv) {
  double        det;

  if(fabs(m[1]) < EPSD) { /*mat diago*/
    minv[0] = 1./m[0];
    minv[1] = 0;
    minv[2] = 1./m[2];
    if(ddebug) printf("Idd : %e %e %e \n",m[0]*minv[0]+m[1]*minv[1],
                      m[1]*minv[0]+m[2]*minv[1],
                      m[1]*minv[1]+m[2]*minv[2]);
  } else {
    det = m[0]*m[2] - m[1]*m[1];
    det = 1. / det;
    minv[0] = det * m[2];
    minv[1] = - det * m[1];
    minv[2] = det * m[0];
    if(ddebug) printf("Id : %e %e %e -- %e\n",m[0]*minv[0]+m[1]*minv[1],
                      m[1]*minv[0]+m[2]*minv[1],
                      m[1]*minv[1]+m[2]*minv[2],det);
  }
  return(1);
}

int interp_ani(double *ma,double *mb,double *mp,double t) {
  double  dma[3],dmb[3],mai[3],mbi[3],mi[3];
  int   i;

  for (i=0; i<3; i++) {
    dma[i] = ma[i];
    dmb[i] = mb[i];
  }
  if(ddebug) {
    printf("met1 : %e %e %e\n",ma[0],ma[1],ma[2]);
    printf("met2 : %e %e %e\n",mb[0],mb[1],mb[2]);
  }
  if ( !MMG2_invmat(dma,mai) || !MMG2_invmat(dmb,mbi) ) {
    fprintf(stderr,"  ## INTERP INVALID METRIC.\n");
    return(0);
  }
  if(ddebug) {
    printf("invmet1 : %e %e %e\n",mai[0],mai[1],mai[2]);
    printf("invmet2 : %e %e %e\n",mbi[0],mbi[1],mbi[2]);
  }

  for (i=0; i<3; i++)
    mi[i] = (1.0-t)*mai[i] + t*mbi[i];
  if(ddebug) {
    printf("meti : %e %e %e\n",mi[0],mi[1],mi[2]);
  }

  if ( !MMG2_invmat(mi,mai) ) {
    fprintf(stderr,"  ## INTERP INVALID METRIC.\n");
    return(0);
  }
  if(ddebug) {
    printf("met : %e %e %e\n",mai[0],mai[1],mai[2]);
  }
  for (i=0; i<3; i++)  mp[i] = mai[i];
  return(1);
}

int interp_iso(double *ma,double *mb,double *mp,double t) {

  *mp = (1.0-t)*(*ma) + t*(*mb);
  return(1);
}

static int cassar(MMG5_pMesh mesh,MMG5_pSol sol,int ia,int ib,double t) {
  MMG5_pPoint   p1,p2;
  //Displ      pd;
  double   c[2],t1,*ma,*mb,*mp;
  int      ip,iadr,memlack;

  memlack = 0;

  p1 = &mesh->point[ia];
  p2 = &mesh->point[ib];
  t1 = 1.0 - t;

  c[0] = t1*p1->c[0] +  t*p2->c[0];
  c[1] = t1*p1->c[1] +  t*p2->c[1];
  ip   = _MMG2D_newPt(mesh,c,0);
  if ( !ip ) {
    /* reallocation of point table */

    _MMG2D_POINT_REALLOC(mesh,sol,ip,mesh->gap,
                         printf("  ## Error: unable to allocate a new point\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         memlack=1;
                         return(-1)
                         ,c,0);
    p1  = &mesh->point[ia];
    p2  = &mesh->point[ib];
  }


  /*interpol metric*/
  iadr = (ia-1)*sol->size + 1;
  ma  = &sol->m[iadr];

  iadr = (ib-1)*sol->size + 1;
  mb  = &sol->m[iadr];

  iadr = (ip-1)*sol->size + 1;
  mp  = &sol->m[iadr];

  if ( sol->size==1 ) {
    if(!interp_iso(ma,mb,mp,t1) ) return(-1);
  }
  else {
    if(!interp_ani(ma,mb,mp,t1) ) return(-1);
  }

  /*interpol dep si option 9*/
  if( mesh->info.lag >=0) {
    printf("comment because of merge needs mmg2d1 option 9\n");
    // pd = mesh->disp;
    //pd.mv[2*(ip-1) + 1 + 0] = t1*pd.mv[2*(ia-1) + 1 + 0] + t*pd.mv[2*(ib-1) + 1 + 0];
    //pd.mv[2*(ip-1) + 1 + 1] = t1*pd.mv[2*(ia-1) + 1 + 1] + t*pd.mv[2*(ib-1) + 1 + 1];
  }
  if ( memlack )  return(-1);
  return(ip);
}

/*compute new vertex on bdry*/
static int cassarbdry(MMG5_pMesh mesh,MMG5_pSol sol,int ied,int ia,int ib,double t,double *tang) {
  MMG5_pPoint   p0,p1,ppt;
  MMG5_pEdge    ped;
  // Displ      pd;
  double   c[2],pc1[2],pc2[2],t0[2],t1[2],t_1,*ma,*mb,*mp;//,dx,dy;
  double   l;
  int      ip,iadr,i,inv,memlack;
  inv = 0;
  memlack = 0;
  p0 = &mesh->point[ia];
  p1 = &mesh->point[ib];

  l = (p0->c[0] - p1->c[0])*(p0->c[0] - p1->c[0]) +
    (p0->c[1] - p1->c[1])*(p0->c[1] - p1->c[1]);
  l = sqrt(l);
  ped = &mesh->edge[ied];
  if(ia != ped->a) {
    assert(ib == ped->a);
    for(i=0 ; i<2 ; i++) {
      t0[i] = l*p0->n[i];
      t1[i] = l*p1->n[i];
    }
    /*if corner, recompute the tangent*/
    if(p0->tag & M_CORNER) {
      for(i=0 ; i<2 ; i++)
        t0[i] = p0->c[i] - p1->c[i];
    }
    if(p1->tag & M_CORNER) {
      for(i=0 ; i<2 ; i++)
        t1[i] = p0->c[i] - p1->c[i];
    }
  } else {
    assert(ib == ped->b);
    for(i=0 ; i<2 ; i++) {
      t0[i] = l*p0->n[i];
      t1[i] = l*p1->n[i];
    }
    /*if corner, recompute the tangent*/
    if(p0->tag & M_CORNER) {
      for(i=0 ; i<2 ; i++)
        t0[i] = p1->c[i] - p0->c[i];
    }
    if(p1->tag & M_CORNER) {
      for(i=0 ; i<2 ; i++)
        t1[i] = p1->c[i] - p0->c[i];
    }

  }


  /*check if t0 has the same sens of vect(P0P1)*/
  if(t0[0]/(p1->c[0]-p0->c[0]) < 0 || t0[1]/(p1->c[1]-p0->c[1])<0) {
    //printf("t0/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
    for(i=0 ; i<2 ; i++) {
      t0[i] *= -1;
    }
    inv = 1;
  }
  /*check if t1 has the opposite sens of vect(P0P1)*/
  if(t1[0]/(p1->c[0]-p0->c[0]) > 0 || t1[1]/(p1->c[1]-p0->c[1])>0) {
    //printf("t1/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
    for(i=0 ; i<2 ; i++) {
      t1[i] *= -1;
    }
  }
  /*control points*/
  for(i=0 ; i<2 ; i++) {
    pc1[i] = (t0[i]+3*p0->c[i])/3.;
    pc2[i] = (t1[i]+3*p1->c[i])/3.;
  }

  /*coor new point*/
  t_1 = 1.0 - t;
  for(i=0 ; i<2 ; i++) {
    c[i] = t_1*t_1*t_1*p0->c[i] + 3*t*t_1*t_1*pc1[i] + 3*t*t*t_1*pc2[i] +  t*t*t*p1->c[i];
  }
  // printf("c %e %e -- mid %e %e\n",c[0],c[1],0.5*(p0->c[0]+p1->c[0]),0.5*(p0->c[1]+p1->c[1]));
  ip   = _MMG2D_newPt(mesh,c,0);
  if ( !ip ) {
    /* reallocation of point table */

    _MMG2D_POINT_REALLOC(mesh,sol,ip,mesh->gap,
                         printf("  ## Error: unable to allocate a new bdry point\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         memlack=1;
                         return(-1);
                         ,c,0);
    p0 = &mesh->point[ia];
    p1 = &mesh->point[ib];
  }

  /*tangent new point*/
  ppt = &mesh->point[ip];
  for(i=0 ; i<2 ; i++) {
    // tang[i] = -(-3*t_1*t_1*p0->c[i] +(3*t_1*(1-3*t))*pc1[i]
    //    + (3*t*(2-3*t))*pc2[i] +  3*t*t*p1->c[i]);
    // printf("tang %e %e diff %e\n",tang[i],(3./8.)*(p1->c[i]+pc2[i]-pc1[i]-p0->c[i]),
    //     fabs(tang[i]-(3./8.)*(p1->c[i]+pc2[i]-pc1[i]-p0->c[i])));
    tang[i] = (3./8.)*(p1->c[i]+pc2[i]-pc1[i]-p0->c[i]);
  }
  l = tang[0]*tang[0] + tang[1]*tang[1];
  l = 1./sqrt(l);
  tang[0] *= l;
  tang[1] *= l;
  /*check if tang has the same sens than P0P*/
  if(inv) {
    if ((tang[0]/(ppt->c[0]-p0->c[0]) > 0 || tang[1]/(ppt->c[1]-p0->c[1])>0)) {
      //printf("tang/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
      for(i=0 ; i<2 ; i++) {
        tang[i] *= -1;
      }
    }
  } else if((tang[0]/(ppt->c[0]-p0->c[0]) < 0 || tang[1]/(ppt->c[1]-p0->c[1])<0)) {
    //printf("tang/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
    for(i=0 ; i<2 ; i++) {
      tang[i] *= -1;
    }
  }

/*   /\*change the two other tangents*\/ */
/*   if(ia != ped->a) { */
/*     assert(ib == ped->a); */
/*     for(i=0 ; i<2 ; i++) { */
/*       /\*t0*\/ */
/*       p0->n[i] = -(3./2.)*(pc1[i]-p0->c[i]); */
/*       /\*t1*\/ */
/*       p1->n[i] = -(3./2.)*(p1->c[i]-pc2[i]); */
/*     } */
/*     /\*check tangent orientation*\/ */
/* #warning remove ? check orientation tangent */
/*     dx = t1[0]/p0->n[0]; */
/*     dy = t1[1]/p0->n[1]; */
/*     if(dy < 0 || dx <0) { */
/*       //printf("1) pbs de colinearite %e\n",fabs(dx-dy)); */
/*     } */
/*     dx = t0[0]/p1->n[0]; */
/*     dy = t0[1]/p1->n[1]; */
/*     if(dx < 0 || dy <0) { */
/*       //printf("3) pbs de colinearite %e %e %e\n",fabs(dx-dy),dx,dy); */
/*     } */
/*   } else { */
/*     assert(ib == ped->b); */
/*     for(i=0 ; i<2 ; i++) { */
/*       /\*t0*\/ */
/*       p0->n[i] = -(3./2.)*(pc1[i]-p0->c[i]); */
/*       /\*t1*\/ */
/*       p1->n[i] = -(3./2.)*(p1->c[i]-pc2[i]); */
/*     } */
/*     /\*check tangent orientation*\/ */
/* #warning remove ? check orientation tangent */
/*     dx = t1[0]/p1->n[0]; */
/*     dy = t1[1]/p1->n[1]; */
/*     if(dx < 0 || dy <0) { */
/*       //printf("2) pbs de colinearite %e %e %e\n",fabs(dx-dy),dx,dy); */
/*     } */
/*     dx = t0[0]/p0->n[0]; */
/*     dy = t0[1]/p0->n[1]; */
/*     if(dx < 0 || dy <0) { */
/*       //printf("4) pbs de colinearite pts %d %e %e %e\n",ia,fabs(dx-dy),dx,dy); */
/*       for(i=0 ; i<2 ; i++) { */
/*         p0->n[i] *= -1; */
/*       } */

/*     } */
/*   } */

  /*interpol metric*/
  iadr = (ia-1)*sol->size + 1;
  ma  = &sol->m[iadr];

  iadr = (ib-1)*sol->size + 1;
  mb  = &sol->m[iadr];

  iadr = (ip-1)*sol->size + 1;
  mp  = &sol->m[iadr];

  if ( sol->size==1 ) {
    if(!interp_iso(ma,mb,mp,t_1) ) return(-1);
  }
  else {
    if(!interp_ani(ma,mb,mp,t_1) ) return(-1);
  }

  /*interpol dep si option 9*/
  if( mesh->info.lag >= 0) {
//#warning option 9
    printf("option 9..........\n");
    exit(EXIT_FAILURE);
    /* pd = mesh->disp; */
    /* pd.mv[2*(ip-1) + 1 + 0] = t_1*pd.mv[2*(ia-1) + 1 + 0] + t*pd.mv[2*(ib-1) + 1 + 0]; */
    /* pd.mv[2*(ip-1) + 1 + 1] = t_1*pd.mv[2*(ia-1) + 1 + 1] + t*pd.mv[2*(ib-1) + 1 + 1]; */
  }
  if ( memlack )  return(-1);
  return(ip);
}
int MMG2_ni,MMG2_nc;
static int analar(MMG5_pMesh mesh,MMG5_pSol sol,pBucket bucket,double declic,int *alert) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppa,ppb;
  double  *ca,*cb,*ma,*mb,tail,t,tang[2];
  int     *adja,voi[3],k,iadr,adj,/*base,*/nbp,npp,ip;
  int     nt;
  int     i,i1,i2;
  int ins,i0,ii0;
//  base  = ++mesh->base;
  MMG2_ni  = 0;
  MMG2_nc  = 0;
  nt    = mesh->nt;
  npp = 0.0;

  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;
    //else if ( /*pt->flag = base-1 || pt->qual < declic*/ )  continue;

    /* base internal edges */
    iadr  = 3*(k-1) + 1;
    adja  = &mesh->adja[iadr];
    voi[0] = adja[0];
    voi[1] = adja[1];
    voi[2] = adja[2];
    i0 = 0;
    if (!voi[1]) i0 = 1;
    if (!voi[2]) i0 = 2;
    for (ii0=i0; ii0<i0+3; ii0++) {
      i = ii0%3;
      adj = voi[i] / 3;

      i1   = pt->v[MMG2_idir[i+1]];
      i2   = pt->v[MMG2_idir[i+2]];

      ppa  = &mesh->point[i1];
      ppb  = &mesh->point[i2];
      //#warning bad test for edge required
      if((ppa->tag & M_REQUIRED) && (ppb->tag & M_REQUIRED)) {
        //printf("edge required %d %d\n",i1,i2);
        continue;
      }
      ca   = &ppa->c[0];
      cb   = &ppb->c[0];
      iadr = (i1-1)*sol->size + 1;
      ma   = &sol->m[iadr];
      iadr = (i2-1)*sol->size + 1;
      mb   = &sol->m[iadr];
      tail = MMG2_length(ca,cb,ma,mb);
      if(ddebug) printf("edge %d -- %d %d: %f tr : %d %d %d\n",
                        i,i1,i2,tail,pt->v[0],pt->v[1],pt->v[2]);
      if ( tail > M_LONG && *alert <= 1 ) {
        if(ddebug) {
          printf("too long %d\n",adj);
        }
        npp++;
        nbp = tail + 0.5;
        if ( nbp*(nbp+1) < 0.99*tail*tail )  nbp++;
        t = 1.0 / (float)nbp;
        if ( nbp < 3 || nbp > 15 )  t = 0.5;
        if ( !adj || pt->ref != mesh->tria[adj].ref )  {
          /*add bdry*/
          if(!pt->edg[i])  {
            if(mesh->info.ddebug) {
              printf("tr %d : %d %d %d mais %d\n",k,pt->edg[0],pt->edg[1],pt->edg[2],i);
              printf("%d %d %d\n",pt->v[0],pt->v[1],pt->v[2]);
            }
            assert(mesh->tria[adj].ref!=pt->ref);
            assert((mesh->tria[adj]).edg[voi[i]%3]);
//#warning find why we have to do that
            pt->edg[i] = (mesh->tria[adj]).edg[voi[i]%3];
          }
          assert(pt->edg[i]);
          ip = cassarbdry(mesh,sol,pt->edg[i],i1,i2,0.5,tang);
        } else {
          ip = cassar(mesh,sol,i1,i2,t);
        }
        if(ddebug) printf("try cut %d %f \n",ip,t);
        if(ddebug) {
          printf("on essaie\n");
          MMG2D_saveMesh(mesh,"in.mesh");
        }
        if(ip < 0) {
          printf("IMPOSSIBLE TO CREATE NEW VERTEX\n");
          //return(0);
          //printf("ahhhhhhhhhhhhhhhh\n");
          *alert = 2;
        } else {
          if ( !adj || pt->ref != mesh->tria[adj].ref )  {
            /*boundary edge*/
            if(!adj) {
              ins = MMG2_splitbdry(mesh,sol,ip,k,i,tang);
              if(!ins) {
                if (ddebug) printf("k= %d on insere pas bdry : %d %d\n",k,i1,i2);
                //mesh->point[ip].tag = M_CORNER;
                //mesh->point[ip].ref = 5;
                //  MMG2D_saveMesh(mesh,"del.mesh");
                //exit(EXIT_FAILURE);
                _MMG2D_delPt(mesh,ip);
                continue;
              }
              mesh->point[ip].tag |= M_BDRY;
              MMG2_ni += 1;
              break;
            } else {
              mesh->point[ip].tag |= M_SD;
              ins = MMG2_split(mesh,sol,ip,k,voi[i]);
              if(!ins) {
                if (ddebug) printf("on insere pas sd\n");
                _MMG2D_delPt(mesh,ip);
                continue;
              }
              MMG2_ni += 1;
              break;
            }
            continue;
          } else {
            if(ddebug) {
              printf("on insere la\n");
              MMG2D_saveMesh(mesh,"cut.o.mesh");
            }
            ins = MMG2_split(mesh,sol,ip,k,voi[i]);
            if(ddebug) printf("cut ? %d\n",ins);
            if(!ins) {
              if (ddebug) printf("on insere pas :::\n");
              _MMG2D_delPt(mesh,ip);
              continue;
            }
            MMG2_ni += 1;
            break;
          }
        }
      }

      else if ( tail < M_SHORT ) {
        if(ddebug) {
          printf("too short %d \n",adj);
          MMG2D_saveMesh(mesh,"ttttt.mesh");
        }
        if ( !adj || pt->ref != mesh->tria[adj].ref )  {
          if(!adj) {
            if(ddebug) printf("edge %d -- %d %d\n",pt->edg[i],mesh->edge[pt->edg[i]].a,mesh->edge[pt->edg[i]].b);
            if (!MMG2_colpoibdry(mesh,sol,k,i,MMG2_iare[i][0],MMG2_iare[i][1],2.75)){
              if (!MMG2_colpoibdry(mesh,sol,k,i,MMG2_iare[i][1],MMG2_iare[i][0],2.75)){
                continue;
              } else {
                MMG2_nc++;
                _MMG2D_delPt(mesh,i1);
                break;
              }
            }
            MMG2_nc++;
            _MMG2D_delPt(mesh,i2);
            break;
          } else {
            if(!MMG2_colpoi(mesh,sol,k,i,MMG2_iare[i][0],MMG2_iare[i][1],2.75)) {
              if(!MMG2_colpoi(mesh,sol,k,i,MMG2_iare[i][1],MMG2_iare[i][0],2.75)) {
                continue;
              } else {
                MMG2_nc++;
                _MMG2D_delPt(mesh,i1);
                break;
              }
            }
            _MMG2D_delPt(mesh,i2);
            MMG2_nc++;
            break;
          }
        } else {
          if(!MMG2_colpoi(mesh,sol,k,i,MMG2_iare[i][0],MMG2_iare[i][1],2.75)) {
            if(!MMG2_colpoi(mesh,sol,k,i,MMG2_iare[i][1],MMG2_iare[i][0],2.75)) {
              continue;
            } else {
              MMG2_nc++;
              _MMG2D_delPt(mesh,i1);
              if(ddebug) {
                printf("del ok\n");
                //MMG2D_saveMesh(mesh,"del.mesh");
                //exit(EXIT_FAILURE);
              }
              break;

            }
          }
          MMG2_nc++;
          _MMG2D_delPt(mesh,i2);
          if(ddebug) {
            printf("del2 ok\n");
            //MMG2D_saveMesh(mesh,"del.mesh");
            //exit(EXIT_FAILURE);
          }

          break;
        }
      }
    }
  }
  if ( mesh->info.imprim > 5 ) {
    fprintf(stdout,"    %8d INSERTED %8d COLLAPSED\n",MMG2_ni,MMG2_nc);
  }
  return(1);
}
static int analargeom(MMG5_pMesh mesh,MMG5_pSol sol,int *alert) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppa,ppb;
  double  tail,tang[2];
  int     *adja,voi[3],k,iadr,adj,/* base ,*/npp,ip;
  int     nt;
  int     i,i1,i2,ni,maxtou,it;
  int     ins,i0,ii0,nitot;
  maxtou = 30;
  it     = 0;
  nitot = 0;
  do {
    ni  = 0;
    // base  = ++mesh->base;
    nt    = mesh->nt;
    npp   = 0;
    for (k=1; k<=nt; k++) {
      pt = &mesh->tria[k];
      if ( !M_EOK(pt) )  continue;
      //else if ( /*pt->flag = base-1 || pt->qual < declic*/ )  continue;

      /* base internal edges */
      iadr  = 3*(k-1) + 1;
      adja  = &mesh->adja[iadr];
      voi[0] = adja[0];
      voi[1] = adja[1];
      voi[2] = adja[2];
      i0 = 0;
      if (!voi[1]) i0 = 1;
      if (!voi[2]) i0 = 2;
      for (ii0=i0; ii0<i0+3; ii0++) {
        i = ii0%3;
        adj = voi[i] / 3;

        i1   = pt->v[MMG2_idir[i+1]];
        i2   = pt->v[MMG2_idir[i+2]];

        ppa  = &mesh->point[i1];
        ppb  = &mesh->point[i2];
        //#warning bad test for edge required
        if((ppa->tag & M_REQUIRED) && (ppb->tag & M_REQUIRED)) {
          //printf("edge required %d %d\n",i1,i2);
          continue;
        }
        if ( adj )  continue;

        tail = MMG2_chkedg(mesh,ppa,ppb);
        if(!tail) continue;
        if ( *alert <= 1 ) {
          npp++;
          if(!pt->edg[i])  {
            printf("tr %d : %d %d %d mais %d\n",k,pt->edg[0],pt->edg[1],pt->edg[2],i);
            printf("%d %d %d\n",pt->v[0],pt->v[1],pt->v[2]);
            printf("adj %d %d %d\n",voi[0]/3,voi[1]/3,voi[2]/3);
            MMG2D_saveMesh(mesh,"titii.mesh");
          }
          assert(pt->edg[i]);
          ip = cassarbdry(mesh,sol,pt->edg[i],i1,i2,0.5,tang);

          if(ip < 0) {
            printf("IMPOSSIBLE TO CREATE NEW VERTEX\n");
            *alert = 2;
          }
          ins = MMG2_splitbdry(mesh,sol,ip,k,i,tang);
          if(!ins) {
            if (ddebug) printf("k= %d on insere pas bdry : %d %d\n",k,i1,i2);
            _MMG2D_delPt(mesh,ip);
            continue;
          }
          mesh->point[ip].tag |= M_BDRY;
          ni += 1;
          break;
        }
      }
    }
    if ( mesh->info.imprim > 5 ) {
      fprintf(stdout,"    %8d INSERTED \n",ni);
    }
    nitot +=ni;
  } while(ni > 0 && it++ < maxtou);
  return(nitot);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return 1 if success.
 *
 * Mesh adaptation.
 *
 **/
int MMG2_mmg2d1(MMG5_pMesh mesh,MMG5_pSol sol) {
  pBucket  bucket;
  double   declic;
  int      ns,base,alert,it,maxtou;
  int      nadd,ndel,nswp,ngeom;

  nadd = ndel = nswp = 0;
  if ( mesh->info.imprim < 0 ) {
    MMG2_outqua(mesh,sol);
    MMG2_prilen(mesh,sol);
  }

  /* 1. Delaunization */
  if ( mesh->info.imprim < -4 )
    fprintf(stdout,"  -- DELAUNIZATION\n");

  /* edge flip */
  if(!mesh->info.noswap) {
    declic = 1.1 / ALPHA;
    base   = mesh->base;
    ns = MMG2_cendel(mesh,sol,declic,base);
    nswp += ns;
    if ( mesh->info.imprim > 5 )
      fprintf(stdout,"  -- %8d SWAPPED\n",ns);
  }
  alert  = 0;

  /* 1. Geometric mesh */
  if ( mesh->info.imprim > 3 )
    fprintf(stdout,"  -- GEOMETRIC MESH\n");
  ngeom = analargeom(mesh,sol,&alert);
  if ( mesh->info.imprim && (abs(mesh->info.imprim) < 6) )
    fprintf(stdout,"     %8d splitted\n",ngeom);

  /* 2. field points */
  bucket = MMG2_newBucket(mesh,M_MAX(mesh->info.bucket,BUCKSIZ));
  assert(bucket);
  declic = 1.5 / ALPHA;
  maxtou = 30;
  it     = 0;
  do {

    analar(mesh,sol,bucket,declic,&alert);
    nadd += MMG2_ni;
    ndel += MMG2_nc;
    if(!mesh->info.noswap) {
      ns = MMG2_cendel(mesh,sol,declic,mesh->base);
      nswp += ns;
      if ( mesh->info.imprim > 5 )
        fprintf(stdout,"  -- %8d SWAPPED\n",ns);
    }
  }
  while ( ++it < maxtou && (MMG2_ni+MMG2_nc > 0));
  //while ( ++it < maxtou && (MMG2_ni+MMG2_nc > 0));//> 0.05*mesh->np));

  if ( mesh->info.imprim && (abs(mesh->info.imprim) < 6) && ( nadd || ndel ) ) {
    fprintf(stdout,"     %8d splitted, %8d collapsed,"
            " %8d swapped.\n",nadd,ndel,nswp);
  }

  if ( mesh->info.imprim < 0 ) {
    MMG2_outqua(mesh,sol);
    MMG2_prilen(mesh,sol);
  }

  /* free memory */
  _MMG5_SAFE_FREE(bucket->head);
  _MMG5_SAFE_FREE(bucket->link);
  _MMG5_SAFE_FREE(bucket);

  return(1);
}
