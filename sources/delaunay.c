#include "eigenv.h"
#include "mmg3d.h"

#define  EPSRAD       1.00005
//For Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2 test case:
//pbs with EPSCON=5e-4 and VOLMIN=1e-15 (MMG3D does not insert enough vertex...)
#define  EPSCON       1e-5//5.0e-4//1.e-4//1.0e-3
#define  VOLMIN       1e-15//1.e-10//1.0e-15  --> vol negatif qd on rejoue
#define LONMAX     4096

int MMG_cas;
extern int MMG_npuiss,MMG_nvol,MMG_npres;
#define KTA     7
#define KTB    11
#define KTC    13

/* hash mesh edge v[0],v[1] (face i of iel) */
int hashEdgeDelone(pMesh mesh,Hash *hash,int iel,int i,int *v) {
  int       *adja,iadr,jel,j,key,mins,maxs;
  hedge     *ha;

  /* compute key */
  if ( v[0] < v[1] ) {
    mins = v[0];
    maxs = v[1];
  }
  else {
    mins = v[1];
    maxs = v[0];
  }
  key = KTA*mins + KTB*maxs;
  key = key % hash->siz;
  ha  = &hash->item[key];

  if ( ha->a ) {
    /* identical face */
    if ( ha->a == mins && ha->b == maxs ) {
      iadr = (iel-1)*4 + 1;
      adja = &mesh->adja[iadr];
      adja[i] = ha->k;

      jel  = ha->k >> 2;
      j    = ha->k % 4;
      iadr = (jel-1)*4 + 1;
      adja = &mesh->adja[iadr];
      adja[j] = iel*4 + i;
      return(1);
    }
    else {
      while ( ha->nxt && ha->nxt < hash->max ) {
        ha = &hash->item[ha->nxt];
        if ( ha->a == mins && ha->b == maxs ) {
          iadr = (iel-1)*4 + 1;
          adja = &mesh->adja[iadr];
          adja[i] = ha->k;

          jel  = ha->k >> 2;
          j    = ha->k % 4;
          iadr = (jel-1)*4 + 1;
          adja = &mesh->adja[iadr];
          adja[j] = iel*4 + i;
          return(1);
        }
      }
    }
    ha->nxt   = hash->nxt;
    ha        = &hash->item[hash->nxt];
    ha->a     = mins;
    ha->b     = maxs;
    ha->k     = iel*4 + i;
    hash->nxt = ha->nxt;
    ha->nxt   = 0;

    if ( hash->nxt >= hash->max ) {
      TAB_RECALLOC(mesh,hash->item,hash->max,0.2,hedge,"face",return(0));
      for (j=hash->nxt; j<hash->max; j++)  hash->item[j].nxt = j+1;
    }
    return(1);
  }

  /* insert */
  ha->a = mins;
  ha->b = maxs;
  ha->k = iel*4 + i;
  ha->nxt = 0;

  return(1);
}
/* cavity -> ball */
int delone(pMesh mesh,pSol sol,int ip,int *list,int ilist) {
  pPoint    ppt;
  pTetra    pt,pt1;
  xTetra    xt;
  pxTetra   pxt0;
  int      *adja,*adjb,i,j,k,l,m,iel,jel,old,v[3],iadr,base,size;
  int       vois[4],iadrold;/*,ii,kk,iare1,iare2;*/
  short     i1;
  char      alert;
  int tref,isused,ixt;
  Hash hedg;
  //obsolete avec la realloc
  // if ( mesh->ne + 2*ilist > mesh->nemax )  {printf("on passe ici boum\n");return(0);}
  base = mesh->mark;
  /* external faces */
  size = 0;
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tetra[old];
    iadr = (old-1)*4 + 1;
    adja = &mesh->adja[iadr];
    vois[0]  = adja[0] >> 2;
    vois[1]  = adja[1] >> 2;
    vois[2]  = adja[2] >> 2;
    vois[3]  = adja[3] >> 2;
    for (i=0; i<4; i++) {
      jel = vois[i];
      if ( !jel || mesh->tetra[jel].mark != base ) {
        for (j=0; j<3; j++) {
          i1  = idir[i][j];
          ppt = &mesh->point[ pt1->v[i1] ];
          ppt->tagdel |= MG_NOM;
        }
        size++;
      }
    }
  }
  /* check isolated vertex */
  alert = 0;
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tetra[old];
    for (i=0; i<4; i++) {
      ppt = &mesh->point[ pt1->v[i] ];
      if ( !(ppt->tagdel & MG_NOM) )  alert = 1;
    }
  }
  /* reset tag */
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tetra[old];
    for (i=0; i<4; i++) {
      ppt = &mesh->point[ pt1->v[i] ];
      ppt->tagdel &= ~MG_NOM;
    }
  }
  if ( alert )  {return(-1);}
  /* hash table params */
  if ( size > 3*LONMAX )  return(0);
  if ( !hashNew(mesh,&hedg,size,3*size) ) { /*3*size suffit */
    fprintf(stdout,"  ## Unable to complete mesh.\n");
    return(-1);
  }
  for (k=0; k<ilist; k++) {
    old  = list[k];

    iadrold = (old-1)*4 + 1;
    adja = &mesh->adja[iadrold];
    vois[0]  = adja[0];
    vois[1]  = adja[1];
    vois[2]  = adja[2];
    vois[3]  = adja[3];

    pt   = &mesh->tetra[old];
    if(pt->xt) {
      pxt0 = &mesh->xtetra[pt->xt];
      memcpy(&xt,pxt0,sizeof(xTetra));
      isused=0;
      ixt = 1;
    } else {
      ixt = 0;
    }
    
    for (i=0; i<4; i++) {
      jel = vois[i] /4;
      j   = vois[i] % 4;

      /* external face */
      if ( !jel || (mesh->tetra[jel].mark != base) ) {
        iel = newElt(mesh);
	
	if ( !iel ) {

	  TETRA_REALLOC(mesh,iel,mesh->gap,
			printf("  ## Warning: unable to allocate a new element.\n");
			return(0));
	  pt   = &mesh->tetra[old];
	  if(pt->xt) 
	    pxt0 = &mesh->xtetra[pt->xt];
	  adja = &mesh->adja[iadrold];
	}        
	pt1 = &mesh->tetra[iel];
        memcpy(pt1,pt,sizeof(Tetra));
        pt1->v[i] = ip;
        pt1->qual = orcal(mesh,iel);
        pt1->ref = mesh->tetra[old].ref;
        if(pt1->qual < 1e-10) {printf("argggg (%d) %d : %e\n",ip,iel,pt1->qual);
          printf("pt1 : %d %d %d %d\n",pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3]);/*exit(0);*/}
        iadr = (iel-1)*4 + 1;
        adjb = &mesh->adja[iadr];
        adjb[i] = adja[i];
        
        if(ixt) {
          if( xt.ref[i] || xt.ftag[i]) {
            // if(!adja[i] || (mesh->tetra[jel].ref != pt1->ref)) {
            /*  if((old==20329) || (iel==20329)) { */
            /*  printf("eh eh on en a un!!!"); */
            /*  printf("adj of %d : %d %d %d %d\n",old,adja[0],adja[1],adja[2],adja[3]); */
            /*  printf("adj of %d : %d %d %d %d\n",iel,adjb[0],adjb[1],adjb[2],adjb[3]); */
            /* printf("on traite %d\n",i); */
            /*  } */
            if(!isused) {
              pt1->xt = pt->xt;
              pt->xt = 0;
              pxt0 = &mesh->xtetra[pt1->xt];
              memset(pxt0,0,sizeof(MMG5_xTetra));
              pxt0->ref[i]   = xt.ref[i] ; pxt0->ftag[i]  = xt.ftag[i];
              pxt0->edg[iarf[i][0]] = xt.edg[iarf[i][0]];
              pxt0->edg[iarf[i][1]] = xt.edg[iarf[i][1]];
              pxt0->edg[iarf[i][2]] = xt.edg[iarf[i][2]];
              pxt0->tag[iarf[i][0]] = xt.tag[iarf[i][0]];
              pxt0->tag[iarf[i][1]] = xt.tag[iarf[i][1]];
              pxt0->tag[iarf[i][2]] = xt.tag[iarf[i][2]];
              pxt0->ori = xt.ori;
              isused=1;
            } else {
              mesh->xt++;
              if ( mesh->xt > mesh->xtmax ) {
                TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             printf("  Exit program.\n");
                             exit(EXIT_FAILURE));
              }
              pt1->xt = mesh->xt;
              pxt0 = &mesh->xtetra[pt1->xt];
              pxt0->ref[i]   = xt.ref[i] ; pxt0->ftag[i]  = xt.ftag[i];
              pxt0->edg[iarf[i][0]] = xt.edg[iarf[i][0]];
              pxt0->edg[iarf[i][1]] = xt.edg[iarf[i][1]];
              pxt0->edg[iarf[i][2]] = xt.edg[iarf[i][2]];
              pxt0->tag[iarf[i][0]] = xt.tag[iarf[i][0]];
              pxt0->tag[iarf[i][1]] = xt.tag[iarf[i][1]];
              pxt0->tag[iarf[i][2]] = xt.tag[iarf[i][2]];
              pxt0->ori = xt.ori;
            }
          }
          else {
            pt1->xt = 0;
          }
        }

        if ( jel ) {
          iadr = (jel-1)*4 + 1;
          adjb = &mesh->adja[iadr];
          adjb[j] = iel*4 + i;
        }

        /* internal faces (p1,p2,ip) */
        for (j=0; j<4; j++) {
          if ( j != i ) {
            m = 0;
            for (l=0; l<3; l++)
              if ( pt1->v[ idir[j][l] ] != ip ) {
                v[m] = pt1->v[ idir[j][l] ];
                m++;
              }
            hashEdgeDelone(mesh,&hedg,iel,j,v);
          }
        }
      }
    }
  }

  /* remove old tetra */
  tref = mesh->tetra[list[0]].ref;
  for (k=0; k<ilist; k++) {
    if(tref!=mesh->tetra[list[k]].ref)
      printf("arg ref ???? %d %d\n",tref,mesh->tetra[list[k]].ref);
    delElt(mesh,list[k]);
  }

  //ppt = &mesh->point[ip];
  //  ppt->flag = mesh->flag;
  DEL_MEM(mesh,hedg.item,(hedg.max+1)*sizeof(hedge));
  return(1);
}


/* /\* clone of delone *\/ */
/* int MMG_delons(pMesh mesh,pSol sol,pQueuenew queue,int ip,pList list,int ilist,double declic) { */
/*   pPoint    ppt; */
/*   pTetra    pt,pt1; */
/*   int      *adja,*adjb,i,j,k,l,m,iel,jel,old,v[3],iadr,base,size; */
/*   int       vois[4],ii,kk,iare1,iare2,tref; */
/*   short     i1; */
/*   char      alert; */

/*   if ( mesh->ne + 2*ilist > mesh->nemax )  return(0); */
/*   base = mesh->mark; */

/*   /\* external faces *\/ */
/*   size = 0; */
/*   for (k=1; k<=ilist; k++) { */
/*     old  = list->tetra[k]; */
/*     pt1  = &mesh->tetra[old]; */
/*     iadr = (old-1)*4 + 1; */
/*     adja = &mesh->adja[iadr]; */
/*     for (i=0; i<4; i++) { */
/*       jel = adja[i] >> 2; */
/*       if ( !jel || mesh->tetra[jel].mark != base ) { */
/* 	for (j=0; j<3; j++) { */
/* 	  i1  = MMG_idir[i][j]; */
/* 	  ppt = &mesh->point[ pt1->v[i1] ]; */
/* 	  ppt->tag |= M_CAVITY; */
/* 	} */
/* 	size++; */
/*       } */
/*     } */
/*   } */

/*   /\* check isolated vertex *\/ */
/*   alert = 0; */
/*   for (k=1; k<=ilist; k++) { */
/*     old  = list->tetra[k]; */
/*     pt1  = &mesh->tetra[old]; */
/*     for (i=0; i<4; i++) { */
/*       ppt = &mesh->point[ pt1->v[i] ]; */
/*       if ( !(ppt->tag & M_CAVITY) )  alert = 1; */
/*     } */
/*   } */
/*   /\* reset tag *\/ */
/*   for (k=1; k<=ilist; k++) { */
/*     old  = list->tetra[k]; */
/*     pt1  = &mesh->tetra[old]; */
/*     for (i=0; i<4; i++) { */
/*       ppt = &mesh->point[ pt1->v[i] ]; */
/*       ppt->tag &= ~M_CAVITY; */
/*     } */
/*   } */
/*   if ( alert )  return(-1); */

/*   /\* hash table params *\/ */
/*   if ( size > 3*LONMAX )  return(0); */
/*   list->hedg.size  = size; */
/*   list->hedg.nhmax = 3*size+1; */
/*   list->hedg.hnxt  = size; */
/*   memset(list->hedg.item,0,list->hedg.nhmax*sizeof(hedge)); */
/*   for (k=size; k<list->hedg.nhmax; k++) */
/*     list->hedg.item[k].nxt = k+1; */

/*   for (k=1; k<=ilist; k++) { */
/*     old  = list->tetra[k]; */
/*     iadr = (old-1)*4 + 1; */
/*     adja = &mesh->adja[iadr]; */
/*     vois[0]  = adja[0]; */
/*     vois[1]  = adja[1]; */
/*     vois[2]  = adja[2]; */
/*     vois[3]  = adja[3]; */
/*     pt   = &mesh->tetra[old]; */

/*     for (i=0; i<4; i++) { */
/*       jel = vois[i] >> 2; */
/*       j   = vois[i] % 4; */

/*       /\* external face *\/ */
/*       if ( !jel || mesh->tetra[jel].mark != base ) { */
/* 	iel = MMG_newElt(mesh); */
/* 	if ( iel < 1 )  return(0); */
/* 	pt1 = &mesh->tetra[iel]; */
/* 	pt1->ref = mesh->tetra[old].ref; */
/* 	mesh->point[ip].tmp = iel; */
/* 	memcpy(pt1,pt,sizeof(Tetra)); */
/* 	pt1->v[i] = ip; */
/* 	pt1->qual = MMG_caltet(mesh,sol,iel); */
/* 	pt1->flag = mesh->flag; */
/* 	pt1->edge = 0; */
/* 	for(ii=0 ; ii<4 ; ii++) */
/* 	  pt1->bdryref[ii] = -1; */
/* 	/\*MAJ bdryinfo *\/ */
/* 	for(ii=0 ; ii<6 ; ii++) */
/* 	  pt1->bdryinfo[ii] = 0; */
/* 	for(ii=0 ; ii<3 ; ii++) { */
/* 	  if(!pt->bdryinfo[MMG_iarf[i][ii]]) continue; */
/* 	  iare1 = pt->v[MMG_iare[MMG_iarf[i][ii]][0]]; */
/* 	  iare2 = pt->v[MMG_iare[MMG_iarf[i][ii]][1]]; */
/* 	  for(kk=0 ; kk<3 ; kk++) { */
/* 	    if(((iare1==pt->v[MMG_iare[MMG_iarf[i][kk]][0]]) && (iare2==pt->v[MMG_iare[MMG_iarf[i][kk]][1]])) */
/* 	       || ((iare2==pt->v[MMG_iare[MMG_iarf[i][kk]][0]]) && (iare1==pt->v[MMG_iare[MMG_iarf[i][kk]][1]])) ) { */
/* 	      pt1->bdryinfo[MMG_iarf[i][kk]] = pt->bdryinfo[MMG_iarf[i][ii]]; */
/* 	      break; */
/* 	    } */
/* 	  } */
/* 	  assert(kk<3); */
/* 	} */
/* 	if(!jel || (mesh->tetra[jel].ref != pt1->ref)) { */
/* 	  pt1->bdryref[i] = mesh->tetra[old].bdryref[i]; */
/* 	  if(pt1->bdryref[i]<0) {puts("delone : pbs sd");exit(0);  } */
/* 	} */
/* 	iadr = (iel-1)*4 + 1; */
/* 	adjb = &mesh->adja[iadr]; */
/* 	adjb[i] = adja[i]; */
/* 	if ( jel ) { */
/* 	  iadr = (jel-1)*4 + 1; */
/* 	  adjb = &mesh->adja[iadr]; */
/* 	  adjb[j] = iel*4 + i; */
/* 	} */
/* 	if ( pt1->qual >= declic ) */
/* 	  MMG_kiuputnew(queue,iel); */

/* 	/\* internal faces (p1,p2,ip) *\/ */
/* 	for (j=0; j<4; j++) { */
/* 	  if ( j != i ) { */
/* 	    m = 0; */
/* 	    for (l=0; l<3; l++) */
/* 	      if ( pt1->v[ MMG_idir[j][l] ] != ip ) { */
/* 		v[m] = pt1->v[ MMG_idir[j][l] ]; */
/* 		m++; */
/* 	      } */
/* 	    MMG_hashEdge(mesh,&list->hedg,iel,j,v); */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   /\* remove old tetra *\/ */
/*   for (k=1; k<=ilist; k++) { */
/*     old = list->tetra[k]; */
/*     MMG_delElt(mesh,old); */
/*     //MMG_kiudel(queue,old); */
/*   } */

/*   ppt = &mesh->point[ip]; */
/*   ppt->flag = mesh->flag; */
/*   return(1); */
/* } */


/* /\* cavity correction for quality *\/ */
/* int MMG_correction_ani(pMesh mesh,pSol sol,int ip,pList list,int ilist,int nedep) { */
/*   pPoint        ppt,p1,p2,p3; */
/*   pTetra        pt; */
/*   double        dd,det,nn,eps,eps2,ux,uy,uz,vx,vy,vz,v1,v2,v3; */
/*   double        *ma,*mb,*mc,*md,mm[6],h1,h2,h3; */
/*   int           *adja,i,j,ipil,iel,lon,iadr,adj,ib,ic,id,base,ncor; */
/*   int           vois[4]; */

/*   ppt  = &mesh->point[ip]; */
/*   if ( ppt->tag & M_UNUSED )  return(ilist); */
/*   base = mesh->mark; */
/*   lon  = ilist; */
/*   eps  = EPSCON; */
/*   eps2 = eps*eps; */

/*   /\* average metric *\/ */
/*   memset(mm,0,6*sizeof(double)); */
/*   iadr = (ip-1)*sol->offset + 1; */
/*   ma   = &sol->met[iadr]; */

/*   do { */
/*     ipil = lon; */
/*     ncor = 0; */

/*     while ( ipil > 0 ) { */
/*       iel  = list->tetra[ipil]; */
/*       iadr = (iel-1)*4 + 1; */
/*       adja = &mesh->adja[iadr]; */
/*       vois[0]  = adja[0] >> 2; */
/*       vois[1]  = adja[1] >> 2; */
/*       vois[2]  = adja[2] >> 2; */
/*       vois[3]  = adja[3] >> 2; */
/*       pt   = &mesh->tetra[iel]; */

/*       MMG_cas=0; */
/*       for (i=0; i<4; i++) { */
/* 	adj = vois[i]; */
/* 	MMG_cas = 0; */
/* 	if ( adj && mesh->tetra[adj].mark == base)  continue; */

/* 	ib = pt->v[ MMG_idir[i][0] ]; */
/* 	ic = pt->v[ MMG_idir[i][1] ]; */
/* 	id = pt->v[ MMG_idir[i][2] ]; */

/* 	p1 = &mesh->point[ib]; */
/* 	p2 = &mesh->point[ic]; */
/* 	p3 = &mesh->point[id]; */

/* 	ux = p2->c[0] - p1->c[0]; */
/* 	uy = p2->c[1] - p1->c[1]; */
/* 	uz = p2->c[2] - p1->c[2]; */

/* 	vx = p3->c[0] - p1->c[0]; */
/* 	vy = p3->c[1] - p1->c[1]; */
/* 	vz = p3->c[2] - p1->c[2]; */

/* 	/\* volume PABC *\/ */
/* 	v1 = uz*vy - uy*vz; */
/* 	v2 = ux*vz - uz*vx; */
/* 	v3 = uy*vx - ux*vy; */
/* 	dd = v1*(ppt->c[0]-p1->c[0]) + v2*(ppt->c[1]-p1->c[1]) \ */
/* 	  + v3*(ppt->c[2]-p1->c[2]); */
/* 	MMG_cas=1; */
/* 	//if ( dd < VOLMIN )  break; */
/* 	/\*test sur le volume avec un eps local*\/ */
/* 	h1 = ux*ux + uy*uy + uz*uz; */
/* 	h2 = vx*vx + vy*vy + vz*vz; */
/* 	h3 = (p2->c[0] - p3->c[0])*(p2->c[0] - p3->c[0]) + (p2->c[1] - p3->c[1])*(p2->c[1] - p3->c[1]) */
/* 	  + (p2->c[2] - p3->c[2])*(p2->c[2] - p3->c[2]); */
/* 	if ( dd < VOLMIN*sqrt(h1*h2*h3) )  break; */

/* 	/\* average metric *\/ */
/* 	iadr = (ib-1)*sol->offset + 1; */
/* 	mb   = &sol->met[iadr]; */
/* 	iadr = (ic-1)*sol->offset + 1; */
/* 	mc   = &sol->met[iadr]; */
/* 	iadr = (id-1)*sol->offset + 1; */
/* 	md   = &sol->met[iadr]; */
/* 	for (j=0; j<6; j++) */
/* 	  mm[j] = 0.25 * (ma[j]+mb[j]+mc[j]+md[j]); */

/* 	det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \ */
/* 	  - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \ */
/* 	  + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]); */
/* 	if ( det < EPSOK )  break; */

/* 	/\* point close to face *\/ */
/* 	/\*nn = (v1*v1 + v2*v2 + v3*v3);*\/ */
/* 	MMG_cas=2; */
/* 	nn = mm[0]*v1*v1 + mm[3]*v2*v2 + mm[5]*v3*v3 \ */
/* 	  + 2.0*(mm[1]*v1*v2 + mm[2]*v1*v3 + mm[4]*v2*v3); */
/* 	/\*if ( det*dd*dd*dd*dd*dd*dd < nn * nn * nn * eps2 * eps2 * eps2 )  break;*\/ */
/* 	/\*\//prendre le min des valeurs propres */
/* 	  eigenv(1,mm,lambda,vv); */
/* 	  det = max(lambda[0],max(lambda[1],lambda[2])); */
/* 	  if ( det*dd*dd < nn * eps2 )  break; */
/* 	*\//\*if ( pow(det,1./3.)*dd*dd < nn * eps2 )  break;*\/ */
/* 	if ( det*dd*dd < nn * eps2 )  break; */
/* 	/\*if ( dd*dd < nn * eps2 ) { */
/* 	  printf("en iso      : %e %e    %e %e\n",dd,nn,dd*dd,nn*eps2); */
/* 	  printf("en iso sqrt : %e %e    %e %e\n",dd,nn,dd/sqrt(nn),(sqrt(mm[0]))*(dd/sqrt(nn))); */

/* 	  dd1 = mm[0]*v1*v1 + mm[3]*v2*v2 + mm[5]*v3*v3 \ */
/* 	  + 2.0*(mm[1]*v1*v2 + mm[2]*v1*v3 + mm[4]*v2*v3); */
/* 	  //len carre = (dd*dd/norm(v1v2v3)^2)*dd1/(norm(v1v2v3)^2 */
/* 	  printf("aniso      : %e %e %e %e %e\n",(dd*dd/nn)*dd1/(nn),sqrt(dd*dd*dd1/(nn*sqrt(nn))),det,det*dd*dd,dd1*eps2); */

/* 	  nn = sqrt(nn); */
/* 	  ph = dd/nn; */
/* 	  v1 /= nn; */
/* 	  v2 /= nn; */
/* 	  v3 /= nn; */
/* 	  xh = ph*v1 + ppt->c[0]; */
/* 	  yh = ph*v2 + ppt->c[1]; */
/* 	  zh = ph*v3 + ppt->c[2]; */

/* 	  //dist PH dans la met/ */
/* 	  ux = xh - ppt->c[0]; */
/* 	  uy = yh - ppt->c[1]; */
/* 	  uz = zh - ppt->c[2]; */
/* 	  dd = ux*ux + uy*uy + uz*uz; */

/* 	  dd2 =      mm[0]*ux*ux + mm[3]*uy*uy + mm[5]*uz*uz \ */
/* 	  + 2.0*(mm[1]*ux*uy + mm[2]*ux*uz + mm[4]*uy*uz); */
/* 	  if ( dd2 <= 0.0 )  dd2 = 0.0; */

/* 	  len = sqrt(dd2); */

/* 	  printf("on trouve len : %e %e %e\n",len,sqrt(eps2)*sqrt(mm[0]),pow(sqrt(eps2)*sqrt(det),1./3.)); */
/* 	  printf("len carre %e %e\n",mm[0]*v1*v1*ph*ph + mm[3]*v2*v2*ph*ph + mm[5]*v3*v3*ph*ph,dd2); */
/* 	  exit(0); */
/* 	  break; */
/* 	  }*\/ */
/* 	MMG_cas=0; */
/*       } */
/*       if ( i < 4 ) { */
/* 	if ( ipil <= nedep )  return(0); */
/* 	/\* remove iel from list *\/ */
/* 	pt->mark = base-1; */
/* 	list->tetra[ipil] = list->tetra[lon]; */
/* 	lon--; */
/* 	ncor = 1; */
/* 	break; */
/*       } */
/*       else */
/* 	ipil--; */
/*     } */
/*   } */
/*   while ( ncor > 0 && lon >= nedep ); */

/*   return(lon); */
/* } */


/* cavity correction for quality */
int correction_iso(pMesh mesh,int ip,int *list,int ilist,int nedep) {
  pPoint   ppt,p1,p2,p3;
  pTetra   pt;
  double   dd,nn,eps,eps2,ux,uy,uz,vx,vy,vz,v1,v2,v3;
  int     *adja,i,ipil,iel,lon,iadr,adj,ib,ic,id,base,ncor;
  int      vois[4];

  ppt  = &mesh->point[ip];
  if ( ppt->tag & MG_NUL )  return(ilist);
  base = mesh->mark;
  lon  = ilist;
  eps  = EPSCON;
  eps2 = eps*eps;
  do {
    ipil = lon-1;
    ncor = 0;

    while ( ipil >= 0 ) {
      iel  = list[ipil];
      iadr = (iel-1)*4 + 1;
      adja = &mesh->adja[iadr];
      vois[0]  = adja[0] >> 2;
      vois[1]  = adja[1] >> 2;
      vois[2]  = adja[2] >> 2;
      vois[3]  = adja[3] >> 2;
      pt   = &mesh->tetra[iel];
      MMG_cas=0;
      for (i=0; i<4; i++) {
        adj = vois[i];
        MMG_cas = 0;
        if ( adj && mesh->tetra[adj].mark == base )  continue;

        ib = pt->v[ idir[i][0] ];
        ic = pt->v[ idir[i][1] ];
        id = pt->v[ idir[i][2] ];

        p1 = &mesh->point[ib];
        p2 = &mesh->point[ic];
        p3 = &mesh->point[id];

        ux = p2->c[0] - p1->c[0];
        uy = p2->c[1] - p1->c[1];
        uz = p2->c[2] - p1->c[2];

        vx = p3->c[0] - p1->c[0];
        vy = p3->c[1] - p1->c[1];
        vz = p3->c[2] - p1->c[2];

        /* volume PABC */
        v1 = uz*vy - uy*vz;
        v2 = ux*vz - uz*vx;
        v3 = uy*vx - ux*vy;
        dd = v1*(ppt->c[0]-p1->c[0]) + v2*(ppt->c[1]-p1->c[1]) \
          + v3*(ppt->c[2]-p1->c[2]);
        MMG_cas=1;
        //printf("on trouve vol %e <? %e\n",dd,VOLMIN);
        if ( dd < VOLMIN )  break;

        /* point close to face */
        nn = (v1*v1 + v2*v2 + v3*v3);
        MMG_cas=2;
        //printf("on trouve close ? %e %e\n",dd*dd,nn*eps2);
        if ( dd*dd < nn * eps2 )  break;
        MMG_cas=0;
      }
      if ( i < 4 ||  pt->tag & MG_REQ ) {
        if ( ipil <= nedep )  {/*printf("on veut tout retirer ? %d %d\n",ipil,nedep);*/return(0);   }
        /* remove iel from list */
        pt->mark = base-1;
        list[ipil] = list[--lon];

        ncor = 1;
        break;
      }
      else
        ipil--;
    }
  }
  while ( ncor > 0 && lon >= nedep );
  
  return(lon);
}


/* /\* mark elements in cavity *\/ */
/* int MMG_cavity_ani(pMesh mesh,pSol sol,int iel,int ip,pList list,int lon) { */
/*   pPoint    ppt; */
/*   pTetra    pt,pt1,ptc; */
/*   double    c[3],eps,dd,ray,ux,uy,uz,crit; */
/*   double    *mj,*mp,ct[12]; */
/*   int       *adja,*adjb,k,adj,adi,voy,i,j,ia,ilist,ipil,jel,iadr,base; */
/*   int       vois[4],l; */

/*   if ( lon < 1 )  return(0); */
/*   ppt = &mesh->point[ip]; */
/*   if ( ppt->tag & M_UNUSED )  return(0); */

/*   for (k=1; k<=lon; k++) */
/*     list->tetra[k] = list->tetra[k] / 6; */

/*   /\* grow cavity by adjacency *\/ */
/*   base  = mesh->mark; */
/*   eps   = EPSRAD * EPSRAD; */
/*   ilist = lon; */
/*   ipil  = 1; */
/*   iadr  = (ip-1)*sol->offset + 1; */
/*   mp    = &sol->met[iadr]; */

/*   do { */
/*     jel  = list->tetra[ipil]; */
/*     iadr = (jel-1)*4 + 1; */
/*     adja = &mesh->adja[iadr]; */
/*     vois[0]  = adja[0]; */
/*     vois[1]  = adja[1]; */
/*     vois[2]  = adja[2]; */
/*     vois[3]  = adja[3]; */
/*     ptc  = &mesh->tetra[jel]; */

/*     for (i=0; i<4; i++) { */
/*       adj = vois[i] >> 2; */
/*       voy = vois[i] % 4; */
/*       if ( !adj )  continue; */
/*       pt  = &mesh->tetra[adj]; */
/*       /\* boundary face *\/ */
/*       if ( pt->mark == base || pt->ref != ptc->ref )  continue; */
/*       for (j=0,l=0; j<4; j++,l+=3) { */
/* 	memcpy(&ct[l],mesh->point[pt->v[j]].c,3*sizeof(double)); */
/*       } */


/*       /\* Delaunay kernel *\/ */
/*       if ( !MMG_cenrad_ani(mesh,ct,mp,c,&ray) )  continue; */

/*       ux = ppt->c[0] - c[0]; */
/*       uy = ppt->c[1] - c[1]; */
/*       uz = ppt->c[2] - c[2]; */
/*       dd =      mp[0]*ux*ux + mp[3]*uy*uy + mp[5]*uz*uz \ */
/* 	+ 2.0*(mp[1]*ux*uy + mp[2]*ux*uz + mp[4]*uy*uz); */
/*       crit = eps * ray; */
/*       if ( dd > crit )  continue; */

/*       /\* mixed metrics *\/ */
/*       crit = sqrt(dd/ray); */
/*       for (j=0; j<4; j++) { */
/* 	ia   = pt->v[j]; */
/* 	iadr = (ia-1)*sol->offset + 1; */
/* 	mj   = &sol->met[iadr]; */
/* 	if ( !MMG_cenrad_ani(mesh,ct,mj,c,&ray) )  continue; */
/* 	ux = ppt->c[0] - c[0]; */
/* 	uy = ppt->c[1] - c[1]; */
/* 	uz = ppt->c[2] - c[2]; */
/* 	dd =      mj[0]*ux*ux + mj[3]*uy*uy + mj[5]*uz*uz \ */
/* 	  + 2.0*(mj[1]*ux*uy + mj[2]*ux*uz + mj[4]*uy*uz); */
/* 	crit += sqrt(dd/ray); */
/*       } */
/*       crit *= EPSRAD; */
/*       if ( crit > 5.0 ) continue; */

/*       /\* lost face(s) *\/ */
/*       iadr = (adj-1)*4 + 1; */
/*       adjb = &mesh->adja[iadr]; */

/*       for (j=0; j<4; j++) { */
/* 	if ( j == voy )  continue; */
/* 	adi = adjb[j] >> 2; */
/* 	if ( !adi )  continue; */
/* 	pt1 = &mesh->tetra[adi]; */
/* 	if ( pt1->mark == base && adi != jel ) { */
/* 	  if ( !adi || pt1->ref != mesh->tetra[adi].ref )  break; */
/* 	} */
/*       } */
/*       /\* store tetra *\/ */
/*       if ( j == 4 ) { */
/* 	pt->mark = base; */
/* 	++ilist; */
/* 	list->tetra[ilist] = adj; */
/*       } */
/*     } */
/*     if ( ilist > LONMAX - 3 )  return(-1); */
/*     ++ipil; */
/*   } */
/*   while ( ipil <= ilist ); */

/*   /\* global overflow *\/ */
/*   if ( mesh->ne + 2*ilist >= mesh->nemax ) */
/*     ilist = -ilist; */
/*   else */
/*     ilist = MMG_correction_ani(mesh,sol,ip,list,ilist,lon); */

/*   if(MMG_cas==1) MMG_nvol++; */
/*   else if(MMG_cas==2 || MMG_cas>20) { */
/*     MMG_npuiss++; */
/*     if(MMG_cas>20) MMG_npres++; */
/*   } */

/*   return(ilist); */
/* } */


/** Return a negative value for ilist if one of the tet of the cavity is required */
int cavity(pMesh mesh,pSol sol,int iel,int ip,int *list,int lon) {
  pPoint    ppt;
  pTetra    pt,pt1,ptc;
  double    c[3],crit,dd,eps,ray,ct[12];
  int      *adja,*adjb,k,adj,adi,voy,i,j,ilist,ipil,jel,iadr,base;
  int       vois[4],l;
  int       tref,isreq;

  if ( lon < 1 )  return(0);
  ppt = &mesh->point[ip];
  if ( ppt->tag & MG_NUL )  return(0);
  base  = ++mesh->mark;

  isreq = 0;

  tref = mesh->tetra[list[0]/6].ref;
  for (k=0; k<lon; k++) {
    mesh->tetra[list[k]/6].mark = base;

    if (tref!=mesh->tetra[list[k]/6].ref) {
      //printf("pbs coquil %d %d tet %d\n",tref,mesh->tetra[list[k]/6].ref,list[k]/6);
      return(0);
    }
  }
  for (k=0; k<lon; k++)
    list[k] = list[k] / 6;

  /* grow cavity by adjacency */
  eps   = EPSRAD*EPSRAD;
  ilist = lon;
  ipil  = 0;

  do {
    jel  = list[ipil];
    iadr = (jel-1)*4 + 1;
    adja = &mesh->adja[iadr];
    vois[0]  = adja[0];
    vois[1]  = adja[1];
    vois[2]  = adja[2];
    vois[3]  = adja[3];
    ptc  = &mesh->tetra[jel];

    for (i=0; i<4; i++) {
      adj = vois[i] >> 2;
      voy = vois[i] % 4;
      if ( !adj )  continue;
      pt  = &mesh->tetra[adj];
      /* boundary face */
      if ( pt->mark == base || pt->ref != ptc->ref )  continue;

      for (j=0,l=0; j<4; j++,l+=3) {
        memcpy(&ct[l],mesh->point[pt->v[j]].c,3*sizeof(double));
      }

      if ( !cenrad_iso(mesh,ct,c,&ray) )  continue;
      crit = eps * ray;

      /* Delaunay criterion */
      dd = (ppt->c[0] - c[0]) * (ppt->c[0] - c[0]) \
        + (ppt->c[1] - c[1]) * (ppt->c[1] - c[1]) \
        + (ppt->c[2] - c[2]) * (ppt->c[2] - c[2]);
      if ( dd > crit )  continue;

      /* lost face(s) */
      iadr = (adj-1)*4 + 1;
      adjb = &mesh->adja[iadr];

      for (j=0; j<4; j++) {
        if ( j == voy )  continue;
        adi = adjb[j] >> 2;
        if ( !adi )  continue;
        pt1 = &mesh->tetra[adi];
        if ( pt1->mark == base && adi != jel ) {
          if ( !adi || pt1->ref != tref )  break;
        }
      }
      /* store tetra */
      if ( j == 4 ) {
	 if ( pt->tag & MG_REQ ) isreq = 1;
        pt->mark = base;
        list[ilist++] = adj;
      }
    }
    if ( ilist > LONMAX - 3 ) return(-1);
    
    ++ipil;
  }
  while ( ipil < ilist );

  /* global overflow: obsolete avec la reallocation */
  //if ( mesh->ne + 2*ilist >= mesh->nemax ) {

  ilist = correction_iso(mesh,ip,list,ilist,lon);
 
  if ( isreq ) ilist = -fabs(ilist);

  if(MMG_cas==1) MMG_nvol++;
  else if(MMG_cas==2 || MMG_cas>20) {
    MMG_npuiss++;
    if(MMG_cas>20) MMG_npres++;
  }
  return(ilist);
}
