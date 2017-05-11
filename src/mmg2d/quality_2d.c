/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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
 * \file mmg2d/quality_2d.c
 * \brief Functions to compute the quality.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg2d.h"

/* Fast calculation of the quality of Tria pt */
double _MMG2_quickcal(MMG5_pMesh mesh, MMG5_pTria pt) {
  MMG5_pPoint        p0,p1,p2;
  double             cal;
  
  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];
  
  cal = MMG2_quickarea(p0->c,p1->c,p2->c);
  return(cal);
}

/* Compute quality of the triangle pt when the supplied metric is isotropic; 
   return 0 in the case that the triangle has inverted orientation */
double _MMG2_caltri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria pt) {
  double    abx,aby,acx,acy,bcx,bcy;
  double    *a,*b,*c,h1,h2,h3,area,hm;

  a  = mesh->point[pt->v[0]].c;
  b  = mesh->point[pt->v[1]].c;
  c  = mesh->point[pt->v[2]].c;

  abx = b[0] - a[0];
  aby = b[1] - a[1];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];

  /* orientation */
  area = abx*acy - aby*acx;
  if ( area <= 0.0 ) return(0.0);

  /* edge lengths */
  h1 = abx*abx + aby*aby;
  h2 = acx*acx + acy*acy;
  h3 = bcx*bcx + bcy*bcy;

  hm = h1 + h2 + h3;
  h1 = sqrt(h1);
  h2 = sqrt(h2);
  h3 = sqrt(h3);
  
  if ( hm > _MMG2_EPSD ) {
    return ( area / hm );
  }
  else {
    return(0.0);
  }
}

/* Compute quality of the triangle pt when the supplied metric is anisotropic;
 return 0 in the case that the triangle has inverted orientation */
double _MMG2_caltri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria pt) {
  double     abx,aby,acx,acy,bcx,bcy;
  double     *a,*b,*c,*ma,*mb,*mc;
  double     area,aream,lab,lac,lbc,l1,l2,hm;
  int        ipa,ipb,ipc;
 
  ipa = pt->v[0];
  ipb = pt->v[1];
  ipc = pt->v[2];
  
  a  = mesh->point[ipa].c;
  b  = mesh->point[ipb].c;
  c  = mesh->point[ipc].c;
  
  ma = &met->m[3*ipa];
  mb = &met->m[3*ipb];
  mc = &met->m[3*ipc];

  /* Isotropic area of pt */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  area = 0.5*(abx*acy - aby*acx);
  if ( area <= 0.0 ) return(0.0);
  
  /* Anisotropic volume of pt */
  aream = 0.0;
  aream += sqrt(ma[0]*ma[2]-ma[1]*ma[1]);
  aream += sqrt(mb[0]*mb[2]-mb[1]*mb[1]);
  aream += sqrt(mc[0]*mc[2]-mc[1]*mc[1]);
  aream = _MMG5_ATHIRD*area*aream;
  
  /* (squared) Anisotropic edge lengths */
  l1 = ma[0]*abx*abx + 2.0*ma[1]*abx*aby + ma[2]*aby*aby;
  l2 = mb[0]*abx*abx + 2.0*mb[1]*abx*aby + mb[2]*aby*aby;
  lab = 0.25*(l1+l2) + 0.5*sqrt(l1*l2);
  
  l1 = ma[0]*acx*acx + 2.0*ma[1]*acx*acy + ma[2]*acy*acy;
  l2 = mc[0]*acx*acx + 2.0*mc[1]*acx*acy + mc[2]*acy*acy;
  lac = 0.25*(l1+l2) + 0.5*sqrt(l1*l2);
  
  l1 = mb[0]*bcx*bcx + 2.0*mb[1]*bcx*bcy + mb[2]*bcy*bcy;
  l2 = mc[0]*bcx*bcx + 2.0*mc[1]*bcx*bcy + mc[2]*bcy*bcy;
  lbc = 0.25*(l1+l2) + 0.5*sqrt(l1*l2);

  /* Quality measure = (Vol_M(T) / (l(ab)^2+l(ac)^2+l(bc)^2)) */
  hm = lab + lac + lbc;
  if ( hm > _MMG2_EPSD ) {
    return ( aream/hm );
  }
  else {
    return (0.0);
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 *
 * Print histogram of mesh qualities.
 *
 */
void MMG2_outqua(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria    pt;
  double   rap,rapmin,rapmax,rapavg,med,good;
  int      i,k,iel,ok,ir,imax,nex,his[5];
  
  /* Compute triangle quality*/
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if( !MG_EOK(pt) )   continue;
    
#warning: add check for null quality
    
    pt->qual = MMG2D_caltri(mesh,met,pt);
  }
  if ( abs(mesh->info.imprim) <= 0 ) return;

  rapmin  = 2.0;
  rapmax  = 0.0;
  rapavg  = med = good = 0.0;
  iel     = 0;

  for (k=0; k<5; k++)  his[k] = 0;

  nex = ok = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if( !MG_EOK(pt) ) {
      nex++;
      continue;
    }
    ok++;
    rap = _MMG2D_ALPHAD * MMG2D_caltri(mesh,met,pt);
    if ( rap < rapmin ) {
      rapmin = rap;
      iel    = ok;
    }
    if ( rap > 0.5 )  med++;
    if ( rap > 0.12 ) good++;
    if ( rap < _MMG2D_BADKAL )  mesh->info.badkal = 1;
    rapavg += rap;
    rapmax  = MG_MAX(rapmax,rap);
    ir = MG_MIN(4,(int)(5.0*rap));
    his[ir] += 1;
  }

#ifndef DEBUG
  fprintf(stdout,"\n  -- MESH QUALITY   %d\n",mesh->nt - nex);
  fprintf(stdout,"     BEST   %8.6f  AVRG.   %8.6f  WRST.   %8.6f (%d)\n",
          rapmax,rapavg / (mesh->nt-nex),rapmin,iel);
#else
  fprintf(stdout,"     BEST   %e  AVRG.   %e  WRST.   %e (%d)\n => %d %d %d\n",
          rapmax,rapavg / (mesh->nt-nex),rapmin,iel,
          _MMG5_indPt(mesh,mesh->tria[iel].v[0]),_MMG5_indPt(mesh,mesh->tria[iel].v[1]),
          _MMG5_indPt(mesh,mesh->tria[iel].v[2]));
#endif

  /* print histo */
  fprintf(stdout,"     HISTOGRAMM:");
  fprintf(stdout,"  %6.2f %% > 0.12\n",100.0*(good/(float)(mesh->nt-nex)));
  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"                  %6.2f %% >  0.5\n",100.0*( med/(float)(mesh->nt-nex)));
    imax = MG_MIN(4,(int)(5.*rapmax));
    for (i=imax; i>=(int)(5*rapmin); i--) {
      fprintf(stdout,"     %5.1f < Q < %5.1f   %7d   %6.2f %%\n",
              i/5.,i/5.+0.2,his[i],100.*(his[i]/(float)(mesh->nt-nex)));
    }
  }
}

/* /\* print mesh quality histo *\/ */
/* void MMG2_outqua(MMG5_pMesh mesh,MMG5_pSol sol) { */
/*   MMG5_pTria     pt; */
/*   double    coef,rap4,rapl,rapmin,rapmax,rapavg; */
/*   int       his10[11],his01[33],rapnum; */
/*   int       k,i,j,imax,iel,ir,nn,nex,ielreal; */

/*   rapmin  =  1.e20; */
/*   rapmax  = -1.e20; */
/*   rapavg  = 0.0; */
/*   rapnum  = 0; */
/*   iel     = 0; */
/*   ielreal = 0; */
/*   nn      = 0; */

/*   for (k=0; k<=33; k++)  his01[k] = 0; */
/*   for (k=0; k<=10; k++)  his10[k] = 0; */

/*   coef = _MMG2D_ALPHA; */
/*   nex  = 0; */
/*   for (k=1; k<=mesh->nt; k++) { */
/*     pt = &mesh->tria[k]; */
/*     if( !MG_EOK(pt) ) { */
/*       nex++; */
/*       continue; */
/*     } */
/*     nn++; */

/*     rap4 = _MMG2D_ALPHA * MMG2_caltri(mesh,sol,pt); */
/*     rap4 = M_MAX(rap4,EPSD); */
/*     ir   = (int)rap4; */
/*     if ( rap4 > rapmax ) { */
/*       rapmax = rap4;  */
/*       iel     = k; */
/*       ielreal = k - nex; */
/*     } */
/*     rapavg += rap4; */
/*     rapnum++; */

/*     if ( rap4 > 1.0 && rap4 < 1e9 ) { */
/*       rapmin = M_MIN(rapmin,rap4); */
/*       if ( rap4 < 10.0 ) { */
/*         his10[ir] += 1; */
/*      his01[0]  += 1; */
/*       } */
/*       else if ( rap4 < 1e9 ) { */
/*         rapl = M_MIN(log10(rap4),32.0); */
/*         his01[(int)rapl] += 1; */
/*      his01[0]  += 1; */
/*       }  */
/*     } else { */
/*        printf("pbs qual %d : %e\n",k,rap4);   */
/*     } */
/*   } */

/*   /\* print histo *\/ */
/*   fprintf(stdout,"\n  -- MESH QUALITY   %d\n",rapnum); */
/*   if ( rapavg / rapnum < 100.0 ) { */
/*     pt = &mesh->tria[iel]; */
/*     fprintf(stdout,"     BEST   %e  AVRG.   %e  WRST.   %e (%d)\n", */
/*             rapmin,rapavg / rapnum,rapmax,ielreal); */
/*   } */


/*   if ( abs(mesh->info.imprim) < 4 )  return; */

/*   fprintf(stdout,"\n     HISTOGRAMM\n"); */
/*   j = 0; */
/*   for (i=1; i<16; i++) */
/*     j += his01[i]; */

/*   for (i=M_MAX((int)rapmin,1); i<=M_MIN((int)rapmax,9); i++) { */
/*     fprintf(stdout,"     %5d < Q < %5d   %7d   %6.2f %%\n", */
/*      i,i+1,his10[i],100.*(his10[i]/(float)his01[0])); */
/*   } */

/*   /\* quality per interval *\/ */
/*   if (j != 0) { */
/*     fprintf(stdout,"\n"); */
/*     imax = (int)(M_MIN(3,log10(rapmax))); */

/*     for (i=1; i<=imax; i++) { */
/*       fprintf(stdout,"     %5.0f < Q < %5.0f   %7d   %6.2f %%\n", */
/*        pow(10.,i),pow(10.,i+1),his01[i],100.*(his01[i]/(float)his01[0])); */
/*     } */
/*     for (i=4; i<=(int)log10(rapmax); i++) { */
/*       fprintf(stdout,"    10**%2d < Q < 10**%2d  %7d   %6.2f %%\n", */
/*        i,i+1,his01[i],100.*(his01[i]/(float)his01[0])); */
/*     } */
/*   } */
/* } */

