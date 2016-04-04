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
 * \file mmg2d/mmg2d9.c
 * \brief Lagrangian meshing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmg2d.h"
//#include "ls_calls.h"
#define _MMG2_DEGTOL  1.e-1

/* Calculate an estimate of the average (isotropic) length of edges in the mesh */
double _MMG2_estavglen(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    p1,p2;
  int            k,na;
  double         len,lent,dna;
  char           i,i1,i2;
  
  na = 0;
  lent = 0.0;
  
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    for (i=0; i<3; i++) {
      i1 = _MMG5_inxt2[i];
      i1 = _MMG5_iprv2[i];
      
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];
      
      len = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]);
      
      lent += sqrt(len);
      na++;
    }
  }
  
  dna = (double)na;
  dna = 1.0 / dna;
  lent *= dna;
  
  return(lent);
}

/** Compute quality of a triangle from the datum of its 3 vertices */
static
inline double _MMG2_caltri_iso_3pt(double *a,double *b,double *c) {
  double        abx,aby,acx,acy,bcx,bcy,area,h1,h2,h3,hm;
  
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

/** Check if moving mesh with disp for a fraction t yields a valid mesh */
int _MMG2_chkmovmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t) {
  MMG5_pTria   pt;
  MMG5_pPoint  ppt;
  double       *v,c[3][2],tau,cal;
  int          k,np;
  char         i,j;
  
  /* Pseudo time-step = fraction of disp to perform */
  tau = (double)t / _MMG2_SHORTMAX;
  
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    
    for (i=0; i<3; i++) {
      np = pt->v[i];
      ppt = &mesh->point[np];
      v = &disp->m[2*(np-1)+1];
      for (j=0; j<2; j++)
        c[i][j] = ppt->c[j]+tau*v[j];
    }
    
    if( _MMG2_caltri_iso_3pt(c[0],c[1],c[2]) < _MMG2_NULKAL) return(0);  //     Other criteria : eg. a rate of degradation, etc... ?
  }
  
  return(1);
}

/* Return the largest fraction t that makes the motion along disp valid */
short _MMG2_dikomv(MMG5_pMesh mesh,MMG5_pSol disp) {
  int     it,maxit;
  short   t,tmin,tmax;
  char    ier;
  
  maxit = 200;
  it    = 0;
  
  tmin  = 0;
  tmax  = _MMG2_SHORTMAX;
  
  /* If full displacement can be achieved */
  if ( _MMG2_chkmovmesh(mesh,disp,tmax) )
    return(tmax);
  
  /* Else, find the largest displacement by dichotomy */
  while( tmin != tmax && it < maxit ) {
    t = (tmin+tmax)/2;
    
    /* Case that tmax = tmin +1 : check move with tmax */
    if ( t == tmin ) {
      ier = _MMG2_chkmovmesh(mesh,disp,tmax);
      if ( ier )
        return(tmax);
      else
        return(tmin);
    }
    
    /* General case: check move with t */
    ier = _MMG2_chkmovmesh(mesh,disp,t);
    if ( ier )
      tmin = t;
    else
      tmax = t;
    
    it++;
  }
  
  return(tmin);
}

/** Perform mesh motion along disp, for a fraction t, and the corresponding updates */
int _MMG2_dispmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t,int itdeg) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  double        *v,tau,ctau,c[3][2],ocal,ncal;
  int           k,np;
  char          i,j;
  
  tau = (double)t /_MMG2_SHORTMAX;
  ctau = 1.0 - tau;
  
  /* Identify elements which are very distorted in the process */
  for (k=1; k<=mesh->nt; k++) {
    pt  = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    
    for (i=0; i<3; i++) {
      np = pt->v[i];
      ppt = &mesh->point[np];
      v = &disp->m[2*(np-1)+1];
      for (j=0; j<2; j++)
        c[i][j] = ppt->c[j];
    }
    
    ocal = _MMG2_caltri_iso_3pt(c[0],c[1],c[2]);
    
    for (i=0; i<3; i++) {
      np = pt->v[i];
      v = &disp->m[2*(np-1)+1];
      for (j=0; j<2; j++)
        c[i][j] += tau*v[j];
    }
    
    ncal = _MMG2_caltri_iso_3pt(c[0],c[1],c[2]);
    
    if ( ncal < _MMG2_DEGTOL*ocal )
      pt->cc = itdeg;
    
  }
  
  /* Perform physical displacement */
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    
    if ( !MG_VOK(ppt) ) continue;
    v = &disp->m[2*(k-1)+1];
    
    for (i=0; i<2; i++) {
      ppt->c[i] = ppt->c[i] + tau*v[i];
      v[i] *= ctau;
    }
  }
  
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param disp pointer toward the displacement structure.
 * \param met pointer toward the metric structure.
 * \param itdeg degraded elements.
 * \param *warn \a warn is set to 1 if not enough memory is available to complete mesh.
 * \return -1 if failed.
 * \return number of new points.
 *
 * Split edges of length bigger than _MMG5_LOPTL, in the Lagrangian mode. 
 * Only affects triangles with cc itdeg
 *
 */
int _MMG2_spllag(MMG5_pMesh mesh,MMG5_pSol disp,MMG5_pSol met,int itdeg,int *warn) {
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param itdeg degraded elements.
 * \return -1 if failed.
 * \return number of collapsed points.
 *
 * Attempt to collapse small internal edges in the Lagrangian mode; only affects tetras with cc itdeg.
 *
 */
static int _MMG2_coleltlag(MMG5_pMesh mesh,MMG5_pSol met,int itdeg) {
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param crit coefficient of quality improvment.
 * \param bucket pointer toward the bucket structure in delaunay mode and
 * toward the \a NULL pointer otherwise.
 * \param itdeg degraded elements.
 *
 * Internal edge flipping in the Lagrangian mode; only affects trias with cc itdeg
 *
 */
int _MMG2_swpmshlag(MMG5_pMesh mesh,MMG5_pSol met,double crit,int itdeg) {
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param itdeg degraded elements.
 * \return -1 if failed, number of moved points otherwise.
 *
 * Analyze trias with cc = itdeg and move internal points so as to make mesh more uniform.
 *
 */
int _MMG2_movtrilag(MMG5_pMesh mesh,MMG5_pSol met,int itdeg) {
  return(1);
}

/* Lagrangian node displacement */
int MMG2_mmg2d9(MMG5_pMesh mesh,MMG5_pSol disp,MMG5_pSol met) {
  double             avlen,tau;
  int                k,itmn,itdc,maxitmn,maxitdc,iit,warn;
  int                nspl,nnspl,nnnspl,nc,nnc,nnnc,ns,nns,nnns,nm,nnm,nnnm;
  short              t;
  char               ier;
  
  maxitmn = 10;
  maxitdc = 100;
  t = 0;
  
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** LAGRANGIAN MOTION\n");
  
  /* Field cc stores information about whether a triangle has been greatly distorted during current step */
  for (k=1; k<=mesh->nt; k++)
    mesh->tria[k].cc = 0;
  
  /* Estimate of the average, maximum and minimum edge lengths */
  avlen = _MMG2_estavglen(mesh);
  mesh->info.hmax = MMG2_LLONG*avlen;
  mesh->info.hmin = MMG2_LSHRT*avlen;
  
  for (itmn=0; itmn<maxitmn; itmn++) {
    
    /* Extension of the displacement field */
    if ( !_MMG2_velextLS(mesh,disp) ) {
      fprintf(stdout,"  ## Problem in func. _MMG2_velextLS. Exit program.\n");
      return(0);
    }
    
    /* Sequence of dichotomy loops to find the largest admissible displacements */
    for (itdc=0; itdc<maxitdc; itdc++) {
      nnspl = nnc = nns = nnm = 0;
      
      t = _MMG2_dikomv(mesh,disp);
      if ( t == 0 ) {
        if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
          printf("   *** Stop: impossible to proceed further\n");
        break;
      }
      
      ier = _MMG2_dispmesh(mesh,disp,t,itdc);
      if ( !ier ) {
        fprintf(stdout,"  ** Impossible motion\n");
        return(0);
      }
      
      tau = tau + ((double)t /_MMG2_SHORTMAX)*(1.0-tau);
      if ( (abs(mesh->info.imprim) > 3 ) || mesh->info.ddebug )
        printf("   ---> Realized displacement: %f\n",tau);
    
      /* Local remeshing depending on the option */
      if ( mesh->info.lag > 0 ) {
        for (iit=0; iit<5; iit++) {
          
          nspl = nc = ns = nm = 0;
          
          if ( mesh->info.lag > 1 ) {
            /* Split of points */
            nspl = _MMG2_spllag(mesh,disp,met,itdc,&warn);
            if ( nspl < 0 ) {
              fprintf(stdout,"  ## Problem in spllag. Exiting.\n");
              return(0);
            }
            
            /* Collapse of points */
            nc = _MMG2_coleltlag(mesh,met,itdc);
            if ( nc < 0 ) {
              fprintf(stdout,"  ## Problem in coltetlag. Exiting.\n");
              return(0);
            }
          }
          
          /* Swap of edges in tetra that have resulted distorted from the process */
          /* I do not know whether it is safe to put NULL in metric here (a
           * priori ok, since there is no vertex creation or suppression) */
          ns = _MMG2_swpmshlag(mesh,met,1.1,itdc);
          if ( ns < 0 ) {
            fprintf(stdout,"  ## Problem in swaptetlag. Exiting.\n");
            return(0);
          }
          
          /* Relocate vertices of tetra which have been distorted in the displacement process */
          nm = _MMG2_movtrilag(mesh,met,itdc);
          if ( nm < 0 ) {
            fprintf(stdout,"  ## Problem in movtetlag. Exiting.\n");
            return(0);
          }
          
          if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && (nspl+nc+ns+nm > 0) )
            printf(" %d edges splitted, %d vertices collapsed, %d elements"
                   " swapped, %d vertices moved.\n",nspl,nc,ns,nm);
          nnspl+= nspl;
          nnm  += nm;
          nnc  += nc;
          nns  += ns;
        }
        
        if ( abs(mesh->info.imprim) > 3 && abs(mesh->info.imprim) < 5 && (nnspl+nnm+nns+nnc > 0) )
          printf(" %d edges splitted, %d vertices collapsed, %d elements"
                 " swapped, %d vertices moved.\n",nnspl,nnc,nns,nnm);
        
        
      }
      
      nnnspl += nnspl;
      nnnm   += nnm;
      nnnc   += nnc;
      nnns   += nns;
    
      if ( t == _MMG2_SHORTMAX ) break;
    }

    if ( abs(mesh->info.imprim) > 2 )
      printf(" %d edges splitted, %d vertices collapsed, %d elements"
               " swapped, %d vertices moved.\n",nnnspl,nnnc,nnns,nnnm);
    
    if ( t == _MMG2_SHORTMAX ) break;
  }
  
  /* Clean memory */
  _MMG5_DEL_MEM(mesh,disp->m,(disp->size*(disp->npmax+1))*sizeof(double));
  
  return(1);
}
