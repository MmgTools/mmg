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
#include "mmg2d.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param severe level of performed check
 * \param base unused argument.
 * \return 0 if fail, 1 if success.
 *
 * Check the mesh validity
 *
 */
int _MMG5_mmg2dChkmsh(MMG5_pMesh mesh, int severe,int base) {
//  MMG5_pPoint ppt;
  MMG5_pTria  pt1,pt2;
  MMG5_pEdge         ped;
  int   *adja,*adja1,adj,adj1,k,i,iadr;
//  int   kk,l,nk,j,ip,lon,len;
  int     *list;
  unsigned char voy,voy1;

  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    if ( !M_EOK(pt1) )  continue;
    iadr = (k-1)*3 + 1;
    adja = &mesh->adja[iadr];

    for (i=0; i<3; i++) {
      adj = adja[i] / 3;
      voy = adja[i] % 3;
      if ( !adj )  continue;

      if ( adj == k ) {
        fprintf(stdout,"  1. Wrong adjacency %d %d\n",k,adj);
        printf("vertices of %d: %d %d %d \n",k,pt1->v[0],pt1->v[1],pt1->v[2]);
        printf("adj of %d: %d %d %d \n",
               k,adja[0]/3,adja[1]/3,adja[2]/3);
        return(0);
      }
      pt2 = &mesh->tria[adj];
      if ( !M_EOK(pt2) ) {
        fprintf(stdout,"  4. Invalid adjacent %d %d\n",adj,k);
        printf("vertices of %d: %d %d %d\n",
               k,pt1->v[0],pt1->v[1],pt1->v[2]);
        printf("vertices adj %d: %d %d %d \n",
               adj,pt2->v[0],pt2->v[1],pt2->v[2]);
        printf("adj of %d: %d %d %d\n",k,adja[0]/3,adja[1]/3,adja[2]/3);
        return(0);
      }
      iadr  = (adj-1)*3 + 1;
      adja1 = &mesh->adja[iadr];
      adj1  = adja1[voy] / 3;
      voy1  = adja1[voy] % 3;
      if ( adj1 != k || voy1 != i ) {
        fprintf(stdout,"  2. Wrong adjacency %d %d\n",k,adj1);
        printf("vertices of %d: %d %d %d \n",k,pt1->v[0],pt1->v[1],pt1->v[2]);
        printf("adj(k) %d: %d %d %d \n",adj,pt2->v[0],pt2->v[1],pt2->v[2]);
        printf("adj(%d): %d %d %d\n",
               k,adja[0]/3,adja[1]/3,adja[2]/3);
        printf("adj(%d): %d %d %d %d\n",
               adj,adja1[0]/3,adja1[1]/3,adja1[2]/3,adja1[3]/3);
        return(0);
      }

      /*chk edge*/
      if(pt1->edg[i]) {
        ped = &mesh->edge[pt1->edg[i]];
        if(!(((ped->a==pt1->v[MMG2_iare[i][0]]) || (ped->a==pt1->v[MMG2_iare[i][1]]))
             || ((ped->b==pt1->v[MMG2_iare[i][0]]) || (ped->b==pt1->v[MMG2_iare[i][1]])))) {
          printf("  3. Wrong edge in triangle %d\n",k);
          printf("vertices of %d: %d %d %d \n",k,pt1->v[0],pt1->v[1],pt1->v[2]);
          printf("edge %d : %d %d\n",i,ped->a,ped->b);
          return(0);
        }
      }

    }
  }

  if ( !severe )  return(1);

  _MMG5_SAFE_CALLOC(list,MMG2D_LMAX,int,0);

  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    if ( !M_EOK(pt1) )  continue;
    iadr = 3*(k-1) + 1;
    adja = &mesh->adja[iadr];

    /*for (i=0; i<3; i++) {
      adj = (adja[i]-1) / 3 + 1;
      voy = (adja[i]-1) % 3;
      if ( !adj )  continue;

      ip  = pt1->v[i];
      ppt = &mesh->point[ip];
      if ( !M_VOK(ppt) ) {
      fprintf(stdout,"  6. Unused vertex %d  %d\n",k,ip);
      printf("%d %d %d\n",pt1->v[0],pt1->v[1],pt1->v[2]);
      return(0);
      }
      lon = boulep(mesh,k,i,list);
      for (l=1; l<=lon; l++) {
      kk  = list[l] / 3;
      nk  = list[l] % 3;
      pt2 = &mesh->tria[kk];
      if ( pt2->v[nk] != ip ) {
      fprintf(stdout,"  5. Wrong ball %d, %d\n",ip,pt2->v[nk]);
      return(0);
      }
      }
      if ( lon < 1 )  continue;
      len = 0;
      for (kk=1; kk<=mesh->nt; kk++) {
      pt2 = &mesh->tria[kk];
      if ( !pt2->v[0] )  continue;
      for (j=0; j<3; j++)
      if ( pt2->v[j] == ip ) {
      len++;
      break;
      }
      }
      if ( len != lon ) {
      fprintf(stdout,"  7. Incorrect ball %d: %d %d\n",pt1->v[i],lon,len);
      return(0);
      }
      } */
  }
  _MMG5_SAFE_FREE(list);
  return(1);
}

/* Check of adjacency relations and edge tags */
int _MMG2_chkmsh(MMG5_pMesh mesh) {
  MMG5_pTria        pt,pt1;
  MMG5_pPoint       p1,p2;
  int               *adja,*adjaj,k,jel;
  char              i,i1,i2,j;

  /* Check adjacencies */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];

    for (i=0; i<3; i++) {
      jel = adja[i] / 3;
      j   = adja[i] % 3;

      if ( !jel ) {
        if ( !(pt->tag[i] & MG_GEO ) ) {
          printf("Wrong tag of edge %d in tria %d \n",i,k);
          return(0);
        }
      }
      else {
        pt1 = &mesh->tria[jel];
        adjaj = &mesh->adja[3*(jel-1)+1];
        if ( adjaj[j] / 3 != k ) {
          printf("Wrong adjacencies %d %d \n",k,jel);
          return(0);
        }
      }
    }
  }

  /* Check consistency between tags of edges and vertices */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      i1 = _MMG5_inxt2[i];
      i2 = _MMG5_iprv2[i];
      if ( pt->tag[i] & MG_GEO ) {
        if ( !(mesh->point[pt->v[i1]].tag & MG_GEO) && !( MG_SIN(mesh->point[pt->v[i1]].tag) )) {
          printf("Tag inconsistency in triangle %d: edge %d, vertex %d\n",k,i,pt->v[i1]);
          return(0);
        }
        if ( !(mesh->point[pt->v[i2]].tag & MG_GEO) && !( MG_SIN(mesh->point[pt->v[i2]].tag) )) {
          printf("Tag inconsistency in triangle %d: edge %d, vertex %d\n",k,i,pt->v[i2]);
          return(0);
        }
      }

      if ( pt->tag[i] & MG_REF ) {
        if ( !(mesh->point[pt->v[i1]].tag & MG_REF) && !( MG_SIN(mesh->point[pt->v[i1]].tag)) ) {
          printf("Tag inconsistency in triangle %d: edge ref %d, vertex %d\n",k,i,pt->v[i1]);
          return(0);
        }
        if ( !(mesh->point[pt->v[i2]].tag & MG_REF) && !( MG_SIN(mesh->point[pt->v[i2]].tag)) ) {
          printf("Tag inconsistency in triangle %d: edge ref %d, vertex %d\n",k,i,pt->v[i2]);
          return(0);
        }
      }

    }
  }

  /* Check consistency between edge tags and triangle refs */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      i1 = _MMG5_inxt2[i];
      i2 = _MMG5_iprv2[i];

      jel = adja[i] / 3;

      if ( ( pt->tag[i] & MG_GEO ) && jel ) {
        printf("edge %d %d is tagged boundary while it has a neighbour\n",pt->v[i1],pt->v[i2]);
        return(0);
      }

      if ( pt->tag[i] & MG_REF ) {
        pt1 = &mesh->tria[jel];
        if ( pt->ref == pt1->ref ) {
          printf("edge %d %d is tagged ref while both corresponding triangles have same ref\n",pt->v[i1],pt->v[i2]);

          {
            printf("Saving mesh...\n");
            if ( !MMG2_hashTria(mesh) ) {
              fprintf(stdout,"  ## Hashing problem. Exit program.\n");
              return(0);
            }

            MMG2_bdryEdge(mesh);
            _MMG2_savemesh_db(mesh,mesh->nameout,0);
            return(0);
          }

          return(0);
        }
      }
    }
  }
  
  /* Check consistency between REF, GEO and BDY tags between edges and points */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    
    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_GEO || pt->tag[i] & MG_REF ) {
        i1 = _MMG5_inxt2[i];
        i2 = _MMG5_iprv2[i];
        
        if ( !(pt->tag[i] & MG_BDY) ) {
          printf("edge %d %d is tagged %d, but not MG_BDY\n",pt->v[i1],pt->v[i2],pt->tag[i]);
          return(0);
        }
        
        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt->v[i2]];
        
        if ( !(p1->tag & MG_BDY) ) {
          printf("edge %d %d is tagged %d, but point %d is not tagged BDY\n",pt->v[i1],pt->v[i2],pt->tag[i],pt->v[i1]);
          return(0);
        }
        
        if ( !(p2->tag & MG_BDY) ) {
          printf("edge %d %d is tagged %d, but point %d is not tagged BDY\n",pt->v[i1],pt->v[i2],pt->tag[i],pt->v[i2]);
          return(0);
        }
      }
    }
  }
  
  return(1);
}

/* Check orientation of elements in the mesh */
int _MMG2_chkor(MMG5_pMesh mesh) {
  MMG5_pTria        pt;
  MMG5_pPoint       p0,p1,p2;
  double            det;
  int               k;
  
  for (k=1; k<=mesh->np; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;
    
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    
    det = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);
    //if( _MMG2_caltri_iso(mesh,NULL,pt) <= 0.001) return(0);
    if ( det <= 0.0 ) return(0);
  }
  
  return(1);
}
