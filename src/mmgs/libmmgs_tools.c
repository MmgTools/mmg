/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
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
 * \file mmgs/libmmgs_tools.c
 * \brief Tools functions for the mmgs library
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmgs.h"
#include "inlined_functions.h"

void MMGS_setfunc(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( met->size < 6 ) {
    MMG5_calelt      = MMG5_caltri_iso;
    MMG5_lenSurfEdg  = MMG5_lenSurfEdg_iso;
    MMG5_compute_meanMetricAtMarkedPoints = MMG5_compute_meanMetricAtMarkedPoints_iso;
    MMGS_defsiz      = MMGS_defsiz_iso;
    MMGS_gradsiz     = MMG5_gradsiz_iso;
    MMGS_gradsizreq  = MMG5_gradsizreq_iso;
    intmet           = intmet_iso;
    movintpt         = movintpt_iso;
    movridpt         = movridpt_iso;
  }
  else {
    if ( (!met->m) && (!mesh->info.optim) && mesh->info.hsiz<=0. ) {
      MMG5_calelt     = MMG5_caltri_iso;
      MMG5_lenSurfEdg = MMG5_lenSurfEdg_iso;
    }
    else {
      MMG5_calelt     = MMG5_caltri_ani;
      MMG5_lenSurfEdg = MMG5_lenSurfEdg_ani;
    }
    MMG5_compute_meanMetricAtMarkedPoints = MMG5_compute_meanMetricAtMarkedPoints_ani;
    MMGS_defsiz      = MMGS_defsiz_ani;
    MMGS_gradsiz     = MMGS_gradsiz_ani;
    MMGS_gradsizreq  = MMG5_gradsizreq_ani;
    intmet        = intmet_ani;
    movintpt      = movintpt_ani;
    movridpt      = movridpt_ani;
  }
}

int MMGS_usage(char *prog) {

  /* Common generic options, file options and mode options */
  MMG5_mmgUsage(prog);

  /* Common parameters (first section) */
  MMG5_paramUsage1( );

  /* Specific options */
  fprintf(stdout,"-keep-ref    preserve initial domain references in level-set mode.\n");

#ifdef USE_SCOTCH
  fprintf(stdout,"-rn [n]      Turn on or off the renumbering using SCOTCH [0/1] \n");
#endif

  fprintf(stdout,"\n");

  /* Common parameters (second section) */
  MMG5_paramUsage2();

  /* Common options for advanced users */
  MMG5_advancedUsage();

  fprintf(stdout,"\n\n");

  return 1;
}

int MMGS_defaultValues(MMG5_pMesh mesh) {

  MMG5_mmgDefaultValues(mesh);
#ifdef USE_SCOTCH
  fprintf(stdout,"SCOTCH renumbering                  : enabled\n");
#else
  fprintf(stdout,"SCOTCH renumbering                  : disabled\n");
#endif
  fprintf(stdout,"\n\n");

  return 1;
}

// In ls mode : metric must be provided using -met option (-sol or default is the ls).
// In adp mode : -sol or -met or default allow to store the metric.
int MMGS_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  MMG5_pSol tmp = NULL;
  int    i;
  char   namein[MMG5_FILESTR_LGTH];

  /* First step: search if user want to see the default parameters values. */
  for ( i=1; i< argc; ++i ) {
    if ( !strcmp(argv[i],"-val") ) {
      MMGS_defaultValues(mesh);
      return 0;
    }
  }

  /* Second step: read all other arguments. */
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        MMGS_usage(argv[0]);
        return 0;
        break;
      case 'a': /* ridge angle */
        if ( !strcmp(argv[i],"-ar") && ++i < argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_angleDetection,
                                    atof(argv[i])) )
            return 0;
        }
        break;
      case 'A': /* anisotropy */
        if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor) )
          return 0;
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hmin,
                                    atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hmax,
                                    atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hsiz") && ++i < argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hsiz,
                                    atof(argv[i])) )
            return 0;

        }
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hausd,
                                    atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hgradreq") && ++i <= argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hgradreq,
                                    atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hgrad,
                                    atof(argv[i])) )
            return 0;
        }
        else {
          MMGS_usage(argv[0]);
          return 0;
        }
        break;
      case 'd':
        if ( !strcmp(argv[i],"-default") ) {
          mesh->mark=1;
        }
        else {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_debug,1) )
            return 0;
        }
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMGS_Set_inputMeshName(mesh, argv[i]) )
              return 0;

            if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,5) )
              return 0;
          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-isoref") && ++i <= argc ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_isoref,
                                    atoi(argv[i])) )
            return 0;
        }
        else {
          MMGS_usage(argv[0]);
          return 0;
        }
        break;
      case 'k':
        if ( !strcmp(argv[i],"-keep-ref") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_keepRef,1) )
            return 0;
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-ls") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_iso,1) )
            return 0;
          if ( ++i < argc && (isdigit(argv[i][0]) ||
                              (argv[i][0]=='-' && isdigit(argv[i][1])) ) ) {
            if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_ls,atof(argv[i])) )
              return 0;
          }
          else i--;
        }
        break;
      case 'm':
        if ( !strcmp(argv[i],"-met") ) {
          if ( !met ) {
            fprintf(stderr,"No metric structure allocated for %c%c%c option\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            return 0;
          }
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMGS_Set_inputSolName(mesh,met,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        else if  (!strcmp(argv[i],"-m") ) {
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_mem,atoi(argv[i])) )
            return 0;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          MMGS_usage(argv[0]);
          return 0;
        }
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_angle,0) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nsd") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_numsubdomain,atoi(argv[i])) )
              return 0;
          }
          else {
            fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_noswap,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_noinsert,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_nomove,1) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nreg") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_nreg,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nosizreq") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_nosizreq,1) ) {
            return 0;
          }
        }
        break;
      case 'o':
        if ( (!strcmp(argv[i],"-out")) || (!strcmp(argv[i],"-o")) ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            if ( !MMGS_Set_outputMeshName(mesh,argv[i]) )
              return 0;
          }else{
            fprintf(stderr,"Missing filname for %c%c%c\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        else if( !strcmp(argv[i],"-optim") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_optim,1) )
            return 0;
        }
        else {
          fprintf(stderr,"Unrecognized option %s\n",argv[i]);
          MMGS_usage(argv[0]);
          return 0;
        }
        break;
#ifdef USE_SCOTCH
      case 'r':
        if ( !strcmp(argv[i],"-rn") ) {
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {
              if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_renum,atoi(argv[i])) )
                return 0;
            }
            else {
              fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
              MMGS_usage(argv[0]);
              return 0;
            }
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        break;
#endif
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          /* For retrocompatibility, store the metric if no sol structure available */
          tmp = sol ? sol : met;
          assert(tmp);
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMGS_Set_inputSolName(mesh,tmp,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) ) {
            if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,atoi(argv[i])) )
              return 0;
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          MMGS_usage(argv[0]);
          return 0;
        }
        break;
      default:
        fprintf(stderr,"Unrecognized option %s\n",argv[i]);
        MMGS_usage(argv[0]);
        return 0;
      }
    }
    else {
      if ( mesh->namein == NULL ) {
        if ( !MMGS_Set_inputMeshName(mesh,argv[i]) )
          return 0;
        if ( mesh->info.imprim == -99 ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,5) )
            return 0;
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMGS_Set_outputMeshName(mesh,argv[i]) )
          return 0;
      }
      else {
        fprintf(stdout,"Argument %s ignored\n",argv[i]);
        MMGS_usage(argv[0]);
        return 0;
      }
    }
    i++;
  }

  /* check file names */
  if ( mesh->info.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    MMG_FSCANF(stdin,"%d",&i);
    if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,i) )
      return 0;
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    MMG_FSCANF(stdin,"%127s",namein);
    if ( !MMGS_Set_inputMeshName(mesh,namein) )
      return 0;
  }

  if ( mesh->nameout == NULL ) {
    if ( !MMGS_Set_outputMeshName(mesh,"") )
      return 0;
  }

  /* adp mode: if the metric name has been stored in sol, move it in met */
  if ( met->namein==NULL && sol && sol->namein && !(mesh->info.iso || mesh->info.lag>=0) ) {
    if ( !MMGS_Set_inputSolName(mesh,met,sol->namein) )
      return 0;
    MMG5_DEL_MEM(mesh,sol->namein);
  }

  /* default : store solution name in iso mode, metric name otherwise */
  tmp = ( mesh->info.iso || mesh->info.lag >=0 ) ? sol : met;
  assert ( tmp );
  if ( tmp->namein == NULL ) {
    if ( !MMGS_Set_inputSolName(mesh,tmp,"") ) { return 0; }
  }
  if ( met->nameout == NULL ) {
    if ( !MMGS_Set_outputSolName(mesh,met,"") )
      return 0;
  }
  return 1;
}

int MMGS_freeLocalPar(MMG5_pMesh mesh) {

  free(mesh->info.par);
  mesh->info.npar = 0;

  return 1;
}

int MMGS_stockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(&mesh->info,info,sizeof(MMG5_Info));
  MMGS_memOption(mesh);
  if( mesh->info.mem > 0) {
    if ( mesh->npmax < mesh->np || mesh->ntmax < mesh->nt ) {
      return 0;
    } else if(mesh->info.mem < 39)
      return 0;
  }
  return 1;
}

void MMGS_destockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(info,&mesh->info,sizeof(MMG5_Info));
  return;
}

int MMGS_Get_numberOfNonBdyEdges(MMG5_pMesh mesh, int* nb_edges) {
  MMG5_pTria pt,pt1;
  MMG5_pEdge ped;
  int        *adja,k,i,j,i1,i2,iel;

  *nb_edges = 0;
  if ( mesh->tria ) {
    /* Create the triangle adjacency if needed */
    if ( !mesh->adja ) {
      if ( !MMGS_hashTria(mesh) ) {
        fprintf(stderr,"\n  ## Error: %s: unable to create "
                "adjacency table.\n",__func__);
        return 0;
      }
    }

    /* Count the number of non boundary edges */
    for ( k=1; k<=mesh->nt; k++ ) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;

      adja = &mesh->adja[3*(k-1)+1];

      for ( i=0; i<3; i++ ) {
        /* Do not treat boundary edges */
        if ( MG_EDG(pt->tag[i]) ) continue;

        iel = adja[i] / 3;
        assert ( iel != k );

        pt1 = &mesh->tria[iel];

        if ( (!iel) || (pt->ref != pt1->ref) ) {
          /* Do not treat boundary edges */
          continue;
        }
        if ( k < iel ) {
          /* Treat edge from the triangle with lowest index */
          ++(*nb_edges);
        }
      }
    }

    /* Append the non boundary edges to the boundary edges array */
    if ( mesh->na ) {
      MMG5_ADD_MEM(mesh,(*nb_edges)*sizeof(MMG5_Edge),"non boundary edges",
                   printf("  Exit program.\n");
                   return 0);
      MMG5_SAFE_RECALLOC(mesh->edge,(mesh->na+1),(mesh->na+(*nb_edges)+1),
                         MMG5_Edge,"non bdy edges arrray",return 0);
    }
    else {
      MMG5_ADD_MEM(mesh,((*nb_edges)+1)*sizeof(MMG5_Edge),"non boundary edges",
                   printf("  Exit program.\n");
                   return 0);
      MMG5_SAFE_RECALLOC(mesh->edge,0,(*nb_edges)+1,
                         MMG5_Edge,"non bdy edges arrray",return 0);
    }

    j = mesh->na+1;
    for ( k=1; k<=mesh->nt; k++ ) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;

      adja = &mesh->adja[3*(k-1)+1];

      for ( i=0; i<3; i++ ) {
        /* Do not treat boundary edges */
        if ( MG_EDG(pt->tag[i]) ) continue;

        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        iel = adja[i] / 3;
        assert ( iel != k );

        pt1 = &mesh->tria[iel];

        if ( (!iel) || (pt->ref != pt1->ref) ) {
          /* Do not treat boundary edges */
          continue;
        }
        if ( k < iel ) {
          /* Treat edge from the triangle with lowest index */
          ped = &mesh->edge[j++];
          assert ( ped );
          ped->a   = pt->v[i1];
          ped->b   = pt->v[i2];
          ped->ref = pt->edg[i];
        }
      }
    }
  }
  return 1;
}

int MMGS_Get_nonBdyEdge(MMG5_pMesh mesh, int* e0, int* e1, int* ref, int idx) {
  MMG5_pEdge ped;
  size_t     na_tot=0;
  char       *ptr_c = (char*)mesh->edge;

  if ( !mesh->edge ) {
    fprintf(stderr,"\n  ## Error: %s: edge array is not allocated.\n"
            " Please, call the MMGS_Get_numberOfNonBdyEdges function"
            " before the %s one.\n",
            __func__,__func__);
    return 0;
  }

  ptr_c = ptr_c-sizeof(size_t);
  na_tot = *((size_t*)ptr_c);

  if ( mesh->namax==na_tot ) {
    fprintf(stderr,"\n  ## Error: %s: no internal edge.\n"
            " Please, call the MMGS_Get_numberOfNonBdyEdges function"
            " before the %s one and check that the number of internal"
            " edges is non null.\n",
            __func__,__func__);
  }

  if ( mesh->namax+idx > na_tot ) {
    fprintf(stderr,"\n  ## Error: %s: Can't get the internal edge of index %d."
            " Index must be between 1 and %zu.\n",
            __func__,idx,na_tot-mesh->namax);
  }

  ped = &mesh->edge[mesh->na+idx];

  *e0  = ped->a;
  *e1  = ped->b;

  if ( ref != NULL ) {
    *ref = mesh->edge[mesh->na+idx].ref;
  }

  return 1;
}

int MMGS_Get_adjaTri(MMG5_pMesh mesh, int kel, int listri[3]) {

  if ( ! mesh->adja ) {
    if (! MMGS_hashTria(mesh))
      return 0;
  }

  listri[0] = mesh->adja[3*(kel-1)+1]/3;
  listri[1] = mesh->adja[3*(kel-1)+2]/3;
  listri[2] = mesh->adja[3*(kel-1)+3]/3;

  return 1;
}

int MMGS_Get_adjaVerticesFast(MMG5_pMesh mesh, int ip,int start, int lispoi[MMGS_LMAX])
{
  MMG5_pTria pt;
  int k,prevk,nbpoi,iploc,i,i1,i2,*adja;

  pt   = &mesh->tria[start];

  for ( iploc=0; iploc<3; ++iploc ) {
    if ( pt->v[iploc] == ip ) break;
  }

  assert(iploc!=3);

  k = start;
  i = iploc;
  nbpoi = 0;
  do {
    if ( nbpoi == MMGS_LMAX ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to compute adjacent"
              " vertices of the vertex %d:\nthe ball of point contain too many"
              " elements.\n",__func__,ip);
      return 0;
    }
    i1 = MMG5_inxt2[i];
    lispoi[nbpoi] = mesh->tria[k].v[i1];
    ++nbpoi;

    adja = &mesh->adja[3*(k-1)+1];
    prevk = k;
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = MMG5_inxt2[i];
  }
  while ( k && k != start );

  if ( k > 0 ) return nbpoi;

  /* store the last point of the boundary triangle */
  if ( nbpoi == MMGS_LMAX ) {
    fprintf(stderr,"\n  ## Warning: %s: unable to compute adjacent vertices of the"
            " vertex %d:\nthe ball of point contain too many elements.\n",
            __func__,ip);
    return 0;
  }
  i1 = MMG5_inxt2[i1];
  lispoi[nbpoi] = mesh->tria[prevk].v[i1];
  ++nbpoi;

  /* check if boundary hit */
  k = start;
  i = iploc;
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i2 = MMG5_iprv2[i];
    k  = adja[i2] / 3;
    if ( k == 0 )  break;

    if ( nbpoi == MMGS_LMAX ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to compute adjacent vertices of the"
              " vertex %d:\nthe ball of point contain too many elements.\n",
              __func__,ip);
      return 0;
    }
    i  = adja[i2] % 3;
    lispoi[nbpoi] = mesh->tria[k].v[i];
    ++nbpoi;

    i  = MMG5_iprv2[i];
  }
  while ( k );

  return nbpoi;
}

int MMGS_doSol(MMG5_pMesh mesh,MMG5_pSol met) {
    MMG5_pTria   ptt;
    MMG5_pPoint  p1,p2;
    double       ux,uy,uz,dd;
    int          i,k,iadr,ipa,ipb,type;
    int          *mark;

    MMG5_SAFE_CALLOC(mark,mesh->np+1,int,return 0);

    /* Memory alloc */
    if ( met->size==1 ) type=1;
    else if ( met->size==6 ) type = 3;
    else {
      fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
              __func__,met->size);
      return 0;
    }

    if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
      return 0;

    /* edges */
    dd = 0.;
    for (k=1; k<=mesh->nt; k++) {
        ptt = &mesh->tria[k];
        if ( !MG_EOK(ptt) )  continue;

        if ( met->size == 1 ) {
          for (i=0; i<3; i++) {
            ipa = ptt->v[i];
            ipb = ptt->v[MMG5_inxt2[i]];
            p1  = &mesh->point[ipa];
            p2  = &mesh->point[ipb];

            ux  = p1->c[0] - p2->c[0];
            uy  = p1->c[1] - p2->c[1];
            uz  = p1->c[2] - p2->c[2];
            dd  = sqrt(ux*ux + uy*uy + uz*uz);

            met->m[ipa] += dd;
            mark[ipa]++;
            met->m[ipb] += dd;
            mark[ipb]++;
          }
        }
        else if ( met->size == 6 ) {
          for (i=0; i<3; i++) {
            ipa = ptt->v[i];
            ipb = ptt->v[MMG5_inxt2[i]];
            p1  = &mesh->point[ipa];
            p2  = &mesh->point[ipb];

            ux  = p1->c[0] - p2->c[0];
            uy  = p1->c[1] - p2->c[1];
            uz  = p1->c[2] - p2->c[2];
            dd  = sqrt(ux*ux + uy*uy + uz*uz);

            iadr = 6*ipa;
            met->m[iadr]   += dd;
            mark[ipa]++;

            iadr = 6*ipb;
            met->m[iadr] += dd;
            mark[ipb]++;
          }
        }
        else {
          MMG5_SAFE_FREE(mark);
          return 0;
        }
    }

    /* if hmax is not specified, compute it from the metric */
    if ( mesh->info.hmax < 0. ) {
      if ( met->size == 1 ) {
        dd = 0.;
        for (k=1; k<=mesh->np; k++) {
          if ( !mark[k] ) continue;
          dd = MG_MAX(dd,met->m[k]);
        }
        assert ( dd );
      }
      else if ( met->size == 6 ) {
        dd = FLT_MAX;
        for (k=1; k<=mesh->np; k++) {
          if ( !mark[k] ) continue;
          iadr = 6*k;
          dd = MG_MIN(dd,met->m[iadr]);
        }
        assert ( dd < FLT_MAX );
        dd = 1./sqrt(dd);
      }
      mesh->info.hmax = 10.*dd;
    }

    /* vertex size */
    if ( met->size == 1 ) {
      for (k=1; k<=mesh->np; k++) {
        if ( !mark[k] ) {
          met->m[k] = mesh->info.hmax;
          continue;
        }
        else
          met->m[k] = met->m[k] / (double)mark[k];
      }
    }
    else if ( met->size == 6 ) {
      for (k=1; k<=mesh->np; k++) {
        iadr = 6*k;
        if ( !mark[k] ) {
          met->m[iadr]   = 1./(mesh->info.hmax*mesh->info.hmax);
          met->m[iadr+3] = met->m[iadr];
          met->m[iadr+5] = met->m[iadr];
          continue;
        }
        else {
          met->m[iadr]   = (double)mark[k]*(double)mark[k]/(met->m[iadr]*met->m[iadr]);
          met->m[iadr+3] = met->m[iadr];
          met->m[iadr+5] = met->m[iadr];
        }
      }
    }

    MMG5_SAFE_FREE(mark);
    return 1;
}

int MMGS_Set_constantSize(MMG5_pMesh mesh,MMG5_pSol met) {
  double      hsiz;
  int         type;

  /* Memory alloc */
  if ( met->size==1 ) type=1;
  else if ( met->size==6 ) type = 3;
  else {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,met->size);
    return 0;
  }
  if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
    return 0;

  if ( !MMG5_Compute_constantSize(mesh,met,&hsiz) )
    return 0;

  mesh->info.hsiz = hsiz;

  MMG5_Set_constantSize(mesh,met,hsiz);

  return 1;
}

int MMGS_Compute_eigenv(double m[6],double lambda[3],double vp[3][3]) {

  return  MMG5_eigenv3d(1,m,lambda,vp);

}

void MMGS_Free_solutions(MMG5_pMesh mesh,MMG5_pSol sol) {

  /* sol */
  if ( !sol ) return;

  if ( sol->m )
    MMG5_DEL_MEM(mesh,sol->m);

  if ( sol->namein ) {
    MMG5_DEL_MEM(mesh,sol->namein);
  }

  if ( sol->nameout ) {
    MMG5_DEL_MEM(mesh,sol->nameout);
  }

  memset ( sol, 0x0, sizeof(MMG5_Sol) );

  /* Reset state to a scalar status */
  sol->dim  = 3;
  sol->ver  = 2;
  sol->size = 1;
  sol->type = 1;

  return;
}
