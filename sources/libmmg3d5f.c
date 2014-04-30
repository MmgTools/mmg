 /**
 * Written by Cecile Dobrzynski (IMB), Charles Dapogny,
 * Pascal Frey (LJLL) and Algiane Froehly
 * Copyright (c) 2004- IMB/LJLL.
 * All rights reserved.
 *
 * Interfaces for fortran call
 *
 * MMG5_mmg3dlib(MMG5_pMesh mesh,MMG5_pSol met [,MMG5_pSingul singul] ):
 *    to use mmg3d via a library
 *
 * Integers parameters:
 *    MMG5_IPARAM_verbose            = [-10..10] , Tune level of verbosity;
 *    MMG5_IPARAM_mem                = [n/-1]    , Set maximal memory size to n Mbytes/keep the default value;
 *    MMG5_IPARAM_debug              = [1/0]     , Turn on/off debug mode;
 *    MMG5_IPARAM_angle              = [1/0]     , Turn on/off angle detection;
 *    MMG5_IPARAM_iso                = [1/0]     , Turn on/off levelset meshing;
 *    MMG5_IPARAM_noinsert           = [1/0]     , avoid/allow point insertion/deletion;
 *    MMG5_IPARAM_noswap             = [1/0]     , avoid/allow edge or face flipping;
 *    MMG5_IPARAM_nomove             = [1/0]     , avoid/allow point relocation;
 *    MMG5_IPARAM_numberOflocalParam = [n]       , number of local parameters;
 *    MMG5_IPARAM_renum              = [1/0]     , Turn on/off the renumbering using SCOTCH;
 *    MMG5_IPARAM_sing               = [1/0]     , Turn on/off the insertion of singularities
 *                                        (need to compile with -DSINGUL flag);
 * Double parameters:
 *    MMG5_DPARAM_dhd   = [val]     , angle detection;
 *    MMG5_DPARAM_hmin  = [val]     , minimal mesh size;
 *    MMG5_DPARAM_hmax  = [val]     , maximal mesh size;
 *    MMG5_DPARAM_hausd = [val]     , control global Hausdorff distance
 *                                    (on all the boundary surfaces of the mesh);
 *    MMG5_DPARAM_hgrad = [val]     , control gradation;
 *    MMG5_DPARAM_ls    = [val]     , level set value;
 **/

#include "libmmg3d5.h"

/** Macro from Scotch **/
#define FORTRAN_NAME(nu,nl,pl,pc) \
  void nu pl;			  \
  void nl pl			  \
  { nu pc; }			  \
  void nl##_ pl			  \
  { nu pc; }			  \
  void nl##__ pl		  \
  { nu pc; }			  \
  void nu pl

/** Deallocations before return */
FORTRAN_NAME(MMG5_FREE_ALL,mmg5_free_all,(MMG5_pMesh *mesh,MMG5_pSol *met
#ifdef SINGUL
				,MMG5_pSingul *singul
#endif
				),(mesh,met
#ifdef SINGUL
				   ,singul
#endif
				   )){

#ifdef SINGUL
  MMG5_Free_all(*mesh,*met,*singul);
#else
  MMG5_Free_all(*mesh,*met);
#endif
  return;
}

FORTRAN_NAME(MMG5_SAVEMESH,mmg5_savemesh,(MMG5_pMesh *mesh, int* retval),
	     (mesh,retval)){
  *retval = MMG5_saveMesh(*mesh);
  return;
}

/** main programm */
FORTRAN_NAME(MMG5_MMG3DLIB,mmg5_mmg3dlib,(MMG5_pMesh *mesh,MMG5_pSol *met
#ifdef SINGUL
				,MMG5_pSingul *sing
#endif
				,int* retval),(mesh,met
#ifdef SINGUL
					       ,sing
#endif
					       ,retval)){

#ifdef SINGUL
  *retval = MMG5_mmg3dlib(*mesh,*met,*sing);
#else
  *retval = MMG5_mmg3dlib(*mesh,*met);
#endif
  return;
}
