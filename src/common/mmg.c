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
 * \file common/mmg.c
 * \brief Common part for functions used in mmgs.c and mmg3d.c files.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 04 2015
 * \copyright GNU Lesser General Public License.
 **/

#include "mmgcommon.h"

/**
 * \param *prog pointer toward the program name.
 *
 * Print help for common options of mmg3d and mmgs.
 *
 */
void _MMG5_mmgUsage(char *prog) {
  fprintf(stdout,"\nUsage: %s [-v [n]] [opts..] filein [fileout]\n",prog);

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-h        Print this message\n");
  fprintf(stdout,"-v [n]    Tune level of verbosity, [-10..10]\n");
  fprintf(stdout,"-m [n]    Set maximal memory size to n Mbytes\n");
  fprintf(stdout,"-d        Turn on debug mode\n");
  fprintf(stdout,"-val      Print the default parameters values\n");
  fprintf(stdout,"-default  Save a local parameters file for default parameters"
          " values\n");

  fprintf(stdout,"\n**  File specifications\n");
  fprintf(stdout,"-in  file  input triangulation\n");
  fprintf(stdout,"-out file  output triangulation\n");
  fprintf(stdout,"-sol file  load solution or metric file\n");

  fprintf(stdout,"\n**  Parameters\n");
  fprintf(stdout,"-ar     val  angle detection\n");
  fprintf(stdout,"-nr          no angle detection\n");
  fprintf(stdout,"-hmin   val  minimal mesh size\n");
  fprintf(stdout,"-hmax   val  maximal mesh size\n");
  fprintf(stdout,"-hausd  val  control Hausdorff distance\n");
  fprintf(stdout,"-hgrad  val  control gradation\n");
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
void _MMG5_mmgDefaultValues(MMG5_pMesh mesh) {
  long long memMax;

  fprintf(stdout,"\nDefault parameters values:\n");

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"verbosity                 (-v)      : %d\n",
          mesh->info.imprim);
  memMax = _MMG5_memSize();
  if ( memMax )
    /* maximal memory = 50% of total physical memory */
    memMax = memMax*50/104857600L;
  else {
    /* default value = 800 Mo */
    memMax = _MMG5_MEMMAX;
  }
  fprintf(stdout,"maximal memory size       (-m)      : %ld MBytes\n",
          _MMG5_safeLL2LCast(memMax));


  fprintf(stdout,"\n**  Parameters\n");
  fprintf(stdout,"angle detection           (-ar)     : %lf\n",
          180/M_PI*acos(mesh->info.dhd) );
  fprintf(stdout,"minimal mesh size         (-hmin)   : 0.01 of "
          "the mesh bounding box if no metric is provided, 0.1 times the "
          "minimum of the metric sizes otherwise.\n");
  fprintf(stdout,"maximal mesh size         (-hmax)   : size of "
          "the mesh bounding box without metric, 10 times the maximum of the "
          "metric sizes otherwise.\n");
  fprintf(stdout,"Hausdorff distance        (-hausd)  : %lf\n",
          mesh->info.hausd);
  fprintf(stdout,"gradation control         (-hgrad)  : %lf\n",
          exp(mesh->info.hgrad));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param out pointer toward the file in which to write.
 * \return 1 if success, 0 otherwise.
 *
 * Write a default file containing the default values of parameters that
 * can be locally defined.
 *
 */
inline
int _MMG5_writeLocalParam( MMG5_pMesh mesh, FILE *out ) {
  _MMG5_Node *bdyRefs, *cur;
  double     dd;
  int        npar,k,ier;

  /** Count the number of different boundary references and list it */
  bdyRefs = NULL;
  npar = 0;

  k = mesh->nt? mesh->tria[1].ref : 0;

  /* Try to alloc the first node */
  ier = _MMG5_Add_node( mesh, &bdyRefs, k );
  if ( ier < 0 ) {
    printf("  ## Error: unable to allocate the first boundary"
           " reference node.\n");
    return(0);
  }
  else {
    assert(ier);
    npar = 1;
  }

  for ( k=1; k<=mesh->nt; ++k ) {
    ier = _MMG5_Add_node( mesh, &bdyRefs, mesh->tria[k].ref );

    if ( ier < 0 ) {
      printf("  ## Warning: unable to list the boundary references.\n"
             "              Uncomplete parameters file.\n" );
      break;
    }
    else if ( ier ) ++npar;
  }

  fprintf(out,"parameters\n %d\n",npar);

  dd = mesh->info.delta;
  cur = bdyRefs;
  while( cur ) {
    fprintf(out,"%d Triangle %e %e %e \n",cur->val,
            mesh->info.hmin, mesh->info.hmax,mesh->info.hausd);
    cur = cur->nxt;
  }

  _MMG5_Free_linkedList(mesh,bdyRefs);

  return(1);
}
