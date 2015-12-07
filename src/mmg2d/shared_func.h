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

/**
 * \file mmg2d/shared_func.h
 * \brief Common functions between MMG3D5 library and executable.
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

/* global variables */
int MMG2_iare[3][2] = {{1,2},{2,0},{0,1}};
int MMG2_iopp[3][2] = {{1,2},{0,2},{0,1}};
unsigned int MMG2_idir[5] = {0,1,2,0,1};

unsigned char _MMG5_iprv2[3] = {2,0,1};
unsigned char _MMG5_inxt2[3] = {1,2,0};

mytime   ctim[TIMEMAX];

/**
 * Set common pointer functions between mmgs and mmg3d to the matching mmg3d
 * functions.
 */
void _MMG2_Set_commonFunc() {
  MMG5_Init_parameters    = _MMG2_Init_parameters;
  _MMG5_chkmsh            = _MMG5_mmg2dChkmsh;
  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * Set function pointers.
 *
 */
int _MMG2_setfunc(int type) {
  if ( type == 3 ) {
    MMG2_length    = long_ani;
    MMG2_caltri    = caltri_ani;
    MMG2_caltri_in = caltri_ani_in;
    MMG2_buckin    = buckin_ani;
    MMG2_lissmet   = lissmet_ani;
    MMG2_optlen    = optlen_ani;
/*    interp     = interp_ani;
 */
  }
  else {
    MMG2_length     = long_iso;
    MMG2_caltri     = caltri_iso;
    MMG2_caltri_in  = caltri_iso_in;
    MMG2_buckin     = buckin_iso;
    MMG2_lissmet    = lissmet_iso;

    MMG2_optlen     = optlen_iso;
/*    interp     = interp_iso;
 */
  }

  return(1);
}

