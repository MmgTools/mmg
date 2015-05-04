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

#ifdef USE_SCOTCH

#ifndef __RENUM__
#define __RENUM__

#include <scotch.h>

#define HASHPRIME 37

#define CHECK_SCOTCH(t,m,e) if(0!=t){perror(m);return(e);}

typedef struct MeshGraphHash_ {
  int vertNum;
  int vertEnd;
} MeshGraphHash;

int    _SCOTCHintSort2asc1(SCOTCH_Num * sortPartTb, int vertNbr);
int    _MMG5_kPartBoxCompute(SCOTCH_Graph, int, int, SCOTCH_Num*,MMG5_pMesh);
void   _MMG5_swapNod(MMG5_pPoint, double*, int*, int, int, int);

#endif /* __RENUM__ */
#endif
