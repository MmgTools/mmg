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
#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <assert.h>

/* prototype (re)definitions */
void  *M_malloc(size_t size,char *call);
void  *M_calloc(size_t nelem,size_t elsize,char *call);
void  *M_realloc(void *ptr, size_t size,char *call);
void   M_free(void *ptr);

/* ptototypes : tools */
int    M_memLeak();
void   M_memDump();
size_t M_memSize();


#ifdef __cplusplus
}
#endif
