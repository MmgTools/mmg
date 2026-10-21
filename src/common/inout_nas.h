/* =============================================================================
**  This file is part of the Mmg software package.
**  It is distributed under the GNU Lesser General Public License, version 3
**  or any later version.
** =============================================================================
*/

/** \file inout_nas.h
 *  \brief Internal Nastran bulk-data mesh parser.
 */

#ifndef MMG_INOUT_NAS_H
#define MMG_INOUT_NAS_H

#include "libmmgtypes.h"

enum MMG5_NastranElementType {
  MMG5_NAS_CTRIA3,
  MMG5_NAS_CQUAD4,
  MMG5_NAS_CTETRA,
  MMG5_NAS_CPENTA,
  MMG5_NAS_CPYRAM,
  MMG5_NAS_CHEXA
};

typedef struct {
  MMG5_int id;
  double   c[3];
} MMG5_NastranPoint;

typedef struct {
  enum MMG5_NastranElementType type;
  MMG5_int id;
  MMG5_int ref;
  MMG5_int vertex[8];
} MMG5_NastranElement;

typedef struct {
  MMG5_NastranPoint   *points;
  MMG5_NastranElement *elements;
  MMG5_int             np;
  MMG5_int             ne;
} MMG5_NastranMesh;

int MMG5_loadNastranMeshData(const char *filename,MMG5_NastranMesh *nas);
void MMG5_freeNastranMeshData(MMG5_NastranMesh *nas);
MMG5_int MMG5_nastranPointIndex(const MMG5_NastranMesh *nas,MMG5_int id);

#endif
