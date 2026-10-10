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
 * Example of use of the mmg3d library progress callback.
 *
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include <stdio.h> /** BEGIN_EXAMPLE (this line is used by Doxygen) */
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32)
#include <io.h>
#define ISATTY(stream) _isatty(_fileno(stream))
#else
#include <unistd.h>
#define ISATTY(stream) isatty(fileno(stream))
#endif

/** Include the mmg3d library header file */
// if the header file is in the "include" directory
// #include "libmmg3d.h"
// if the header file is in "include/mmg/mmg3d"
#include "mmg/mmg3d/libmmg3d.h"

#define PROGRESS_PHASE_COUNT 4

typedef struct {
  int lastPhase;
  int barIsOpen;
  int isInteractive;
  int seen[PROGRESS_PHASE_COUNT];
  int iteration[PROGRESS_PHASE_COUNT];
  int maxIterations[PROGRESS_PHASE_COUNT];
  int percent[PROGRESS_PHASE_COUNT];
  int64_t nSplit[PROGRESS_PHASE_COUNT];
  int64_t nCollapse[PROGRESS_PHASE_COUNT];
  int64_t nSwap[PROGRESS_PHASE_COUNT];
  int64_t nMove[PROGRESS_PHASE_COUNT];
} ProgressData;

static const char *phaseName(int phase) {
  switch ( phase ) {
  case MMG5_PHASE_GEOMETRIC_MESH:
    return "Geometric Mesh";
  case MMG5_PHASE_COMPUTATIONAL_MESH:
    return "Computational Mesh";
  case MMG5_PHASE_ADAPTATION:
    return "Adaptation";
  case MMG5_PHASE_OPTIMIZATION:
    return "Optimization";
  default:
    return "Unknown";
  }
}

static int phaseIndex(int phase) {
  switch ( phase ) {
  case MMG5_PHASE_GEOMETRIC_MESH:
  case MMG5_PHASE_COMPUTATIONAL_MESH:
  case MMG5_PHASE_ADAPTATION:
  case MMG5_PHASE_OPTIMIZATION:
    return phase;
  default:
    return -1;
  }
}

static void printProgressBarLine(const char *phaseName,int iteration,
                                 int maxIterations,int percent,
                                 int64_t nSplit,int64_t nCollapse,
                                 int64_t nSwap,int64_t nMove) {
  const int width = 30;
  int currentIteration, filled, i;

  currentIteration = iteration + 1;
  filled = width * percent / 100;

  fprintf(stdout,"  %-18s [",phaseName);
  for ( i=0; i<width; ++i ) {
    fputc(i < filled ? '=' : ' ',stdout);
  }
  fprintf(stdout,"] %3d/%-3d split=%" PRId64 " collapse=%" PRId64
          " swap=%" PRId64 " move=%" PRId64,
          currentIteration,maxIterations,nSplit,nCollapse,nSwap,nMove);
}

static void clearProgressLine(void) {
  int i;

  fprintf(stdout,"\r");
  for ( i=0; i<120; ++i ) {
    fputc(' ',stdout);
  }
  fprintf(stdout,"\r");
}

static void saveProgress(ProgressData *data,int phase,int iteration,
                         int maxIterations,int percent,
                         int64_t nSplit,int64_t nCollapse,
                         int64_t nSwap,int64_t nMove) {
  int idx;

  idx = phaseIndex(phase);
  if ( idx < 0 ) {
    return;
  }

  data->seen[idx] = 1;
  data->iteration[idx] = iteration;
  data->maxIterations[idx] = maxIterations;
  data->percent[idx] = percent;
  data->nSplit[idx] = nSplit;
  data->nCollapse[idx] = nCollapse;
  data->nSwap[idx] = nSwap;
  data->nMove[idx] = nMove;
}

static void printProgressSummary(const ProgressData *data) {
  int phase;

  for ( phase=0; phase<PROGRESS_PHASE_COUNT; ++phase ) {
    if ( !data->seen[phase] ) {
      continue;
    }
    printProgressBarLine(phaseName(phase),
                         data->iteration[phase],
                         data->maxIterations[phase],
                         100,
                         data->nSplit[phase],
                         data->nCollapse[phase],
                         data->nSwap[phase],
                         data->nMove[phase]);
    fprintf(stdout,"\n");
  }
}

static int progressCallback(void *mesh,
                            int phase,
                            int iteration,
                            int maxIterations,
                            int64_t nSplit,
                            int64_t nCollapse,
                            int64_t nSwap,
                            int64_t nMove,
                            void *userData) {
  ProgressData *data;
  int percent;

  (void)mesh;

  data = (ProgressData *)userData;
  if ( maxIterations > 0 ) {
    percent = 100 * (iteration + 1) / maxIterations;
    if ( percent > 100 ) {
      percent = 100;
    }
  }
  else {
    maxIterations = iteration + 1;
    percent = 100;
  }

  if ( !data ) {
    return 1;
  }

  saveProgress(data,phase,iteration,maxIterations,percent,
               nSplit,nCollapse,nSwap,nMove);

  if ( !data->isInteractive ) {
    return 1;
  }

  if ( data->barIsOpen && data->lastPhase != phase ) {
    int lastPhaseIndex;

    lastPhaseIndex = phaseIndex(data->lastPhase);
    if ( lastPhaseIndex >= 0 ) {
      clearProgressLine();
      printProgressBarLine(phaseName(data->lastPhase),
                           data->iteration[lastPhaseIndex],
                           data->maxIterations[lastPhaseIndex],
                           100,
                           data->nSplit[lastPhaseIndex],
                           data->nCollapse[lastPhaseIndex],
                           data->nSwap[lastPhaseIndex],
                           data->nMove[lastPhaseIndex]);
    }
    fprintf(stdout,"\n");
    data->barIsOpen = 0;
  }

  data->lastPhase = phase;
  data->barIsOpen = 1;

  clearProgressLine();
  printProgressBarLine(phaseName(phase),iteration,maxIterations,percent,
                       nSplit,nCollapse,nSwap,nMove);
  fflush(stdout);

  if ( percent == 100 ) {
    fprintf(stdout,"\n");
    data->lastPhase = -1;
    data->barIsOpen = 0;
  }

  return 1;
}

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  ProgressData    progressData;
  int             ier;
  char            *filename, *fileout;

  fprintf(stdout,"  -- TEST MMG3DLIB PROGRESS CALLBACK\n");

  if ( argc != 3 ) {
    printf(" Usage: %s filein fileout \n",argv[0]);
    return(1);
  }

  filename = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filename,argv[1]);

  fileout = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[2]);

  mmgMesh = NULL;
  mmgSol  = NULL;

  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_verbose,-1) != 1 )
    exit(EXIT_FAILURE);

  if ( MMG3D_loadMesh(mmgMesh,filename) != 1 )  exit(EXIT_FAILURE);

  if ( MMG3D_loadSol(mmgMesh,mmgSol,filename) != 1 )
    exit(EXIT_FAILURE);

  if ( MMG3D_Chk_meshData(mmgMesh,mmgSol) != 1 ) exit(EXIT_FAILURE);

  memset(&progressData,0,sizeof(ProgressData));
  progressData.lastPhase = -1;
  progressData.isInteractive = ISATTY(stdout);
  if ( MMG3D_Set_progressCallback(mmgMesh,progressCallback,&progressData) != 1 )
    exit(EXIT_FAILURE);

  ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);

  if ( progressData.isInteractive && progressData.barIsOpen ) {
    fprintf(stdout,"\n");
  }
  else if ( !progressData.isInteractive ) {
    printProgressSummary(&progressData);
  }

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  if ( MMG3D_saveMesh(mmgMesh,fileout) != 1 ) {
    fprintf(stdout,"UNABLE TO SAVE MESH\n");
    return(MMG5_STRONGFAILURE);
  }

  if ( MMG3D_saveSol(mmgMesh,mmgSol,fileout) != 1 ) {
    fprintf(stdout,"UNABLE TO SAVE SOL\n");
    return(MMG5_LOWFAILURE);
  }

  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  free(filename);
  filename = NULL;

  free(fileout);
  fileout = NULL;

  return(ier);
}   /** END_EXAMPLE (this line is used by Doxygen) */
