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
 * Print help for common options of the 3 codes (first section).
 *
 */
void MMG5_mmgUsage(char *prog) {
  fprintf(stdout,"\nUsage: %s [-v [n]] [opts..] filein [fileout]\n",prog);

  fprintf(stdout,"\n** Generic options\n");
  fprintf(stdout,"-h        Print this message\n");
  fprintf(stdout,"-v [n]    Tune level of verbosity, [-1..10]\n");
  fprintf(stdout,"-m [n]    Set maximal memory size to n Mbytes\n");
  fprintf(stdout,"-d        Turn on debug mode\n");
  fprintf(stdout,"-val      Print the default parameters values\n");
  fprintf(stdout,"-default  Save a local parameters file for default parameters"
          " values\n");

  fprintf(stdout,"\n**  File specifications\n");
  fprintf(stdout,"-in  file  input triangulation\n");
  fprintf(stdout,"-out file  output triangulation\n");
  fprintf(stdout,"-sol file  load solution or metric file\n");
  fprintf(stdout,"-met file  load metric file\n");

  fprintf(stdout,"\n**  Mode specifications (mesh adaptation by default)\n");
  fprintf(stdout,"-ls  val     create mesh of isovalue val (0 if no argument provided)\n");

}

/**
 *
 * Print help for common parameters options of the 3 codes (first section).
 *
 */
void MMG5_paramUsage1(void) {
  fprintf(stdout,"\n**  Parameters\n");
  fprintf(stdout,"-A           enable anisotropy (without metric file).\n");
  fprintf(stdout,"-ar     val  angle detection\n");
  fprintf(stdout,"-nr          no angle detection\n");
  fprintf(stdout,"-hausd  val  control Hausdorff distance\n");
  fprintf(stdout,"-hgrad  val  control gradation\n");
  fprintf(stdout,"-hmax   val  maximal mesh size\n");
  fprintf(stdout,"-hmin   val  minimal mesh size\n");
  fprintf(stdout,"-hsiz   val  constant mesh size\n");
}

/**
 *
 * Print help for common options of the 3 codes (second section).
 *
 */
void MMG5_paramUsage2(void) {

  fprintf(stdout,"-noinsert    no point insertion/deletion \n");
  fprintf(stdout,"-nomove      no point relocation\n");
  fprintf(stdout,"-noswap      no edge or face flipping\n");
  fprintf(stdout,"-nreg        normal regul.\n");
  fprintf(stdout,"-nsd    val  save the subdomain number val (0==all subdomain)\n");
  fprintf(stdout,"-optim       mesh optimization\n");

}

/**
 *
 * Print help for lagrangian motion option.
 *
 */
void MMG5_lagUsage(void) {

#ifdef USE_ELAS
  fprintf(stdout,"-lag [n]     lagrangian mesh displacement according to mode [0/1/2]\n");
  fprintf(stdout,"               0: displacement\n");
  fprintf(stdout,"               1: displacement + remeshing (swap and move)\n");
  fprintf(stdout,"               2: displacement + remeshing (split, collapse,"
          " swap and move)\n");
#endif
}

/**
 *
 * Print help for common options between 2D and 3D.
 *
 */
void MMG5_2d3dUsage(void) {

  fprintf(stdout,"-opnbdy      preserve input triangles at the interface of"
          " two domains of the same reference.\n");

  fprintf(stdout,"-rmc   [val] enable the removal of componants whose volume fraction is less than\n"
          "             val (1e-5 if not given) of the mesh volume (ls mode).\n");
}

/**
 *
 * Print help for advanced users of mmg.
 *
 */
void MMG5_advancedUsage(void) {

  fprintf(stdout,"\n**  Parameters for advanced users\n");
  fprintf(stdout,"-nosizreq       disable setting of required edge sizes over required vertices.\n");
  fprintf(stdout,"-hgradreq  val  control gradation from required entities toward others\n");

}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
void MMG5_mmgDefaultValues(MMG5_pMesh mesh) {

  fprintf(stdout,"\nDefault parameters values:\n");

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"verbosity                 (-v)      : %d\n",
          mesh->info.imprim);
  fprintf(stdout,"maximal memory size       (-m)      : %zu MB\n",
          mesh->memMax/MMG5_MILLION);


  fprintf(stdout,"\n**  Parameters\n");
  fprintf(stdout,"angle detection           (-ar)     : %lf\n",
          180/M_PI*acos(mesh->info.dhd) );
  fprintf(stdout,"minimal mesh size         (-hmin)   : %lf\n"
          "If not yet computed: 0.001 of "
          "the mesh bounding box if no metric is provided, 0.1 times the "
          "minimum of the metric sizes otherwise.\n",mesh->info.hmin);
  fprintf(stdout,"maximal mesh size         (-hmax)   : %lf\n"
          " If not yet computed: size of "
          "the mesh bounding box without metric, 10 times the maximum of the "
          "metric sizes otherwise.\n",mesh->info.hmax);
  fprintf(stdout,"Hausdorff distance        (-hausd)  : %lf\n",
          mesh->info.hausd);

  fprintf(stdout,"gradation control         (-hgrad)  : %lf\n",
          (mesh->info.hgrad < 0) ? mesh->info.hgrad : exp(mesh->info.hgrad) );

  fprintf(stdout,"gradation control for required entities (-hgradreq)  : %lf\n",
          (mesh->info.hgradreq < 0) ? mesh->info.hgradreq : exp(mesh->info.hgradreq) );
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param bdryRefs pointer toward the list of the boundary references.
 * \return npar, the number of local parameters at triangles if success,
 * 0 otherwise.
 *
 * Count the local default values at triangles and fill the list of the boundary
 * references.
 *
 */
inline
int MMG5_countLocalParamAtTri( MMG5_pMesh mesh,MMG5_iNode **bdryRefs) {
  int         npar,k,ier;

  /** Count the number of different boundary references and list it */
  (*bdryRefs) = NULL;

  k = mesh->nt? mesh->tria[1].ref : 0;

  /* Try to alloc the first node */
  ier = MMG5_Add_inode( mesh, bdryRefs, k );
  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate the first boundary"
           " reference node.\n",__func__);
    return 0;
  }
  else {
    assert(ier);
    npar = 1;
  }

  for ( k=1; k<=mesh->nt; ++k ) {
    ier = MMG5_Add_inode( mesh, bdryRefs, mesh->tria[k].ref );

    if ( ier < 0 ) {
      printf("  ## Warning: %s: unable to list the tria references."
             " Uncomplete parameters file.\n",__func__ );
      break;
    }
    else if ( ier ) ++npar;
  }

  return npar;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param bdryRefs pointer toward the list of the boundary references.
 * \param out pointer toward the file in which to write.
 * \return 1 if success, 0 otherwise.
 *
 * Write the local default values at triangles in the parameter file.
 *
 */
inline
int MMG5_writeLocalParamAtTri( MMG5_pMesh mesh, MMG5_iNode *bdryRefs,
                                FILE *out ) {
  MMG5_iNode *cur;

  cur = bdryRefs;
  while( cur ) {
    fprintf(out,"%d Triangle %e %e %e \n",cur->val,
            mesh->info.hmin, mesh->info.hmax,mesh->info.hausd);
    cur = cur->nxt;
  }

  MMG5_Free_ilinkedList(mesh,bdryRefs);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the solution structure.
 *
 * Truncate the metric computed by the DoSol function by hmax and hmin values
 * (if setted by the user). Set hmin and hmax if they are not setted.
 *
 * \warning works only for a metric computed by the DoSol function because we
 * suppose that we have a diagonal tensor in aniso.
 *
 */
void MMG5_solTruncatureForOptim(MMG5_pMesh mesh, MMG5_pSol met) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  double      isqhmin, isqhmax;
  int         i,k,iadr,sethmin,sethmax;

  assert ( mesh->info.optim );

  /* Detect the point used only by prisms */
  if ( mesh->nprism ) {
    for (k=1; k<=mesh->np; k++) {
      mesh->point[k].flag = 1;
    }
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;

      for (i=0; i<4; i++) {
        mesh->point[pt->v[i]].flag = 0;
      }
    }
  }

  /* If not provided by the user, compute hmin/hmax from the metric computed by
   * the DoSol function. */
  sethmin = sethmax = 1;
  if ( mesh->info.hmin < 0 ) {
    sethmin = 0;
    if ( met->size == 1 ) {
      mesh->info.hmin = FLT_MAX;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        mesh->info.hmin = MG_MIN(mesh->info.hmin,met->m[k]);
      }
    }
    else if ( met->size == 6 ){
      mesh->info.hmin = 0.;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        iadr = met->size*k;
        mesh->info.hmin = MG_MAX(mesh->info.hmin,met->m[iadr]);
        mesh->info.hmin = MG_MAX(mesh->info.hmin,met->m[iadr+3]);
        mesh->info.hmin = MG_MAX(mesh->info.hmin,met->m[iadr+5]);
      }
      mesh->info.hmin = 1./sqrt(mesh->info.hmin);
    }
  }

  if ( mesh->info.hmax < 0 ) {
    sethmax = 0;
    if ( met->size == 1 ) {
      mesh->info.hmax = 0.;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        mesh->info.hmax = MG_MAX(mesh->info.hmax,met->m[k]);
      }
    }
    else if ( met->size == 6 ){
      mesh->info.hmax = FLT_MAX;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        iadr = met->size*k;
        mesh->info.hmax = MG_MIN(mesh->info.hmax,met->m[iadr]);
        mesh->info.hmax = MG_MIN(mesh->info.hmax,met->m[iadr+3]);
        mesh->info.hmax = MG_MIN(mesh->info.hmax,met->m[iadr+5]);
      }
      mesh->info.hmax = 1./sqrt(mesh->info.hmax);
    }
  }

  /* Check the compatibility between the user settings and the automatically
   * computed values */
  MMG5_check_hminhmax(mesh,sethmin,sethmax);

  /* vertex size */
  if ( met->size == 1 ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      met->m[k] = MG_MIN(mesh->info.hmax,MG_MAX(mesh->info.hmin,met->m[k]));
    }
  }
  else if ( met->size == 6 ) {
    isqhmin = 1./(mesh->info.hmin*mesh->info.hmin);
    isqhmax = 1./(mesh->info.hmax*mesh->info.hmax);

    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      iadr = 6*k;
      met->m[iadr]   = MG_MAX(isqhmax,MG_MIN(isqhmin,met->m[iadr]));
      met->m[iadr+3] = met->m[iadr];
      met->m[iadr+5] = met->m[iadr];
    }
  }

  return;
}

/**
 * \param filename string containing a filename
 *
 * \return pointer toward the filename extension or toward the end of the string
 * if no extension have been founded
 *
 * Get the extension of the filename string. Do not consider '.o' as an extension.
 *
 */
char *MMG5_Get_filenameExt( char *filename ) {
  const char pathsep='/';
  char       *dot,*lastpath;

  if ( !filename ) {
    return NULL;
  }

  dot = strrchr(filename, '.');
  lastpath = (pathsep == 0) ? NULL : strrchr (filename, pathsep);

  if ( (!dot) || dot == filename || (lastpath>dot) || (!strcmp(dot,".o")) ) {
    /* No extension */
    return filename + strlen(filename);
  }

  return dot;
}

/**
 * \param path string containing a filename and its path
 *
 * \return a pointer toward the allocated string that contains the file basename.
 *
 * Extract basename from a path (allocate a string to store it).
 *
 */
char *MMG5_Get_basename(char *path) {
  char *s = strrchr(path, '/');

  if (!s)
    return strdup(path);
  else
    return strdup(s + 1);
}

/**
 * \param path string containing a filename and its path
 *
 * \return a pointer toward the path allocated here
 *
 * Remove filename from a path and return the path in a newly allocated string.
 *
 */
char *MMG5_Get_path(char *path) {
  char *lastpath,*retpath;
  int len;

  if ( path == NULL) return NULL;

  lastpath = (MMG5_PATHSEP == 0) ? NULL : strrchr (path, MMG5_PATHSEP);

  if ( !lastpath ) {
    return NULL;
  }


  len = 0;
  while ( path+len != lastpath ) {
    ++len;
  }

  MMG5_SAFE_MALLOC(retpath,len+1,char,return NULL);

  /* Copy the string without the extension and add \0 */
  strncpy ( retpath, path, len );
  retpath[len] = '\0';

  return retpath;
}

/**
 * \param path path from which we want to remove the extension.
 *
 * \return allocated string or NULL if the allocation fail.
 *
 * Allocate a new string and copy \a path without extension in it.
 *
 */
char *MMG5_Remove_ext (char* path,char *ext) {
  int        len;
  char       *retpath, *lastext, *lastpath;
  char       *extloc;

  /* Default extension if not provided */
  if ( (!ext) || !*ext ) {
    extloc = ".";
  }
  else {
    extloc = ext;
  }

  /* Error checks and string allocation. */
  if ( path == NULL) return NULL;

  /* Find the relevant characters and the length of the string without
   * extension */
  lastext = strstr (path, extloc);
  lastpath = (MMG5_PATHSEP == 0) ? NULL : strrchr (path, MMG5_PATHSEP);

  if ( lastext == NULL || (lastpath != NULL && lastpath > lastext) ) {
    /* No extension or the extension is left from a separator (i.e. it is not an
     * extension) */
    len = strlen(path);
  }
  else {
    /* An extension is found */
    len = 0;
    while ( path+len != lastext ) {
      ++len;
    }
  }

  MMG5_SAFE_MALLOC(retpath,len+1,char,return NULL);

  /* Copy the string without the extension and add \0 */
  strncpy ( retpath, path, len );
  retpath[len] = '\0';

  return retpath;
}

/**
 * \param ptr pointer toward the file extension (dot included)
 * \param fmt default file format.
 *
 * \return and index associated to the file format detected from the extension.
 *
 * Get the wanted file format from the mesh extension. If \a fmt is provided, it
 * is used as default file format (\a ptr==NULL), otherwise, the default file
 * format is the medit one.
 *
 */
int MMG5_Get_format( char *ptr, int fmt ) {
  /* Default is the Medit file format or a format given as input */
  int defFmt = fmt;

  if ( !ptr ) return defFmt;

  if ( !strncmp ( ptr,".meshb",strlen(".meshb") ) ) {
    return MMG5_FMT_MeditBinary;
  }
  else if ( !strncmp( ptr,".mesh",strlen(".mesh") ) ) {
    return MMG5_FMT_MeditASCII;
  }
  else if ( !strncmp( ptr,".mshb",strlen(".mshb") ) ) {
    return MMG5_FMT_GmshBinary;
  }
  else if ( !strncmp( ptr,".msh",strlen(".msh") ) ) {
    return MMG5_FMT_GmshASCII;
  }
  else if ( !strncmp ( ptr,".pvtu",strlen(".pvtu") ) ) {
    return MMG5_FMT_VtkPvtu;
  }
  else if ( !strncmp ( ptr,".vtu",strlen(".vtu") ) ) {
    return MMG5_FMT_VtkVtu;
  }
  else if ( !strncmp ( ptr,".pvtp",strlen(".pvtu") ) ) {
    return MMG5_FMT_VtkPvtp;
  }
  else if ( !strncmp ( ptr,".vtp",strlen(".vtp") ) ) {
    return MMG5_FMT_VtkVtp;
  }
  else if ( !strncmp ( ptr,".vtk",strlen(".vtk") ) ) {
    return MMG5_FMT_VtkVtk;
  }
  else if ( !strncmp ( ptr,".node",strlen(".node") ) ) {
    return MMG5_FMT_Tetgen;
  }

  return defFmt;
}

/**
 * \param fmt file format.
 *
 * \return The name of the file format in a string.
 *
 * Print the name of the file format associated to \a fmt.
 *
 */
const char* MMG5_Get_formatName(enum MMG5_Format fmt)
{
  switch (fmt)
  {
  case MMG5_FMT_MeditASCII:
    return "MMG5_FMT_MeditASCII";
    break;
  case MMG5_FMT_MeditBinary:
    return "MMG5_FMT_MeditBinary";
    break;
  case MMG5_FMT_VtkVtu:
    return "MMG5_FMT_VtkVtu";
    break;
  case MMG5_FMT_VtkVtp:
    return "MMG5_FMT_VtkVtp";
    break;
  case MMG5_FMT_VtkPvtu:
    return "MMG5_FMT_VtkPvtu";
    break;
  case MMG5_FMT_VtkPvtp:
    return "MMG5_FMT_VtkPvtp";
    break;
  case MMG5_FMT_VtkVtk:
    return "MMG5_FMT_VtkVtk";
    break;
  case MMG5_FMT_GmshASCII:
    return "MMG5_FMT_GmshASCII";
    break;
  case MMG5_FMT_GmshBinary:
    return "MMG5_FMT_GmshBinary";
    break;
  case MMG5_FMT_Tetgen:
    return "MMG5_FMT_Tetgen";
    break;
  default:
    return "MMG5_Unknown";
  }
}
