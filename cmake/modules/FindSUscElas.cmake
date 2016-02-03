## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
##
##  mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the mmg distribution only if you accept them.
## =============================================================================

IF ((NOT WIN32) AND (NOT WIN64))
  SET ( SUSCELAS_INCLUDE_DIR SUSCELAS_INCLUDE_DIR-NOTFOUND )
  SET ( SUSCELAS_LIBRARY SUSCELAS_LIBRARY-NOTFOUND )
ENDIF()

FIND_PATH(SUSCELAS_INCLUDE_DIR
  NAMES elastic.h
  HINTS ${SUSCELAS_INCLUDE_DIR}
  $ENV{SUSCELAS_INCLUDE_DIR}
  $ENV{HOME}/include/SUscElas
  ${SUSCELAS_DIR}/include/SUscElas
  $ENV{SUSCELAS_DIR}/include/SUscElas
  PATH_SUFFIXES SUscElas
  DOC "Directory of SUscElas Header")

# Check for elas library (search in common locations in order to prevent
# modifications in SUscElas installation)
FIND_LIBRARY(SUSCELAS_LIBRARY
  NAMES Elas Elas${SUSCELAS_LIB_SUFFIX}
  HINTS ${SUSCELAS_LIBRARY}
  $ENV{SUSCELAS_LIBRARY}
  $ENV{HOME}/lib
  ${SUSCELAS_DIR}/lib
  $ENV{SUSCELAS_DIR}/lib
  DOC "The SUscElas library"
  )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SUSCELAS DEFAULT_MSG
  SUSCELAS_INCLUDE_DIR SUSCELAS_LIBRARY)
IF ((NOT WIN32) AND (NOT WIN64))
  MARK_AS_ADVANCED(SUSCELAS_INCLUDE_DIR SUSCELAS_LIBRARY)
ENDIF()
