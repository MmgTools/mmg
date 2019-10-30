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
  SET ( TETGEN_EXEC TETGEN_EXEC-NOTFOUND )
ENDIF()

FIND_PROGRAM(TETGEN_EXEC
  NAMES tetgen
  HINTS ${TETGEN_EXEC} ${TETGEN} ${TETGEN}/build
  $ENV{TETGEN}  $ENV{TETGEN}/build
  $ENV{HOME}/bin
  DOC "The tetgen executable"
  )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(TETGEN DEFAULT_MSG
  TETGEN_EXEC)
