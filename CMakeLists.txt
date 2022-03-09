## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

# if PROJECT_NAME is defined, mmg is a subproject
if(DEFINED PROJECT_NAME)
  set(MMG_SUBPROJECT ON)
endif()

INCLUDE(CMakeDependentOption)

LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/modules
  ${CMAKE_SOURCE_DIR}/cmake/testing )

INCLUDE(macros)

PROJECT (mmg)

# Executable path
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(COMMON_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/common)
SET(COMMON_BINARY_DIR ${PROJECT_BINARY_DIR}/src/common)

FILE(MAKE_DIRECTORY ${COMMON_BINARY_DIR})

# Path toward the continuous integration directory
SET ( CI_DIR ${PROJECT_SOURCE_DIR}/ci_tests
  CACHE PATH "Path toward the continuous integration directory" )
MARK_AS_ADVANCED ( CI_DIR )


include(CMakePackageConfigHelpers)
# various variables to help standard install directories
include(GNUInstallDirs)


###############################################################################
#####
#####         Release version and date
#####
###############################################################################

SET (CMAKE_RELEASE_VERSION_MAJOR "5" )
SET (CMAKE_RELEASE_VERSION_MINOR "6" )
SET (CMAKE_RELEASE_VERSION_PATCH "0" )
SET (CMAKE_RELEASE_DATE "Nov. 05, 2021" )

SET (CMAKE_RELEASE_VERSION
  "${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}")

###############################################################################
#####
#####  Set preprocessor flags and manage automatically generated headers data
#####
###############################################################################
# Architecture detection
IF (NOT WIN32)
  # Set preprocessor flags to say that we are posix and gnu compatible
  SET (DEF_POSIX "#define MMG_POSIX")
  SET (DEF_GNU "#define MMG_GNU")
ELSEIF(MINGW)
  SET (DEF_GNU "#define MMG_GNU")
ENDIF ()

# Create mmgcmakedefines.h file with the architecture infos and
# the definition ot the USE_POINTMAP flag.
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/common/mmgcmakedefines.h.in
  ${PROJECT_BINARY_DIR}/src/common/mmgcmakedefines.h @ONLY)

# Create mmgversion.h file with the release infos
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/common/mmgversion.h.in
  ${PROJECT_BINARY_DIR}/src/common/mmgversion.h @ONLY)

# Print git hash to source file
IF (NOT WIN32 OR MINGW)

  ADD_CUSTOM_TARGET(GenerateGitHash
    COMMAND ./git_log_mmg.sh ${PROJECT_SOURCE_DIR} ${COMMON_BINARY_DIR}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/
    COMMENT "Getting git commit hash"
  )
ENDIF ()

###############################################################################
#####
#####         CMake flags / variables
#####
###############################################################################
# Hide some options
MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES
  CMAKE_OSX_DEPLOYMENT_TARGET
  CMAKE_OSX_SYSROOT)

include(add_build_types)

# List CMake build types
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  MESSAGE(STATUS "Setting build type to 'Release' as none was specified.")
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo" "Maintainer" "RelWithAssert")
ELSEIF( CMAKE_CONFIGURATION_TYPES )
  LIST(APPEND CMAKE_CONFIGURATION_TYPES Maintainer RelWithAssert )
ENDIF()

# Explicitly set the DNDEBUG flag in case the user or a parent project overrides
# it.
if (NOT CMAKE_BUILD_TYPE MATCHES Debug AND NOT CMAKE_BUILD_TYPE MATCHES RelWithAssert)
    add_definitions(-DNDEBUG)
endif()

# Do not print warnings about the use of secure function with MSVC and the use
# of POSIX strdup
if ( MSVC )
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-Dstrdup=_strdup)
endif()

###############################################################################
#####
#####         Choose executable target to compile
#####
###############################################################################
SET ( BUILD "MMG" CACHE STRING
  "\"MMG3D\", \"MMGS\",\"MMG2D\", or \"MMG\" to compile both MMGS, MMG3D and MMG2D" )

IF ( BUILD MATCHES "MMGS" )
  SET(BUILD_MMG        OFF )
  SET(BUILD_MMG3D      OFF )
  SET(BUILD_MMGS       ON  )
  SET(BUILD_MMG2D      OFF  )
ELSEIF ( BUILD MATCHES "MMG3D" )
  SET(BUILD_MMG3D      ON  )
  SET(BUILD_MMGS       OFF )
  SET(BUILD_MMG        OFF )
  SET(BUILD_MMG2D      OFF  )
ELSEIF ( BUILD MATCHES "MMG2D" )
  SET(BUILD_MMG        OFF )
  SET(BUILD_MMG3D      OFF )
  SET(BUILD_MMGS       OFF  )
  SET(BUILD_MMG2D      ON  )
ELSEIF ( BUILD MATCHES "MMG" )
  SET(BUILD_MMG        ON  )
  SET(BUILD_MMGS       ON )
  SET(BUILD_MMG3D      ON )
  SET(BUILD_MMG2D      ON  )
ELSE ( )
  MESSAGE ( ERROR " Nothing to build.\n"
    "Please choose the target to build (set the BUILD CMake's"
    " variable to \"MMG\", \"MMG2D\", \"MMGS\" OR"
    " \"MMG3D\")" )

ENDIF()

IF ( BUILD_MMG3D OR BUILD_MMGS )
 SET(BUILD_MMGS3D ON)
ELSE ()
  SET(BUILD_MMGS3D OFF)
ENDIF()

############################################################################
#####
#####         Fortran header: libmmgtypesf.h
#####
############################################################################
ADD_EXECUTABLE(genheader ${PROJECT_SOURCE_DIR}/scripts/genheader.c)
TARGET_LINK_LIBRARIES ( genheader PRIVATE ${LIBRARIES}  )
GENERATE_FORTRAN_HEADER ( mmg
  ${COMMON_SOURCE_DIR} libmmgtypes.h
  ${COMMON_BINARY_DIR}
  ${COMMON_BINARY_DIR} libmmgtypesf.h
  )

###############################################################################
#####
#####         libraries
#####
###############################################################################

# Libraries
IF ( NOT (WIN32 OR MINGW) )

  IF(NOT DEFINED M_LIB)
    MESSAGE(STATUS "M_LIB not defined. Searching it")
    FIND_LIBRARY(M_LIB m)
  ENDIF(NOT DEFINED M_LIB)

  MARK_AS_ADVANCED(FORCE M_LIB)
  SET( LIBRARIES ${M_LIB} ${LIBRARIES})
ENDIF()

############################################################################
#####
#####         Point map
#####
############################################################################

# add point map?
OPTION ( USE_POINTMAP "Use map for point tracking" OFF )

############################################################################
#####
#####         Scotch
#####
############################################################################
# Find SCOTCH library?
SET(SCOTCH_DIR "" CACHE PATH "Installation directory for scotch")

# add Scotch library?
OPTION ( USE_SCOTCH "Use SCOTCH TOOL for renumbering" ON )

IF ( USE_SCOTCH )
  FIND_PACKAGE(SCOTCH)

  IF ( NOT SCOTCH_FOUND )
    MESSAGE ( WARNING "Scotch library not found:"
      "Using scotch reduce the execution time of mmg3d "
      "(see https://gforge.inria.fr/frs/?group_id=248 to download it)."
      "If you have already installed Scotch and want to use it, "
      "please set the CMake variable or environment variable SCOTCH_DIR "
      "to your scotch directory.")
  ENDIF ( )
ENDIF ( )

# add LinearElasticity library?
SET(ELAS_DIR "" CACHE PATH "Installation directory for Elas")
OPTION (
  USE_ELAS "Use the Elas library for lagrangian motion option" ON )

IF ( USE_ELAS )
  INCLUDE(FindElas)
  IF ( NOT ELAS_FOUND )
    MESSAGE ( WARNING "Elas is a library to solve the linear elasticity "
      "problem (see https://github.com/ISCDtoolbox/LinearElasticity to"
      " download it). "
      "This library is needed to use the lagrangian motion option. "
      "If you have already installed Elas and want to use it, "
      "please set the CMake variable or environment variable ELAS_DIR "
      "to your Elas directory.")
  ENDIF ( )
ENDIF ( )

############################################################################
#####
#####         VTK (to parse (p)vtp/(p)vtu files )
#####
############################################################################
# add the VTK library ?
OPTION ( USE_VTK "Use VTK I/O" ON )
IF ( USE_VTK )
  FIND_PACKAGE(VTK QUIET)
  IF ( NOT VTK_FOUND )
    MESSAGE ( WARNING "VTK library not found: vtk I/O will not be available.")
  ENDIF ( )
ENDIF ( )

############################################################################
#####
##### RPATH for MacOSX
#####
############################################################################

SET(CMAKE_MACOSX_RPATH 1)

###############################################################################
#####
#####         Add dependent options
#####
###############################################################################

FUNCTION (INVERT_BOOL OUTVAR INVAR)
  IF(${INVAR})
    SET(${OUTVAR} OFF PARENT_SCOPE)
  ELSE()
    SET(${OUTVAR} ON PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

OPTION(BUILD_SHARED_LIBS "Build shared libraries" OFF)
INVERT_BOOL("BUILD_STATIC_LIBS" ${BUILD_SHARED_LIBS})
IF ( (${BUILD_STATIC_LIBS} EQUAL ON) AND NOT CMAKE_POSITION_INDEPENDENT_CODE)
  SET(CMAKE_POSITION_INDEPENDENT_CODE ON)
ENDIF()

#------------------------------- mmg2d
CMAKE_DEPENDENT_OPTION ( LIBMMG2D_STATIC "Compile the mmg2d static library" ${BUILD_STATIC_LIBS}
  "BUILD_MMG2D" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMG2D_SHARED "Compile the mmg2d dynamic library" ${BUILD_SHARED_LIBS}
  "BUILD_MMG2D" OFF )
IF (LIBMMG2D_STATIC OR LIBMMG2D_SHARED)
  SET(LIBMMG2D_INTERNAL ON )
ELSE()
  SET(LIBMMG2D_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMG2D "Compile tests for the mmg2d library" OFF
  "BUILD_MMG2D;LIBMMG2D_INTERNAL" OFF
  )

CMAKE_DEPENDENT_OPTION(
  MMG2D_CI "Enable/Disable continuous integration for mmg2d" ON
  "BUILD_TESTING;BUILD_MMG2D" OFF)

#------------------------------- mmgs
CMAKE_DEPENDENT_OPTION ( LIBMMGS_STATIC "Compile the mmgs static library" ${BUILD_STATIC_LIBS}
  "BUILD_MMGS" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMGS_SHARED "Compile the mmgs dynamic library" ${BUILD_SHARED_LIBS}
  "BUILD_MMGS" OFF)
IF (LIBMMGS_STATIC OR LIBMMGS_SHARED)
  SET(LIBMMGS_INTERNAL ON )
ELSE()
  SET(LIBMMGS_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMGS "Compile tests for the mmgs library" OFF
  "LIBMMGS_INTERNAL;BUILD_MMGS" OFF
  )

CMAKE_DEPENDENT_OPTION(
  MMGS_CI "Enable/Disable continuous integration for mmgs" ON
  "BUILD_TESTING;BUILD_MMGS" OFF )

#------------------------------- mmg3d
CMAKE_DEPENDENT_OPTION ( LIBMMG3D_STATIC "Compile the mmg3d static library" ${BUILD_STATIC_LIBS}
  "BUILD_MMG3D" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMG3D_SHARED "Compile the mmg3d dynamic library" ${BUILD_SHARED_LIBS}
  "BUILD_MMG3D" OFF )
IF (LIBMMG3D_STATIC OR LIBMMG3D_SHARED)
  SET(LIBMMG3D_INTERNAL ON )
ELSE()
  SET(LIBMMG3D_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMG3D "Compile tests for the mmg3d library" OFF
  "BUILD_MMG3D;LIBMMG3D_INTERNAL" OFF
  )
CMAKE_DEPENDENT_OPTION(
  MMG3D_CI "Enable/Disable continuous integration for mmg3d" ON
  "BUILD_TESTING;BUILD_MMG3D" OFF)

CMAKE_DEPENDENT_OPTION ( PATTERN
  "If ON, insertion by patterns, otherwise insertion by delaunay kernel"
  OFF
  "BUILD_MMG3D" OFF)
#-- Remove the next line to have this option visible in basic cmake mode
MARK_AS_ADVANCED(PATTERN)


IF( SCOTCH_FOUND )
  add_definitions(-DUSE_SCOTCH)

  MESSAGE(STATUS
    "Compilation with scotch: ${SCOTCH_LIBRARIES}")
  SET( LIBRARIES ${SCOTCH_LIBRARIES} ${LIBRARIES})
ENDIF()


IF ( VTK_FOUND )
  add_definitions(-DUSE_VTK)

  MESSAGE ( STATUS "Compilation with VTK: add vtk, vtp and vtu I/O." )

  IF( "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}" LESS 8.90 )
    INCLUDE ( ${VTK_USE_FILE} )
  ENDIF()

  SET( LIBRARIES ${VTK_LIBRARIES} ${LIBRARIES} )
ENDIF ( )

#---------------------------- mmg library gathering mmg2d, mmgs and mmg3d libs
CMAKE_DEPENDENT_OPTION ( LIBMMG_STATIC
  "Compile the mmg static library (mmg2d + mmgs + mmg3d)" ${BUILD_STATIC_LIBS}
  "BUILD_MMG" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMG_SHARED
  "Compile the mmg dynamic library (mmg2d + mmgs + mmg3d)" ${BUILD_SHARED_LIBS}
  "BUILD_MMG" OFF )
IF (LIBMMG_STATIC OR LIBMMG_SHARED)
  SET(LIBMMG_INTERNAL ON )
ELSE()
  SET(LIBMMG_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMG "Compile tests for the mmg library" OFF
  "BUILD_MMG;LIBMMG_INTERNAL" OFF
  )
CMAKE_DEPENDENT_OPTION(
  MMG_CI "Enable/Disable continuous integration for mmg" ON
  "BUILD_TESTING;BUILD_MMG" OFF)


###############################################################################
#####
#####         Build Targets and libmmg if needed
#####
###############################################################################
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# CTEST_OUTPUT_DIR MUST BE SETTED HERE TO AVOID UNINITIALIZATION FOR LIBRARY
# TESTS DEFS
OPTION ( BUILD_TESTING "Enable/Disable continuous integration" OFF )

IF( BUILD_TESTING )
  SET ( CTEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/TEST_OUTPUTS
    CACHE PATH "Path toward the tests outputs" )
  MARK_AS_ADVANCED ( CTEST_OUTPUT_DIR )
  FILE ( MAKE_DIRECTORY  ${CTEST_OUTPUT_DIR} )
ENDIF ( )

IF ( BUILD_MMG2D )
  INCLUDE(mmg2d)
ENDIF ( )

IF ( BUILD_MMGS )
  INCLUDE(mmgs)
ENDIF ( )

IF ( BUILD_MMG3D )
  INCLUDE(mmg3d)
ENDIF ( )

IF ( BUILD_MMG )
  INCLUDE(mmg)
ENDIF ( )

IF( SCOTCH_FOUND )
  # Include Scotch Dir here to ensure that Mmg doesn't search it's own headers in /usr/local
  INCLUDE_DIRECTORIES(${SCOTCH_INCLUDE_DIRS})
ENDIF()

IF ( MmgTargetsExported )
  set( MMG_CMAKE_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/mmg )

  configure_package_config_file(cmake/config/mmgConfig.cmake.in
    ${PROJECT_BINARY_DIR}/mmgConfig.cmake
    INSTALL_DESTINATION ${MMG_CMAKE_INSTALL_DIR} )

  install(EXPORT MmgTargets
    FILE MmgTargets.cmake
    NAMESPACE Mmg::
    DESTINATION ${MMG_CMAKE_INSTALL_DIR}
    )
  install(FILES ${PROJECT_BINARY_DIR}/mmgConfig.cmake
    DESTINATION ${MMG_CMAKE_INSTALL_DIR} )

  install(FILES
    ${PROJECT_SOURCE_DIR}/cmake/modules/FindSCOTCH.cmake
    ${PROJECT_SOURCE_DIR}/cmake/modules/FindElas.cmake
    DESTINATION ${MMG_CMAKE_INSTALL_DIR} )
ENDIF()

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################
CMAKE_DEPENDENT_OPTION (
  ONLY_VERY_SHORT_TESTS "Enable/Disable very short tests" OFF
  "BUILD_TESTING;NOT LONG_TESTS" OFF )
MARK_AS_ADVANCED(ONLY_VERY_SHORT_TESTS)

CMAKE_DEPENDENT_OPTION (
  LONG_TESTS "Enable/Disable time consuming tests for the project" OFF
  "BUILD_TESTING" OFF )
MARK_AS_ADVANCED(LONG_TESTS)

IF( BUILD_TESTING )

  IF ( NOT ONLY_VERY_SHORT_TESTS )
    ##-------------------------------------------------------------------##
    ##------- Get the continuous integration tests ----------------------##
    ##-------------------------------------------------------------------##
    INCLUDE(LoadCiTests)
  ENDIF ( )

  ##-------------------------------------------------------------------##
  ##------- Set the continuous integration options --------------------##
  ##-------------------------------------------------------------------##
  SET ( MMG_CI_TESTS ${CI_DIR}/mmg )

  # Copy config file and custom file in the good directories
  CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/cmake/config/CTestCustom.cmake
    ./CTestCustom.cmake)

  ENABLE_TESTING()
  INCLUDE(CTest)

  IF ( ( MMG3D_CI OR MMGS_CI ) AND NOT ONLY_VERY_SHORT_TESTS )
    # Add common tests for mmgs/3d appli
    INCLUDE( mmg_tests )

  ENDIF ( )

ENDIF ( BUILD_TESTING )

###############################################################################
#####
#####         Create API Documentation
#####
###############################################################################
include(doxygen)

###############################################################################
#####
#####         Create Package
#####
###############################################################################
include(packages)
