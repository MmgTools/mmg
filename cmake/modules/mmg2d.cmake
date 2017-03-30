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

## =============================================================================
##
## Compilation of mmg2d executable, libraries and tests
##
## =============================================================================

SET(MMG2D_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/mmg2d)
SET(MMG2D_BINARY_DIR ${CMAKE_BINARY_DIR}/src/mmg2d)

FILE(MAKE_DIRECTORY ${MMG2D_BINARY_DIR})

############################################################################
#####
#####         Fortran header: libmmg2df.h
#####
############################################################################

GENERATE_FORTRAN_HEADER ( mmg2d
  ${MMG2D_SOURCE_DIR} libmmg2d.h ${MMG2D_BINARY_DIR} libmmg2df.h
  )

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################

# Source files
FILE(
  GLOB
  mmg2d_library_files
  ${MMG2D_SOURCE_DIR}/*.c
  ${COMMON_SOURCE_DIR}/*.c
  )
LIST(REMOVE_ITEM mmg2d_library_files
  ${MMG2D_SOURCE_DIR}/mmg2d.c
  ${REMOVE_FILE} )
FILE(
  GLOB
  mmg2d_main_file
  ${MMG2D_SOURCE_DIR}/mmg2d.c
  )

############################################################################
#####
#####         Elastic
#####
############################################################################

IF( USE_ELAS )
# Set flags for building test program
INCLUDE_DIRECTORIES(${ELAS_INCLUDE_DIR})

SET(CMAKE_REQUIRED_INCLUDES ${ELAS_INCLUDE_DIR})
SET(CMAKE_REQUIRED_LIBRARIES ${ELAS_LIBRARY})

SET(CMAKE_C_FLAGS "-DUSE_ELAS ${CMAKE_C_FLAGS}")
MESSAGE(STATUS
"Compilation with the Elas library: ${ELAS_LIBRARY} ")
SET( LIBRARIES ${ELAS_LINK_FLAGS} ${LIBRARIES})
SET( LIBRARIES ${ELAS_LIBRARY} ${LIBRARIES})
ENDIF()

IF (ELAS_NOTFOUND)
MESSAGE ( WARNING "Elas is a library to solve the linear elasticity "
    "problem (see https://github.com/SUscTools/Elas to download it). "
"This library is needed to use the lagrangian motion option. "
    "If you have already installed Elas and want to use it, "
"please set the CMake variable or environment variable ELAS_DIR "
"to your Elas directory.")
ENDIF ( )

############################################################################
#####
#####         Compile mmg2d libraries
#####
############################################################################
# Compile static library
IF ( LIBMMG2D_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}2d_a STATIC
    "${mmg2d_library_files}" ${PROJECT_NAME}2d )
ENDIF()

# Compile shared library
IF ( LIBMMG2D_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}2d_so SHARED
    "${mmg2d_library_files}" ${PROJECT_NAME}2d )
ENDIF()

IF ( LIBMMG2D_STATIC OR LIBMMG2D_SHARED )
  # mmg2d header files needed for library
  SET( mmg2d_headers
    ${MMG2D_SOURCE_DIR}/libmmg2d.h
    ${MMG2D_BINARY_DIR}/libmmg2df.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${COMMON_BINARY_DIR}/libmmgtypesf.h
    )
  SET(MMG2D_INCLUDE ${CMAKE_SOURCE_DIR}/include/mmg/mmg2d )

  # Install header files in /usr/local or equivalent
  INSTALL(FILES ${mmg2d_headers} DESTINATION include/mmg/mmg2d)

  COPY_FORTRAN_HEADER_AND_CREATE_TARGET ( ${MMG2D_BINARY_DIR} ${MMG2D_INCLUDE} 2d )

  # Copy header files in project directory at configuration step
  # (generated file don't exists yet or are outdated)
  FILE(INSTALL  ${mmg2d_headers} DESTINATION ${MMG2D_INCLUDE}
    PATTERN "libmmg*f.h"  EXCLUDE)

ENDIF()

############################################################################
#####
#####         Compile program to test library
#####
############################################################################

IF ( TEST_LIBMMG2D )
  INCLUDE(cmake/testing/libmmg2d_tests.cmake)
ENDIF ( )

###############################################################################
#####
#####         Compile MMG2D executable
#####
###############################################################################

ADD_AND_INSTALL_EXECUTABLE ( ${PROJECT_NAME}2d
  "${mmg2d_library_files}" ${mmg2d_main_file} )

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

IF ( BUILD_TESTING )
  ##-------------------------------------------------------------------##
  ##------- Set the continuous integration options --------------------##
  ##-------------------------------------------------------------------##
  SET(MMG2D_CI_TESTS ${CI_DIR}/mmg2d )

  ##-------------------------------------------------------------------##
  ##--------------------------- Add tests and configure it ------------##
  ##-------------------------------------------------------------------##
  # Add runtime that we want to test for mmg2d
  IF ( MMG2D_CI )

    SET ( CTEST_OUTPUT_DIR ${CMAKE_BINARY_DIR}/TEST_OUTPUTS )
    FILE ( MAKE_DIRECTORY  ${CTEST_OUTPUT_DIR} )

    ADD_EXEC_TO_CI_TESTS ( ${PROJECT_NAME}2d EXECUT_MMG2D )

    IF ( TEST_LIBMMG2D )
      SET(LIBMMG2D_EXEC0_a ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_example0_a )

      SET(LIBMMG2D_EXEC0_b ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_example0_b )
      SET(LIBMMG2D_EXEC1 ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_example1 )

      ADD_TEST(NAME libmmg2d_example0_a   COMMAND ${LIBMMG2D_EXEC0_a}
        "${CMAKE_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0/example0_a/init.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_0_a-init.o"
        )
      ADD_TEST(NAME libmmg2d_example0_b   COMMAND ${LIBMMG2D_EXEC0_b}
        "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_0_b.o.mesh"
        )
      ADD_TEST(NAME libmmg2d_example1   COMMAND ${LIBMMG2D_EXEC1}
        "${CMAKE_SOURCE_DIR}/libexamples/mmg2d/adaptation_example1/dom.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_1-dom.o"
       )

      IF ( CMAKE_Fortran_COMPILER)
        SET(LIBMMG2D_EXECFORTRAN_a ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_a )

        SET(LIBMMG2D_EXECFORTRAN_b ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_b )
        ADD_TEST(NAME libmmg2d_fortran_a   COMMAND ${LIBMMG2D_EXECFORTRAN_a}
          "${CMAKE_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0_fortran/example0_a/init.mesh"
          "${CTEST_OUTPUT_DIR}/libmmg2d-Adaptation_Fortran_0_a-init.o"
         )
        ADD_TEST(NAME libmmg2d_fortran_b   COMMAND ${LIBMMG2D_EXECFORTRAN_b}
          "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_Fortran_0_b.o"
         )
      ENDIF()

    ENDIF()
    # Add mmg2d tests
    INCLUDE( ${CMAKE_SOURCE_DIR}/cmake/testing/mmg2d_tests.cmake )

  ENDIF( MMG2D_CI )

ENDIF ( BUILD_TESTING )
