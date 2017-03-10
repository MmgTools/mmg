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
## Compilation of mmgs executable, libraries and tests
##
## =============================================================================

SET(MMGS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/mmgs)
SET(MMGS_BINARY_DIR ${CMAKE_BINARY_DIR}/src/mmgs)

FILE(MAKE_DIRECTORY ${MMGS_BINARY_DIR})

############################################################################
#####
#####         Fortran header: libmmgsf.h
#####
############################################################################

GENERATE_FORTRAN_HEADER ( mmgs
  ${MMGS_SOURCE_DIR} libmmgs.h ${MMGS_BINARY_DIR} libmmgsf.h
  )

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################

# Source files
FILE(
  GLOB
  mmgs_library_files
  ${MMGS_SOURCE_DIR}/*.c
  ${COMMON_SOURCE_DIR}/*.c
  )
LIST(REMOVE_ITEM mmgs_library_files
  ${MMGS_SOURCE_DIR}/mmgs.c
  ${REMOVE_FILE} )
FILE(
  GLOB
  mmgs_main_file
  ${MMGS_SOURCE_DIR}/mmgs.c
  )

############################################################################
#####
#####         Compile mmgs libraries
#####
############################################################################

# Compile static library
IF ( LIBMMGS_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}s_a STATIC
    "${mmgs_library_files}" ${PROJECT_NAME}s )
ENDIF()

# Compile shared library
IF ( LIBMMGS_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}s_so SHARED
    "${mmgs_library_files}" ${PROJECT_NAME}s )
ENDIF()

IF ( LIBMMGS_STATIC OR LIBMMGS_SHARED )
  # mmgs header files needed for library
  SET( mmgs_headers
    ${MMGS_SOURCE_DIR}/libmmgs.h
    ${MMGS_BINARY_DIR}/libmmgsf.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${COMMON_BINARY_DIR}/libmmgtypesf.h
    )
  SET(MMGS_INCLUDE ${CMAKE_SOURCE_DIR}/include/mmg/mmgs )

  # Install header files in /usr/local or equivalent
  INSTALL(FILES ${mmgs_headers} DESTINATION include/mmg/mmgs)

  COPY_FORTRAN_HEADER_AND_CREATE_TARGET ( ${MMGS_BINARY_DIR} ${MMGS_INCLUDE} s )

  # Copy header files in project directory at configuration step
  # (generated file don't exists yet or are outdated)
  FILE(INSTALL  ${mmgs_headers} DESTINATION ${MMGS_INCLUDE}
    PATTERN "libmmg*f.h"  EXCLUDE)

ENDIF()

############################################################################
#####
#####         Compile program to test library
#####
############################################################################

IF ( TEST_LIBMMGS )
  INCLUDE(cmake/testing/libmmgs_tests.cmake)
ENDIF()

###############################################################################
#####
#####         Compile MMGS executable
#####
###############################################################################
ADD_AND_INSTALL_EXECUTABLE ( ${PROJECT_NAME}s
  "${mmgs_library_files}" ${mmgs_main_file} )

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

IF ( BUILD_TESTING )
  ##-------------------------------------------------------------------##
  ##------- Set the continuous integration options --------------------##
  ##-------------------------------------------------------------------##
  SET(MMGS_CI_TESTS ${CI_DIR}/mmgs )
  SET(MMG_CI_TESTS ${CI_DIR}/mmg )

  ##-------------------------------------------------------------------##
  ##--------------------------- Add tests and configure it ------------##
  ##-------------------------------------------------------------------##
  # Add runtime that we want to test for mmgs
  IF( MMGS_CI )

    ADD_EXEC_TO_CI_TESTS ( ${PROJECT_NAME}s EXECUT_MMGS )

    IF ( TEST_LIBMMGS )
      SET(LIBMMGS_EXEC0_a ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example0_a)
      SET(LIBMMGS_EXEC0_b ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example0_b)
      SET(LIBMMGS_EXEC1   ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example1)
      SET(LIBMMGS_EXEC2   ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example2)

      ADD_TEST(NAME libmmgs_example0_a   COMMAND ${LIBMMGS_EXEC0_a})
      ADD_TEST(NAME libmmgs_example0_b   COMMAND ${LIBMMGS_EXEC0_b})
      ADD_TEST(NAME libmmgs_example1   COMMAND ${LIBMMGS_EXEC1})
      ADD_TEST(NAME libmmgs_example2   COMMAND ${LIBMMGS_EXEC2})

      IF ( CMAKE_Fortran_COMPILER)
        SET(LIBMMGS_EXECFORTRAN_a ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_a)
        SET(LIBMMGS_EXECFORTRAN_b ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_b)
        ADD_TEST(NAME libmmgs_fortran_a   COMMAND ${LIBMMGS_EXECFORTRAN_a})
        ADD_TEST(NAME libmmgs_fortran_b   COMMAND ${LIBMMGS_EXECFORTRAN_b})
      ENDIF()

    ENDIF()
    # Add mmgs tests
    INCLUDE( ${CMAKE_SOURCE_DIR}/cmake/testing/mmgs_tests.cmake )

    IF ( RUN_AGAIN )
      INCLUDE( ${CMAKE_SOURCE_DIR}/cmake/testing/mmgs_rerun_tests.cmake )
    ENDIF()


  ENDIF ( MMGS_CI )

ENDIF ( BUILD_TESTING )
