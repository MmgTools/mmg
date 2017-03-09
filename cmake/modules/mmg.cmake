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
## Compilation of mmg library and tests
##
## =============================================================================

# Note: mmg3d.cmake, mmg2d.cmake and mmgs.cmake must have been called before
# (thus the MMG*_SOURCE_DIR variables are already setted)

############################################################################
#####
#####         Compile mmg libraries
#####
############################################################################

# Compile static library
IF ( LIBMMG_STATIC )
  ADD_LIBRARY(lib${PROJECT_NAME}_a  STATIC
    ${mmg2d_library_files} ${mmgs_library_files} ${mmg3d_library_files}
    )
  SET_TARGET_PROPERTIES(lib${PROJECT_NAME}_a PROPERTIES OUTPUT_NAME
    ${PROJECT_NAME})
  TARGET_LINK_LIBRARIES(lib${PROJECT_NAME}_a ${LIBRARIES})
  INSTALL(TARGETS lib${PROJECT_NAME}_a
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
ENDIF()

# Compile shared library
IF ( LIBMMG_SHARED )
  ADD_LIBRARY(lib${PROJECT_NAME}_so SHARED
    ${mmg2d_library_files} ${mmgs_library_files} ${mmg3d_library_files}
    )
  SET_TARGET_PROPERTIES(lib${PROJECT_NAME}_so PROPERTIES
    VERSION ${CMAKE_RELEASE_VERSION} SOVERSION 5)
  SET_TARGET_PROPERTIES(lib${PROJECT_NAME}_so PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME})
  TARGET_LINK_LIBRARIES(lib${PROJECT_NAME}_so ${LIBRARIES})
  INSTALL(TARGETS lib${PROJECT_NAME}_so
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
ENDIF()

IF ( LIBMMG_STATIC OR LIBMMG_SHARED )
  # mmg header files needed for library
  SET( mmg2d_headers
    ${MMG2D_SOURCE_DIR}/libmmg2d.h
    ${MMG2D_BINARY_DIR}/libmmg2df.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${COMMON_BINARY_DIR}/libmmgtypesf.h
    )
  SET( mmg3d_headers
    ${MMG3D_SOURCE_DIR}/libmmg3d.h
    ${MMG3D_BINARY_DIR}/libmmg3df.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${COMMON_BINARY_DIR}/libmmgtypesf.h
    )
  SET( mmgs_headers
    ${MMGS_SOURCE_DIR}/libmmgs.h
    ${MMGS_BINARY_DIR}/libmmgsf.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${COMMON_BINARY_DIR}/libmmgtypesf.h
    )
  SET( mmg_headers
    ${CMAKE_SOURCE_DIR}/src/mmg/libmmg.h
    ${CMAKE_SOURCE_DIR}/src/mmg/libmmgf.h
    )
  SET(MMG2D_INCLUDE ${CMAKE_SOURCE_DIR}/include/mmg/mmg2d )
  SET(MMGS_INCLUDE ${CMAKE_SOURCE_DIR}/include/mmg/mmgs )
  SET(MMG3D_INCLUDE ${CMAKE_SOURCE_DIR}/include/mmg/mmg3d )
  SET(MMG_INCLUDE ${CMAKE_SOURCE_DIR}/include/mmg )

  # Install header files in /usr/local or equivalent
  INSTALL(FILES ${mmg2d_headers} DESTINATION include/mmg/mmg2d)
  INSTALL(FILES ${mmgs_headers} DESTINATION include/mmg/mmgs)
  INSTALL(FILES ${mmg3d_headers} DESTINATION include/mmg/mmg3d)
  INSTALL(FILES ${mmg_headers} DESTINATION include/mmg)

  IF ( (NOT LIBMMG2D_STATIC) AND (NOT LIBMMG2D_SHARED) )
    COPY_FORTRAN_HEADER (
      ${COMMON_BINARY_DIR} libmmgtypesf.h ${MMG2D_INCLUDE} libmmgtypesf.h
      mmg_fortran_header copy2d_libmmgtypesf )

    COPY_FORTRAN_HEADER (
      ${MMG2D_BINARY_DIR} libmmg2df.h ${MMG2D_INCLUDE} libmmg2df.h
      mmg2d_fortran_header copy_libmmg2df
      )
  ENDIF ()

  IF ( (NOT LIBMMGS_STATIC) AND (NOT LIBMMGS_SHARED) )
    COPY_FORTRAN_HEADER (
      ${COMMON_BINARY_DIR} libmmgtypesf.h ${MMGS_INCLUDE} libmmgtypesf.h
      mmg_fortran_header copys_libmmgtypesf )

    COPY_FORTRAN_HEADER (
      ${MMGS_BINARY_DIR} libmmgsf.h ${MMGS_INCLUDE} libmmgsf.h
      mmgs_fortran_header copy_libmmgsf
    )
  ENDIF()

  IF ( (NOT LIBMMG3D_STATIC) AND (NOT LIBMMG3D_SHARED) )
    COPY_FORTRAN_HEADER (
      ${COMMON_BINARY_DIR} libmmgtypesf.h ${MMG3D_INCLUDE} libmmgtypesf.h
      mmg_fortran_header copy3d_libmmgtypesf )

    COPY_FORTRAN_HEADER (
      ${MMG3D_BINARY_DIR} libmmg3df.h ${MMG3D_INCLUDE} libmmg3df.h
      mmg3d_fortran_header copy_libmmg3df
      )
  ENDIF()

  FILE(INSTALL ${CMAKE_SOURCE_DIR}/src/mmg/libmmgf.h DESTINATION  ${CMAKE_SOURCE_DIR}/include/mmg/)


  # Install header files in project directory
  FILE(INSTALL  ${mmg2d_headers}
    DESTINATION ${CMAKE_SOURCE_DIR}/include/mmg/mmg2d/
    PATTERN "libmmg*f.h"  EXCLUDE)
  FILE(INSTALL  ${mmgs_headers}
    DESTINATION ${CMAKE_SOURCE_DIR}/include/mmg/mmgs/
    PATTERN "libmmg*f.h"  EXCLUDE)
  FILE(INSTALL  ${mmg3d_headers}
    DESTINATION ${CMAKE_SOURCE_DIR}/include/mmg/mmg3d/
    PATTERN "libmmg*f.h"  EXCLUDE)
  FILE(INSTALL  ${mmg_headers}
    DESTINATION ${CMAKE_SOURCE_DIR}/include/mmg/
    PATTERN "libmmg*f.h"  EXCLUDE)


  ADD_CUSTOM_TARGET(copy_mmg_headers ALL
    DEPENDS
    copy_2d_headers copy_s_headers copy_3d_headers
    ${CMAKE_SOURCE_DIR}/include/mmg/libmmgf.h
    ${CMAKE_SOURCE_DIR}/include/mmg/libmmg.h
    ${CMAKE_SOURCE_DIR}/include/mmg/mmg3d/libmmgtypes.h )

ENDIF()

############################################################################
#####
#####         Compile program to test library
#####
############################################################################

IF ( TEST_LIBMMG )
  INCLUDE(cmake/testing/libmmg_tests.cmake)
ENDIF()

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

IF ( BUILD_TESTING )
  ##-------------------------------------------------------------------##
  ##--------------------------- Add tests and configure it ------------##
  ##-------------------------------------------------------------------##
  # Add runtime that we want to test for mmg
  IF( MMG_CI )
    # Add libmmg tests
    IF ( TEST_LIBMMG )
      SET(LIBMMG_EXEC0_a ${EXECUTABLE_OUTPUT_PATH}/libmmg_example0_a)
      SET(LIBMMG_CPP_a   ${EXECUTABLE_OUTPUT_PATH}/libmmg_cpp_a)

     ADD_TEST(NAME libmmg_example0_a   COMMAND ${LIBMMG_EXEC0_a})
     ADD_TEST(NAME libmmg_cpp_a        COMMAND ${LIBMMG_CPP_a})

      IF ( CMAKE_Fortran_COMPILER)
        SET(LIBMMG_FORTRAN_a ${EXECUTABLE_OUTPUT_PATH}/libmmg_fortran_a)
        ADD_TEST(NAME libmmg_fortran   COMMAND ${LIBMMG_FORTRAN_a})
      ENDIF()
    ENDIF ()

  ENDIF( MMG_CI )

ENDIF ( BUILD_TESTING )
