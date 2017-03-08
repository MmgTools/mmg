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
    ${sourcemmg2d_files} ${libmmg2d_file}
    ${sourcemmgs_files} ${libmmgs_file}
    ${source_files} ${lib_file}
    ${CMAKE_SOURCE_DIR}/src/mmg/libmmg.h
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
    ${sourcemmg2d_files} ${libmmg2d_file}
    ${sourcemmgs_files} ${libmmgs_file}
    ${source_files} ${lib_file}
    ${CMAKE_SOURCE_DIR}/src/mmg/libmmg.h
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
  SET( mmg2d_includes
    ${MMG2D_INCLUDE}/libmmg2d.h
    ${MMG2D_INCLUDE}/libmmg2df.h
    )
  SET( mmgs_includes
    ${MMGS_INCLUDE}/libmmgs.h
    ${MMGS_INCLUDE}/libmmgsf.h
    )
  SET( mmg3d_includes
    ${MMG3D_INCLUDE}/libmmg3d.h
    ${MMG3D_INCLUDE}/libmmg3df.h
    )
  SET( mmg_includes
    ${MMG_INCLUDE}/libmmg.h
    ${MMG_INCLUDE}/libmmgf.h
    ${MMG_INCLUDE}/mmg3d/libmmgtypes.h
    ${MMG_INCLUDE}/mmg3d/libmmgtypesf.h
    )

  # Install header files in /usr/local or equivalent
  INSTALL(FILES ${mmg2d_headers} DESTINATION include/mmg/mmg2d)
  INSTALL(FILES ${mmgs_headers} DESTINATION include/mmg/mmgs)
  INSTALL(FILES ${mmg3d_headers} DESTINATION include/mmg/mmg3d)
  INSTALL(FILES ${mmg_headers} DESTINATION include/mmg)

  ADD_CUSTOM_COMMAND(OUTPUT ${MMG_INCLUDE}/libmmgtypesf.h
    COMMAND ${CMAKE_COMMAND} -E copy ${COMMON_BINARY_DIR}/libmmgtypesf.h ${MMG_INCLUDE}/libmmgtypesf.h
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS ${COMMON_BINARY_DIR}/libmmgtypesf.h)
  ADD_CUSTOM_TARGET ( copy_libmmgtypesf ALL
    DEPENDS ${MMG_INCLUDE}/libmmgtypesf.h )
  ADD_DEPENDENCIES ( copy_libmmgtypesf mmg_fortran_header )

  IF ( (NOT LIBMMG2D_STATIC) AND (NOT LIBMMG2D_SHARED) )
    ADD_CUSTOM_COMMAND(OUTPUT ${MMG2D_INCLUDE}/libmmg2df.h
      COMMAND ${CMAKE_COMMAND} -E copy  ${MMG2D_BINARY_DIR}/libmmg2df.h ${MMG2D_INCLUDE}/libmmg2df.h
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )
    ADD_CUSTOM_TARGET ( copy_libmmg2df ALL
      DEPENDS ${MMG2D_INCLUDE}/libmmg2df.h )
    ADD_DEPENDENCIES ( copy_libmmg2df mmg2d_fortran_header )

    ADD_CUSTOM_TARGET(copy_2d_headers ALL
      DEPENDS
      copy_libmmg2df copy_libmmgtypesf
      ${CMAKE_SOURCE_DIR}/include/mmg/mmg2d/libmmg2d.h
      ${CMAKE_SOURCE_DIR}/include/mmg/mmg2d/libmmgtypes.h )
  ENDIF ()

  IF ( (NOT LIBMMGS_STATIC) AND (NOT LIBMMGS_SHARED) )
    ADD_CUSTOM_COMMAND(OUTPUT ${MMGS_INCLUDE}/libmmgsf.h
      COMMAND ${CMAKE_COMMAND} -E copy  ${MMGS_BINARY_DIR}/libmmgsf.h ${MMGS_INCLUDE}/libmmgsf.h
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      DEPENDS ${MMGS_BINARY_DIR}/libmmgsf.h)
    ADD_CUSTOM_TARGET ( copy_libmmgsf ALL
      DEPENDS ${MMGS_INCLUDE}/libmmgsf.h )
    ADD_DEPENDENCIES ( copy_libmmgsf mmgs_fortran_header )

    ADD_CUSTOM_TARGET(copy_s_headers ALL
      DEPENDS
      copy_libmmgsf copy_libmmgtypesf
      ${CMAKE_SOURCE_DIR}/include/mmg/mmgs/libmmgs.h
      ${CMAKE_SOURCE_DIR}/include/mmg/mmgs/libmmgtypes.h )
  ENDIF()
  IF ( (NOT LIBMMG3D_STATIC) AND (NOT LIBMMG3D_SHARED) )
    ADD_CUSTOM_COMMAND(OUTPUT ${MMG3D_INCLUDE}/libmmg3df.h
      COMMAND ${CMAKE_COMMAND} -E copy  ${MMG3D_BINARY_DIR}/libmmg3df.h ${MMG3D_INCLUDE}/libmmg3df.h
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      DEPENDS ${MMG3D_BINARY_DIR}/libmmg3df.h)
    ADD_CUSTOM_TARGET ( copy_libmmg3df ALL
      DEPENDS ${MMG3D_INCLUDE}/libmmg3df.h )
    ADD_DEPENDENCIES ( copy_libmmg3df mmg3d_fortran_header )

    ADD_CUSTOM_TARGET(copy_3d_headers ALL
      DEPENDS
      copy_libmmg3df copy_libmmgtypesf
      ${CMAKE_SOURCE_DIR}/include/mmg/mmg3d/libmmg3d.h
      ${CMAKE_SOURCE_DIR}/include/mmg/mmg3d/libmmgtypes.h )

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
    copy_2d_headers copy_s_headers copy_3d_headers copy_libmmgtypesf
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
