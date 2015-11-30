## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
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
  ADD_LIBRARY(${PROJECT_NAME}_a  STATIC
    ${MMG2D_SOURCE_DIR}/lib${PROJECT_NAME}2df.h
    ${sourcemmg2d_files} ${libmmg2d_file}
    ${MMGS_SOURCE_DIR}/lib${PROJECT_NAME}sf.h
    ${sourcemmgs_files} ${libmmgs_file}
    ${MMG3D_SOURCE_DIR}/lib${PROJECT_NAME}3df.h
    ${source_files} ${lib_file}
    ${CMAKE_SOURCE_DIR}/src/mmg/libmmg.h
    ${CMAKE_SOURCE_DIR}/src/mmg/libmmgf.h
    )
  SET_TARGET_PROPERTIES(${PROJECT_NAME}_a PROPERTIES OUTPUT_NAME
    ${PROJECT_NAME})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}_a ${LIBRARIES})
  INSTALL(TARGETS ${PROJECT_NAME}_a
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
ENDIF()

# Compile shared library
IF ( LIBMMG_SHARED )
  ADD_LIBRARY(${PROJECT_NAME}_so SHARED
    ${MMG2D_SOURCE_DIR}/lib${PROJECT_NAME}2df.h
    ${sourcemmg2d_files} ${libmmg2d_file}
    ${MMGS_SOURCE_DIR}/lib${PROJECT_NAME}sf.h
    ${sourcemmgs_files} ${libmmgs_file}
    ${MMG3D_SOURCE_DIR}/lib${PROJECT_NAME}3df.h
    ${source_files} ${lib_file}
    ${CMAKE_SOURCE_DIR}/src/mmg/libmmg.h
    ${CMAKE_SOURCE_DIR}/src/mmg/libmmgf.h
    )
  SET_TARGET_PROPERTIES(${PROJECT_NAME}_so PROPERTIES
    VERSION ${CMAKE_RELEASE_VERSION} SOVERSION 5)
  SET_TARGET_PROPERTIES(${PROJECT_NAME}_so PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}_so ${LIBRARIES})
  INSTALL(TARGETS ${PROJECT_NAME}_so
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
ENDIF()

IF ( LIBMMG_STATIC OR LIBMMG_SHARED )
  # mmg header files needed for library
  SET( mmg2d_headers
    ${MMG2D_SOURCE_DIR}/libmmg2d.h
    ${MMG2D_SOURCE_DIR}/libmmg2df.h
    )
  SET( mmg3d_headers
    ${MMG3D_SOURCE_DIR}/libmmg3d.h
    ${MMG3D_SOURCE_DIR}/libmmg3df.h
    )
  SET( mmgs_headers
    ${MMGS_SOURCE_DIR}/libmmgs.h
    ${MMGS_SOURCE_DIR}/libmmgsf.h
    )
  SET( mmg_headers
    ${COMMON_SOURCE_DIR}/mmgcommon.h
    ${COMMON_SOURCE_DIR}/eigenv.h
    ${COMMON_SOURCE_DIR}/libmmgcommon.h
    ${COMMON_SOURCE_DIR}/libmmgcommonf.h
    ${COMMON_SOURCE_DIR}/chrono.h
    ${CMAKE_SOURCE_DIR}/src/mmg/libmmg.h
    ${CMAKE_SOURCE_DIR}/src/mmg/libmmgf.h
    )
  # Install header files in /usr/local or equivalent
  INSTALL(FILES ${mmg2d_headers} DESTINATION include/mmg/mmg2d)
  INSTALL(FILES ${mmgs_headers} DESTINATION include/mmg/mmgs)
  INSTALL(FILES ${mmg3d_headers} DESTINATION include/mmg/mmg3d)
  INSTALL(FILES ${mmg_headers} DESTINATION include/mmg)
  # Install header files in project directory
  INSTALL(FILES ${mmg2d_headers} DESTINATION ${CMAKE_SOURCE_DIR}/include/mmg/mmg2d)
  INSTALL(FILES ${mmgs_headers} DESTINATION ${CMAKE_SOURCE_DIR}/include/mmg/mmgs)
  INSTALL(FILES ${mmg3d_headers} DESTINATION ${CMAKE_SOURCE_DIR}/include/mmg/mmg3d)
  INSTALL(FILES ${mmg_headers} DESTINATION ${CMAKE_SOURCE_DIR}/include/mmg)
ENDIF()

############################################################################
#####
#####         Compile program to test library
#####
############################################################################

IF ( TEST_LIBMMG )
  INCLUDE(cmake/testing/libmmg_tests.cmake)
ENDIF()
