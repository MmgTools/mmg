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
## Compilation of mmg3d executable, libraries and tests
##
## =============================================================================

SET(MMG3D_SOURCE_DIR      ${PROJECT_SOURCE_DIR}/src/mmg3d)
SET(MMG3D_BINARY_DIR      ${PROJECT_BINARY_DIR}/src/mmg3d)
SET(MMG3D_SHRT_INCLUDE    mmg/mmg3d )
SET(MMG3D_INCLUDE         ${PROJECT_BINARY_DIR}/include/${MMG3D_SHRT_INCLUDE} )

FILE(MAKE_DIRECTORY ${MMG3D_BINARY_DIR})

############################################################################
#####
#####         Fortran header: libmmg3df.h
#####
############################################################################


GENERATE_FORTRAN_HEADER ( mmg3d
  ${MMG3D_SOURCE_DIR} libmmg3d.h
  ${MMG3D_SHRT_INCLUDE}
  ${MMG3D_BINARY_DIR} libmmg3df.h
  )

############################################################################
#####
#####         Choose type of insertion: delaunay kernel or patterns
#####
############################################################################
# Should we use patterns for insertion istead of delaunay kernel
IF ( NOT PATTERN )
  MESSAGE(STATUS "Vertex insertion by delaunay kernel")
ELSE()
  MESSAGE(STATUS "Vertex insertion by patterns")
  SET(CMAKE_C_FLAGS "-DPATTERN ${CMAKE_C_FLAGS}")
ENDIF()

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################

# Library files
FILE(
  GLOB
  mmg3d_library_files
  ${MMG3D_SOURCE_DIR}/*.c
  ${COMMON_SOURCE_DIR}/*.c
  ${MMG3D_SOURCE_DIR}/inoutcpp_3d.cpp
  )
LIST(REMOVE_ITEM mmg3d_library_files
  ${MMG3D_SOURCE_DIR}/${PROJECT_NAME}3d.c
  )

IF ( VTK_FOUND )
  LIST(APPEND  mmg3d_library_files
    ${COMMON_SOURCE_DIR}/vtkparser.cpp)
ENDIF ( )

FILE(
  GLOB
  mmg3d_main_file
  ${MMG3D_SOURCE_DIR}/mmg3d.c
  )

############################################################################
#####
#####         Elastic
#####
############################################################################

IF( ELAS_FOUND )
  # Set flags for building test program
  INCLUDE_DIRECTORIES(PUBLIC ${ELAS_INCLUDE_DIR})

  SET(CMAKE_REQUIRED_INCLUDES ${ELAS_INCLUDE_DIR})
  SET(CMAKE_REQUIRED_LIBRARIES ${ELAS_LIBRARY})

  SET(CMAKE_C_FLAGS "-DUSE_ELAS ${CMAKE_C_FLAGS}")
  MESSAGE(STATUS
    "Compilation with the Elas library: ${ELAS_LIBRARY} ")
  SET( LIBRARIES ${ELAS_LINK_FLAGS} ${LIBRARIES})
  SET( LIBRARIES ${ELAS_LIBRARY} ${LIBRARIES})
ENDIF()


############################################################################
#####
#####         Compile mmg3d libraries
#####
############################################################################

# Compile static library
IF ( LIBMMG3D_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}3d_a STATIC copy_3d_headers
    "${mmg3d_library_files}" ${PROJECT_NAME}3d )
ENDIF()

# Compile shared library
IF ( LIBMMG3D_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}3d_so SHARED copy_3d_headers
    "${mmg3d_library_files}" ${PROJECT_NAME}3d )
ENDIF()

# mmg3d header files needed for library
SET( mmg3d_headers
  ${MMG3D_SOURCE_DIR}/libmmg3d.h
  ${MMG3D_BINARY_DIR}/libmmg3df.h
  ${COMMON_SOURCE_DIR}/libmmgtypes.h
  ${COMMON_BINARY_DIR}/libmmgtypesf.h
  ${COMMON_BINARY_DIR}/mmgcmakedefines.h
  ${COMMON_BINARY_DIR}/mmgversion.h
  )
IF (NOT WIN32 OR MINGW)
  LIST(APPEND mmg3d_headers  ${COMMON_BINARY_DIR}/git_log_mmg.h )
ENDIF()

# install man pages
INSTALL(FILES ${PROJECT_SOURCE_DIR}/doc/man/mmg3d.1.gz DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# Install header files in /usr/local or equivalent
INSTALL(FILES ${mmg3d_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmg3d COMPONENT headers)

# Copy header files in project directory at build step
COPY_HEADERS_AND_CREATE_TARGET ( ${MMG3D_SOURCE_DIR} ${MMG3D_BINARY_DIR} ${MMG3D_INCLUDE} 3d )

############################################################################
#####
#####         Compile program to test library
#####
############################################################################
SET(MMG3D_CI_TESTS ${CI_DIR}/mmg3d )
SET(MMG_CI_TESTS ${CI_DIR}/mmg )

IF ( TEST_LIBMMG3D )
  INCLUDE(cmake/testing/libmmg3d_tests.cmake)
ENDIF()

###############################################################################
#####
#####         Compile MMG3D executable
#####
###############################################################################
ADD_AND_INSTALL_EXECUTABLE ( ${PROJECT_NAME}3d copy_3d_headers
  "${mmg3d_library_files}" ${mmg3d_main_file} )

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

IF ( BUILD_TESTING )
  ##-------------------------------------------------------------------##
  ##------- Set the continuous integration options --------------------##
  ##-------------------------------------------------------------------##

  ##-------------------------------------------------------------------##
  ##--------------------------- Add tests and configure it ------------##
  ##-------------------------------------------------------------------##
  # Add runtime that we want to test for mmg3d
  IF ( MMG3D_CI )

    IF ( LONG_TESTS )
      # Run some tests twice with the output of the previous test as input
      OPTION ( RUN_AGAIN "Enable/Disable second run of some tests" ON )
      MARK_AS_ADVANCED(RUN_AGAIN)
    ELSE ( )
      SET ( RUN_AGAIN OFF )
    ENDIF ( )

    ADD_EXEC_TO_CI_TESTS ( ${PROJECT_NAME}3d EXECUT_MMG3D )
    SET ( LISTEXEC_MMG ${EXECUT_MMG3D} )

    IF ( TEST_LIBMMG3D )
      SET(LIBMMG3D_EXEC0_a ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example0_a)
      SET(LIBMMG3D_EXEC0_b ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example0_b)
      SET(LIBMMG3D_EXEC1   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example1)
      SET(LIBMMG3D_EXEC2   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example2)
      SET(LIBMMG3D_EXEC4   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example4)
      SET(LIBMMG3D_EXEC5   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example5)
      SET(LIBMMG3D_EXEC6   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example6_io)
      SET(LIBMMG3D_GENERICIO ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_generic_io)
      SET(LIBMMG3D_LSONLY ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_lsOnly )
      SET(LIBMMG3D_LSANDMETRIC ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_lsAndMetric )
      SET(TEST_API3D_EXEC0 ${EXECUTABLE_OUTPUT_PATH}/test_api3d_0)
      SET(TEST_API3D_DOMSEL ${EXECUTABLE_OUTPUT_PATH}/test_api3d_domain-selection)
      SET(TEST_API3D_VTK2MESH ${EXECUTABLE_OUTPUT_PATH}/test_api3d_vtk2mesh)
      SET(TEST_MET3D ${EXECUTABLE_OUTPUT_PATH}/test_met3d)
      SET(TEST_COMPARE_PARA_TRIA ${EXECUTABLE_OUTPUT_PATH}/test_compare-para-tria)

      ADD_TEST(NAME libmmg3d_example0_a COMMAND ${LIBMMG3D_EXEC0_a}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0/example0_a/cube.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_0_a-cube.o"
        )
      ADD_TEST(NAME libmmg3d_example0_b COMMAND ${LIBMMG3D_EXEC0_b}
       "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_0_b.o.mesh"
        )
      ADD_TEST(NAME libmmg3d_example1   COMMAND ${LIBMMG3D_EXEC1}
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_1.o.mesh"
        )
      ADD_TEST(NAME libmmg3d_example2   COMMAND ${LIBMMG3D_EXEC2}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/adaptation_example2/2spheres.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_1-2spheres_1.o"
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_1-2spheres_2.o"
        )
      IF ( ELAS_FOUND )
        ADD_TEST(NAME libmmg3d_example4   COMMAND ${LIBMMG3D_EXEC4}
          "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/LagrangianMotion_example0/tinyBoxt"
          "${CTEST_OUTPUT_DIR}/libmmg3d_LagrangianMotion_0-tinyBoxt.o"
          )
      ENDIF ()
      ADD_TEST(NAME libmmg3d_example6_io_0   COMMAND ${LIBMMG3D_EXEC6}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_multisols_example6/torus.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_io_6-naca.o" "0"
       )
      ADD_TEST(NAME libmmg3d_example6_io_1   COMMAND ${LIBMMG3D_EXEC6}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_multisols_example6/torus.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_io_6-naca.o" "1"
       )
      ADD_TEST(NAME libmmg3d_lsOnly   COMMAND ${LIBMMG3D_LSONLY}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/IsosurfDiscretization_lsOnly/plane.mesh"
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/IsosurfDiscretization_lsOnly/m.sol"
        "${CTEST_OUTPUT_DIR}/libmmg3d_lsOnly_multimat.o"
        )
      ADD_TEST(NAME libmmg3d_lsAndMetric   COMMAND ${LIBMMG3D_LSANDMETRIC}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/IsosurfDiscretization_lsOnly/plane.mesh"
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/IsosurfDiscretization_lsOnly/m.sol"
        "${CTEST_OUTPUT_DIR}/libmmg3d_lsAndMetric_multimat.o"
        )
      ADD_TEST(NAME test_api3d_0   COMMAND ${TEST_API3D_EXEC0}
        "${MMG3D_CI_TESTS}/API_tests/2dom.mesh"
        "${CTEST_OUTPUT_DIR}/test_API3d.o"
       )
      ADD_TEST(NAME test_api3d_domain-selection   COMMAND ${TEST_API3D_DOMSEL}
        "${MMG3D_CI_TESTS}/OptLs_plane/plane.mesh"
        "${MMG3D_CI_TESTS}/OptLs_plane/p.sol"
        "${CTEST_OUTPUT_DIR}/test_API3d-domsel-whole.o"
        "${CTEST_OUTPUT_DIR}/test_API3d-domsel-dom2.o"
        )
      ADD_TEST(NAME test_api3d_vtk2mesh   COMMAND ${TEST_API3D_VTK2MESH}
        "${MMG3D_CI_TESTS}/API_tests/cellsAndNode-data.vtk"
        "${CTEST_OUTPUT_DIR}/test_API3d-vtk2mesh.o"
        )
      ADD_TEST(NAME test_met3d   COMMAND ${TEST_MET3D}
        )
      ADD_TEST(NAME libmmg3d_generic_io_msh   COMMAND ${LIBMMG3D_GENERICIO}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_generic_and_get_adja/cube.msh"
        "${CTEST_OUTPUT_DIR}/cube.o.msh" "1"
        )
      ADD_TEST(NAME libmmg3d_generic_io_mesh   COMMAND ${LIBMMG3D_GENERICIO}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_generic_and_get_adja/cube.mesh"
        "${CTEST_OUTPUT_DIR}/cube.o.mesh" "1"
        )
      ADD_TEST(NAME libmmg3d_generic_io_vtk   COMMAND ${LIBMMG3D_GENERICIO}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_generic_and_get_adja/cube.vtk"
        "${CTEST_OUTPUT_DIR}/cube.o.vtk" "1"
        )
      ADD_TEST(NAME libmmg3d_generic_io_vtu   COMMAND ${LIBMMG3D_GENERICIO}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_generic_and_get_adja/cube.vtu"
        "${CTEST_OUTPUT_DIR}/cube.o.vtu" "1"
        )

      IF ( NOT VTK_FOUND )
        SET(expr "VTK library not founded")
        SET_PROPERTY(TEST test_api3d_vtk2mesh
          PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
        SET_PROPERTY(TEST libmmg3d_generic_io_vtk
          PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
        SET_PROPERTY(TEST libmmg3d_generic_io_vtu
          PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
      ENDIF ( )

      ADD_TEST(NAME test_para_tria
        COMMAND ${EXECUT_MMG3D}
        -ar 0.02 -nofem -nosizreq -hgradreq -1 -hgrad -1
        ${MMG3D_CI_TESTS}/test_para_tria/proc0.mesh
        -sol ${MMG3D_CI_TESTS}/test_para_tria/proc0.sol
        ${CTEST_OUTPUT_DIR}/proc0.o.mesh
        )

      SET_TESTS_PROPERTIES ( test_para_tria
        PROPERTIES FIXTURES_SETUP test_para_tria )

      ADD_TEST(NAME test_compare_para_tria
        COMMAND ${TEST_COMPARE_PARA_TRIA}
        ${MMG3D_CI_TESTS}/test_para_tria/proc0.mesh
        ${CTEST_OUTPUT_DIR}/proc0.o.mesh
        )
      SET_TESTS_PROPERTIES ( test_compare_para_tria
        PROPERTIES FIXTURES_REQUIRED test_para_tria )

      IF ( CMAKE_Fortran_COMPILER)
        SET(LIBMMG3D_EXECFORTRAN_a  ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_fortran_a)
        SET(LIBMMG3D_EXECFORTRAN_b  ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_fortran_b)
        SET(LIBMMG3D_EXECFORTRAN_IO ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_fortran_io)
        SET(LIBMMG3D_EXECFORTRAN_LSONLY ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_fortran_lsOnly )
        SET(LIBMMG3D_EXECFORTRAN_LSANDMETRIC ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_fortran_lsAndMetric )
        SET(TEST_API3D_FORTRAN_EXEC0 ${EXECUTABLE_OUTPUT_PATH}/test_api3d_fortran_0)


        ADD_TEST(NAME libmmg3d_fortran_a  COMMAND ${LIBMMG3D_EXECFORTRAN_a}
          "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0_fortran/example0_a/cube.mesh"
          "${CTEST_OUTPUT_DIR}/libmmg3d-Adaptation_Fortran_0_a-cube.o"
          )
        ADD_TEST(NAME libmmg3d_fortran_b  COMMAND ${LIBMMG3D_EXECFORTRAN_b}
          "${CTEST_OUTPUT_DIR}/libmmg3d-Adaptation_Fortran_0_b-cube.o"
          )
        ADD_TEST(NAME libmmg3d_fortran_io_0   COMMAND ${LIBMMG3D_EXECFORTRAN_IO}
          "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_multisols_example6/torus.mesh"
          "${CTEST_OUTPUT_DIR}/libmmg3d_Fortran_io-torus.o" "0"
          )
        ADD_TEST(NAME libmmg3d_fortran_io_1   COMMAND ${LIBMMG3D_EXECFORTRAN_IO}
          "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_multisols_example6/torus.mesh"
          "${CTEST_OUTPUT_DIR}/libmmg3d_Fortran_io-torus.o" "1"
          )
       ADD_TEST(NAME libmmg3d_fortran_lsOnly3d   COMMAND ${LIBMMG3D_EXECFORTRAN_LSONLY}
         "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/IsosurfDiscretization_lsOnly/plane.mesh"
         "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/IsosurfDiscretization_lsOnly/m.sol"
         "${CTEST_OUTPUT_DIR}/libmmg3d_lsOnly_multimat.o" )

       ADD_TEST(NAME libmmg3d_fortran_lsAndMetric3d   COMMAND ${LIBMMG3D_EXECFORTRAN_LSANDMETRIC}
          "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/IsosurfDiscretization_lsOnly/plane.mesh"
          "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/IsosurfDiscretization_lsOnly/m.sol"
          "${CTEST_OUTPUT_DIR}/libmmg3d_lsAndMetric_multimat.o" )
       ADD_TEST(NAME test_api3d_fortran_0   COMMAND ${TEST_API3D_FORTRAN_EXEC0}
         "${MMG3D_CI_TESTS}/API_tests/2dom.mesh"
         "${CTEST_OUTPUT_DIR}/test_API3d.o"
         )

      ENDIF()

    ENDIF ( TEST_LIBMMG3D )

    IF ( ONLY_VERY_SHORT_TESTS )
      SET ( CTEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/TEST_OUTPUTS )

      ADD_TEST(NAME mmg3d_very_short COMMAND ${EXECUT_MMG3D}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0/example0_a/cube.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_0_a-cube.o"
        )
    ELSE ( )

      # Add more tests
      INCLUDE( ${PROJECT_SOURCE_DIR}/cmake/testing/mmg3d_tests.cmake )
      INCLUDE( ${PROJECT_SOURCE_DIR}/cmake/testing/mmg_tests.cmake )
    ENDIF ( )

  ENDIF ( MMG3D_CI )

ENDIF ( BUILD_TESTING )
