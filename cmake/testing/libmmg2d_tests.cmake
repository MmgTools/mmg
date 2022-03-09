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


###############################################################################
#####
#####        Build Mmg2d library executables and add tests if needed
#####
###############################################################################

SET ( MMG2D_LIB_TESTS
  libmmg2d_adp_example0_a
  libmmg2d_adp_example0_b
  libmmg2d_adp_example1
  libmmg2d_adp_example2
  libmmg2d_gene_example0
  libmmg2d_ls_example0
  libmmg2d_lsOnly
  libmmg2d_lsAndMetric
  )

# Additional tests that needs to download ci meshes
IF ( MMG2D_CI )
  LIST ( APPEND MMG2D_LIB_TESTS test_api2d_0 )
ENDIF ( )

SET ( MMG2D_LIB_TESTS_MAIN_PATH
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0/example0_a/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0/example0_b/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example1/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example2/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/squareGeneration_example0/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsAndMetric/main.c
  ${PROJECT_SOURCE_DIR}/cmake/testing/code/test_met2d.c
  )

# Additional tests that needs to download ci meshes
IF ( MMG2D_CI )
  LIST ( APPEND MMG2D_LIB_TESTS_MAIN_PATH
    ${MMG2D_CI_TESTS}/API_tests/2d.c
    )
ENDIF( )

IF ( LIBMMG2D_STATIC )
  SET ( lib_name lib${PROJECT_NAME}2d_a )
ELSEIF ( LIBMMG2D_SHARED )
  SET ( lib_name lib${PROJECT_NAME}2d_so )
ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ( )

#####         Fortran Tests
IF ( CMAKE_Fortran_COMPILER )
  ENABLE_LANGUAGE ( Fortran )

  SET ( MMG2D_LIB_TESTS ${MMG2D_LIB_TESTS}
    libmmg2d_fortran_a
    libmmg2d_fortran_b
    libmmg2d_fortran_io
    libmmg2d_fortran_lsOnly
    libmmg2d_fortran_lsAndMetric
    test_io2d_fortran
    )
  # Additional tests that needs to download ci meshes
  IF ( MMG2D_CI )
    LIST ( APPEND MMG2D_LIB_TESTS test_api2d_fortran_0 )
  ENDIF( )

  SET ( MMG2D_LIB_TESTS_MAIN_PATH ${MMG2D_LIB_TESTS_MAIN_PATH}
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0_fortran/example0_a/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0_fortran/example0_b/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsAndMetric/main.F90
    ${PROJECT_SOURCE_DIR}/cmake/testing/code/mmg2d_io.F90
    )
  # Additional tests that needs to download ci meshes
  IF ( MMG2D_CI )
    LIST ( APPEND MMG2D_LIB_TESTS_MAIN_PATH
      ${MMG2D_CI_TESTS}/API_tests/2d.F90
      )
  ENDIF( )


ENDIF ( CMAKE_Fortran_COMPILER )

LIST(LENGTH MMG2D_LIB_TESTS nbTests_tmp)
MATH(EXPR nbTests "${nbTests_tmp} - 1")

FOREACH ( test_idx RANGE ${nbTests} )
  LIST ( GET MMG2D_LIB_TESTS           ${test_idx} test_name )
  LIST ( GET MMG2D_LIB_TESTS_MAIN_PATH ${test_idx} main_path )

  ADD_LIBRARY_TEST ( ${test_name} ${main_path} copy_2d_headers ${lib_name} )

ENDFOREACH ( )

IF ( BUILD_TESTING )

  SET(LIBMMG2D_ADP0_a ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_adp_example0_a )
  SET(LIBMMG2D_ADP0_b ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_adp_example0_b )
  SET(LIBMMG2D_ADP1 ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_adp_example1 )
  SET(LIBMMG2D_ADP2 ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_adp_example2 )
  SET(LIBMMG2D_GENE0 ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_gene_example0 )
  SET(LIBMMG2D_LS0 ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_ls_example0 )
  SET(LIBMMG2D_LSONLY ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_lsOnly )
  SET(LIBMMG2D_LSANDMETRIC ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_lsAndMetric )
  SET(TEST_API2D_EXEC0 ${EXECUTABLE_OUTPUT_PATH}/test_api2d_0)
  SET(TEST_MET2D ${EXECUTABLE_OUTPUT_PATH}/test_met2d)


  ADD_TEST(NAME libmmg2d_adp_example0_a   COMMAND ${LIBMMG2D_ADP0_a}
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0/example0_a/init.mesh"
    "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_0_a-init.o"
    )
  ADD_TEST(NAME libmmg2d_adp_example0_b   COMMAND ${LIBMMG2D_ADP0_b}
    "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_0_b.o.mesh"
    )
  ADD_TEST(NAME libmmg2d_adp_example1   COMMAND ${LIBMMG2D_ADP1}
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example1/dom.mesh"
    "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_1-dom.o"
    )
  ADD_TEST(NAME libmmg2d_adp_example2   COMMAND ${LIBMMG2D_ADP2}
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example1/dom.mesh"
    "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_2-dom.o"
    "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_2-dom-end.o"
    )
  ADD_TEST(NAME libmmg2d_gene_example0   COMMAND ${LIBMMG2D_GENE0}
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/squareGeneration_example0/carretest.mesh"
    "${CTEST_OUTPUT_DIR}/libmmg2d_Generation_0-carre.o"
    )
  ADD_TEST(NAME libmmg2d_ls0_io_0   COMMAND ${LIBMMG2D_LS0}
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/naca-multiSols.mesh"
    "${CTEST_OUTPUT_DIR}/libmmg2d_io_0-naca.o" "0"
    )
  ADD_TEST(NAME libmmg2d_ls0_io_1   COMMAND ${LIBMMG2D_LS0}
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/naca-multiSols.mesh"
    "${CTEST_OUTPUT_DIR}/libmmg2d_io_0-naca.o" "1"
    )
  ADD_TEST(NAME libmmg2d_lsOnly   COMMAND ${LIBMMG2D_LSONLY}
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat.mesh"
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
    "${CTEST_OUTPUT_DIR}/libmmg2d_lsOnly_multimat.o"
    )
  ADD_TEST(NAME libmmg2d_lsAndMetric   COMMAND ${LIBMMG2D_LSANDMETRIC}
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat.mesh"
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
    "${CTEST_OUTPUT_DIR}/libmmg2d_lsAndMetric_multimat.o"
    )
  ADD_TEST(NAME test_api2d_0   COMMAND ${TEST_API2D_EXEC0}
    "${MMG2D_CI_TESTS}/API_tests/2dom.mesh"
    "${CTEST_OUTPUT_DIR}/test_API2d.o"
    )
  ADD_TEST(NAME test_met2d   COMMAND ${TEST_MET2D}
    )


  IF ( CMAKE_Fortran_COMPILER)
    SET(LIBMMG2D_EXECFORTRAN_a ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_a )
    SET(LIBMMG2D_EXECFORTRAN_b ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_b )
    SET(LIBMMG2D_EXECFORTRAN_IO ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_io )
    SET(LIBMMG2D_EXECFORTRAN_LSONLY ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_lsOnly )
    SET(LIBMMG2D_EXECFORTRAN_LSANDMETRIC ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_lsAndMetric )
    SET(TEST_API2D_FORTRAN_EXEC0 ${EXECUTABLE_OUTPUT_PATH}/test_api2d_fortran_0)
    SET(TEST_IO2D_FORTRAN_EXEC ${EXECUTABLE_OUTPUT_PATH}/test_io2d_fortran)


    ADD_TEST(NAME libmmg2d_fortran_a   COMMAND ${LIBMMG2D_EXECFORTRAN_a}
      "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0_fortran/example0_a/init.mesh"
      "${CTEST_OUTPUT_DIR}/libmmg2d-Adaptation_Fortran_0_a-init.o"
      )
    ADD_TEST(NAME libmmg2d_fortran_b   COMMAND ${LIBMMG2D_EXECFORTRAN_b}
      "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_Fortran_0_b.o"
      )
    ADD_TEST(NAME libmmg2d_fortran_io_0   COMMAND ${LIBMMG2D_EXECFORTRAN_IO}
      "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/naca-multiSols.mesh"
      "${CTEST_OUTPUT_DIR}/libmmg2d_Fortran_io-naca.o" "0"
      )
    ADD_TEST(NAME libmmg2d_fortran_io_1   COMMAND ${LIBMMG2D_EXECFORTRAN_IO}
      "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/naca-multiSols.mesh"
      "${CTEST_OUTPUT_DIR}/libmmg2d_Fortran_io-naca.o" "1"
      )
    ADD_TEST(NAME libmmg2d_fortran_lsOnly   COMMAND ${LIBMMG2D_EXECFORTRAN_LSONLY}
      "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat.mesh"
      "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
      "${CTEST_OUTPUT_DIR}/libmmg2d_lsOnly_multimat.o"
      )
    ADD_TEST(NAME libmmg2d_fortran_lsAndMetric   COMMAND ${LIBMMG2D_EXECFORTRAN_LSANDMETRIC}
      "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat.mesh"
      "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
      "${CTEST_OUTPUT_DIR}/libmmg2d_lsAndMetric_multimat.o"
      )
    ADD_TEST(NAME test_api2d_fortran_0   COMMAND ${TEST_API2D_FORTRAN_EXEC0}
      "${MMG2D_CI_TESTS}/API_tests/2dom.mesh"
      "${CTEST_OUTPUT_DIR}/test_API2d.o"
      )
    ADD_TEST(NAME test_io2d_fortran_scalar   COMMAND ${TEST_IO2D_FORTRAN_EXEC}
      "${MMG2D_CI_TESTS}/Hybrid/hybrid.mesh"
      "${CTEST_OUTPUT_DIR}/hybrid-2d-scal.o" 0
      )
    ADD_TEST(NAME test_io2d_fortran_array   COMMAND ${TEST_IO2D_FORTRAN_EXEC}
      "${MMG2D_CI_TESTS}/Hybrid/hybrid.mesh"
      "${CTEST_OUTPUT_DIR}/hybrid-2d-array.o" 1
      )

  ENDIF()
ENDIF ( )
