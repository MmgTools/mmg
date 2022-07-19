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

# Execut_mmg contains the list of executables to test (for now the mmg3d
# executable is in EXECUT_MMG3D (filled by mmg3d.cmake) and the mmgs executable
# is in EXECUTE_MMGS (filled by mmgs.cmake))
SET(EXECUT_MMG ${EXECUT_MMGS} ${EXECUT_MMG3D})

# Make some files not openable
IF ( EXISTS ${CTEST_OUTPUT_DIR}/unwrittable7.meshb
    AND NOT IS_DIRECTORY ${CTEST_OUTPUT_DIR}/unwrittable7.meshb )
  FILE ( REMOVE ${CTEST_OUTPUT_DIR}/unwrittable7.meshb )
ENDIF ()
IF ( NOT EXISTS ${CTEST_OUTPUT_DIR}/unwrittable7.meshb)
  FILE(MAKE_DIRECTORY ${CTEST_OUTPUT_DIR}/unwrittable7.meshb)
ENDIF()

IF ( EXISTS ${CTEST_OUTPUT_DIR}/unwrittable8.sol
    AND NOT IS_DIRECTORY ${CTEST_OUTPUT_DIR}/unwrittable8.sol )
  FILE ( REMOVE ${CTEST_OUTPUT_DIR}/unwrittable8.sol )
ENDIF ()
IF ( NOT EXISTS ${CTEST_OUTPUT_DIR}/unwrittable8.sol)
  FILE(MAKE_DIRECTORY ${CTEST_OUTPUT_DIR}/unwrittable8.sol)
ENDIF()

# Lists of tests that are common to mmgs and mmg3d
FOREACH(EXEC ${EXECUT_MMG})
  GET_FILENAME_COMPONENT ( SHRT_EXEC ${EXEC} NAME )

  ###############################################################################
  #####
  #####         Input/Output
  #####
  ###############################################################################

  # Gmsh without metric: see mmg3d_tests.cmake and mmgs_tests.cmake

  # Binary gmsh iso metric
  MMG_ADD_TEST(mmg_binary_gmsh_iso_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/GmshInout" "iso.mshb"
    )

  # Ascii gmsh iso metric
  MMG_ADD_TEST(mmg_ascii_gmsh_iso_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/GmshInout" "iso.msh"
    )

  # Binary gmsh iso metric
  MMG_ADD_TEST(mmg_binary_gmsh_ani_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/GmshInout" "ani.mshb"
    )

  # Ascii gmsh iso metric
  MMG_ADD_TEST(mmg_ascii_gmsh_ani_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/GmshInout" "ani.msh"
    )

  # VTK .vtk no metric
  MMG_ADD_TEST(mmg_vtkvtk_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/VtkInout" "c1.vtk"
    )

  # VTK .vtu no metric
  MMG_ADD_TEST(mmg_vtkvtu_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/VtkInout" "c1.vtu"
    )

  # VTK .vtk with iso metric
  MMG_ADD_TEST(mmg_vtkvtk_iso_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/VtkInout" "iso.vtk"
    )

  # VTK .vtu with iso metric
  MMG_ADD_TEST(mmg_vtkvtu_iso_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/VtkInout" "iso.vtu"
    )

  # VTK .vtk with aniso metric
  MMG_ADD_TEST(mmg_vtkvtk_ani_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/VtkInout" "ani.vtk"
    )

  # VTK .vtu with aniso metric
  MMG_ADD_TEST(mmg_vtkvtu_ani_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/VtkInout" "ani.vtu"
    )

  IF ( (NOT VTK_FOUND) OR USE_VTK MATCHES OFF )
    SET(expr "VTK library not founded")
    SET_PROPERTY(TEST mmg_vtkvtk_${SHRT_EXEC}
      PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
    SET_PROPERTY(TEST mmg_vtkvtu_${SHRT_EXEC}
      PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
    SET_PROPERTY(TEST mmg_vtkvtk_iso_${SHRT_EXEC}
      PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
    SET_PROPERTY(TEST mmg_vtkvtu_iso_${SHRT_EXEC}
      PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
    SET_PROPERTY(TEST mmg_vtkvtk_ani_${SHRT_EXEC}
      PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
    SET_PROPERTY(TEST mmg_vtkvtu_ani_${SHRT_EXEC}
      PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  ENDIF ( )

  ##############################################################################
  #####
  #####         Aniso test case
  #####
  ##############################################################################
  #####

  SET ( test_names
    mmg_TorusholesAni_${SHRT_EXEC}
    mmg_TorusholesAni_chocCyl_${SHRT_EXEC}
    )

  SET ( input_paths
    ${MMG_CI_TESTS}/TorusholesAni
    ${MMG_CI_TESTS}/TorusholesAni_chocCyl
    )

  SET ( input_files
    torusholes
    torusholesTiny
    )

  SET ( args
    "-v 5 -hgrad 1.25"
    "-v 5 -hgrad 1.25"
    )

  ADD_RUN_AGAIN_TESTS ( ${EXEC} "${test_names}" "${args}"
    "${input_paths}" "${input_files}" )

  MMG_ADD_TEST(mmg_CubeVolAni_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/CubeVolAni" "cube"
    )

  MMG_ADD_TEST(mmg_CubeVolAni2_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/CubeVolAni2" "cube"
    )

  MMG_ADD_TEST(mmg_SphereVolAni_${SHRT_EXEC}
    "${EXEC} -v 5 ${common_args}"
    "${MMG_CI_TESTS}/SphereVolAni" "sphere"
    )

  ##############################################################################
  #####
  #####         Check Memory Leak
  #####
  ##############################################################################
  #####
  ADD_TEST(NAME mmg_LeakCheck_AbnormalEnd2_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 ${common_args}
    ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd2/d)
  SET(passRegex "unable to scale mesh:")
  SET_PROPERTY(TEST mmg_LeakCheck_AbnormalEnd2_${SHRT_EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  ADD_TEST(NAME mmg_LeakCheck_AbnormalEnd7_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 ${common_args}
    ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd7/d
    -out ${CTEST_OUTPUT_DIR}/unwrittable7.meshb)
  SET(passRegex "\\*\\* UNABLE TO OPEN.*")
  SET_PROPERTY(TEST mmg_LeakCheck_AbnormalEnd7_${SHRT_EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  ADD_TEST(NAME mmg_LeakCheck_AbnormalEnd8_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 ${common_args}
    ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd8/d
    -out ${CTEST_OUTPUT_DIR}/unwrittable8.meshb)
  SET(passRegex "\\*\\* UNABLE TO OPEN.*.sol")
  SET_PROPERTY(TEST mmg_LeakCheck_AbnormalEnd8_${SHRT_EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  #####
  ADD_TEST(NAME mmg_LeakCheck_args0_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 ${common_args}
    ${MMG_CI_TESTS}/LeakCheck_args0/d
    ${CTEST_OUTPUT_DIR}/mmg_LeakCheck_args0_${SHRT_EXEC}-d.o)
  #####
  ADD_TEST(NAME mmg_LeakCheck_args1_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 ${common_args}
    -in ${MMG_CI_TESTS}/LeakCheck_args1/d -sol
    ${MMG_CI_TESTS}/LeakCheck_args1/dsol.sol
    -out  ${CTEST_OUTPUT_DIR}/mmg_LeakCheck_args1_${SHRT_EXEC}-dout.meshb)

  ##############################################################################
  #####
  #####         Check Local parameters
  #####
  ##############################################################################
  #####
  MMG_ADD_TEST(mmg_HausdLoc_2Spheres${SHRT_EXEC}
    "${EXEC} -v 5 -hgrad 2 ${common_args}"
    "${MMG_CI_TESTS}/HausdLoc_2Spheres" "2spheres"
    )

  #####
  MMG_ADD_TEST(mmg_hminmaxLoc_2Spheres${SHRT_EXEC}
    "${EXEC} -v 5 -hgrad 2 ${common_args}"
    "${MMG_CI_TESTS}/HausdLoc_2Spheres" "2spheres"
    )

  MMG_ADD_TEST(mmg_HausdLoc_2SpheresAni${SHRT_EXEC}
    "${EXEC} -v 5 -hgrad 2 -A ${common_args}"
    "${MMG_CI_TESTS}/HausdLoc_2Spheres" "2spheres"
    )
  #####
  MMG_ADD_TEST(mmg_hminmaxLoc_2SpheresAni${SHRT_EXEC}
    "${EXEC} -v 5 -hgrad 2 -A ${common_args}"
    "${MMG_CI_TESTS}/HausdLoc_2Spheres" "2spheres"
    )

  ##############################################################################
  #####
  #####         Check Precision
  #####
  ##############################################################################
  #####
  ADD_TEST(NAME mmg_MeshVersionFormatted1_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 ${common_args}
    -in ${MMG_CI_TESTS}/MeshVersionFormatted1/d
    -sol ${MMG_CI_TESTS}/MeshVersionFormatted1/dsol.sol
    ${CTEST_OUTPUT_DIR}/mmg_MeshVersionFormatted1_${SHRT_EXEC}-d.o
    )
  #####
  ADD_TEST(NAME mmg_MeshVersionFormatted2_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 ${common_args}
    -in ${MMG_CI_TESTS}/MeshVersionFormatted2/d
    -sol ${MMG_CI_TESTS}/MeshVersionFormatted2/dsol.sol
    ${CTEST_OUTPUT_DIR}/mmg_MeshVersionFormatted2_${SHRT_EXEC}-d.o
    )

  ###############################################################################
  #####
  #####         Options
  #####
  ###############################################################################
  ADD_TEST(NAME mmg_help_${SHRT_EXEC}
    COMMAND ${EXEC} -h
    )
  SET_PROPERTY(TEST mmg_help_${SHRT_EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "File specifications")

  MMG_ADD_TEST(mmg_hsizOption_${SHRT_EXEC}
    "${EXEC} -v 5 -hsiz 0.1 ${common_args}"
    "${MMG_CI_TESTS}/Cube" "cube"
    )

  # hsiz Aniso
  MMG_ADD_TEST(mmg_hsizAni_${SHRT_EXEC}
    "${EXEC} -v 5 -hsiz 0.1 -sol 2 -A ${common_args}"
    "${MMG_CI_TESTS}/TorusholesAni_chocCyl" "torusholesTiny"
    )

  ADD_TEST(NAME mmg_hsizHmax_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 -hsiz 0.1 -hmax 0.05 ${common_args}
    ${MMG_CI_TESTS}/Cube/cube
    -out ${CTEST_OUTPUT_DIR}/mmg_hsizHmax_${SHRT_EXEC}.o.meshb)
  SET(passRegex "Mismatched options")
  SET_PROPERTY(TEST mmg_hsizHmax_${SHRT_EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

  ADD_TEST(NAME mmg_hsizHmin_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 -hsiz 0.1 -hmin 0.2 ${common_args}
    ${MMG_CI_TESTS}/Cube/cube
    -out ${CTEST_OUTPUT_DIR}/mmg_hsizHmin_${SHRT_EXEC}.o.meshb)
  SET_PROPERTY(TEST mmg_hsizHmin_${SHRT_EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

  # Required entities
  MMG_ADD_TEST(mmg_MultiDom_Ellipse_ReqEntities_${SHRT_EXEC}
    "${EXEC} -v 5 -hausd 0.002 ${common_args}"
    "${MMG_CI_TESTS}/MultiDom_Ellipse_ReqEntities" "c.d"
    )

  MMG_ADD_TEST(mmg_MultiDom_Cube_ReqEntities_${SHRT_EXEC}
    "${EXEC} -v 5 -hsiz 0.02 ${common_args}"
    "${MMG_CI_TESTS}/MultiDom_Cube_ReqEntities" "c"
    )

  MMG_ADD_TEST(mmg_MultiDom_Cube_ReqEntities_nosizreq_${SHRT_EXEC}
    "${EXEC} -v 5 -nosizreq -hgradreq -1 ${common_args}"
    "${MMG_CI_TESTS}/MultiDom_Cube_ReqEntities" "c"
    )

  # Too long for github actions on windows (more than 1 hour). Lasts 15mins on winnie
  # (physical windows), 4mins on ubuntu/mac. See gprof-profile.out file of
  # ${MMG_CI_TESTS}/MultiDom_Ellipse_ReqEntities/ (lot of time is passed in gradsiz_ani)
  #
  # MMG_ADD_TEST(mmg_MultiDom_Ellipse_ReqEntitiesAni_${SHRT_EXEC}
  # "${EXEC} -v 5 -hausd 0.002 -A ${common_args}"
  #  "${MMG_CI_TESTS}/MultiDom_Ellipse_ReqEntities/c" "d"
  # )


  # -A
  MMG_ADD_TEST(mmg_CommandLineAni_${SHRT_EXEC}
    "${EXEC} -v 5 -hausd 0.005 -sol 2 -A ${common_args}"
    "${MMG_CI_TESTS}/TorusholesAni_chocCyl" "torusholesTiny"
    )

  # -Optim
  MMG_ADD_TEST(mmg_optimOption_${SHRT_EXEC}
    "${EXEC} -v 5 -hausd 1 -optim ${common_args}"
    "${MMG_CI_TESTS}/MecaPart" "geom_1_before.mesh"
    )

  MMG_ADD_TEST(mmg_optimHmax_${SHRT_EXEC}
    "${EXEC} -v 5 -hausd 1 -optim -hmax 1 ${common_args}"
    "${MMG_CI_TESTS}/MecaPart" "geom_1_before.mesh"
    )

  # -nreg
  MMG_ADD_TEST(mmg_nreg_${SHRT_EXEC}
    "${EXEC} -v 5 -nreg ${common_args}"
    "${MMG_CI_TESTS}/c1" "c1.meshb"
    )

  ##############################################################################
  #####
  #####         Various test cases
  #####
  ##############################################################################
  #####

  # Lot of reference edges, ridges, corners and singularities
  MMG_ADD_TEST(mmg_SurfEdges_${SHRT_EXEC}
    "${EXEC} -v 5 -hausd 0.1 -A ${common_args}"
    "${MMG_CI_TESTS}/SurfEdges_house" "housebad.meshb"
    )

  MMG_ADD_TEST(mmg_SurfEdgesAni_${SHRT_EXEC}
    "${EXEC} -v 5 -hausd 0.1 ${common_args}"
    "${MMG_CI_TESTS}/SurfEdges_house" "housebad.meshb"
    )


ENDFOREACH(EXEC)
