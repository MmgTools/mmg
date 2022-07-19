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

GET_FILENAME_COMPONENT ( SHRT_EXECUT_MMG3D ${EXECUT_MMG3D} NAME )

##############################################################################
#####
#####         Tests that may be run twice
#####
##############################################################################

SET ( test_names
  # Simple test: must already pass
  mmg3d_SimpleCube_fast
  # MultiDomain
  mmg3d_MultiDom_Ellipse_fast
  # Non-manifold test case
  mmg3d_NM_Cube_fast
  mmg3d_NM_Complex_fast
  # test case with non-manifold, ridges, ref edges and a curve surface
  mmg3d_NM_cone_fast
  # mmg3d_NM_cone_ani_fast #Fail because at second run a tetra we have a tet with 4 ridge vertices
  )

SET ( input_paths
  ${MMG3D_CI_TESTS}/Cube
  ### Multidom
  ${MMG3D_CI_TESTS}/MultiDom_Ellipse
   ### non-manifold
  ${MMG3D_CI_TESTS}/NM_Cube
  ${MMG3D_CI_TESTS}/NM_Complex
  ${MMG3D_CI_TESTS}
  #${MMG3D_CI_TESTS}/cone-nm.mesh
  )

SET ( input_files
  cube
  ### Multidom
  c.d
   ### non-manifold
  nm
  nm4
  cone-nm.mesh
  #cone-nm.mesh
  )

SET ( args
  "-v 5"
  ### MultiDomain
  "-v 5 -hausd 0.002"
  ### non-manifold
  "-v 5 -hmax 0.1"
  "-v 5"
  "-v 5"
  #"-v 5 -A"
  )

IF ( LONG_TESTS )
  SET ( test_names ${test_names}
    # Check what happend when we refine an isotropic cube of size h with a
    # constant metric (h, h/2, h/4, h/8 and h/16) ---First with hmin=hmax
    mmg3d_CubeIso_h_hminMax
    mmg3d_CubeIso_0.5h_hminMax
    mmg3d_CubeIso_0.25h_hminMax
    #---Second with sol file
    mmg3d_CubeIso_h_met
    mmg3d_CubeIso_0.5h_met
    mmg3d_CubeIso_0.25h_met
    mmg3d_CubeIso_0.125h_met
    mmg3d_CubeAniIso_0.125h_met
    #####
    mmg3d_SphereIso_h_met
    mmg3d_SphereIso_0.5h_met
    mmg3d_SphereIso_0.25h_met
    mmg3d_SphereIso_0.125h_met
    mmg3d_SphereIso_0.020_met
    # mmg3d_SphereIso_0.020-0.015_met # not enough mem on windows 4G
    mmg3d_SphereAni_0.02
    # Check what happend when we unrefine a sphere of size smallh with a
    # constant metric (2*smallh, 4*smallh and 8*smallh)
    mmg3d_SphereIso_2smallh_met
    mmg3d_SphereIso_4smallh_met
    mmg3d_SphereIso_8smallh_met
    # Check what happend when we use hausdorff number to refine the skin and a
    # big hgrad to have an inside of the initial size (0.5)
    mmg3d_SphereIso_h_hausd0.001
    mmg3d_SphereIso_h_hausd0.005
    # Check what happend when we refine a cube whose skin has already the good size
    mmg3d_CubeSkin0.05_Inside0.4
    mmg3d_CubeSkin0.1_Inside0.4
    mmg3d_CubeSkin0.2_Inside0.4
    mmg3d_CubeSkin0.0125_Inside0.125
    # mmg3d_CubeSkin0.0125_Inside0.25 # too long on OSX
    # mmg3d_CubeSkin0.0125_Inside0.5 # too long on all machine
    # Check results on various meshes
    # First: Meshes that we want unrefined
    mmg3d_Various_unref_Linkrods_met0.2
    mmg3d_Various_unref_Linkrods_met0.2_hausd0.01
    # Second: Meshes that we want refined
    mmg3d_Various_ref_Linkrods_met0.05
    mmg3d_Various_ref_Linkrods_met0.05_hausd0.01
    mmg3d_Various_ref_Linkrods_met0.05_hausd0.001
    # Third: We refine some parts and unrefined others
    mmg3d_Various_refunref_Santa_met0.05_hausd0.001_ar90
    mmg3d_Various_refunref_Santa_met0.05_hausd0.0001_ar90
    # 5: MultiDomain
    mmg3d_MultiDom_Cube
    mmg3d_MultiDom_Ellipse
    # Non-manifold test case
    mmg3d_NM_Cube
    mmg3d_NM_Complex
    )

  SET ( input_paths  ${input_paths}
    ### Cube
    ${MMG3D_CI_TESTS}/CubeIso_h_hminMax
    ${MMG3D_CI_TESTS}/CubeIso_0.5h_hminMax
    ${MMG3D_CI_TESTS}/CubeIso_0.25h_hminMax
    ###
    ${MMG3D_CI_TESTS}/CubeIso_h_met
    ${MMG3D_CI_TESTS}/CubeIso_0.5h_met
    ${MMG3D_CI_TESTS}/CubeIso_0.25h_met
    ${MMG3D_CI_TESTS}/CubeIso_0.125h_met
    ${MMG3D_CI_TESTS}/CubeAniIso_0.125h_met
    ### Sphere
    ${MMG3D_CI_TESTS}/SphereIso_h_met
    ${MMG3D_CI_TESTS}/SphereIso_0.5h_met
    ${MMG3D_CI_TESTS}/SphereIso_0.25h_met
    ${MMG3D_CI_TESTS}/SphereIso_0.125h_met
    ${MMG3D_CI_TESTS}/SphereIso_0.020_met
    # ${MMG3D_CI_TESTS}/SphereIso_0.020-0.015_met
    ${MMG3D_CI_TESTS}/SphereAni_0.02
    ###
    ${MMG3D_CI_TESTS}/SphereIso_2smallh_met
    ${MMG3D_CI_TESTS}/SphereIso_4smallh_met
    ${MMG3D_CI_TESTS}/SphereIso_8smallh_met
    ###
    ${MMG3D_CI_TESTS}/SphereIso_h_hausd0.001
    ${MMG3D_CI_TESTS}/SphereIso_h_hausd0.005
    ### CubeSkin
    ${MMG3D_CI_TESTS}/CubeSkin0.05_Inside0.4
    ${MMG3D_CI_TESTS}/CubeSkin0.1_Inside0.4
    ${MMG3D_CI_TESTS}/CubeSkin0.2_Inside0.4
    ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.125
    # ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.25
    # ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.5
    ### Linkrods
    ${MMG3D_CI_TESTS}/Various_unref_Linkrods_met0.2
    ${MMG3D_CI_TESTS}/Various_unref_Linkrods_met0.2_hausd0.01
    ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05
    ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05_hausd0.01
    ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05_hausd0.001
    ### Santa
    ${MMG3D_CI_TESTS}/Various_refunref_Santa_met0.05_hausd0.001_ar90
    ${MMG3D_CI_TESTS}/Various_refunref_Santa_met0.05_hausd0.0001_ar90
    ### MultiDomain
    ${MMG3D_CI_TESTS}/MultiDom_Cube
    ${MMG3D_CI_TESTS}/MultiDom_Ellipse
    ${MMG3D_CI_TESTS}/NM_Cube
    ${MMG3D_CI_TESTS}/NM_Complex
    )

  SET ( input_files  ${input_files}
    ### Cube
    CubeIso0.1
    CubeIso0.1
    CubeIso0.1
    ###
    CubeIso0.1
    CubeIso0.1
    CubeIso0.1
    CubeIso0.1
    CubeIso0.1
    ### Sphere
    SphereIso0.5
    SphereIso0.5
    SphereIso0.5
    SphereIso0.5
    SphereIso0.5
    #SphereIso0.020
    sphere
    ###
    SphereIso0.0625
    SphereIso0.0625
    SphereIso0.0625
    ###
    SphereIso0.5
    SphereIso0.5
    ### CubeSkin
    CubeSkin0.05
    CubeSkin0.1
    CubeSkin0.2
    CubeSkin0.125
    #CubeSkin0.25
    #CubeSkin0.5
    ### Linkrods
    linkrods
    linkrods
    linkrods
    linkrods
    linkrods
    ### Santa
    santa
    santa
    ### MultiDomain
    c
    c.d
    nm
    nm4
    )

  SET ( args ${args}
    ### Cube
    "-v 5 -hmax 0.1 -hmin 0.1"
    "-v 5 -hmax 0.05 -hmin 0.05"
    "-v 5 -hmax 0.025 -hmin 0.025"
    ###
    "-v 5"
    "-v 5"
    "-v 5"
    "-v 5"
    "-v 5"
    ### Sphere
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    # "-v 5 -hausd 0.1"
    "-v 5"
    ###
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    ###
    "-v 5 -hausd 0.001 -hgrad -1"
    "-v 5 -hausd 0.005 -hgrad -1"
    ### CubeSkin
    "-v 5"
    "-v 5"
    "-v 5"
    "-v 5"
    # "-v 5"
    # "-v 5"
    ### Linkrods
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.01"
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.01"
    "-v 5 -hausd 0.001"
    ### Santa
    "-v 5 -hausd 0.001 -ar 90"
    "-v 5 -hausd 0.0001 -ar 90"
    ### MultiDomain
    "-v 5 -hmax 0.02"
    "-v 5 -hausd 0.0003"
    "-v 5 -hmax 0.05"
    "-v 5"
    )

ENDIF ( )

ADD_RUN_AGAIN_TESTS ( ${EXECUT_MMG3D} "${test_names}" "${args}" "{input_paths}" "${input_files}" )

IF ( LONG_TESTS )
  ### M6
  SET ( test_name
    # 4: Refinment on a solution
    mmg3d_Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2
    )
  MMG_ADD_TEST(${test_name}
    "${EXECUT_MMG3D}
    -v 5 -sol ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/metM6.sol -hausd 0.1 -ar 60"
    ### M6
    "${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2" "M6"
    )

  SET_TESTS_PROPERTIES ( ${test_name}
    PROPERTIES FIXTURES_SETUP ${test_name} )

  IF ( RUN_AGAIN )
    MMG_ADD_TEST(${test_name}_2
      "${EXECUT_MMG3D}
      -v 5 -hausd 0.1 -ar 60 -hgrad -1"
      "${CTEST_OUTPUT_DIR}" "${test_name}-out.o.meshb"
      )

    SET_TESTS_PROPERTIES ( ${test_name}_2
      PROPERTIES FIXTURES_REQUIRED ${test_name} )

  ENDIF ( RUN_AGAIN )

  SET ( test_name
      # 4: Refinment on a solution
      mmg3d_Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3
      )

  MMG_ADD_TEST(${test_name}
    "${EXECUT_MMG3D}
    -v 5 -sol ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/metM6.sol -hausd 0.1 -ar 60"
    ### M6
    "${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3" "M6"
    )

  SET_TESTS_PROPERTIES ( ${test_name}
    PROPERTIES FIXTURES_SETUP ${test_name} )

    IF ( RUN_AGAIN )
      MMG_ADD_TEST(${test_name}_2
        "${EXECUT_MMG3D}
        -v 5 -hausd 0.1 -ar 60 -hgrad -1"
        "${CTEST_OUTPUT_DIR}" "${test_name}-out.o.meshb"
        )

    SET_TESTS_PROPERTIES ( ${test_name}_2
      PROPERTIES FIXTURES_REQUIRED ${test_name} )

  ENDIF ( )
ENDIF ( LONG_TESTS )


###############################################################################
#####
#####         Input/Output
#####
###############################################################################

# Binary gmsh
MMG_ADD_TEST(mmg3d_binary_gmsh_3d
  "${EXECUT_MMG3D} -v 5"
  "${MMG3D_CI_TESTS}/GmshInout" "cube.mshb"
  )

# Ascii gmsh
MMG_ADD_TEST(mmg3d_ascii_gmsh_3d
  "${EXECUT_MMG3D} -v 5"
  "${MMG3D_CI_TESTS}/GmshInout" "cube.msh"
  )

# Tetgen
# Default Tetgen behaviour saves only boundary tria (resp. edges) in
# .face (resp. .edge) file.
MMG_ADD_TEST ( mmg3d_cube-tetgen
  "${EXECUT_MMG3D} -v 5"
  "${MMG3D_CI_TESTS}/Cube" "cube"
  )

##############################################################################
#####
#####         Check Memory Leak
#####
##############################################################################
#####
ADD_TEST ( mmg3d_LeakCheck_AbnormalEnd3
  COMMAND ${EXECUT_MMG3D} -v 5
  ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd3/d -sol
  ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd3/dsol.sol -ls
  -out ${CTEST_OUTPUT_DIR}/mmg3d_LeakCheck_AbnormalEnd3-d.o.meshb)
SET(passRegex "## ERROR: UNABLE TO LOAD SOLUTION")
SET_PROPERTY(TEST mmg3d_LeakCheck_AbnormalEnd3
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
MMG_ADD_TEST(mmg3d_LeakCheck_optLevelSet
  "${EXECUT_MMG3D} -v 5  -ls -hgrad 1.5"
  "${MMG3D_CI_TESTS}/LeakCheck_optLevelSet" "rect03d"
  )

##############################################################################
#####
#####         Check Options
#####
##############################################################################
#####
MMG_ADD_TEST(mmg3d_memOption
  "${EXECUT_MMG3D} -v 5 -m 100"
  "${MMG3D_CI_TESTS}/Cube" "cube"
  )

MMG_ADD_TEST(mmg3d_hsizAndNosurfOption
  "${EXECUT_MMG3D} -v 5 -hsiz 0.1 -nosurf"
  "${MMG3D_CI_TESTS}/Cube" "cube"
  )

MMG_ADD_TEST(mmg3d_hsizAndNosurfAni
  "${EXECUT_MMG3D} -v 5 -hsiz 0.1 -nosurf -A"
  "${MMG3D_CI_TESTS}/Cube" "cube"
  )

MMG_ADD_TEST(mmg3d_val
  "${EXECUT_MMG3D} -v 5 -val"
  "${MMG3D_CI_TESTS}/Cube" "cube"
  )

#MMG_ADD_TEST(mmg3d_default
#  "${EXECUT_MMG3D} -v 5 -default"
#  ${MMG3D_CI_TESTS}/Cube/cube
#  )

SET_PROPERTY(TEST mmg3d_val #mmg3d_default
  PROPERTY WILL_FAIL TRUE)

# default hybrid
MMG_ADD_TEST(mmg3d_hybrid_3d
  "${EXECUT_MMG3D} -v 5"
  "${MMG3D_CI_TESTS}/Hybrid" "prism.mesh"
  )

# nsd + hybrid
MMG_ADD_TEST(mmg3d_hybrid-nsd1
  "${EXECUT_MMG3D} -v 5 -nsd 1"
  "${MMG3D_CI_TESTS}/Hybrid" "prism.mesh"
  )

###############################################################################
#####
#####         Check Boundaries
#####
###############################################################################
#####
MMG_ADD_TEST(mmg3d_ChkBdry_optls_test4
  "${EXECUT_MMG3D} -v 5  -ls -hgrad 1.5
  -sol ${MMG3D_CI_TESTS}/ChkBdry_optls_test4/test4.sol -in"
  "${MMG3D_CI_TESTS}/ChkBdry_optls_test4" "test4"
  )

#####
MMG_ADD_TEST(mmg3d_ChkBdry_optls_temp
  "${EXECUT_MMG3D} -v 5 -ls -hmin 5 -hmax 6
  -nr -hausd 0.5 -hgrad 1.2
  -sol ${MMG3D_CI_TESTS}/ChkBdry_optls_temp/temp.sol
  -in " "${MMG3D_CI_TESTS}/ChkBdry_optls_temp" "temp"
  )
####
MMG_ADD_TEST(mmg3d_ChkBdry_optls_temp2
  "${EXECUT_MMG3D} -v 5  -ls -hmin 5 -hmax 6
  -nr -hausd 0.5 -hgrad 1.2
  -sol ${MMG3D_CI_TESTS}/ChkBdry_optls_temp/temp.sol
  -in " "${MMG3D_CI_TESTS}/ChkBdry_optls_temp" "temp"
  )
#####
MMG_ADD_TEST(mmg3d_ChkBdry_cube
  "${EXECUT_MMG3D} -v 5"
  "${MMG3D_CI_TESTS}/ChkBdry_cube" "cube"
  )
#####
MMG_ADD_TEST(mmg3d_ChkBdry_multidomCube
  "${EXECUT_MMG3D} -v 5 -hmax 0.1"
  "${MMG3D_CI_TESTS}/ChkBdry_multidomCube" "c"
  )
#####
MMG_ADD_TEST(mmg3d_ChkBdry_multidomCube2
  "${EXECUT_MMG3D} -v 5 -hmax 0.1"
  "${MMG3D_CI_TESTS}/ChkBdry_multidomCube2" "c"
  )
#####
MMG_ADD_TEST(mmg3d_ChkBdry_multidomCube3
  "${EXECUT_MMG3D} -v 5 -hmax 0.1"
  "${MMG3D_CI_TESTS}/ChkBdry_multidomCube3" "c"
  )

MMG_ADD_TEST(mmg3d_opnbdy_unref_peninsula
  "${EXECUT_MMG3D} -v 5 -opnbdy -in"
  "${MMG3D_CI_TESTS}/OpnBdy_peninsula" "peninsula"
  )

MMG_ADD_TEST(mmg3d_opnbdy_ls_peninsula
  "${EXECUT_MMG3D} -v 5 -opnbdy -ls
  -sol  ${MMG3D_CI_TESTS}/OpnBdy_peninsula/ls.sol
  -in" "${MMG3D_CI_TESTS}/OpnBdy_peninsula" "peninsula"
  )

# ls + nsd
MMG_ADD_TEST(mmg3d_opnbdy_ls_peninsula-nsd3
  "${EXECUT_MMG3D} -v 5 -opnbdy -ls -nsd 3
  -sol  ${MMG3D_CI_TESTS}/OpnBdy_peninsula/ls.sol
  -in" "${MMG3D_CI_TESTS}/OpnBdy_peninsula" "peninsula"
  )

MMG_ADD_TEST(mmg3d_opnbdy_ref_peninsula
  "${EXECUT_MMG3D} -v 5 -hmax 0.06 -opnbdy
  -in" "${MMG3D_CI_TESTS}/OpnBdy_peninsula" "peninsula"
  )

MMG_ADD_TEST(mmg3d_opnbdy_unref_island
  "${EXECUT_MMG3D} -v 5 -opnbdy
  -in " "${MMG3D_CI_TESTS}/OpnBdy_island" "island"
  )

MMG_ADD_TEST(mmg3d_opnbdy_ref_island
  "${EXECUT_MMG3D} -v 5 -hmax 0.06 -opnbdy
  -in " "${MMG3D_CI_TESTS}/OpnBdy_island" "island"
  )

###############################################################################
#####
#####         Check Lagrangian motion option
#####
###############################################################################
#####
IF ( ELAS_FOUND AND NOT USE_ELAS MATCHES OFF )
  MMG_ADD_TEST(mmg3d_LagMotion0_tinyBoxt
    "${EXECUT_MMG3D} -v 5  -lag 0
    -sol ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt.sol
    -in " "${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt" "tinyBoxt"
    )
  MMG_ADD_TEST(mmg3d_LagMotion1_tinyBoxt
    "${EXECUT_MMG3D} -v 5  -lag 1
    -sol ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt.sol
    -in " "${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt" "tinyBoxt"
    )
  MMG_ADD_TEST(mmg3d_LagMotion2_tinyBoxt
    "${EXECUT_MMG3D} -v 5  -lag 2
    -sol ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt.sol
    -in " "${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt" "tinyBoxt"
    )
  # nsd
  MMG_ADD_TEST(mmg3d_LagMotion2_tinyBoxt-nsd3
    "${EXECUT_MMG3D} -v 5  -lag 2 -nsd 3
    -sol ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt.sol
    -in " "${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt" "tinyBoxt"
    )

ENDIF()

##############################################################################
#####
#####         Check Local parameters at tetra
#####
##############################################################################
#####
MMG_ADD_TEST(mmg3d_TetLoc_Ellipse
  "${EXECUT_MMG3D} -v 5 -hgrad -1"
  "${MMG3D_CI_TESTS}/TetLoc_Ellipse" "c"
  )

##############################################################################
#####
#####         Check optim + aniso option
#####
##############################################################################
#####
MMG_ADD_TEST(mmg3d_OptimAni_Sphere
  "${EXECUT_MMG3D} -v 5 -optim -A -sol 2"
  "${MMG3D_CI_TESTS}/SphereIso_h_met" "SphereIso0.5.meshb"
  )

##############################################################################
#####
#####         Check optimLES
#####
##############################################################################
#####
MMG_ADD_TEST(mmg3d_OptimLES_sphere
  "${EXECUT_MMG3D} -v 5 -optimLES"
  "${MMG3D_CI_TESTS}/SphereIso_0.25h_met" "SphereIso0.5"
  )

###############################################################################
#####
#####         Check Results
#####
###############################################################################
#####
MMG_ADD_TEST(mmg3d_LSMultiMat
  "${EXECUT_MMG3D} -v 5 -ls -nr
  -sol ${MMG3D_CI_TESTS}/LSMultiMat/step.0.phi.sol"
  "${MMG3D_CI_TESTS}/LSMultiMat" "step.0.mesh"
  )

#multi-mat + opnbdy + non-manifold check
MMG_ADD_TEST(mmg3d_LSMultiMat_nm
  "${EXECUT_MMG3D} -v 5 -ls -0.1 -hausd 0.05 -hgrad 1.8 -nr -opnbdy"
  "${MMG3D_CI_TESTS}/LSMultiMat" "3d-opn"
  )

MMG_ADD_TEST(mmg3d_OptLs_plane_val
  "${EXECUT_MMG3D} -v 5 -ls -val
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
  -met ${MMG3D_CI_TESTS}/OptLs_plane/met.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

#MMG_ADD_TEST(mmg3d_OptLs_plane_default
#  COMMAND ${EXECUT_MMG3D} -v 5 -ls -default
#  ${MMG3D_CI_TESTS}/OptLs_plane/plane
#  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
#  -met ${MMG3D_CI_TESTS}/OptLs_plane/met.sol
#  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-nonzero.o.meshb)

SET_PROPERTY(TEST  mmg3d_OptLs_plane_val #mmg3d_OptLs_plane_default
  PROPERTY WILL_FAIL TRUE)

# ls oritentation
MMG_ADD_TEST(mmg3d_OptLs_plane_p
  "${EXECUT_MMG3D} -v 5 -ls
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/p.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

MMG_ADD_TEST(mmg3d_OptLs_plane_m
  "${EXECUT_MMG3D} -v 5 -ls
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

# non-zero ls
MMG_ADD_TEST(mmg3d_OptLs_plane_nonzero
  "${EXECUT_MMG3D} -v 5 -ls 0.1
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

# ls discretization + optim
MMG_ADD_TEST(mmg3d_OptLs_plane_optim
  "${EXECUT_MMG3D} -v 5 -ls -optim
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

# ls discretization + optim + aniso
MMG_ADD_TEST(mmg3d_OptLs_plane_optimAni
  "${EXECUT_MMG3D} -v 5 -ls -optim -A
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

# ls discretization + hsiz
MMG_ADD_TEST(mmg3d_OptLs_plane_hsiz
  "${EXECUT_MMG3D} -v 5 -ls -hsiz 0.2
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

# ls discretization + hsiz
MMG_ADD_TEST(mmg3d_OptLs_plane_hsizAni
  "${EXECUT_MMG3D} -v 5 -ls -hsiz 0.2 -A
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

# ls discretization + metric
MMG_ADD_TEST(mmg3d_OptLs_plane_withMetAndLs
  "${EXECUT_MMG3D} -v 5 -ls
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
  -met ${MMG3D_CI_TESTS}/OptLs_plane/met.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

# ls + rmc + LSBaseReference
MMG_ADD_TEST(mmg3d_OptLs_LSBaseReferences-rmc
  "${EXECUT_MMG3D} -v 5 -ls -rmc -nr
  -sol ${MMG3D_CI_TESTS}/LSBaseReferences/box.sol"
  "${MMG3D_CI_TESTS}/LSBaseReferences" "box"
  )

MMG_ADD_TEST(mmg3d_OptLs_LSBaseReferences-normc
  "${EXECUT_MMG3D} -v 5 -ls -nr
  -sol ${MMG3D_CI_TESTS}/LSBaseReferences/box.sol"
  "${MMG3D_CI_TESTS}/LSBaseReferences" "box"
  )

# ls + rmc
MMG_ADD_TEST(mmg3d_OptLs_plane_withbub
  "${EXECUT_MMG3D} -v 5 -ls
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/bub.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

# ls + rmc: max pile bug
MMG_ADD_TEST(mmg3d_OptLs_plane_rmcmaxpile
  "${EXECUT_MMG3D} -v 5 -ls -rmc
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/whole.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

MMG_ADD_TEST(mmg3d_OptLs_plane_rembub
  "${EXECUT_MMG3D} -v 5 -ls -rmc
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/bub.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

MMG_ADD_TEST(mmg3d_OptLs_plane_rembub2
  "${EXECUT_MMG3D} -v 5 -ls -rmc 0.1
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/bub.sol"
  "${MMG3D_CI_TESTS}/OptLs_plane" "plane"
  )

# Preservation of orphan points
MMG_ADD_TEST(mmg3d_OptLs_temp_orphan
  "${EXECUT_MMG3D} -v 5 -ls
  -sol ${MMG3D_CI_TESTS}/OptLs_temp_hminMax_hgrad1.2_hausd0.1/temp.sol
  -hausd 0.5 -nr -hgrad -1 -nsd 3"
  "${MMG3D_CI_TESTS}/OptLs_temp_hminMax_hgrad1.2_hausd0.1" "temp"
  )

# OptLs and isoref option: compare the result of ls discretization with ref 10
# and results of the same case with ref 5
#include(FindUnixCommands)

MMG_ADD_TEST(
  mmg3d_OptLs_isoref_defaut
  "${EXECUT_MMG3D} -v 5 -ls
  -sol ${MMG3D_CI_TESTS}/OptLs_isoref/3d-mesh.sol"
  "${MMG3D_CI_TESTS}/OptLs_isoref" "3d-mesh.mesh"
  )
MMG_ADD_TEST(
  mmg3d_OptLs_isoref_5
  "${EXECUT_MMG3D} -v 5 -isoref 5 -ls
  -sol ${MMG3D_CI_TESTS}/OptLs_isoref/3d-mesh.sol"
  "${MMG3D_CI_TESTS}/OptLs_isoref" "3d-mesh-isoref5.mesh"
  )

if (BASH)
  add_test(
    NAME mmg3d_optLs_isoref
    COMMAND ${BASH} -c "diff <(wc -wl ${CTEST_OUTPUT_DIR}/mmg3d_isoref.o.mesh  | awk '{print $1 $2}') <(wc -wl ${CTEST_OUTPUT_DIR}/mmg3d_isoref5.o.mesh | awk '{print $1 $2}')"
    )
endif()

IF ( TEST_LIBMMG3D )
  MMG_ADD_TEST(test_para_tria
    "${EXECUT_MMG3D}
    -ar 0.02 -nofem -nosizreq -hgradreq -1 -hgrad -1
    -sol ${MMG3D_CI_TESTS}/test_para_tria/proc0.sol"
    "${MMG3D_CI_TESTS}/test_para_tria" "proc0.mesh"
    )

  SET_TESTS_PROPERTIES ( test_para_tria
    PROPERTIES FIXTURES_SETUP test_para_tria )

  ADD_TEST(test_compare_para_tria
    COMMAND ${TEST_COMPARE_PARA_TRIA}
    ${MMG3D_CI_TESTS}/test_para_tria/proc0.mesh
    ${CTEST_OUTPUT_DIR}/test_para_tria.o.mesh
    )
  SET_TESTS_PROPERTIES ( test_compare_para_tria
    PROPERTIES FIXTURES_REQUIRED test_para_tria )
ENDIF()

IF ( LONG_TESTS )
  # Test the Ls option
  MMG_ADD_TEST(mmg3d_OptLs_cube303d_hminMax_hgrad1.2_hausd0.005
    "${EXECUT_MMG3D} -v 5 -ls
    -sol ${MMG3D_CI_TESTS}/OptLs_cube303d_hminMax_hgrad1.2_hausd0.005/cube303d.sol
    -hausd 0.005 -nr -hgrad 1.2 -hmin 0.001 -hmax 0.1"
    "${MMG3D_CI_TESTS}/OptLs_cube303d_hminMax_hgrad1.2_hausd0.005" "cube303d"
    )
  MMG_ADD_TEST(mmg3d_OptLs_temp_hminMax_hgrad1.2_hausd0.1
    "${EXECUT_MMG3D} -v 5 -ls
    -sol ${MMG3D_CI_TESTS}/OptLs_temp_hminMax_hgrad1.2_hausd0.1/temp.sol
    -hausd 0.1 -nr -hgrad 1.2 -hmin 3 -hmax 4"
    "${MMG3D_CI_TESTS}/OptLs_temp_hminMax_hgrad1.2_hausd0.1" "temp"
    )

  ###############################################################################
  #####
  #####         Check Lagrangian motion option
  #####
  ###############################################################################
  #####
  IF ( ELAS_FOUND AND NOT USE_ELAS MATCHES OFF )
    MMG_ADD_TEST(mmg3d_LagMotion0_boxt
      "${EXECUT_MMG3D} -v 5  -lag 0
      -sol ${MMG3D_CI_TESTS}/LagMotion1_boxt/boxt.sol
      -in" "${MMG3D_CI_TESTS}/LagMotion1_boxt" "boxt"
      )
    MMG_ADD_TEST(mmg3d_LagMotion1_boxt
      "${EXECUT_MMG3D} -v 5  -lag 1
      -sol ${MMG3D_CI_TESTS}/LagMotion1_boxt/boxt.sol
      -in" "${MMG3D_CI_TESTS}/LagMotion1_boxt" "boxt"
      )
    MMG_ADD_TEST(mmg3d_LagMotion2_boxt
      "${EXECUT_MMG3D} -v 5  -lag 2
      -sol ${MMG3D_CI_TESTS}/LagMotion1_boxt/boxt.sol
      -in" "${MMG3D_CI_TESTS}/LagMotion1_boxt" "boxt"
      )
  ENDIF()

ENDIF()

###############################################################################
#####
#####         Bug Fix
#####
###############################################################################
#####
#MMG_ADD_TEST(mmg3d_BUG_OptLsSingularities
# "${EXECUT_MMG3D} -v 5  -ls"
# "${MMG3D_CI_TESTS}/BUG_OptLsSingularities" "test4"
# )
#
#MMG_ADD_TEST(mmg3d_TestDoSol_1
# "${EXECUT_MMG3D} -v 5  -hgrad -1 -hausd 1 -m 100"
# "${MMG3D_CI_TESTS}/TestDoSol_1" "66_shaver3.mesh"
# )
