###############################################################################
#####
#####         Continuous Integration
#####
###############################################################################
SET(CTEST_TIMEOUT           "7200")
SET(REG_TESTS_PATH ${CMAKE_SOURCE_DIR}/../RegTests)

# simple test: must already pass
ADD_TEST(NAME SimpleCube
  COMMAND ${EXECUT} -v 6 -d
  ${REG_TESTS_PATH}/Cube/cube
  -out ${REG_TESTS_PATH}/Cube/cube.o.meshb)

###############################################################################
#####
#####         Check Build
#####
###############################################################################
# on agni we can test gcc and icc
#  ADD_CUSTOM_COMMAND(
#    COMMAND ${CMAKE_COMMAND} -E  make clean
#    COMMENT "make clean")


###############################################################################
#####
#####         Check Memory Leak
#####
###############################################################################
#####
#ADD_TEST(NAME LeakCheck_AbnormalEnd0
#  COMMAND ${EXEC} -v 5
#  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd0/d
#  -DCOMMAND1="chmod -r" -P "chmod.cmake")
#SET(passRegex "LeakCheck_AbnormalEnd0/d  NOT FOUND.")
#SET_PROPERTY(TEST LeakCheck_AbnormalEnd0
#  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#
#####
#ADD_TEST(NAME LeakCheck_AbnormalEnd1
#  COMMAND ${EXEC} -v 5
#  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd1) ## no possible for now
#SET(passRegex " ")
#SET_PROPERTY(TEST LeakCheck_AbnormalEnd1
#  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
#SET ( LISTEXEC ${EXECUT} )
#IF ( FIND_LIBEXEC2 )
#  SET( LISTEXEC ${LISTEXEC} ${LIBEXEC2} )
#ENDIF ()

FOREACH(EXEC ${LISTEXEC})
  ADD_TEST(NAME LeakCheck_AbnormalEnd2_${EXEC}
    COMMAND ${EXEC} -v 5
    ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd2/d
    -out ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd2/d.o.meshb)
  SET(passRegex "## Unable to scale mesh.")
  SET_PROPERTY(TEST LeakCheck_AbnormalEnd2_${EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  ADD_TEST(NAME LeakCheck_AbnormalEnd3_${EXEC}
    COMMAND ${EXEC} -v 5
    ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd3/d -sol
    ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd3/dsol.sol -ls 2
    -out ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd3/d.o.meshb)
  SET(passRegex "## ERROR: A VALID SOLUTION FILE IS NEEDED")
  SET_PROPERTY(TEST LeakCheck_AbnormalEnd3_${EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  #ADD_TEST(NAME LeakCheck_AbnormalEnd4_${EXEC}
  #  COMMAND ${EXEC} -v 5
  #  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd4) ## no possible for now
  #SET(passRegex " ")
  #SET_PROPERTY(TEST LeakCheck_AbnormalEnd4_${EXEC}
  #  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  #ADD_TEST(NAME LeakCheck_AbnormalEnd5.1_${EXEC}
  #  COMMAND ${EXEC} -v 5
  #  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd5.1) ## todo mmg3d2
  #SET(passRegex " ")
  #SET_PROPERTY(TEST LeakCheck_AbnormalEnd5.1_${EXEC}
  #  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  #ADD_TEST(NAME LeakCheck_AbnormalEnd5.2_${EXEC}
  #  COMMAND ${EXEC} -v 5
  #  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd5.2) ## todo mmg3d1
  #SET(passRegex " ")
  #SET_PROPERTY(TEST LeakCheck_AbnormalEnd5.2_${EXEC}
  #  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  #ADD_TEST(NAME LeakCheck_AbnormalEnd6_${EXEC}
  #  COMMAND ${EXEC} -v 5
  #  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd6) ## todo analys.c
  #SET(passRegex " ")
  #SET_PROPERTY(TEST LeakCheck_AbnormalEnd6_${EXEC}
  #  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  ADD_TEST(NAME LeakCheck_AbnormalEnd7_${EXEC}
    COMMAND ${EXEC} -v 5
    ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd7/d
    -out ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd7/d.o)
  SET(passRegex "\\*\\* UNABLE TO OPEN.*d.o")
  SET_PROPERTY(TEST LeakCheck_AbnormalEnd7_${EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  ADD_TEST(NAME LeakCheck_AbnormalEnd8_${EXEC}
    COMMAND ${EXEC} -v 5
    ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd8/d
    -out ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd8/d.o.meshb)
  SET(passRegex "\\*\\* UNABLE TO OPEN.*d.o.sol")
  SET_PROPERTY(TEST LeakCheck_AbnormalEnd8_${EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  #####
  ADD_TEST(NAME LeakCheck_args0_${EXEC}
    COMMAND ${EXEC} -v 5
    ${REG_TESTS_PATH}/LeakCheck_args0/d
    -out ${REG_TESTS_PATH}/LeakCheck_args0/d.o.meshb)
  #####
  ADD_TEST(NAME LeakCheck_args1_${EXEC}
    COMMAND ${EXEC} -v 5
    -in ${REG_TESTS_PATH}/LeakCheck_args1/d -sol
    ${REG_TESTS_PATH}/LeakCheck_args1/dsol.sol
    -out ${REG_TESTS_PATH}/LeakCheck_args1/dout.meshb)
  #####
  ADD_TEST(NAME LeakCheck_optLevelSet_${EXEC}
    COMMAND ${EXEC} -v 5 -ls 0
    ${REG_TESTS_PATH}/LeakCheck_optLevelSet/rect03d
    -out ${REG_TESTS_PATH}/LeakCheck_optLevelSet/rect03d.o.meshb)

  ###############################################################################
  #####
  #####         Check Precision
  #####
  ###############################################################################
  #####
  ADD_TEST(NAME MeshVersionFormatted1_${EXEC}
    COMMAND ${EXEC} -v 5
    -in ${REG_TESTS_PATH}/MeshVersionFormatted1/d
    -sol ${REG_TESTS_PATH}/MeshVersionFormatted1/dsol.sol
    -out ${REG_TESTS_PATH}/MeshVersionFormatted1/d.o.meshb)
  #####
  ADD_TEST(NAME MeshVersionFormatted2_${EXEC}
    COMMAND ${EXEC} -v 5
    -in ${REG_TESTS_PATH}/MeshVersionFormatted2/d
    -sol ${REG_TESTS_PATH}/MeshVersionFormatted2/dsol.sol
    -out ${REG_TESTS_PATH}/MeshVersionFormatted2/d.o.meshb)

  ###############################################################################
  #####
  #####         Check Boundaries
  #####
  ###############################################################################
  #####
  ADD_TEST(NAME ChkBdry_optls_temp_${EXEC}
    COMMAND ${EXEC} -v 5 -ls -hmin 5 -hmax 6
    -nr -hausd 0.5 -hgrad 1.2
    -in ${REG_TESTS_PATH}/ChkBdry_optls_temp/temp
    -sol ${REG_TESTS_PATH}/ChkBdry_optls_temp/temp.sol
    -out ${REG_TESTS_PATH}/ChkBdry_optls_temp/temp.o.meshb)
  ####
  ADD_TEST(NAME ChkBdry_optls_temp2_${EXEC}
    COMMAND ${EXEC} -v 5 -ls -hmin 5 -hmax 6
    -nr -hausd 0.5 -hgrad 1.2
    -in ${REG_TESTS_PATH}/ChkBdry_optls_temp/temp
    -sol ${REG_TESTS_PATH}/ChkBdry_optls_temp/temp.sol
    -out ${REG_TESTS_PATH}/ChkBdry_optls_temp/temp.o.meshb)
  #####
  ADD_TEST(NAME ChkBdry_cube_${EXEC}
    COMMAND ${EXEC}
    ${REG_TESTS_PATH}/ChkBdry_cube/cube)
  #####
  ADD_TEST(NAME ChkBdry_multidomCube_${EXEC}
    COMMAND ${EXEC} -hmax 0.1
    ${REG_TESTS_PATH}/ChkBdry_multidomCube/c)
  #####
  ADD_TEST(NAME ChkBdry_multidomCube2_${EXEC}
    COMMAND ${EXEC} -hmax 0.1
    ${REG_TESTS_PATH}/ChkBdry_multidomCube2/c)
 #####
  ADD_TEST(NAME ChkBdry_multidomCube3_${EXEC}
    COMMAND ${EXEC} -hmax 0.1
    ${REG_TESTS_PATH}/ChkBdry_multidomCube3/c)
ENDFOREACH(EXEC)



###############################################################################
#####
#####         Check Results
#####
###############################################################################
#####

# Check what happend when we refine an isotropic cube of size h with a constant
# metric (h, h/2, h/4, h/8 and h/16)
#---First with hmin=hmax
ADD_TEST(NAME CubeIso_h_hminMax
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeIso_h_hminMax/CubeIso0.1 -hmax 0.1 -hmin 0.1
  -out ${REG_TESTS_PATH}/CubeIso_h_hminMax/CubeIso0.1.o.meshb)
ADD_TEST(NAME CubeIso_0.5h_hminMax
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeIso_0.5h_hminMax/CubeIso0.1 -hmax 0.05 -hmin 0.05
  -out ${REG_TESTS_PATH}/CubeIso_0.5h_hminMax/CubeIso0.1.o.meshb)
ADD_TEST(NAME CubeIso_0.25h_hminMax
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeIso_0.25h_hminMax/CubeIso0.1 -hmax 0.025 -hmin 0.025
  -out ${REG_TESTS_PATH}/CubeIso_0.25h_hminMax/CubeIso0.1.o.meshb)
#ADD_TEST(NAME CubeIso_0.125h_hminMax
#  COMMAND ${EXECUT} -v 5
#  ${REG_TESTS_PATH}/CubeIso_0.125h_hminMax/CubeIso0.1 -hmax 0.0125 -hmin 0.0125)

#---Second with sol file
ADD_TEST(NAME CubeIso_h_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeIso_h_met/CubeIso0.1
  -out ${REG_TESTS_PATH}/CubeIso_h_met/CubeIso0.1.o.meshb)
ADD_TEST(NAME CubeIso_0.5h_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeIso_0.5h_met/CubeIso0.1
  -out ${REG_TESTS_PATH}/CubeIso_0.5h_met/CubeIso0.1.o.meshb)
ADD_TEST(NAME CubeIso_0.25h_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeIso_0.25h_met/CubeIso0.1
  -out ${REG_TESTS_PATH}/CubeIso_0.25h_met/CubeIso0.1.o.meshb)
ADD_TEST(NAME CubeIso_0.125h_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeIso_0.125h_met/CubeIso0.1 -m 600
  -out ${REG_TESTS_PATH}/CubeIso_0.125h_met/CubeIso0.1.o.meshb)

#####

# Check what happend when we refine a sphere of size h with a constant metric
# (h, h/2, h/4 and h/8)
##---First with hmin=hmax
#ADD_TEST(NAME SphereIso_h_hminMax
#  COMMAND ${EXECUT} -v 5
#  ${REG_TESTS_PATH}/SphereIso_h_hminMax/SphereIso0.5
#  -hmax 0.5 -hmin 0.5 -hausd 0.1)
#ADD_TEST(NAME SphereIso_0.5h_hminMax
#  COMMAND ${EXECUT} -v 5
#  ${REG_TESTS_PATH}/SphereIso_0.5h_hminMax/SphereIso0.5
#  -hmax 0.25 -hmin 0.25 -hausd 0.1)
#ADD_TEST(NAME SphereIso_0.25h_hminMax
#  COMMAND ${EXECUT} -v 5
#  ${REG_TESTS_PATH}/SphereIso_0.25h_hminMax/SphereIso0.5
#  -hmax 0.125 -hmin 0.125 -hausd 0.1)
#ADD_TEST(NAME SphereIso_0.125h_hminMax
#  COMMAND ${EXECUT} -v 5
#  ${REG_TESTS_PATH}/SphereIso_0.125h_hminMax/SphereIso0.5
#  -hmax 0.0625 -hmin 0.0625 -hausd 0.1)
#---Second with sol file
ADD_TEST(NAME SphereIso_h_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/SphereIso_h_met/SphereIso0.5 -hausd 0.1
  -out ${REG_TESTS_PATH}/SphereIso_h_met/SphereIso0.5.o.meshb)
ADD_TEST(NAME SphereIso_0.5h_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/SphereIso_0.5h_met/SphereIso0.5 -hausd 0.1
  -out ${REG_TESTS_PATH}/SphereIso_0.5h_met/SphereIso0.5.o.meshb)
ADD_TEST(NAME SphereIso_0.25h_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/SphereIso_0.25h_met/SphereIso0.5 -hausd 0.1
  -out ${REG_TESTS_PATH}/SphereIso_0.25h_met/SphereIso0.5.o.meshb)
ADD_TEST(NAME SphereIso_0.125h_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/SphereIso_0.125h_met/SphereIso0.5 -hausd 0.1
  -out ${REG_TESTS_PATH}/SphereIso_0.125h_met/SphereIso0.5.o.meshb)

# Check what happend when we unrefine a sphere of size smallh with a constant metric
# (2*smallh, 4*smallh and 8*smallh)
ADD_TEST(NAME SphereIso_2smallh_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/SphereIso_2smallh_met/SphereIso0.0625 -hausd 0.1
  -out ${REG_TESTS_PATH}/SphereIso_2smallh_met/SphereIso0.0625.o.meshb)
ADD_TEST(NAME SphereIso_4smallh_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/SphereIso_4smallh_met/SphereIso0.0625 -hausd 0.1
  -out ${REG_TESTS_PATH}/SphereIso_4smallh_met/SphereIso0.0625.o.meshb)
ADD_TEST(NAME SphereIso_8smallh_met
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/SphereIso_8smallh_met/SphereIso0.0625 -hausd 0.1
  -out ${REG_TESTS_PATH}/SphereIso_8smallh_met/SphereIso0.0625.o.meshb)

# Check what happend when we use hausdorff number to refine the skin and a big hgrad
# to have an inside of the initial size (0.5)
ADD_TEST(NAME SphereIso_h_hausd0.001
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/SphereIso_h_hausd0.001/SphereIso0.5 -hausd 0.001 -hgrad 500
  -out ${REG_TESTS_PATH}/SphereIso_h_hausd0.001/SphereIso0.5.o.meshb)
ADD_TEST(NAME SphereIso_h_hausd0.005
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/SphereIso_h_hausd0.005/SphereIso0.5 -hausd 0.005 -hgrad 100
  -out ${REG_TESTS_PATH}/SphereIso_h_hausd0.005/SphereIso0.5.o.meshb)

# Check what happend when we refine a cube whose skin has already the good size
ADD_TEST(NAME CubeSkin0.05_Inside0.4
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeSkin0.05_Inside0.4/CubeSkin0.05
  ${REG_TESTS_PATH}/CubeSkin0.05_Inside0.4/CubeSkin0.05.o.meshb)
ADD_TEST(NAME CubeSkin0.1_Inside0.4
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeSkin0.1_Inside0.4/CubeSkin0.1
  ${REG_TESTS_PATH}/CubeSkin0.1_Inside0.4/CubeSkin0.1.o.meshb)
ADD_TEST(NAME CubeSkin0.2_Inside0.4
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeSkin0.2_Inside0.4/CubeSkin0.2
  ${REG_TESTS_PATH}/CubeSkin0.2_Inside0.4/CubeSkin0.2.o.meshb)
ADD_TEST(NAME CubeSkin0.0125_Inside0.125
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeSkin0.0125_Inside0.125/CubeSkin0.125 -m 600
  -out ${REG_TESTS_PATH}/CubeSkin0.0125_Inside0.125/CubeSkin0.125.o.meshb)
ADD_TEST(NAME CubeSkin0.0125_Inside0.25
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeSkin0.0125_Inside0.25/CubeSkin0.25 -m 600
  -out ${REG_TESTS_PATH}/CubeSkin0.0125_Inside0.25/CubeSkin0.25.o.meshb)
ADD_TEST(NAME CubeSkin0.0125_Inside0.5
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/CubeSkin0.0125_Inside0.5/CubeSkin0.5 -m 600
  ${REG_TESTS_PATH}/CubeSkin0.0125_Inside0.5/CubeSkin0.5.o.meshb)


# Check results on various meshes
# First: Meshes that we want unrefined
ADD_TEST(NAME Various_unref_Linkrods_met0.2
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/Various_unref_Linkrods_met0.2/linkrods -hausd 0.1
  ${REG_TESTS_PATH}/Various_unref_Linkrods_met0.2/linkrods.o.meshb)
ADD_TEST(NAME Various_unref_Linkrods_met0.2_hausd0.01
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/Various_unref_Linkrods_met0.2_hausd0.01/linkrods
  -hausd 0.01
  ${REG_TESTS_PATH}/Various_unref_Linkrods_met0.2_hausd0.01/linkrods.o.meshb)



# Second: Meshes that we want refined
ADD_TEST(NAME Various_ref_Linkrods_met0.05
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/Various_ref_Linkrods_met0.05/linkrods -hausd 0.1
  ${REG_TESTS_PATH}/Various_ref_Linkrods_met0.05/linkrods.o.meshb)
ADD_TEST(NAME Various_ref_Linkrods_met0.05_hausd0.01
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/Various_ref_Linkrods_met0.05_hausd0.01/linkrods
  -hausd 0.01
  ${REG_TESTS_PATH}/Various_ref_Linkrods_met0.05_hausd0.01/linkrods.o.meshb)
ADD_TEST(NAME Various_ref_Linkrods_met0.05_hausd0.001
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/Various_ref_Linkrods_met0.05_hausd0.001/linkrods
  -hausd 0.001
  ${REG_TESTS_PATH}/Various_ref_Linkrods_met0.05_hausd0.001/linkrods.o.meshb)

# Third: We refine some parts and unrefined others
ADD_TEST(NAME Various_refunref_Santa_met0.05_hausd0.001_ar90
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/Various_refunref_Santa_met0.05_hausd0.001_ar90/santa
  -hausd 0.001 -ar 90
  ${REG_TESTS_PATH}/Various_refunref_Santa_met0.05_hausd0.001_ar90/santa.o.meshb)
ADD_TEST(NAME Various_refunref_Santa_met0.05_hausd0.0001_ar90
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/Various_refunref_Santa_met0.05_hausd0.0001_ar90/santa
  -hausd 0.0001 -ar 90
  ${REG_TESTS_PATH}/Various_refunref_Santa_met0.05_hausd0.0001_ar90/santa.o.meshb)

# 4: Refinment on a solution
ADD_TEST(NAME Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/M6
  -sol ${REG_TESTS_PATH}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/metM6.sol
  -hausd 0.1 -ar 60 -hgrad 1 -m 700
  ${REG_TESTS_PATH}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/M6.o.meshb)
ADD_TEST(NAME Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3
  COMMAND ${EXECUT} -v 5
  ${REG_TESTS_PATH}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/M6
  -sol
  ${REG_TESTS_PATH}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/metM6.sol
  -hausd 0.1 -ar 60 -hgrad 1 -m 1000
  ${REG_TESTS_PATH}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/M6.o.meshb)


# Test the Ls option
ADD_TEST(NAME OptLs_cube303d_hminMax_hgrad1.2_hausd0.005
  COMMAND ${EXECUT} -ls
  ${REG_TESTS_PATH}/OptLs_cube303d_hminMax_hgrad1.2_hausd0.005/cube303d
  -sol ${REG_TESTS_PATH}/OptLs_cube303d_hminMax_hgrad1.2_hausd0.005/cube303d.sol
  -hausd 0.005 -nr -hgrad 1.2 -hmin 0.001 -hmax 0.1
  ${REG_TESTS_PATH}/OptLs_cube303d_hminMax_hgrad1.2_hausd0.005/cube303d.o.meshb)
ADD_TEST(NAME OptLs_temp_hminMax_hgrad1.2_hausd0.1
  COMMAND ${EXECUT} -v 6 -d -ls
  ${REG_TESTS_PATH}/OptLs_temp_hminMax_hgrad1.2_hausd0.1/temp
  -sol ${REG_TESTS_PATH}/OptLs_temp_hminMax_hgrad1.2_hausd0.1/temp.sol
  -hausd 0.1 -nr -hgrad 1.2 -hmin 3 -hmax 4
  ${REG_TESTS_PATH}/OptLs_temp_hminMax_hgrad1.2_hausd0.1/temp.o.meshb)


# Test multi-domain remeshing
ADD_TEST(NAME MultiDom_Cube
  COMMAND ${EXECUT} -v 6 -hmax 0.02 ${REG_TESTS_PATH}/MultiDom_Cube/c
  -out ${REG_TESTS_PATH}/MultiDom_Cube/c.o.meshb)

ADD_TEST(NAME MultiDom_Ellipse
  COMMAND ${EXECUT} -v 6 -m 700 -hausd 0.0003 ${REG_TESTS_PATH}/MultiDom_Ellipse/c.d
  -out ${REG_TESTS_PATH}/MultiDom_Ellipse/c.d.o.meshb)

# Non-manifold test case
ADD_TEST(NAME NM_Cube
  COMMAND ${EXECUT} -v 6 -d -hmax 0.05 ${REG_TESTS_PATH}/NM_Cube/nm
  -out ${REG_TESTS_PATH}/NM_Cube/nm.o.meshb)


# Compare with a reference result when we run
#ADD_TEST(NAME RefCube
#  COMMAND ${EXECUT} -v 5
#  ${REG_TESTS_PATH}/RefCube/cube) marre... a finir




###############################################################################
#####
#####         Bug Fix
#####
###############################################################################
#####
ADD_TEST(NAME BUG_OptLsSingularities
 COMMAND ${EXECUT} -v 5 -ls
 ${REG_TESTS_PATH}/BUG_OptLsSingularities/test4
 ${REG_TESTS_PATH}/BUG_OptLsSingularities/test4.o.meshb)
