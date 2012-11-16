###############################################################################
#####
#####         Continuous Integration
#####
###############################################################################

SET(REG_TESTS_PATH ${CMAKE_SOURCE_DIR}/../../RegTests)

# simple test: must already pass
ADD_TEST(NAME SimpleCube
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/Cube/cube.mesh)

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
#  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
#  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd0/d.mesh
#  -DCOMMAND1="chmod -r" -P "chmod.cmake")
#SET(passRegex "LeakCheck_AbnormalEnd0/d.mesh  NOT FOUND.")
#SET_PROPERTY(TEST LeakCheck_AbnormalEnd0
#  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#
#####
#ADD_TEST(NAME LeakCheck_AbnormalEnd1
#  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
#  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd1) ## no possible for now
#SET(passRegex " ")
#SET_PROPERTY(TEST LeakCheck_AbnormalEnd1
#  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
ADD_TEST(NAME LeakCheck_AbnormalEnd2
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd2/d.mesh)
SET(passRegex "## Unable to scale mesh.")
SET_PROPERTY(TEST LeakCheck_AbnormalEnd2
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
ADD_TEST(NAME LeakCheck_AbnormalEnd3
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd3/d.mesh -sol
  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd3/dsol.sol -ls 2)
SET(passRegex "## ERROR : A VALID SOLUTION FILE IS NEEDED")
SET_PROPERTY(TEST LeakCheck_AbnormalEnd3
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
#ADD_TEST(NAME LeakCheck_AbnormalEnd4
#  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
#  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd4) ## no possible for now
#SET(passRegex " ")
#SET_PROPERTY(TEST LeakCheck_AbnormalEnd4
#  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
#ADD_TEST(NAME LeakCheck_AbnormalEnd5.1
#  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
#  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd5.1) ## todo mmg3d2
#SET(passRegex " ")
#SET_PROPERTY(TEST LeakCheck_AbnormalEnd5.1
#  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
#ADD_TEST(NAME LeakCheck_AbnormalEnd5.2
#  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
#  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd5.2) ## todo mmg3d1
#SET(passRegex " ")
#SET_PROPERTY(TEST LeakCheck_AbnormalEnd5.2
#  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
#ADD_TEST(NAME LeakCheck_AbnormalEnd6
#  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
#  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd6) ## todo analys.c
#SET(passRegex " ")
#SET_PROPERTY(TEST LeakCheck_AbnormalEnd6
#  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
ADD_TEST(NAME LeakCheck_AbnormalEnd7
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd7/d.mesh
  -out ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd7/dout.mesh)
SET(passRegex "\\*\\* UNABLE TO OPEN.*dout.mesh")
SET_PROPERTY(TEST LeakCheck_AbnormalEnd7
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
ADD_TEST(NAME LeakCheck_AbnormalEnd8
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd8/d.mesh)
SET(passRegex "\\*\\* UNABLE TO OPEN.*d.o.sol")
SET_PROPERTY(TEST LeakCheck_AbnormalEnd8
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
#####
ADD_TEST(NAME LeakCheck_args0
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/LeakCheck_args0/d.mesh)
#####
ADD_TEST(NAME LeakCheck_args1
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  -in ${REG_TESTS_PATH}/LeakCheck_args1/d.mesh -sol
  ${REG_TESTS_PATH}/LeakCheck_args1/dsol.sol
  -out ${REG_TESTS_PATH}/LeakCheck_args1/dout.mesh)

###############################################################################
#####
#####         Check Precision
#####
###############################################################################
#####
ADD_TEST(NAME MeshVersionFormatted1
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  -in ${REG_TESTS_PATH}/MeshVersionFormatted1/d.mesh
  -sol ${REG_TESTS_PATH}/MeshVersionFormatted1/dsol.sol)
#####
if(WITH_MEMCHECK)
  MESSAGE(STATUS "coucou")
else()
  ADD_TEST(NAME MeshVersionFormatted2
    COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
    -in ${REG_TESTS_PATH}/MeshVersionFormatted2/d.mesh
    -sol ${REG_TESTS_PATH}/MeshVersionFormatted2/dsol.sol)
endif()
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
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/CubeIso_h_hminMax/CubeIso0.1.mesh -hmax 0.1 -hmin 0.1)
ADD_TEST(NAME CubeIso_0.5h_hminMax
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/CubeIso_0.5h_hminMax/CubeIso0.1.mesh -hmax 0.05 -hmin 0.05)
ADD_TEST(NAME CubeIso_0.25h_hminMax
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/CubeIso_0.25h_hminMax/CubeIso0.1.mesh -hmax 0.025 -hmin 0.025)
#---Second with sol file
ADD_TEST(NAME CubeIso_h_met
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/CubeIso_h_met/CubeIso0.1.mesh)
ADD_TEST(NAME CubeIso_0.5h_met
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/CubeIso_0.5h_met/CubeIso0.1.mesh)
ADD_TEST(NAME CubeIso_0.25h_met
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/CubeIso_0.25h_met/CubeIso0.1.mesh)

#####

# Check what happend when we refine a sphere of size h with a constant metric
# (h, h/2, h/4 and h/8)
#---First with hmin=hmax
ADD_TEST(NAME SphereIso_h_hminMax
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/SphereIso_h_hminMax/SphereIso0.5.mesh -hmax 0.5 -hmin 0.5)
ADD_TEST(NAME SphereIso_0.5h_hminMax
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/SphereIso_0.5h_hminMax/SphereIso0.5.mesh -hmax 0.25 -hmin 0.25)
ADD_TEST(NAME SphereIso_0.25h_hminMax
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/SphereIso_0.25h_hminMax/SphereIso0.5.mesh -hmax 0.125 -hmin 0.125)
ADD_TEST(NAME SphereIso_0.125h_hminMax
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/SphereIso_0.125h_hminMax/SphereIso0.5.mesh -hmax 0.0625 -hmin 0.0625)
#---Second with sol file
ADD_TEST(NAME SphereIso_h_met
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/SphereIso_h_met/SphereIso0.5.mesh)
ADD_TEST(NAME SphereIso_0.5h_met
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/SphereIso_0.5h_met/SphereIso0.5.mesh)
ADD_TEST(NAME SphereIso_0.25h_met
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/SphereIso_0.25h_met/SphereIso0.5.mesh)
ADD_TEST(NAME SphereIso_0.125h_met
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/SphereIso_0.125h_met/SphereIso0.5.mesh)

# Check what happend when we refine a cube whose skin has already the good size
ADD_TEST(NAME CubeSkin_0.05
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/CubeSkin_0.05/CubeSkin_0.05.mesh -hmax 0.05 -hmin 0.05)
ADD_TEST(NAME CubeSkin_0.1
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/CubeSkin_0.1/CubeSkin_0.1.mesh -hmax 0.1 -hmin 0.1)
ADD_TEST(NAME CubeSkin_0.2
  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
  ${REG_TESTS_PATH}/CubeSkin_0.2/CubeSkin_0.2.mesh -hmax 0.2 -hmin 0.2)


# Check results on various meshes
# First: Meshes that we want unrefined
#ADD_TEST(NAME VariousMeshes_unrefined_santa
#  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
#  ${REG_TESTS_PATH}/VariousMeshes_unrefined_santa/santa.meshb -hmax 0.5 -hmin 0.5)

# Second: Meshes that we want refined




# Compare with a reference result when we run
#ADD_TEST(NAME RefCube
#  COMMAND $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -v 5
#  ${REG_TESTS_PATH}/RefCube/cube.mesh) marre... a finir



