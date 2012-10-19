###############################################################################
#####
#####         Continuous Integration
#####
###############################################################################

SET(REG_TESTS_PATH ${CMAKELISTS_PATH}/../../RegTests)

# simple test: must already pass
ADD_TEST(SimpleCube $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/Cube/cube.mesh)

###############################################################################
#####
#####         Check Build
#####
###############################################################################



###############################################################################
#####
#####         Check Memory Leak
#####
###############################################################################

ADD_TEST(LeakCheck_AbnormalEnd0 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd0/d.mesh)
#ADD_TEST(LeakCheck_AbnormalEnd1 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd1) ## no possible for now
ADD_TEST(LeakCheck_AbnormalEnd2 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd2/d.mesh)
ADD_TEST(LeakCheck_AbnormalEnd3 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd3/d.mesh -sol ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd3/dsol.sol -ls 2)
#ADD_TEST(LeakCheck_AbnormalEnd4 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd4) ## no possible for now
#ADD_TEST(LeakCheck_AbnormalEnd5.1 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd5.1) ## todo mmg3d2
#ADD_TEST(LeakCheck_AbnormalEnd5.2 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd5.2) ## todo mmg3d1
#ADD_TEST(LeakCheck_AbnormalEnd6 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd6) ## todo analys.c
ADD_TEST(LeakCheck_AbnormalEnd7 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd7/d.mesh -out dout.mesh)
ADD_TEST(LeakCheck_AbnormalEnd8 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_AbnormalEnd8/d.mesh)

ADD_TEST(LeakCheck_args0 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/LeakCheck_args0/d.mesh)
ADD_TEST(LeakCheck_args1 $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 -in ${REG_TESTS_PATH}/LeakCheck_args1/d.mesh -sol ${REG_TESTS_PATH}/LeakCheck_args1/dsol.sol -out dout.mesh)