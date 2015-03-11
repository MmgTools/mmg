## =============================================================================
##  This file is part of the Mmg software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
##
##  Mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  Mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with Mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the Mmg distribution only if you accept them.
## =============================================================================

###############################################################################
#####
#####         Check Memory Leak
#####
###############################################################################
#####

FOREACH(EXEC ${LISTEXEC})
  ###############################################################################
  #####
  #####         Check Boundaries
  #####
  ###############################################################################
  #####
  ADD_TEST(NAME ChkBdry_optSing_insertFile_${EXEC}
    COMMAND ${EXEC}
    -sing ${REG_TESTS_PATH}/ChkBdry_optSing/cube
    -sol  ${REG_TESTS_PATH}/ChkBdry_optSing/sol.sol
    -sf   ${REG_TESTS_PATH}/ChkBdry_optSing/singularities.mesh
    -out  ${REG_TESTS_PATH}/ChkBdry_optSing/sing.o.meshb)
 #####
  ADD_TEST(NAME ChkBdry_optSing_noInsertFile_${EXEC}
    COMMAND ${EXEC}
    -sing ${REG_TESTS_PATH}/ChkBdry_optSing2/cube
    -sol  ${REG_TESTS_PATH}/ChkBdry_optSing2/sol.sol
    -out  ${REG_TESTS_PATH}/ChkBdry_optSing2/sing.o.meshb)
 #####
  ADD_TEST(NAME ChkBdry_optSing_optLs_${EXEC}
    COMMAND ${EXEC}
    -ls -sing ${REG_TESTS_PATH}/ChkBdry_optSing2/cube
    -sol  ${REG_TESTS_PATH}/ChkBdry_optSing2/ls.sol
    -out  ${REG_TESTS_PATH}/ChkBdry_optSing2/ls_sing.o.meshb)

  ###############################################################################
  #####
  #####         Check Insertion of singularities mode
  #####
  ###############################################################################
  #####
  ADD_TEST(NAME LeakCheck_optSing_insertFile_${EXEC}
    COMMAND ${EXEC}
    -sing ${REG_TESTS_PATH}/LeakCheck_optSing/cube
    -sol ${REG_TESTS_PATH}/LeakCheck_optSing/sol.sol
    -sf ${REG_TESTS_PATH}/LeakCheck_optSing/singularities.mesh
    -out ${REG_TESTS_PATH}/LeakCheck_optSing/cube.o.meshb)
  #####
  ADD_TEST(NAME LeakCheck_optSing_noInsertFile_${EXEC}
    COMMAND ${EXEC}
    -sing ${REG_TESTS_PATH}/LeakCheck_optSing2/cube
    -sol ${REG_TESTS_PATH}/LeakCheck_optSing2/sol.sol
    -out ${REG_TESTS_PATH}/LeakCheck_optSing2/sing.o.meshb)
  #####
  ADD_TEST(NAME LeakCheck_noSing_insertFile_${EXEC}
    COMMAND ${EXEC}
    ${REG_TESTS_PATH}/LeakCheck_optSing2/cube
    -sf ${REG_TESTS_PATH}/LeakCheck_optSing2/singularities.mesh
    -out ${REG_TESTS_PATH}/LeakCheck_optSing2/insert.o.meshb)
  #####
  ADD_TEST(NAME LeakCheck_noSing_noInsertFile_${EXEC}
    COMMAND ${EXEC}
    ${REG_TESTS_PATH}/LeakCheck_optSing2/cube
    -out ${REG_TESTS_PATH}/LeakCheck_optSing2/cube.o.meshb)
  #####
  ADD_TEST(NAME LeakCheck_optSing_optLs_noInsertFile_${EXEC}
    COMMAND ${EXEC}
    -ls -sing ${REG_TESTS_PATH}/LeakCheck_optSing2/cube
    -sol ${REG_TESTS_PATH}/LeakCheck_optSing2/ls.sol
    -out ${REG_TESTS_PATH}/LeakCheck_optSing2/sing_ls.o.meshb)
  #####
  ADD_TEST(NAME LeakCheck_optSing_optLs_insertFile_${EXEC}
    COMMAND ${EXEC}
    -ls -sing ${REG_TESTS_PATH}/LeakCheck_optSing2/cube
    -sf ${REG_TESTS_PATH}/LeakCheck_optSing2/singularities.mesh
    -sol ${REG_TESTS_PATH}/LeakCheck_optSing2/ls.sol
    -out ${REG_TESTS_PATH}/LeakCheck_optSing2/sing_ls_insert.o.meshb)
  #####
  ADD_TEST(NAME LeakCheck_noSing_optLs_noInsertFile_${EXEC}
    COMMAND ${EXEC}
    -ls ${REG_TESTS_PATH}/LeakCheck_optSing2/cube
    -sol ${REG_TESTS_PATH}/LeakCheck_optSing2/ls.sol
    -out ${REG_TESTS_PATH}/LeakCheck_optSing2/ls.o.meshb)
  #####
  ADD_TEST(NAME LeakCheck_noSing_optLs_insertFile_${EXEC}
    COMMAND ${EXEC}
    ${REG_TESTS_PATH}/LeakCheck_optSing2/cube -ls
    -sf ${REG_TESTS_PATH}/LeakCheck_optSing2/singularities.mesh
    -sol ${REG_TESTS_PATH}/LeakCheck_optSing2/ls.sol
    -out ${REG_TESTS_PATH}/LeakCheck_optSing2/ls_insert.o.meshb)
  #####

ENDFOREACH(EXEC)
