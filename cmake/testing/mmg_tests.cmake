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

FOREACH(EXEC ${LISTEXEC_MMG})
  ##############################################################################
  #####
  #####         Check Memory Leak
  #####
  ##############################################################################
  #####
  ADD_TEST(NAME LeakCheck_AbnormalEnd2_${EXEC}
    COMMAND ${EXEC} -v 5
    ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd2/d)
  SET(passRegex "## Unable to scale mesh.")
  SET_PROPERTY(TEST LeakCheck_AbnormalEnd2_${EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  ADD_TEST(NAME LeakCheck_AbnormalEnd7_${EXEC}
    COMMAND ${EXEC} -v 5
    ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd7/d
    -out ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd7/unwrittable.meshb)
  SET(passRegex "\\*\\* UNABLE TO OPEN.*")
  SET_PROPERTY(TEST LeakCheck_AbnormalEnd7_${EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  ADD_TEST(NAME LeakCheck_AbnormalEnd8_${EXEC}
    COMMAND ${EXEC} -v 5
    ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd8/d
    -out ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd8/unwrittable.meshb)
  SET(passRegex "\\*\\* UNABLE TO OPEN.*.sol")
  SET_PROPERTY(TEST LeakCheck_AbnormalEnd8_${EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  #####
  ADD_TEST(NAME LeakCheck_args0_${EXEC}
    COMMAND ${EXEC} -v 5
    ${MMG_CI_TESTS}/LeakCheck_args0/d)
  #####
  ADD_TEST(NAME LeakCheck_args1_${EXEC}
    COMMAND ${EXEC} -v 5
    -in ${MMG_CI_TESTS}/LeakCheck_args1/d -sol
    ${MMG_CI_TESTS}/LeakCheck_args1/dsol.sol
    -out ${MMG_CI_TESTS}/LeakCheck_args1/dout.meshb)

  ##############################################################################
  #####
  #####         Check Precision
  #####
  ##############################################################################
  #####
  ADD_TEST(NAME MeshVersionFormatted1_${EXEC}
    COMMAND ${EXEC} -v 5
    -in ${MMG_CI_TESTS}/MeshVersionFormatted1/d
    -sol ${MMG_CI_TESTS}/MeshVersionFormatted1/dsol.sol)
  #####
  ADD_TEST(NAME MeshVersionFormatted2_${EXEC}
    COMMAND ${EXEC} -v 5
    -in ${MMG_CI_TESTS}/MeshVersionFormatted2/d
    -sol ${MMG_CI_TESTS}/MeshVersionFormatted2/dsol.sol)

ENDFOREACH(EXEC)