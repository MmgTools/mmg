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
#####         Continuous Integration
#####
###############################################################################

# Simple tests: must already pass
ADD_TEST(NAME SimpleTeapot
  COMMAND ${EXECUT_MMGS} -v 6 -d
  ${MMGS_CI_TESTS}/Teapot/teapot
  -out ${MMGS_CI_TESTS}/Teapot/teapot.d.meshb)

ADD_TEST(NAME CubeAni
  COMMAND ${EXECUT_MMGS} -v 6 -d
  ${MMGS_CI_TESTS}/CubeAni/cube
  -out ${MMGS_CI_TESTS}/CubeAni/cube.d.meshb)

ADD_TEST(NAME CubeVolAni
  COMMAND ${EXECUT_MMGS} -v 6 -d
  ${MMGS_CI_TESTS}/CubeVolAni/cube
  -out ${MMGS_CI_TESTS}/CubeVolAni/cube.d.meshb)

ADD_TEST(NAME SphereAni
  COMMAND ${EXECUT_MMGS} -v 6 -d
  ${MMGS_CI_TESTS}/SphereAni/sphere
  -out ${MMGS_CI_TESTS}/SphereAni/sphere.d.meshb)

###############################################################################
#####
#####         Check Memory Leak
#####
###############################################################################
#####
