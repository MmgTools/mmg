## =============================================================================
##  This file is part of the Mmg software package for the tetrahedral
##  mesh modification.
##**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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
SET ( test_names mmgs_SimpleTeapot )
SET ( input_paths ${MMGS_CI_TESTS}/Teapot )
SET ( input_files teapot )
SET ( args  "-v 5" )
SET ( common_args "" )

ADD_RUN_AGAIN_TESTS ( ${EXECUT_MMGS} "${test_names}" "${args}" "${input_paths}" "${input_files}" )

MMG_ADD_TEST(mmgs_CubeAni
  "${EXECUT_MMGS}"
  "${MMGS_CI_TESTS}/CubeAni" "cube"
  )

MMG_ADD_TEST(mmgs_SphereAni
  "${EXECUT_MMGS}"
  "${MMGS_CI_TESTS}/SphereAni" "sphere"
  )

###############################################################################
#####
#####         Options
#####
###############################################################################

MMG_ADD_TEST(mmgs_memOption
  "${EXECUT_MMGS} -v 5 -m 100 ${common_args}"
  "${MMGS_CI_TESTS}/Teapot" "teapot"
  )

MMG_ADD_TEST(mmgs_val
  "${EXECUT_MMGS} -val"
  "${MMGS_CI_TESTS}/Teapot" "teapot"
  )

# nsd
MMG_ADD_TEST(mmgs_nsd24
  "${EXECUT_MMGS} -v 5 -nsd 24 ${common_args}"
  "${MMGS_CI_TESTS}/Teapot" "teapot"
  )

# MMG_ADD_TEST(mmgs_default
# "${EXECUT_MMGS} -default"
# "${MMGS_CI_TESTS}/Teapot" "teapot"
# )

SET_PROPERTY(TEST mmgs_val #mmgs_default
  PROPERTY WILL_FAIL TRUE)

###############################################################################
#####
#####         Input/Output
#####
###############################################################################

# Binary gmsh
MMG_ADD_TEST(mmgs_binary_gmsh_s
  "${EXECUT_MMGS} -v 5 ${common_args}"
  "${MMGS_CI_TESTS}/GmshInout" "cube.mshb"
  )

# Ascii gmsh
MMG_ADD_TEST(mmgs_ascii_gmsh_s
  "${EXECUT_MMGS} -v 5 ${common_args}"
  "${MMGS_CI_TESTS}/GmshInout" "cube.msh"
  )

# VTK .vtp no metric
MMG_ADD_TEST(mmgs_vtkvtp
  "${EXECUT_MMGS} -v 5"
  "${MMGS_CI_TESTS}/VtkInout" "c1.vtp"
  )

# VTK .vtp with iso metric
MMG_ADD_TEST(mmgs_vtkvtp_iso
  "${EXECUT_MMGS} -v 5"
  "${MMGS_CI_TESTS}/VtkInout" "iso.vtp"
  )

# VTK .vtp with aniso metric
MMG_ADD_TEST(mmgs_vtkvtp_ani
  "${EXECUT_MMGS} -v 5"
  "${MMGS_CI_TESTS}/VtkInout" "ani.vtp"
  )

IF ( (NOT VTK_FOUND) OR USE_VTK MATCHES OFF )
  SET(expr "VTK library not founded")
  SET_PROPERTY(TEST mmgs_vtkvtp
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  SET_PROPERTY(TEST mmgs_vtkvtp_iso
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  SET_PROPERTY(TEST mmgs_vtkvtp_ani
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
ENDIF ( )

###############################################################################
#####
#####         Check Memory Leaks
#####
###############################################################################


###############################################################################
#####
#####         Manifold cases
#####
###############################################################################
MMG_ADD_TEST(mmgs_Rhino_M
  "${EXECUT_MMGS} -v 5 ${common_args} -hausd 1"
  "${MMGS_CI_TESTS}/Rhino_M" "rhino"
  )

MMG_ADD_TEST(mmgs_moebius
  "${EXECUT_MMGS} -v 5 ${common_args} -d -nr"
  "${MMGS_CI_TESTS}/Moebius-strip" "moebius-strip.mesh"
  )

###############################################################################
#####
#####         Non manifold cases
#####
###############################################################################
MMG_ADD_TEST ( mmgs_Cow_NM_hausd10
  "${EXECUT_MMGS} -v 5 ${common_args} -hausd 10"
  "${MMGS_CI_TESTS}/Cow_NM" "cow"
  )

###############################################################################
#####
#####         Test results
#####
###############################################################################
# Test the Ls option
MMG_ADD_TEST(mmgs_OptLs_val
  "${EXECUT_MMGS} -v 5 -ls -val"
  "${MMGS_CI_TESTS}/OptLs_teapot" "teapot"
  )

# MMG_ADD_TEST(mmgs_OptLs_default
#  "${EXECUT_MMGS} -v 5 -ls -default"
#  "${MMGS_CI_TESTS}/OptLs_teapot" "teapot"
# )

SET_PROPERTY(TEST mmgs_OptLs_val #mmgs_OptLs_default
  PROPERTY WILL_FAIL TRUE)

MMG_ADD_TEST(mmgs_OptLs_teapot
  "${EXECUT_MMGS} -v 5 -ls ${common_args}"
  "${MMGS_CI_TESTS}/OptLs_teapot" "teapot"
  )

MMG_ADD_TEST(mmgs_OptLs_teapot_keepRef
  "${EXECUT_MMGS} -v 5 -ls -keep-ref ${common_args}"
  "${MMGS_CI_TESTS}/OptLs_teapot" "teapot"
  )

MMG_ADD_TEST(mmgs_OptLs_teapot_0.5_keepRef
  "${EXECUT_MMGS} -v 5 -ls 0.5 -keep-ref ${common_args}"
  "${MMGS_CI_TESTS}/OptLs_teapot" "teapot"
  )

MMG_ADD_TEST(mmgs_OptLs_teapot2
  "${EXECUT_MMGS} -v 5 -ls -nr ${common_args}"
  "${MMGS_CI_TESTS}/OptLs_teapot" "teapot"
  )

MMG_ADD_TEST(
  mmgs_OptLs_isoref_defaut
  "${EXECUT_MMGS} -v 5 -ls
  -sol ${MMGS_CI_TESTS}/OptLs_isoref/surf-mesh.sol"
  "${MMGS_CI_TESTS}/OptLs_isoref" "surf-mesh.mesh"
  )
MMG_ADD_TEST(
  mmgs_OptLs_isoref_5
  "${EXECUT_MMGS} -v 5 -isoref 5 -ls
  -sol ${MMGS_CI_TESTS}/OptLs_isoref/surf-mesh.sol"
  "${MMGS_CI_TESTS}/OptLs_isoref" "surf-mesh-isoref5.mesh"
  )

if (BASH)
  add_test(
    NAME mmgs_optLs_isoref
    COMMAND ${BASH} -c "diff <(wc -wl ${CTEST_OUTPUT_DIR}/mmgs_isoref.o.mesh  | awk '{print $1 $2}') <(wc -wl ${CTEST_OUTPUT_DIR}/mmgs_isoref5.o.mesh | awk '{print $1 $2}')"
    )
endif()


####### -met option
MMG_ADD_TEST(mmgs_2squares-withMet
  "${EXECUT_MMGS} -v 5  -met ${MMG2D_CI_TESTS}/2squares/2s.sol"
  "${MMG2D_CI_TESTS}/2squares" "2squares"
  )

####### -sol option
MMG_ADD_TEST(mmgs_2squares-withSol
  "${EXECUT_MMGS} -v 5  -sol ${MMG2D_CI_TESTS}/2squares/2s.sol"
  "${MMG2D_CI_TESTS}/2squares" "2squares"
  )

####### orphan points
MMG_ADD_TEST(mmgs_2squares-orphan
  "${EXECUT_MMGS} -v 5 -nsd 10"
  "${MMG2D_CI_TESTS}/2squares" "2squares"
  )

# nsd + ls
MMG_ADD_TEST(mmgs_OptLs_teapot-nsd3
  "${EXECUT_MMGS} -v 5 -ls -nsd 3 ${common_args}"
  "${MMGS_CI_TESTS}/OptLs_teapot" "teapot"
  )


###############################################################################
#####
#####         Detected Bugs
#####
###############################################################################
MMG_ADD_TEST(mmgs_Car_NM
  "${EXECUT_MMGS} -v 5 ${common_args}"
  "${MMGS_CI_TESTS}/Car_NM" "car"
  )

MMG_ADD_TEST(mmgs_Cow_NM_hausd20
  "${EXECUT_MMGS} -v 5 ${common_args}  -hausd 20"
  "${MMGS_CI_TESTS}/Cow_NM" "cow"
  )

###############################################################################
#####
#####         Implicit domain discretization
#####
###############################################################################
MMG_ADD_TEST(mmgs_LSMultiMat_val
  "${EXECUT_MMGS} -v 5 -ls -hmin 0.005 -hmax 0.1 -hausd 0.001 -hgrad 1.3 -val
  -met ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMGS_CI_TESTS}/LSMultiMat" "multi-mat"
  )

#MMG_ADD_TEST(mmgs_LSMultiMat_default
#  "${EXECUT_MMGS} -v 5 -ls -hmin 0.005 -hmax 0.1 -hausd 0.001 -hgrad 1.3 -default
#  -met ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-met.sol
#  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
#  "${MMGS_CI_TESTS}/LSMultiMat" "multi-mat"
#  )

SET_PROPERTY(TEST mmgs_LSMultiMat_val #mmgs_LSMultiMat_default
  PROPERTY WILL_FAIL TRUE)


MMG_ADD_TEST(mmgs_LSMultiMat
  "${EXECUT_MMGS} -v 5 -ls -hmin 0.005 -hmax 0.1 -hausd 0.001 -hgrad 1.3
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMGS_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# non 0 ls
MMG_ADD_TEST(mmgs_LSMultiMat_nonzero
  "${EXECUT_MMGS} -v 5 -ls 0.01 -hausd 0.001
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMGS_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + optim option
MMG_ADD_TEST(mmgs_LSMultiMat_optim
  "${EXECUT_MMGS} -v 5 -ls -optim -hausd 0.001
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMGS_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + optim + aniso option
MMG_ADD_TEST(mmgs_LSMultiMat_optimAni
  "${EXECUT_MMGS} -v 5 -ls -optim -A -hausd 0.001
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMGS_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + hsiz option
MMG_ADD_TEST(mmgs_LSMultiMat_hsiz
  "${EXECUT_MMGS} -v 5 -ls -hsiz 0.05 -hausd 0.001
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMGS_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + hsiz Ani option
MMG_ADD_TEST(mmgs_LSMultiMat_hsizAni
  "${EXECUT_MMGS} -v 5 -ls -hsiz 0.05 -A -hausd 0.001
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMGS_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + metric
MMG_ADD_TEST(mmgs_LSMultiMat_withMet
  "${EXECUT_MMGS} -v 5 -ls -hausd 0.001
  -met ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMGS_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + metric + ls
MMG_ADD_TEST(mmgs_LSMultiMat_withMetAndLs
  "${EXECUT_MMGS} -v 5 -ls -hausd 0.001
  -met ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMGS_CI_TESTS}/LSMultiMat" "multi-mat"
  )
