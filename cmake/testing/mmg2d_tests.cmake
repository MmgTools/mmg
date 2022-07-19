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

###############################################################################
#####
#####         Continuous Integration
#####
###############################################################################
# Simple test: must already pass (-d option allows to cover chkmsh function)

SET ( MMG_SCRIPTS_DIR ${PROJECT_BINARY_DIR}/cmake_scripts )
FILE ( MAKE_DIRECTORY  ${MMG_SCRIPTS_DIR} )

MMG_ADD_TEST ( mmg2d_SimpleCircle
  "${EXECUT_MMG2D} -v 5 -d"
  "${MMG2D_CI_TESTS}/Circle" "cercle" )

###############################################################################
#####
#####         Options
#####
###############################################################################
ADD_TEST(NAME mmg2d_help
  COMMAND ${EXECUT_MMG2D} -h
  )
SET_PROPERTY(TEST mmg2d_help
  PROPERTY PASS_REGULAR_EXPRESSION "File specifications")

ADD_TEST(NAME mmg2d_memOption
  COMMAND ${EXECUT_MMG2D} -v 5 -m 100
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_memOption.o.meshb)

ADD_TEST(NAME mmg2d_val
  COMMAND ${EXECUT_MMG2D} -v val
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_val.o.meshb)

#ADD_TEST(NAME mmg2d_default
#  COMMAND ${EXECUT_MMG2D} -default
#  ${MMG2D_CI_TESTS}/Circle/cercle
#  -out ${CTEST_OUTPUT_DIR}/mmg2d_default.o.meshb)

SET_PROPERTY(TEST mmg2d_val #mmg2d_default
  PROPERTY WILL_FAIL TRUE)

MMG_ADD_TEST(mmg2d_hsizOption
  "${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2"
  "${MMG2D_CI_TESTS}/Circle" "cercle"
  )

MMG_ADD_TEST(mmg2d_hsizAni
  "${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -A"
  "${MMG2D_CI_TESTS}/Circle" "cercle"
  )

MMG_ADD_TEST(mmg2d_hsizAndNosurfOption
  "${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -nosurf"
  "${MMG2D_CI_TESTS}/Circle" "cercle"
  )

MMG_ADD_TEST(mmg2d_hsizAndNosurfAni
  "${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -nosurf -A"
  "${MMG2D_CI_TESTS}/Circle" "cercle"
  )

MMG_ADD_TEST(mmg2d_hsizAndNosurfOption2
  "${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -nosurf -3dMedit 2"
  "${MMG2D_CI_TESTS}/2squares" "2squares"
  )

ADD_TEST(NAME mmg2d_hsizHmax
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -hmax 0.05
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_hsizHmax-circle.o.meshb)
SET(passRegex "Mismatched options")
SET_PROPERTY(TEST mmg2d_hsizHmax
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

ADD_TEST(NAME mmg2d_hsizHmin
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -hmin 0.2
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_hsizHmin-circle.o.meshb)
SET(passRegex "Mismatched options")
SET_PROPERTY(TEST mmg2d_hsizHmin
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

MMG_ADD_TEST(mmg2d_reqEntities-ref
  "${EXECUT_MMG2D} -v 5 -hsiz 0.02"
  "${MMG2D_CI_TESTS}/Disk_ReqEntities" "disk.mesh"
  )

MMG_ADD_TEST(mmg2d_orphanPoint
  "${EXECUT_MMG2D} -v 5 -hausd 10 -hgradreq -1 -nosizreq"
  "${MMG2D_CI_TESTS}/Disk_ReqEntities" "disk.mesh"
  )

MMG_ADD_TEST(mmg2d_reqEntitiesAni-ref
  "${EXECUT_MMG2D} -v 5 -hsiz 0.02 -A"
  "${MMG2D_CI_TESTS}/Disk_ReqEntities" "disk.mesh"
  )

MMG_ADD_TEST(mmg2d_reqEntities-unref
  "${EXECUT_MMG2D} -v 5 -hsiz 0.1"
  "${MMG2D_CI_TESTS}/Disk_ReqEntities" "disk-tiny.mesh"
  )

MMG_ADD_TEST(mmg2d_reqEntitiesAni-unref
  "${EXECUT_MMG2D} -v 5 -hsiz 0.1 -A"
  "${MMG2D_CI_TESTS}/Disk_ReqEntities" "disk-tiny.mesh"
  )

MMG_ADD_TEST(mmg2d_locParam
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/LocParams" "circle2refs.mesh"
  )

MMG_ADD_TEST(mmg2d_locParam_ani
  "${EXECUT_MMG2D} -v 5 -A"
  "${MMG2D_CI_TESTS}/LocParams" "circle2refs.mesh"
  )

MMG_ADD_TEST(mmg2d_opnbdy_yes
  "${EXECUT_MMG2D} -v 5 -opnbdy -hausd 0.001"
  "${MMG2D_CI_TESTS}/Opnbdy" "opnbdy-mesh.msh"
  )

MMG_ADD_TEST(mmg2d_opnbdy_no
  "${EXECUT_MMG2D} -v 5 -hausd 0.001"
  "${MMG2D_CI_TESTS}/Opnbdy" "opnbdy-mesh.msh"
  )

MMG_ADD_TEST(mmg2d_opnbdy_ls
  "${EXECUT_MMG2D} -v 5 -opnbdy -ls 3.4 -hausd 0.001
   -sol  ${MMG2D_CI_TESTS}/Opnbdy/ls.sol -in "
  "${MMG2D_CI_TESTS}/Opnbdy" "opnbdy.mesh"
  )

MMG_ADD_TEST(mmg2d_opnbdy_yes_ani
  "${EXECUT_MMG2D} -v 5 -hausd 0.001 -A -opnbdy"
  "${MMG2D_CI_TESTS}/Opnbdy" "opnbdy-mesh.msh"
  )

# default hybrid
MMG_ADD_TEST(mmg2d_hybrid_2d
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/Hybrid" "hybrid.mesh"
  )

# hybrid opnbdy
MMG_ADD_TEST(mmg2d_hybrid_opnbdy_2d
  "${EXECUT_MMG2D} -v 5 -opnbdy"
  "${MMG2D_CI_TESTS}/Hybrid" "hybrid.mesh"
  )

# hybrid hsiz
MMG_ADD_TEST(mmg2d_hybrid_hsiz_2d
  "${EXECUT_MMG2D} -v 5 -hsiz 0.05 -hgradreq -1"
  "${MMG2D_CI_TESTS}/Hybrid" "hybrid.mesh"
  )

MMG_ADD_TEST(mmg2d_hybrid_nosizreq_nohgradreq_2d
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/Hybrid" "hybrid.mesh -nosizreq -hgradreq -1"
  )

# hybrid nsd: remove the triangular domain as it is of ref 1001
MMG_ADD_TEST(mmg2d_hybrid-nsd1
  "${EXECUT_MMG2D} -v 5 -nsd 1"
  "${MMG2D_CI_TESTS}/Hybrid" "hybrid.mesh"
  )

###############################################################################
#####
#####         Input/Output
#####
###############################################################################

# Binary gmsh no metric
MMG_ADD_TEST(mmg2d_binary_gmsh_2d
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/GmshInout" "cercle1.mshb"
  )

# Ascii gmsh no metric
MMG_ADD_TEST(mmg2d_ascii_gmsh_2d
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/GmshInout" "cercle1.msh"
  )

# Ascii gmsh no metric hybrid
MMG_ADD_TEST(mmg2d_gmsh_hybrid_2d
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/Hybrid" "hybrid.msh"
  )

# Binary gmsh iso metric
MMG_ADD_TEST(mmg2d_binary_gmsh_iso
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/GmshInout" "iso.mshb"
  )

# Ascii gmsh iso metric
MMG_ADD_TEST(mmg2d_ascii_gmsh_iso
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/GmshInout" "iso.msh"
  )

# Binary gmsh iso metric
MMG_ADD_TEST(mmg2d_binary_gmsh_ani
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/GmshInout" "ani.mshb"
  )

# Ascii gmsh iso metric
MMG_ADD_TEST(mmg2d_ascii_gmsh_ani
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/GmshInout" "ani.msh"
  )

# VTK .vtk no metric
MMG_ADD_TEST(mmg2d_vtkvtk
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/VtkInout" "cercle.vtk"
  )

# VTK .vtp no metric
MMG_ADD_TEST(mmg2d_vtkvtp
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/VtkInout" "cercle.vtp"
  )

# VTK .vtu no metric
MMG_ADD_TEST(mmg2d_vtkvtu
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/VtkInout" "cercle.vtu"
  )

# VTK .vtk with iso metric
MMG_ADD_TEST(mmg2d_vtkvtk_iso
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/VtkInout" "iso.vtk"
  )

# VTK .vtp with iso metric
MMG_ADD_TEST(mmg2d_vtkvtp_iso
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/VtkInout" "iso.vtp"
  )

# VTK .vtu with iso metric
MMG_ADD_TEST(mmg2d_vtkvtu_iso
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/VtkInout" "iso.vtu"
  )

# VTK .vtk with aniso metric
MMG_ADD_TEST(mmg2d_vtkvtk_ani
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/VtkInout" "ani.vtk"
  )

# VTK .vtp with aniso metric
MMG_ADD_TEST(mmg2d_vtkvtp_ani
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/VtkInout" "ani.vtp"
  )

# VTK .vtu with aniso metric
MMG_ADD_TEST(mmg2d_vtkvtu_ani
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/VtkInout" "ani.vtu"
  )

IF ( (NOT VTK_FOUND) OR USE_VTK MATCHES OFF )
  SET(expr "VTK library not founded")
  SET_PROPERTY(TEST mmg2d_vtkvtk
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  SET_PROPERTY(TEST mmg2d_vtkvtp
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  SET_PROPERTY(TEST mmg2d_vtkvtu
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  SET_PROPERTY(TEST mmg2d_vtkvtk_iso
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  SET_PROPERTY(TEST mmg2d_vtkvtp_iso
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  SET_PROPERTY(TEST mmg2d_vtkvtu_iso
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  SET_PROPERTY(TEST mmg2d_vtkvtk_ani
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  SET_PROPERTY(TEST mmg2d_vtkvtp_ani
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  SET_PROPERTY(TEST mmg2d_vtkvtu_ani
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
ENDIF()

# Triangle output
#
# Respect the default Tetgen behaviour: saves only boundary edges in
# .edge file.
MMG_ADD_TEST(mmg2d_Circle-triangle
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/Circle" "cercle"
  )

###############################################################################
#####
#####         Isotropic cases
#####
###############################################################################
MMG_ADD_TEST(mmg2d_SquareIso
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/SquareIso" "carretest"
  )

MMG_ADD_TEST(mmg2d_SquareIso_nonConstant
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/SquareIso" "non-constant"
  )

MMG_ADD_TEST(mmg2d_SquareIso_nonConstant2
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/SquareIso" "non-constant-2"
  )

####### -nosurf option
MMG_ADD_TEST(mmg2d_2squares
  "${EXECUT_MMG2D} -3dMedit 2 -hmax 1 -nosurf -v 5"
  "${MMG2D_CI_TESTS}/2squares" "2squares"
  )

####### -nsd
MMG_ADD_TEST(mmg2d_2squares-nsd16
  "${EXECUT_MMG2D} -3dMedit 2 -v 5 -nsd 16"
  "${MMG2D_CI_TESTS}/2squares" "2squares"
  )

####### orphan
MMG_ADD_TEST(mmg2d_2squares-orphan
  "${EXECUT_MMG2D} -3dMedit 2 -v 5 -nsd 10"
  "${MMG2D_CI_TESTS}/2squares" "2squares"
  )

####### -met option
MMG_ADD_TEST(mmg2d_2squares-withMet
  "${EXECUT_MMG2D} -3dMedit 2  -v 5"
  "${MMG2D_CI_TESTS}/2squares" "2squares -met ${MMG2D_CI_TESTS}/2squares/2s.sol"
  )

####### -sol option
MMG_ADD_TEST(mmg2d_2squares-withSol
  "${EXECUT_MMG2D} -3dMedit 2  -v 5"
  "${MMG2D_CI_TESTS}/2squares" "2squares -sol ${MMG2D_CI_TESTS}/2squares/2s.sol"
  )

# -nreg
MMG_ADD_TEST(mmg2d_nreg
  "${EXECUT_MMG2D} -v 5 -nreg"
  "${MMG2D_CI_TESTS}/SquareIso" "carretest"
  )

###############################################################################
#####
#####         Anisotropic cases
#####
###############################################################################
MMG_ADD_TEST(mmg2d_SquareAniso
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/SquareAniso" "adap1"
  )

MMG_ADD_TEST(mmg2d_Circle-optimAni
  "${EXECUT_MMG2D} -v 5 -optim -A -sol 2"
  "${MMG2D_CI_TESTS}/Circle" "cercle"
  )

MMG_ADD_TEST(mmg2d_Circle-hsizAni
  "${EXECUT_MMG2D} -v 5 -hsiz 0.01 -A -sol 2"
  "${MMG2D_CI_TESTS}/Circle" "cercle"
  )

###############################################################################
#####
#####         Mesh generation
#####
###############################################################################
MMG_ADD_TEST(mmg2d_SquareGeneration
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/SquareGeneration" "carretest"
  )

MMG_ADD_TEST(mmg2d_NacaGeneration
  "${EXECUT_MMG2D} -v 5 -hausd 0.001"
  "${MMG2D_CI_TESTS}/NacaGeneration" "naca"
  )

MMG_ADD_TEST(mmg2d_NacaGenerationAni
  "${EXECUT_MMG2D} -v 5 -hausd 0.001 -A"
  "${MMG2D_CI_TESTS}/NacaGeneration" "naca"
  )

# optim
MMG_ADD_TEST(mmg2d_NacaGeneration-optim
  "${EXECUT_MMG2D} -v 5 -hausd 0.001 -optim"
  "${MMG2D_CI_TESTS}/NacaGeneration" "naca"
  )

# hsiz
MMG_ADD_TEST(mmg2d_NacaGeneration-hsiz
  "${EXECUT_MMG2D} -v 5 -hausd 0.001 -hsiz 0.01"
  "${MMG2D_CI_TESTS}/NacaGeneration" "naca"
  )

# hsiz + ani
MMG_ADD_TEST(mmg2d_NacaGeneration-hsizAni
  "${EXECUT_MMG2D} -v 5 -hausd 0.001 -hsiz 0.01 -A"
  "${MMG2D_CI_TESTS}/NacaGeneration" "naca"
  )

# non convex test cases
MMG_ADD_TEST(mmg2d_ACDCGeneration
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/ACDCGeneration" "acdcBdy.mesh"
  )

# nsd option: keep only domain of ref 2
MMG_ADD_TEST(mmg2d_ACDCGeneration-nsd2
  "${EXECUT_MMG2D} -v 5 -nsd 2"
  "${MMG2D_CI_TESTS}/ACDCGeneration" "acdcBdy.mesh"
  )

MMG_ADD_TEST(mmg2d_GaronneGeneration
  "${EXECUT_MMG2D} -v 5"
  "${MMG2D_CI_TESTS}/GaronneGeneration" "garonneEdges.mesh"
  )

###############################################################################
#####
#####         Implicit domain discretization
#####
###############################################################################
MMG_ADD_TEST(mmg2d_LSMultiMat_val
  "${EXECUT_MMG2D} -val -v 5 -ls -hausd 0.001
  -met ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  -sol ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMG2D_CI_TESTS}/LSMultiMat" "multi-mat"
  )

MMG_ADD_TEST(mmg2d_OptLs_Bridge
  "${EXECUT_MMG2D} -v 5 -ls
  -sol ${MMG2D_CI_TESTS}/OptLs_bridge/bridge.sol"
  "${MMG2D_CI_TESTS}/OptLs_bridge" "bridge"
  )

#multi-mat + opnbdy + non-manifold check
MMG_ADD_TEST(mmg2d_LSMultiMat_nm
  "${EXECUT_MMG2D} -v 5 -ls 3 -opnbdy -nr"
  "${MMG2D_CI_TESTS}/LSMultiMat" "2d-opn.mesh"
  )

####### -nsd
MMG_ADD_TEST(mmg2d_LSMultiMat-nsd22
  "${EXECUT_MMG2D} -nsd 22 -v 5 -ls"
  "${MMG2D_CI_TESTS}/LSMultiMat" "multi-mat"
  )

#ADD_TEST(NAME mmg2d_LSMultiMat_default
#  COMMAND ${EXECUT_MMG2D} -val -v 5 -ls -hausd 0.001
#  -met ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-met.sol
#  -sol ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
#  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
#  ${CTEST_OUTPUT_DIR}/mmg2d_multi-mat-default.o.meshb
#  )
SET_PROPERTY(TEST mmg2d_LSMultiMat_val #mmg2d_LSMultiMat_default
  PROPERTY WILL_FAIL TRUE)

MMG_ADD_TEST(mmg2d_LSDiscretization
  "${EXECUT_MMG2D} -v 5 -ls"
  "${MMG2D_CI_TESTS}/LSDiscretization" "dom"
  )

MMG_ADD_TEST(mmg2d_LSDiscretization2
  "${EXECUT_MMG2D} -v 5 -ls"
  "${MMG2D_CI_TESTS}/LSDiscretization" "nacai"
  )

MMG_ADD_TEST(mmg2d_LSMultiMat
  "${EXECUT_MMG2D} -v 5 -ls -hmin 0.005 -hmax 0.1 -hausd 0.001 -hgrad 1.3"
  "${MMG2D_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# non 0 ls
MMG_ADD_TEST(mmg2d_LSMultiMat_nonzero
  "${EXECUT_MMG2D} -v 5 -ls 0.01 -hausd 0.001"
  "${MMG2D_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls + rmc
MMG_ADD_TEST(mmg2d_OptLs_dom_withbub
  "${EXECUT_MMG2D} -v 5 -ls
  -sol ${MMG2D_CI_TESTS}/LSDiscretization/bub.sol"
  "${MMG2D_CI_TESTS}/LSDiscretization" "dom"
  )

# ls + rmc + LSBaseReference
MMG_ADD_TEST(mmg2d_OptLs_LSBaseReferences-rmc
  "${EXECUT_MMG2D} -v 5 -ls -rmc
  -sol ${MMG2D_CI_TESTS}/LSBaseReferences/box.sol"
  "${MMG2D_CI_TESTS}/LSBaseReferences" "box"
  )

MMG_ADD_TEST(mmg2d_OptLs_LSBaseReferences-normc
  "${EXECUT_MMG2D} -v 5 -ls
  -sol ${MMG2D_CI_TESTS}/LSBaseReferences/box.sol"
  "${MMG2D_CI_TESTS}/LSBaseReferences" "box"
  )

# ls + rmc: max pile size bug
MMG_ADD_TEST(mmg2d_OptLs_dom_rmcmaxpile
  "${EXECUT_MMG2D} -v 5 -ls -rmc
  -sol ${MMG2D_CI_TESTS}/LSDiscretization/whole.sol"
  "${MMG2D_CI_TESTS}/LSDiscretization" "dom"
  )

MMG_ADD_TEST(mmg2d_OptLs_dom_rembub
  "${EXECUT_MMG2D} -v 5 -ls
  -sol ${MMG2D_CI_TESTS}/LSDiscretization/bub.sol"
  "${MMG2D_CI_TESTS}/LSDiscretization" "dom"
  )

MMG_ADD_TEST(mmg2d_OptLs_dom_rembub2
  "${EXECUT_MMG2D} -v 5 -ls -rmc 0.1
  -sol ${MMG2D_CI_TESTS}/LSDiscretization/bub.sol"
  "${MMG2D_CI_TESTS}/LSDiscretization" "dom"
  )

MMG_ADD_TEST(
  mmg2d_OptLs_isoref_defaut
  "${EXECUT_MMG2D} -v 5 -ls
  -sol ${MMG2D_CI_TESTS}/OptLs_isoref/2d-mesh.sol"
  "${MMG2D_CI_TESTS}/OptLs_isoref" "2d-mesh.mesh"
  )
MMG_ADD_TEST(
  mmg2d_OptLs_isoref_5
  "${EXECUT_MMG2D} -v 5 -isoref 5 -ls
  -sol ${MMG2D_CI_TESTS}/OptLs_isoref/2d-mesh.sol"
  "${MMG2D_CI_TESTS}/OptLs_isoref" "2d-mesh-isoref5.mesh"
  )

if (BASH)
  add_test(
    NAME mmg2d_optLs_isoref
    COMMAND ${BASH} -c "diff <(wc -wl ${CTEST_OUTPUT_DIR}/mmg2d_isoref.o.mesh  | awk '{print $1 $2}') <(wc -wl ${CTEST_OUTPUT_DIR}/mmg2d_isoref5.o.mesh | awk '{print $1 $2}')"
    )
endif()

# ls discretisation + optim option
MMG_ADD_TEST(mmg2d_LSMultiMat_optim
  "${EXECUT_MMG2D} -v 5 -ls -optim -hausd 0.001"
  "${MMG2D_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + optim + aniso option
MMG_ADD_TEST(mmg2d_LSMultiMat_optimAni
  "${EXECUT_MMG2D} -v 5 -ls -optim -A -hausd 0.001"
  "${MMG2D_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + hsiz option
MMG_ADD_TEST(mmg2d_LSMultiMat_hsiz
  "${EXECUT_MMG2D} -v 5 -ls -hsiz 0.05 -hausd 0.001"
  "${MMG2D_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + hsiz Ani option
MMG_ADD_TEST(mmg2d_LSMultiMat_hsizAni
  "${EXECUT_MMG2D} -v 5 -ls -hsiz 0.05 -A -hausd 0.001"
  "${MMG2D_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + metric
MMG_ADD_TEST(mmg2d_LSMultiMat_withMet
  "${EXECUT_MMG2D} -v 5 -ls -hausd 0.001
  -met ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-met.sol"
  "${MMG2D_CI_TESTS}/LSMultiMat" "multi-mat"
  )

# ls discretisation + metric + ls
MMG_ADD_TEST(mmg2d_LSMultiMat_withMetAndLs
  "${EXECUT_MMG2D} -v 5 -ls -hausd 0.001
  -met ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  -sol ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-sol.sol"
  "${MMG2D_CI_TESTS}/LSMultiMat" "multi-mat"
  )

###############################################################################
#####
#####         Check Lagrangian motion option
#####
###############################################################################
#####
IF ( ELAS_FOUND AND NOT USE_ELAS MATCHES OFF )
  MMG_ADD_TEST(mmg2d_LagMotion0_circle
    "${EXECUT_MMG2D} -v 5  -lag 0"
    "${MMG2D_CI_TESTS}/LagMotion_circle" "circle"
    )
  MMG_ADD_TEST(mmg2d_LagMotion1_circle
    "${EXECUT_MMG2D} -v 5  -lag 1"
    "${MMG2D_CI_TESTS}/LagMotion_circle" "circle"
    )
  MMG_ADD_TEST(mmg2d_LagMotion2_circle
    "${EXECUT_MMG2D} -v 5  -lag 2"
    "${MMG2D_CI_TESTS}/LagMotion_circle" "circle"
    )

  # nsd
  MMG_ADD_TEST(mmg2d_LagMotion2_circle-nsd3
    "${EXECUT_MMG2D} -v 5  -lag 2 -nsd 3"
    "${MMG2D_CI_TESTS}/LagMotion_circle" "circle"
    )

ENDIF()
