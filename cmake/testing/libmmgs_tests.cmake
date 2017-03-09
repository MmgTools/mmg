## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

###############################################################################
#####
#####         C Tests
#####
###############################################################################

ADD_EXECUTABLE(libmmgs_example0_a
  ${CMAKE_SOURCE_DIR}/libexamples/mmgs/adaptation_example0/example0_a/main.c ${mmgs_includes})

ADD_EXECUTABLE(libmmgs_example0_b
  ${CMAKE_SOURCE_DIR}/libexamples/mmgs/adaptation_example0/example0_b/main.c ${mmgs_includes})

ADD_EXECUTABLE(libmmgs_example1
  ${CMAKE_SOURCE_DIR}/libexamples/mmgs/adaptation_example1/main.c ${mmgs_includes})

ADD_EXECUTABLE(libmmgs_example2
  ${CMAKE_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_example0/main.c ${mmgs_includes})


 IF ( WIN32 AND ((NOT MINGW) AND USE_SCOTCH) )
    my_add_link_flags(libmmg3d_example0_a "/SAFESEH:NO")
    my_add_link_flags(libmmg3d_example0_b "/SAFESEH:NO")
    my_add_link_flags(libmmgs_example1 "/SAFESEH:NO")
    my_add_link_flags(libmmgs_example2 "/SAFESEH:NO")
 ENDIF ( )

IF ( LIBMMGS_STATIC )

  TARGET_LINK_LIBRARIES(libmmgs_example0_a lib${PROJECT_NAME}s_a)
  TARGET_LINK_LIBRARIES(libmmgs_example0_b lib${PROJECT_NAME}s_a)
  TARGET_LINK_LIBRARIES(libmmgs_example1   lib${PROJECT_NAME}s_a)
  TARGET_LINK_LIBRARIES(libmmgs_example2   lib${PROJECT_NAME}s_a)

ELSEIF ( LIBMMGS_SHARED )

  TARGET_LINK_LIBRARIES(libmmgs_example0_a lib${PROJECT_NAME}s_so)
  TARGET_LINK_LIBRARIES(libmmgs_example0_b lib${PROJECT_NAME}s_so)
  TARGET_LINK_LIBRARIES(libmmgs_example1   lib${PROJECT_NAME}s_so)
  TARGET_LINK_LIBRARIES(libmmgs_example2   lib${PROJECT_NAME}s_so)

ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ()

INSTALL(TARGETS libmmgs_example0_a  RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmgs_example0_b  RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmgs_example1  RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmgs_example2  RUNTIME DESTINATION bin )

###############################################################################
#####
#####         Fortran Tests
#####
###############################################################################

IF (CMAKE_Fortran_COMPILER)
  ENABLE_LANGUAGE (Fortran)

  ADD_EXECUTABLE(libmmgs_fortran_a
    ${CMAKE_SOURCE_DIR}/libexamples/mmgs/adaptation_example0_fortran/example0_a/main.F90)
  ADD_DEPENDENCIES(libmmgs_fortran_a copy_s_headers)

  ADD_EXECUTABLE(libmmgs_fortran_b
    ${CMAKE_SOURCE_DIR}/libexamples/mmgs/adaptation_example0_fortran/example0_b/main.F90)
  ADD_DEPENDENCIES(libmmgs_fortran_b copy_s_headers)

  IF ( WIN32 AND ((NOT MINGW) AND USE_SCOTCH) )
    my_add_link_flags(libmmgs_fortran "/SAFESEH:NO")
  ENDIF ( )

  IF ( LIBMMGS_STATIC )

    TARGET_LINK_LIBRARIES(libmmgs_fortran_a lib${PROJECT_NAME}s_a)
    TARGET_LINK_LIBRARIES(libmmgs_fortran_b lib${PROJECT_NAME}s_a)

  ELSEIF ( LIBMMGS_SHARED )

    TARGET_LINK_LIBRARIES(libmmgs_fortran_a lib${PROJECT_NAME}s_so)
    TARGET_LINK_LIBRARIES(libmmgs_fortran_b lib${PROJECT_NAME}s_so)

  ELSE ()
    MESSAGE(WARNING "You must activate the compilation of the static or"
      " shared ${PROJECT_NAME} library to compile this tests." )
  ENDIF ( )

  INSTALL(TARGETS libmmgs_fortran_a RUNTIME DESTINATION bin )
  INSTALL(TARGETS libmmgs_fortran_b RUNTIME DESTINATION bin )

ENDIF ( CMAKE_Fortran_COMPILER )
