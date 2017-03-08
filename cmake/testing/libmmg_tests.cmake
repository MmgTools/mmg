INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

###############################################################################
#####
#####         C Tests
#####
###############################################################################

ADD_EXECUTABLE(libmmg_example0_a
  ${CMAKE_SOURCE_DIR}/libexamples/mmg/adaptation_example0/main.c)
ADD_EXECUTABLE(libmmg_cpp_a
  ${CMAKE_SOURCE_DIR}/libexamples/mmg/adaptation_example0_cpp/main.cpp)

IF ( WIN32 AND NOT MINGW AND USE_SCOTCH )
  my_add_link_flags(libmmg_example0_a "/SAFESEH:NO")
  my_add_link_flags(libmmg_cpp_a "/SAFESEH:NO")
ENDIF ( )

IF ( LIBMMG_STATIC )

  TARGET_LINK_LIBRARIES(libmmg_example0_a lib${PROJECT_NAME}_a)
  TARGET_LINK_LIBRARIES(libmmg_cpp_a lib${PROJECT_NAME}_a)

ELSEIF ( LIBMMG_SHARED )

  TARGET_LINK_LIBRARIES(libmmg_example0_a lib${PROJECT_NAME}_so)
  TARGET_LINK_LIBRARIES(libmmg_cpp_a lib${PROJECT_NAME}_so)

ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ()

INSTALL(TARGETS libmmg_example0_a RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmg_cpp_a RUNTIME DESTINATION bin )

###############################################################################
#####
#####         Fortran Tests
#####
###############################################################################
IF (CMAKE_Fortran_COMPILER)
  ENABLE_LANGUAGE (Fortran)

  ADD_EXECUTABLE(libmmg_fortran_a
    ${CMAKE_SOURCE_DIR}/libexamples/mmg/adaptation_example0_fortran/main.F90)
  ADD_DEPENDENCIES(libmmg_fortran_a copy_mmg_headers )

  IF ( WIN32 AND NOT MINGW AND USE_SCOTCH )
    my_add_link_flags(libmmg_fortran_a "/SAFESEH:NO")
  ENDIF ( )

  IF ( LIBMMG_STATIC )

    TARGET_LINK_LIBRARIES(libmmg_fortran_a  lib${PROJECT_NAME}_a)

  ELSEIF ( LIBMMG_SHARED )

    TARGET_LINK_LIBRARIES(libmmg_fortran_a  lib${PROJECT_NAME}_so)

  ELSE ()
    MESSAGE(WARNING "You must activate the compilation of the static or"
      " shared ${PROJECT_NAME} library to compile this tests." )
  ENDIF ()

  INSTALL(TARGETS libmmg_fortran_a  RUNTIME DESTINATION bin )

ENDIF (CMAKE_Fortran_COMPILER)
