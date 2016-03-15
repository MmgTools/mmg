INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

###############################################################################
#####
#####         C Tests
#####
###############################################################################

ADD_EXECUTABLE(libmmgs_example0
  ${CMAKE_SOURCE_DIR}/libexamples/mmgs/adaptation_example0/main.c ${mmgs_includes})

ADD_EXECUTABLE(libmmgs_example1
  ${CMAKE_SOURCE_DIR}/libexamples/mmgs/adaptation_example1/main.c ${mmgs_includes})

ADD_EXECUTABLE(libmmgs_example2
  ${CMAKE_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_example0/main.c ${mmgs_includes})


 IF ( WIN32 AND ((NOT MINGW) AND USE_SCOTCH) )
    my_add_link_flags(libmmgs_example0 "/SAFESEH:NO")
    my_add_link_flags(libmmgs_example1 "/SAFESEH:NO")
    my_add_link_flags(libmmgs_example2 "/SAFESEH:NO")
 ENDIF ( )

IF ( LIBMMGS_STATIC )

  TARGET_LINK_LIBRARIES(libmmgs_example0 ${PROJECT_NAME}s_a)
  TARGET_LINK_LIBRARIES(libmmgs_example1 ${PROJECT_NAME}s_a)
  TARGET_LINK_LIBRARIES(libmmgs_example2 ${PROJECT_NAME}s_a)

ELSEIF ( LIBMMGS_SHARED )

  TARGET_LINK_LIBRARIES(libmmgs_example0 ${PROJECT_NAME}s_so)
  TARGET_LINK_LIBRARIES(libmmgs_example1 ${PROJECT_NAME}s_so)
  TARGET_LINK_LIBRARIES(libmmgs_example2 ${PROJECT_NAME}s_so)

ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ()

INSTALL(TARGETS libmmgs_example0  RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmgs_example1  RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmgs_example2  RUNTIME DESTINATION bin )

###############################################################################
#####
#####         Fortran Tests
#####
###############################################################################

IF (CMAKE_Fortran_COMPILER)
  ENABLE_LANGUAGE (Fortran)

  ADD_EXECUTABLE(libmmgs_fortran
    ${CMAKE_SOURCE_DIR}/libexamples/mmgs/adaptation_example0_fortran/main.F90 ${mmgs_includes})

  IF ( WIN32 AND ((NOT MINGW) AND USE_SCOTCH) )
    my_add_link_flags(libmmgs_fortran "/SAFESEH:NO")
  ENDIF ( )

  IF ( LIBMMGS_STATIC )
    TARGET_LINK_LIBRARIES(libmmgs_fortran ${PROJECT_NAME}s_a)
  ELSEIF ( LIBMMGS_SHARED )
    TARGET_LINK_LIBRARIES(libmmgs_fortran ${PROJECT_NAME}s_so)
  ELSE ()
    MESSAGE(WARNING "You must activate the compilation of the static or"
      " shared ${PROJECT_NAME} library to compile this tests." )
  ENDIF ( )

  INSTALL(TARGETS libmmgs_fortran RUNTIME DESTINATION bin )

ENDIF ( CMAKE_Fortran_COMPILER )
