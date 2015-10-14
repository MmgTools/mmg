ADD_EXECUTABLE(libmmg2d_example0
  ${CMAKE_SOURCE_DIR}/libexamples/mmg2d/example0/main.c)

ADD_EXECUTABLE(libmmg2d_example1
  ${CMAKE_SOURCE_DIR}/libexamples/mmg2d/example1/main.c)

IF ( LIBMMG2D_STATIC )
  ENABLE_LANGUAGE (Fortran)
  ADD_EXECUTABLE(libmmg2d_fortran_a
    ${CMAKE_SOURCE_DIR}/libexamples/mmg2d/example0_fortran/main.F90)

  TARGET_LINK_LIBRARIES(libmmg2d_example0 ${PROJECT_NAME}2d_a)
  TARGET_LINK_LIBRARIES(libmmg2d_example1 ${PROJECT_NAME}2d_a)
  TARGET_LINK_LIBRARIES(libmmg2d_fortran_a ${PROJECT_NAME}2d_a)

  INSTALL(TARGETS libmmg2d_fortran_a RUNTIME DESTINATION bin )

ELSEIF ( LIBMMG2D_SHARED )
  TARGET_LINK_LIBRARIES(libmmg2d_example0 ${PROJECT_NAME}2d_so)
  TARGET_LINK_LIBRARIES(libmmg2d_example1 ${PROJECT_NAME}2d_so)
ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ()

INSTALL(TARGETS libmmg2d_example0  RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmg2d_example1  RUNTIME DESTINATION bin )
