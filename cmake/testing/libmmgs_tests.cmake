ADD_EXECUTABLE(libmmgs_example0
  ${CMAKE_SOURCE_DIR}/libexamples/mmgs/example0/main.c)

ADD_EXECUTABLE(libmmgs_example1
  ${CMAKE_SOURCE_DIR}/libexamples/mmgs/example1/main.c)

IF ( LIBMMGS_STATIC )
  ENABLE_LANGUAGE (Fortran)
  ADD_EXECUTABLE(libmmgs_fortran
    ${CMAKE_SOURCE_DIR}/libexamples/mmgs/example0_fortran/main.F90)

  TARGET_LINK_LIBRARIES(libmmgs_example0 ${PROJECT_NAME}s_a)
  TARGET_LINK_LIBRARIES(libmmgs_example1 ${PROJECT_NAME}s_a)
  TARGET_LINK_LIBRARIES(libmmgs_fortran ${PROJECT_NAME}s_a)

  INSTALL(TARGETS libmmgs_fortran RUNTIME DESTINATION bin )

ELSEIF ( LIBMMGS_SHARED )
  TARGET_LINK_LIBRARIES(libmmgs_example0 ${PROJECT_NAME}s_so)
  TARGET_LINK_LIBRARIES(libmmgs_example1 ${PROJECT_NAME}s_so)
ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ()

INSTALL(TARGETS libmmgs_example0  RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmgs_example1  RUNTIME DESTINATION bin )
