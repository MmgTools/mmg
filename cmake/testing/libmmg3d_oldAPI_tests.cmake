ADD_EXECUTABLE(libmmg3d_example0_a_oldAPI
  ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/example0/example0_a/main.c)

ADD_EXECUTABLE(libmmg3d_example0_b_oldAPI
  ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/example0/example0_b/main.c)

ADD_EXECUTABLE(libmmg3d_example1_oldAPI
  ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/example1/main.c)

ADD_EXECUTABLE(libmmg3d_example2_oldAPI
  ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/example2/main.c)

ADD_EXECUTABLE(libmmg3d_example4_oldAPI
  ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/example4/main.c)

IF ( LIBMMG3D_STATIC )
  ENABLE_LANGUAGE (Fortran)
  ADD_EXECUTABLE(libmmg3d_fortran_a_oldAPI
    ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/example0_fortran/example0_a/main.F90)

  ADD_EXECUTABLE(libmmg3d_fortran_b_oldAPI
    ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/example0_fortran/example0_b/main.F90)

  TARGET_LINK_LIBRARIES(libmmg3d_example0_a_oldAPI ${PROJECT_NAME}3d_a)
  TARGET_LINK_LIBRARIES(libmmg3d_example0_b_oldAPI ${PROJECT_NAME}3d_a)
  TARGET_LINK_LIBRARIES(libmmg3d_example1_oldAPI   ${PROJECT_NAME}3d_a)
  TARGET_LINK_LIBRARIES(libmmg3d_example2_oldAPI   ${PROJECT_NAME}3d_a)
  TARGET_LINK_LIBRARIES(libmmg3d_fortran_a_oldAPI  ${PROJECT_NAME}3d_a)
  TARGET_LINK_LIBRARIES(libmmg3d_fortran_b_oldAPI  ${PROJECT_NAME}3d_a)
  TARGET_LINK_LIBRARIES(libmmg3d_example4_oldAPI   ${PROJECT_NAME}3d_a)

  INSTALL(TARGETS libmmg3d_fortran_b_oldAPI  RUNTIME DESTINATION bin )
  INSTALL(TARGETS libmmg3d_fortran_a_oldAPI  RUNTIME DESTINATION bin )

ELSEIF ( LIBMMG3D_SHARED )
  TARGET_LINK_LIBRARIES(libmmg3d_example0_a_oldAPI ${PROJECT_NAME}3d_so)
  TARGET_LINK_LIBRARIES(libmmg3d_example0_b_oldAPI ${PROJECT_NAME}3d_so)
  TARGET_LINK_LIBRARIES(libmmg3d_example1_oldAPI   ${PROJECT_NAME}3d_so)
  TARGET_LINK_LIBRARIES(libmmg3d_example2_oldAPI   ${PROJECT_NAME}3d_so)
  TARGET_LINK_LIBRARIES(libmmg3d_example4_oldAPI   ${PROJECT_NAME}3d_so)

ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ()

INSTALL(TARGETS libmmg3d_example0_a_oldAPI RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmg3d_example0_b_oldAPI RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmg3d_example1_oldAPI   RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmg3d_example2_oldAPI   RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmg3d_example4_oldAPI   RUNTIME DESTINATION bin )
