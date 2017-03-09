###############################################################################
#####
#####         Generation of a fortran header file
#####
###############################################################################

MACRO ( GENERATE_FORTRAN_HEADER name
    in_dir in_file out_dir out_file
    )
  # Wrap add_custom_command into add_custom target to remove dpendencies from
  # the custom command and thus allow parallel build.
  ADD_CUSTOM_COMMAND (
    OUTPUT ${out_dir}/${out_file}
    COMMAND genheader ${out_dir}/${out_file} ${in_dir}/${in_file}  ${CMAKE_SOURCE_DIR}/scripts/genfort.pl
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS genheader ${in_dir}/${in_file}
    ${CMAKE_SOURCE_DIR}/scripts/genfort.pl
    COMMENT "Generating Fortran header for ${name}"
    )

  ADD_CUSTOM_TARGET (
    ${name}_fortran_header
    ALL
    DEPENDS ${out_dir}/${out_file}
    )

ENDMACRO ( )

###############################################################################
#####
#####         Copy an automatically generated header file to another place
#####
###############################################################################
MACRO ( COPY_FORTRAN_HEADER
    in_dir in_file out_dir out_file
    file_dependencies
    target_name
    )
  # Wrap add_custom_command into add_custom target to remove dpendencies from
  # the custom command and thus allow parallel build.
  ADD_CUSTOM_COMMAND (
    OUTPUT  ${out_dir}/${out_file}
    COMMAND ${CMAKE_COMMAND} -E copy  ${in_dir}/${in_file} ${out_dir}/${out_file}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS ${file_dependencies} )
  ADD_CUSTOM_TARGET ( ${target_name} ALL
    DEPENDS ${out_dir}/${out_file} )
  ADD_DEPENDENCIES ( ${target_name} ${file_dependencies} )

ENDMACRO ( )

###############################################################################
#####
#####         Add a target postfix depending on the build type
#####
###############################################################################

MACRO ( ADD_TARGET_POSTFIX target_name )
  IF ( CMAKE_BUILD_TYPE MATCHES "Debug" )
    # in debug mode we name the executable mmgs_debug
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES DEBUG_POSTFIX _debug)
  ELSEIF ( CMAKE_BUILD_TYPE MATCHES "Release" )
    # in Release mode we name the executable mmgs_O3
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES RELEASE_POSTFIX _O3)
  ELSEIF ( CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo" )
    # in RelWithDebInfo mode we name the executable mmgs_O3d
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES RELWITHDEBINFO_POSTFIX _O3d)
  ELSEIF ( CMAKE_BUILD_TYPE MATCHES "MinSizeRel" )
    # in MinSizeRel mode we name the executable mmgs_O3
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES MINSIZEREL_POSTFIX _Os)
  ENDIF ( )
ENDMACRO ( )

###############################################################################
#####
#####         Add Executable that must be tested by ci
#####
###############################################################################

MACRO ( ADD_EXEC_TO_CI_TESTS exec_name )

  IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}_debug)
    SET(BUILDNAME ${BUILDNAME}_debug CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Release")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}_O3)
    SET(BUILDNAME ${BUILDNAME}_O3 CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "RelWithDebInfo")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}3d_O3d)
    SET(BUILDNAME ${BUILDNAME}_O3d CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "MinSizeRel")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}3d_Os)
    SET(BUILDNAME ${BUILDNAME}_Os CACHE STRING "build name variable")
  ELSE()
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}3d)
    SET(BUILDNAME ${BUILDNAME} CACHE STRING "build name variable")
  ENDIF()

  SET ( LISTEXEC_${exec_name} ${EXECUT} )

ENDMACRO ( )


###############################################################################
#####
#####         Add a library test
#####
###############################################################################

MACRO ( ADD_LIBRARY_TEST target_name main_path target_dependency lib_name )
  ADD_EXECUTABLE ( ${target_name} ${main_path} )
  ADD_DEPENDENCIES( ${target_name} ${target_dependency} )

  TARGET_INCLUDE_DIRECTORIES ( ${target_name} PUBLIC ${CMAKE_SOURCE_DIR}/include )

  IF ( WIN32 AND ((NOT MINGW) AND USE_SCOTCH) )
    MY_ADD_LINK_FLAGS ( ${target_name} "/SAFESEH:NO" )
  ENDIF ( )

  TARGET_LINK_LIBRARIES ( ${target_name} ${lib_name} )
  INSTALL(TARGETS ${target_name} RUNTIME DESTINATION bin )

ENDMACRO ( )
