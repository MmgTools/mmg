CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT (mmg3d5)

# To see flags and options of compilation
#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

# Executable path
SET(EXECUTABLE_OUTPUT_PATH $ENV{HOME}/bin/$ENV{ARCHI}/)
SET(MMG_SOURCE_DIR ${CMAKE_SOURCE_DIR}/sources/)
# Find "exotic" compilers
IF (${CMAKE_C_COMPILER} MATCHES Clang OR ${CMAKE_C_COMPILER} MATCHES clang)
  # using clang
  SET(CMAKE_COMPILER_IS_CLANG TRUE)
ELSEIF(${CMAKE_C_COMPILER} MATCHES Icc OR ${CMAKE_C_COMPILER} MATCHES icc)
  # using icc
  SET(CMAKE_COMPILER_IS_INTEL TRUE)
ENDIF()


IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "-Wno-char-subscripts ${CMAKE_C_FLAGS}")
  IF(APPLE)
    # Hide some options
    MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES
      CMAKE_OSX_DEPLOYMENT_TARGET
      CMAKE_OSX_SYSROOT)
    # Add flags to the compiler to work on old mac
    ADD_DEFINITIONS( -mmacosx-version-min=10.4 -arch x86_64)

    # To avoid pbs with binary files...
    SET(CMAKE_EXE_LINKER_FLAGS "-arch x86_64 ${CMAKE_EXE_LINKER_FLAGS}")

    # Determine if the processor supports 64bit execution
    EXECUTE_PROCESS(
      COMMAND sysctl hw.cpu64bit_capable
      ERROR_QUIET
      OUTPUT_VARIABLE 64_CMD
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    STRING(REGEX REPLACE "^hw.cpu64bit_capable: (.*)" "\\1" 64_BIT "${64_CMD}")
    #  ELSEIF(UNIX)# UNIX must be after APPLE becaus APPLE is UNIX too
  ENDIF()
ENDIF()


IF(NOT CMAKE_COMPILER_IS_CLANG)
  # Compiler options for profiling... but not possible with clang
  OPTION ( PROFILING "Enable/Disable PROFILING" OFF )
  IF(PROFILING)
    ADD_DEFINITIONS(-pg)
    SET(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")
  ENDIF(PROFILING)
ENDIF(NOT CMAKE_COMPILER_IS_CLANG)

###############################################################################
#####
#####         Redefinitions: mmg3d5_redefine.h
#####
###############################################################################
EXECUTE_PROCESS(COMMAND ctags  -x --c-kinds=fpx  ${MMG_SOURCE_DIR}/mmg3d.h
  ERROR_QUIET
  OUTPUT_VARIABLE CTAG)
IF ( CTAG )
  # We want to create mmg3d5_redefine.h but we don't want it to be removed by make clean so:
  # we use a .dummy output that is created/removed by make/make clean. And the true command create
  # the wanted file. We can't use ADD_CUSTOM_TARGET 'cause mmg3d5_redefine.h would be rebuild each
  # time.
  ADD_CUSTOM_COMMAND(OUTPUT .dummy
    COMMAND touch .dummy
    COMMAND ctags -x  --c-kinds=fpx ${MMG_SOURCE_DIR}/mmg3d.h|awk '{printf \"\#define %-30s\t _MMG5_%s\\n\", $$1, $$1}' > ${MMG_SOURCE_DIR}/${PROJECT_NAME}_redefine.h
    COMMAND ctags -x --c-kinds=fpst ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}.h|awk '{print $$1}'|sed -e 's/MMG5_//g'|awk '{printf \"\#define %-30s\\t MMG5_%s\\n\", $$1, $$1}' >> ${MMG_SOURCE_DIR}/${PROJECT_NAME}_redefine.h
    DEPENDS ${MMG_SOURCE_DIR}/mmg3d.h ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}.h
    COMMENT "update of ${PROJECT_NAME}_redefine.h"
    )
ELSE ()
  ADD_CUSTOM_COMMAND(OUTPUT .dummy
    COMMAND touch .dummy
    DEPENDS ${MMG_SOURCE_DIR}/mmg3d.h ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}.h
    )
ENDIF ()

###############################################################################
#####
#####         Choose type of insertion: delaunay kernel or patterns
#####
###############################################################################
# Should we use patterns for insertion istead of delaunay kernel
OPTION ( PATTERN
  "If ON, insertion by patterns, otherwise insertion by delaunay kernel"
  ON)
IF ( NOT PATTERN )
  MESSAGE(STATUS "Vertex insertion by delaunay kernel")
ELSE()
  SET(CMAKE_C_FLAGS "-DPATTERN ${CMAKE_C_FLAGS}")
ENDIF()
# Remove the next line to have this option visible in basic cmake mode
MARK_AS_ADVANCED(PATTERN)

###############################################################################
#####
#####         Remove compilation of singul.c if not needed
#####
###############################################################################
SET(REMOVE_FILE "${MMG_SOURCE_DIR}/singul.c")
# Should we add insertion of singularities part
OPTION ( INSER_SINGULARITIES "Compile the part of code that allow to insert singularties" OFF)
IF ( INSER_SINGULARITIES )
  SET(CMAKE_C_FLAGS "-DSINGUL ${CMAKE_C_FLAGS}")
  SET(REMOVE_FILE "")
  MESSAGE(STATUS "Compilation with insertion of singularities available")
ENDIF()

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################
# Header files
INCLUDE_DIRECTORIES(${MMG_SOURCE_DIR})
# Source files
FILE(
  GLOB
  source_files
  ${MMG_SOURCE_DIR}/*.c   ${MMG_SOURCE_DIR}/*.h
  )
LIST(REMOVE_ITEM source_files
  ${MMG_SOURCE_DIR}/mmg3d.c ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}.c ${REMOVE_FILE})
FILE(
  GLOB
  main_file
  ${MMG_SOURCE_DIR}/mmg3d.c
  )
FILE(
  GLOB
  lib_file
  ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}.c
  )

IF ((NOT WIN32) AND (NOT WIN64))
  FIND_LIBRARY(M_LIB m)
  SET( LIBRARIES ${M_LIB} ${LIBRARIES})
ENDIF()

###############################################################################
#####
#####         Scotch
#####
###############################################################################
# add SCOTCH library?
SET(SCOTCH_DIR "" CACHE PATH "Installation directory for scotch")
INCLUDE(FindScotch.cmake)

IF(SCOTCH_FOUND)
  OPTION ( USE_SCOTCH "Use SCOTCH TOOL for renumbering" ON)
  IF( USE_SCOTCH )
    # Set flags for building test program
    INCLUDE_DIRECTORIES(${SCOTCH_INCLUDE_DIR})

    SET(CMAKE_REQUIRED_INCLUDES ${SCOTCH_INCLUDE_DIR})
    SET(CMAKE_REQUIRED_LIBRARIES ${SCOTCH_LIBRARY})
    SET(CMAKE_REQUIRED_LIBRARIES ${SCOTCHERR_LIBRARY})

    SET(CMAKE_C_FLAGS "-DUSE_SCOTCH ${CMAKE_C_FLAGS}")
    SET(SCOTCH_LINK_FLAGS "-lpthread")
    MESSAGE(STATUS "Compilation with scotch: ${SCOTCH_DIR}$ENV{SCOTCH_DIR}")
    SET( LIBRARIES ${SCOTCH_LINK_FLAGS} ${LIBRARIES})
    SET( LIBRARIES ${SCOTCHERR_LIBRARY} ${LIBRARIES})
    SET( LIBRARIES ${SCOTCH_LIBRARY} ${LIBRARIES})
  ENDIF()
ELSE ( )
  MESSAGE ( WARNING "Scotch library not found:
 Using scotch reduce the execution time of mmg3d (see https://gforge.inria.fr/frs/?group_id=248 to download it).
 If you have already installed Scotch and want to use it, please set the CMake variable or environment variable SCOTCH_DIR to your scotch directory.")
  OPTION ( USE_SCOTCH "Use SCOTCH TOOL for renumbering" OFF)
ENDIF ( )

###############################################################################
#####
#####         Compile mmg3d5 libraries
#####
###############################################################################
SET(LIBRARY_OUTPUT_PATH ${MMG_SOURCE_DIR}/../lib)

OPTION ( LIBMMG3D5_STATIC "Compile static library" ON)
IF ( LIBMMG3D5_STATIC )
  ADD_LIBRARY(${PROJECT_NAME}_a  STATIC
    .dummy ${MMG_SOURCE_DIR}/${PROJECT_NAME}_redefine.h
    ${source_files} ${lib_file} )
  SET_TARGET_PROPERTIES(${PROJECT_NAME}_a PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}_a ${LIBRARIES})
  INSTALL(FILES ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}.h DESTINATION include)
  INSTALL(FILES ${MMG_SOURCE_DIR}/chrono.h DESTINATION include)
  INSTALL(TARGETS ${PROJECT_NAME}_a
   ARCHIVE DESTINATION lib
   LIBRARY DESTINATION lib)
ENDIF()

# Compile shared library
OPTION ( LIBMMG3D5_SHARED "Compile dynamic library" OFF)
IF ( LIBMMG3D5_SHARED )
  ADD_LIBRARY(${PROJECT_NAME}_so SHARED
    .dummy ${MMG_SOURCE_DIR}/${PROJECT_NAME}_redefine.h
    ${source_files} ${lib_file})
  SET_TARGET_PROPERTIES(${PROJECT_NAME}_so PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}_so ${LIBRARIES})
  INSTALL(FILES ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}.h DESTINATION include)
  INSTALL(FILES ${MMG_SOURCE_DIR}/chrono.h DESTINATION include)
  INSTALL(TARGETS ${PROJECT_NAME}_so
   ARCHIVE DESTINATION
   LIBRARY DESTINATION lib)
ENDIF()

###############################################################################
#####
#####         Compile program to test library
#####
###############################################################################
OPTION ( TEST_LIBMMG3D5 "Compile tests of library" OFF)
IF ( TEST_LIBMMG3D5 )
  IF ( LIBMMG3D5_STATIC )
    IF ( INSER_SINGULARITIES AND TEST_LIBMMG3D5 )
      MESSAGE(WARNING "\"libmmg5_example3\" is the only library test available with the flag of insertion of singularities set to \"ON\"")
    ELSE()
      ADD_EXECUTABLE(libmmg5_example0_a ${MMG_SOURCE_DIR}/../libexamples/example0/example0_a/main.c)
      TARGET_LINK_LIBRARIES(libmmg5_example0_a ${PROJECT_NAME}_a)
      INSTALL(TARGETS libmmg5_example0_a RUNTIME DESTINATION bin )

      ADD_EXECUTABLE(libmmg5_example0_b ${MMG_SOURCE_DIR}/../libexamples/example0/example0_b/main.c)
      TARGET_LINK_LIBRARIES(libmmg5_example0_b ${PROJECT_NAME}_a)
      INSTALL(TARGETS libmmg5_example0_b RUNTIME DESTINATION bin )

      ADD_EXECUTABLE(libmmg5_example1 ${MMG_SOURCE_DIR}/../libexamples/example1/main.c)
      TARGET_LINK_LIBRARIES(libmmg5_example1 ${PROJECT_NAME}_a)
      INSTALL(TARGETS libmmg5_example1 RUNTIME DESTINATION bin )
      ADD_EXECUTABLE(libmmg5_example2 ${MMG_SOURCE_DIR}/../libexamples/example2/main.c)
      TARGET_LINK_LIBRARIES(libmmg5_example2 ${PROJECT_NAME}_a)
      INSTALL(TARGETS libmmg5_example2 RUNTIME DESTINATION bin )

    ENDIF()
    ADD_EXECUTABLE(libmmg5_example3  ${MMG_SOURCE_DIR}/../libexamples/example3/mmg3d.c)
    TARGET_LINK_LIBRARIES(libmmg5_example3 ${PROJECT_NAME}_a)
    INSTALL(TARGETS libmmg5_example3 RUNTIME DESTINATION bin )
  ELSEIF ( LIBMMG3D5_SHARED )
    IF ( INSER_SINGULARITIES AND TEST_LIBMMG3D5 )
      MESSAGE(WARNING "\"libmmg5_example2\" is the only library test available with the flag of insertion of singularities set to \"ON\"")
    ELSE()
      ADD_EXECUTABLE(libmmg5_example0_a ${MMG_SOURCE_DIR}/../libexamples/example0/example0_a/main.c)
      TARGET_LINK_LIBRARIES(libmmg5_example0_a ${PROJECT_NAME}_a)
      INSTALL(TARGETS libmmg5_example0_a RUNTIME DESTINATION bin )

      ADD_EXECUTABLE(libmmg5_example0_b ${MMG_SOURCE_DIR}/../libexamples/example0/example0_b/main.c)
      TARGET_LINK_LIBRARIES(libmmg5_example0_b ${PROJECT_NAME}_a)
      INSTALL(TARGETS libmmg5_example0_b RUNTIME DESTINATION bin )

      ADD_EXECUTABLE(libmmg5_example1 ${MMG_SOURCE_DIR}/../libexamples/example1/main.c)
      TARGET_LINK_LIBRARIES(libmmg5_example1 ${PROJECT_NAME}_so)
      INSTALL(TARGETS libmmg5_example1 RUNTIME DESTINATION bin )
      ADD_EXECUTABLE(libmmg5_example2  ${MMG_SOURCE_DIR}/../libexamples/example2/mmg3d.c)
      TARGET_LINK_LIBRARIES(libmmg5_example2 ${PROJECT_NAME}_so)
      INSTALL(TARGETS libmmg5_example2 RUNTIME DESTINATION bin )
    ENDIF()
    ADD_EXECUTABLE(libmmg5_example3  ${MMG_SOURCE_DIR}/../libexamples/example3/mmg3d.c)
    TARGET_LINK_LIBRARIES(libmmg5_example3 ${PROJECT_NAME}_so)
    INSTALL(TARGETS libmmg5_example3 RUNTIME DESTINATION bin )
  ELSE ()
    MESSAGE(WARNING "You must activate the compilation of the static or shared ${PROJECT_NAME} library to compile this tests." )
  ENDIF ()
ENDIF()

###############################################################################
#####
#####         Compile executable
#####
###############################################################################
ADD_EXECUTABLE(${PROJECT_NAME}
  .dummy  ${MMG_SOURCE_DIR}/${PROJECT_NAME}_redefine.h
  ${source_files} ${main_file})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${LIBRARIES})
INSTALL(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)


SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX _debug) # in debug mode we name the executable mmg3d5_debug
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES RELEASE_POSTFIX _O3) # in Release mode we name the executable mmg3d5_O3
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES MINSIZEREL_POSTFIX _Os)

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################
OPTION ( BUILD_TESTING "Enable/Disable continuous integration" ON )
IF(BUILD_TESTING)
  # Run some tests twice with the output of the previous test as input
  OPTION ( RUN_AGAIN "Enable/Disable second run of some tests" ON )

  # Ignore errors coming from the libmesh5 library
  SET(MEMORYCHECK_SUPPRESSIONS_FILE "${CMAKE_SOURCE_DIR}/libmesh5.supp")
  IF(APPLE)
    # and errors from the libSystem.B library
    SET(MEMORYCHECK_SUPPRESSIONS_FILE "${CMAKE_SOURCE_DIR}/libSystem+libmesh5.supp")
  ENDIF(APPLE)
  # Copy custom file in the good directory
  IF(NOT EXISTS "CTestCustom.cmake")
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/CTestCustom.cmake ./CTestCustom.cmake)
  ENDIF(NOT EXISTS "CTestCustom.cmake")

  ENABLE_TESTING()
  INCLUDE(CTest)

  IF ( CMAKE_COMPILER_IS_GNUCC AND NOT APPLE )
    SET(CMAKE_C_FLAGS_DEBUG "-g -Wall -fprofile-arcs -ftest-coverage --coverage ${CMAKE_C_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "-g -fprofile-arcs -ftest-coverage --coverage ${CMAKE_EXE_LINKER_FLAGS}")
  ENDIF ( )

  # Add runtime that we want to test
  GET_TARGET_PROPERTY(EXECUT ${PROJECT_NAME} LOCATION)
  IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    SET(EXECUT ${EXECUT}_debug)
    SET(BUILDNAME ${BUILDNAME}_debug CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Release")
    SET(EXECUT ${EXECUT}_O3)
    SET(BUILDNAME ${BUILDNAME}_O3 CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Minsizerel")
    SET(EXECUT ${EXECUT}_Os)
    SET(BUILDNAME ${BUILDNAME}_Os CACHE STRING "build name variable")
  ENDIF()
  SET ( LISTEXEC ${EXECUT} )
  IF ( TEST_LIBMMG3D5 )
    IF ( LIBMMG3D5_STATIC OR LIBMMG3D5_SHARED )
      GET_TARGET_PROPERTY(LIBEXEC0_a libmmg5_example0_a LOCATION)
      GET_TARGET_PROPERTY(LIBEXEC0_b libmmg5_example0_b LOCATION)
      GET_TARGET_PROPERTY(LIBEXEC1   libmmg5_example1 LOCATION)
      GET_TARGET_PROPERTY(LIBEXEC2   libmmg5_example2 LOCATION)

      IF ( NOT INSER_SINGULARITIES )
	ADD_TEST(NAME libexample0_a COMMAND ${LIBEXEC0_a})
	ADD_TEST(NAME libexample0_b COMMAND ${LIBEXEC0_b})
	ADD_TEST(NAME libexample1   COMMAND ${LIBEXEC1})
	ADD_TEST(NAME libexample2   COMMAND ${LIBEXEC2})
      ENDIF()
      GET_TARGET_PROPERTY(LIBEXEC3 libmmg5_example3 LOCATION)
      SET( LISTEXEC ${LISTEXEC} ${LIBEXEC3} )
    ENDIF()
  ENDIF()

  # Add more tests
  INCLUDE( ${CMAKE_SOURCE_DIR}/testing.cmake )
  IF ( INSER_SINGULARITIES )
    INCLUDE( ${CMAKE_SOURCE_DIR}/testing_insert_singul.cmake )
  ENDIF ()
  IF ( RUN_AGAIN )
    INCLUDE( ${CMAKE_SOURCE_DIR}/testing_rerun.cmake )
  ENDIF()
ENDIF(BUILD_TESTING)

###############################################################################
#####
#####         Create API Documentation
#####
###############################################################################
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxygen/Doxyfile.in ${CMAKE_CURRENT_SOURCE_DIR}/Doxygen/Doxyfile @ONLY)
  ADD_CUSTOM_TARGET(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxygen/Doxyfile
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/Doxygen/latex make
    # COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/Doxygen/latex/refman.pdf ${CMAKE_CURRENT_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Doxygen
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
ENDIF(DOXYGEN_FOUND)

###############################################################################
#####
#####         Create Package
#####
###############################################################################
INCLUDE(InstallRequiredSystemLibraries)
SET( CPACK_SOURCE_GENERATOR TGZ )
SET( CPACK_GENERATOR TGZ )

SET(CPACK_PACKAGE_NAME ${PROJECT_NAME})
SET(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MMG3D: ISOTROPIC TETRAHEDRAL REMESHER")
SET(CPACK_PACKAGE_VENDOR "Cécile Dobrzynski, Pascal Frey, Charles Dapogny, Algiane Froehly")
SET(CPACK_PACKAGE_CONTACT "cecile.dobrzynski@inria.fr")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/Copyright.txt")
SET(CPACK_PACKAGE_VERSION "5.0.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "5")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}")
SET(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_INSTALL_DIRECTORY}-${CMAKE_SYSTEM})
SET(CPACK_PACKAGE_EXECUTABLES "${PROJECT_NAME}")
SET(CPACK_OUTPUT_FILE_PREFIX "../archives")

SET(CPACK_IGNORE_FILES ";.\\\\\\\\.#;/\\\\\\\\.svn;.~;.#" ${CPACK_IGNORE_FILES})
SET(CPACK_IGNORE_FILES ";libmmg5.*" ${CPACK_IGNORE_FILES})

SET(CPACK_SOURCE_IGNORE_FILES ";.\\\\\\\\.#;/\\\\\\\\.svn;.~;.#" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";lib/;bin;archives;def" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";emacs-format-file" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/example0/example0_a/mesh.o.*" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/example0/example0_b/mesh.o.*" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/example1/result0.*" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/example2/result1.*" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/example3/SphereIso.*" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";Scripts" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";README_developper.txt" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";CMakeFiles" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libmesh5.supp" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libSystem+libmesh5.supp" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";testing.cmake" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";.#$"      ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";build*" ${CPACK_SOURCE_IGNORE_FILES}  )
SET(CPACK_SOURCE_IGNORE_FILES ";Doxygen/html;Doxygen/latex;Doxygen/Doxyfile"  ${CPACK_SOURCE_IGNORE_FILES}  )

IF ((WIN32 OR WIN64) AND NOT UNIX)
  # Warning compilation under windows not tested
  #SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/logo.pdf")
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\${PROJECT_NAME}.exe")
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} MMG3D: ISOTROPIC TETRAHEDRAL REMESHER")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.math.u-bordeaux1.fr/~dobrzyns/.../mmg3d.php‎")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.math.u-bordeaux1.fr/~dobrzyns/.../mmg3d.php‎")
  SET(CPACK_NSIS_CONTACT ${CPACK_PACKAGE_CONTACT})
  SET(CPACK_NSIS_MODIFY_PATH ON) # Choice to add prog dir to the system PATH
  SET(CPACK_PACKAGING_INSTALL_PREFIX  "${PROJECT_NAME}" CACHE PATH "Where you want to install your package")
  MARK_AS_ADVANCED(CPACK_PACKAGING_INSTALL_PREFIX)
ELSEIF (APPLE )
  SET(CPACK_BINARY_OSX11 "ON")
ELSEIF ( UNIX )
ENDIF()

INCLUDE(CPack)
