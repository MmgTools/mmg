CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT (MMG3D5)

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
#####         Redefinitions: mmg3d_redefine.h
#####
###############################################################################
EXECUTE_PROCESS(COMMAND ctags  -x --c-kinds=fp  ${MMG_SOURCE_DIR}/mmg3d.h
		ERROR_QUIET
		OUTPUT_VARIABLE CTAG)
IF ( CTAG )
# We want to create mmg3d_redefine.h but we don't want it to be removed by make clean so:
# we use a .dummy output that is created/removed by make/make clean. And the true command create
# the wanted file. We can't use ADD_CUSTOM_TARGET 'cause mmg3d_redefine.h would be rebuild each
# time.
  ADD_CUSTOM_COMMAND(OUTPUT .dummy
    COMMAND touch .dummy
    COMMAND ctags -x  --c-kinds=fp ${MMG_SOURCE_DIR}/mmg3d.h|awk '{printf \"\#define %-30s\t _MMG5_%s\\n\", $$1, $$1}' > ${MMG_SOURCE_DIR}/mmg3d_redefine.h
    COMMAND ctags -x --c-kinds=fpst ${MMG_SOURCE_DIR}/libmmg3d5.h|awk '{print $$1}'|sed -e 's/MMG5_//g'|awk '{printf \"\#define %-30s\\t MMG5_%s\\n\", $$1, $$1}' >> ${MMG_SOURCE_DIR}/mmg3d_redefine.h
    DEPENDS ${MMG_SOURCE_DIR}/mmg3d.h ${MMG_SOURCE_DIR}/libmmg3d5.h
    COMMENT "update of mmg3d_redefine.h"
    )
ELSE ()
  ADD_CUSTOM_COMMAND(OUTPUT .dummy
    COMMAND touch .dummy
    DEPENDS ${MMG_SOURCE_DIR}/mmg3d.h ${MMG_SOURCE_DIR}/libmmg3d5.h
    )
ENDIF ()

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
  ${MMG_SOURCE_DIR}/mmg3d.c ${MMG_SOURCE_DIR}/libmmg3d5.c ${REMOVE_FILE})
FILE(
  GLOB
  main_file
  ${MMG_SOURCE_DIR}/mmg3d.c
  )
FILE(
  GLOB
  lib_file
  ${MMG_SOURCE_DIR}/libmmg3d5.c
  )

FIND_LIBRARY(M_LIB m)


###############################################################################
#####
#####         Scotch
#####
###############################################################################
OPTION ( USE_SCOTCH "Use SCOTCH TOOL for renumbering" OFF)
# add SCOTCH library?
IF ( USE_SCOTCH )
  # Check for SCOTCH Path
  IF ( $ENV{SCOTCH_DIR} MATCHES scotch )

    # Check for header file
    INCLUDE_DIRECTORIES(${SCOTCH_INCLUDE_DIR})
    FIND_PATH(SCOTCH_INCLUDE_DIR scotch.h
      HINTS ${SCOTCH_DIR}/include $ENV{SCOTCH_DIR}/include
      ${SCOTCH_INCLUDE_DIR} $ENV{SCOTCH_INCLUDE_DIR}
      PATH_SUFFIXES scotch
      DOC "Directory of SCOTCH Header")

    # Check for scotch
    FIND_LIBRARY(SCOTCH_LIBRARY
      NAMES scotch scotch${SCOTCH_LIB_SUFFIX}
      HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
      ${SCOTCH_LIB_DIR} $ENV{SCOTCH_LIB_DIR}
      DOC "The SCOTCH library"
      )

    FIND_LIBRARY(SCOTCHERR_LIBRARY
      NAMES scotcherr scotcherr${SCOTCH_LIB_SUFFIX}
      HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
      ${SCOTCH_LIB_DIR} $ENV{SCOTCH_LIB_DIR}
      DOC "The SCOTCH-ERROR library"
      )

    # Set flags for building test program
    SET(CMAKE_REQUIRED_INCLUDES ${SCOTCH_INCLUDE_DIR})
    SET(CMAKE_REQUIRED_LIBRARIES ${SCOTCH_LIBRARY})
    SET(CMAKE_REQUIRED_LIBRARIES ${SCOTCHERR_LIBRARY})

    SET(CMAKE_C_FLAGS "-DUSE_SCOTCH ${CMAKE_C_FLAGS}")
    SET(SCOTCH_LINK_FLAGS "-lpthread")
    MESSAGE(STATUS "Compilation with scotch")
  ELSE ()
    MESSAGE(WARNING "Please, disable the USE_SCOTCH option or set the environment variable SCOTCH_DIR to your SCOTCH directory whose name should contain the \"scotch\" string" )
  ENDIF()

ENDIF ( USE_SCOTCH )

###############################################################################
#####
#####         Update date of compilation: compil.date
#####
###############################################################################
ADD_CUSTOM_TARGET( date COMMAND rm -f ${MMG_SOURCE_DIR}/compil.date )
IF (NOT EXISTS ${MMG_SOURCE_DIR}/compil.date)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${MMG_SOURCE_DIR}/compil.date
    COMMAND ${CMAKE_COMMAND} -E echo \"\#define COMPIL \" '\"' `date` '\"' > ${MMG_SOURCE_DIR}/compil.date
    DEPENDS ${source_file} ${main_file} ${lib_file}
    COMMENT "date modification"
    )
ELSEIF(EXISTS ${MMG_SOURCE_DIR}/compil.date)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${MMG_SOURCE_DIR}/compil.date
    COMMAND ${CMAKE_COMMAND} -E echo \"\#define COMPIL \" '\"' `date` '\"' > ${MMG_SOURCE_DIR}/compil.date
    DEPENDS ${source_files} ${main_file} ${lib_file}
    COMMENT "date modification"
    )
ENDIF(NOT EXISTS ${MMG_SOURCE_DIR}/compil.date)

###############################################################################
#####
#####         Compile mmg3d5 libraries
#####
###############################################################################
SET(LIBRARY_OUTPUT_PATH ${MMG_SOURCE_DIR}/../lib)

OPTION ( LIBMMG3D5_STATIC "Compile static library" OFF)
IF ( LIBMMG3D5_STATIC )
  ADD_LIBRARY(mmg3d5_a  STATIC
    ${MMG_SOURCE_DIR}/compil.date ${MMG_SOURCE_DIR}/mmg3d_redefine.h
    ${source_files} ${lib_file} )
  SET_TARGET_PROPERTIES(mmg3d5_a PROPERTIES OUTPUT_NAME mmg3d5)
  ADD_DEPENDENCIES(mmg3d5_a date)
  TARGET_LINK_LIBRARIES(mmg3d5_a ${M_LIB} ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY} ${SCOTCH_LINK_FLAGS})
  INSTALL(TARGETS mmg3d5_a ARCHIVE DESTINATION lib/static)
ENDIF()

# Compile shared library
OPTION ( LIBMMG3D5_SHARED "Compile dynamic library" OFF)
IF ( LIBMMG3D5_SHARED )
  ADD_LIBRARY(mmg3d5_so SHARED
    ${MMG_SOURCE_DIR}/compil.date ${MMG_SOURCE_DIR}/mmg3d_redefine.h
    ${source_files} ${lib_file})
  SET_TARGET_PROPERTIES(mmg3d5_so PROPERTIES OUTPUT_NAME mmg3d5)
  ADD_DEPENDENCIES(mmg3d5_so date)
  TARGET_LINK_LIBRARIES(mmg3d5_so ${M_LIB} ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY} ${SCOTCH_LINK_FLAGS})
  INSTALL(TARGETS mmg3d5_so LIBRARY DESTINATION lib)
ENDIF()

INSTALL(FILES ${MMG_SOURCE_DIR}/libmmg3d5.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
INSTALL(FILES ${MMG_SOURCE_DIR}/mmg3d_redefine.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

###############################################################################
#####
#####         Compile program to test library
#####
###############################################################################
OPTION ( TEST_LIBMMG3D5 "Compile tests of library" OFF)
IF ( TEST_LIBMMG3D5 )
  IF ( LIBMMG3D5_STATIC )
    IF ( INSER_SINGULARITIES AND TEST_LIBMMG3D5 )
      MESSAGE(WARNING "Test \"libmmg5_example1\" not available with the insertion of singularities flag set to \"ON\"")
    ELSE()
      ADD_EXECUTABLE(libmmg5_example1 ${MMG_SOURCE_DIR}/../libexamples/example1/main.c)
      TARGET_LINK_LIBRARIES(libmmg5_example1 mmg3d5_a)
      INSTALL(TARGETS libmmg5_example1 RUNTIME DESTINATION bin )
    ENDIF()
    ADD_EXECUTABLE(libmmg5_example2  ${MMG_SOURCE_DIR}/../libexamples/example2/mmg3d.c)
    TARGET_LINK_LIBRARIES(libmmg5_example2 mmg3d5_a)
    INSTALL(TARGETS libmmg5_example2 RUNTIME DESTINATION bin )
  ELSEIF ( LIBMMG3D5_SHARED )
    IF ( INSER_SINGULARITIES AND TEST_LIBMMG3D5 )
      MESSAGE(WARNING "Test \"libmmg5_example1\" not available with the insertion of singularities flag set to \"ON\"")
    ELSE()
      ADD_EXECUTABLE(libmmg5_example1 ${MMG_SOURCE_DIR}/../libexamples/example1/main.c)
      TARGET_LINK_LIBRARIES(libmmg5_example1 mmg3d5_so)
      INSTALL(TARGETS libmmg5_example1 RUNTIME DESTINATION bin )
    ENDIF()
    ADD_EXECUTABLE(libmmg5_example2  ${MMG_SOURCE_DIR}/../libexamples/example2/mmg3d.c)
    TARGET_LINK_LIBRARIES(libmmg5_example2 mmg3d5_so)
    INSTALL(TARGETS libmmg5_example2 RUNTIME DESTINATION bin )
  ELSE ()
    MESSAGE(WARNING "You must activate the compilation of the static or shared mmg3d5 library to compile this tests." )
  ENDIF ()
ENDIF()

###############################################################################
#####
#####         Compile executable
#####
###############################################################################
ADD_EXECUTABLE(mmg3d5
  ${MMG_SOURCE_DIR}/compil.date .dummy  ${MMG_SOURCE_DIR}/mmg3d_redefine.h
  ${source_files} ${main_file})
ADD_DEPENDENCIES(mmg3d5 date)
TARGET_LINK_LIBRARIES(mmg3d5 ${M_LIB} ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY} ${SCOTCH_LINK_FLAGS})
INSTALL(TARGETS mmg3d5 RUNTIME DESTINATION bin)


SET_TARGET_PROPERTIES(mmg3d5 PROPERTIES DEBUG_POSTFIX _debug) # in debug mode we name the executable mmg3d5_debug
SET_TARGET_PROPERTIES(mmg3d5 PROPERTIES RELEASE_POSTFIX _O3) # in Release mode we name the executable mmg3d5_O3
SET_TARGET_PROPERTIES(mmg3d5 PROPERTIES MINSIZEREL_POSTFIX _Os)

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################
OPTION ( BUILD_TESTING "Enable/Disable continuous integration" ON )
IF(BUILD_TESTING)
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

  IF (CMAKE_COMPILER_IS_GNUCC)
    SET(CMAKE_C_FLAGS_DEBUG "-g -Wall -fprofile-arcs -ftest-coverage --coverage ${CMAKE_C_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "-g -fprofile-arcs -ftest-coverage --coverage ${CMAKE_EXE_LINKER_FLAGS}")
  ENDIF(CMAKE_COMPILER_IS_GNUCC)

  # Add runtime that we want to test
  GET_TARGET_PROPERTY(EXECUT mmg3d5 LOCATION)
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
      GET_TARGET_PROPERTY(LIBEXEC1 libmmg5_example1 LOCATION)
      ADD_TEST(NAME libexample1	COMMAND ${LIBEXEC1})

      GET_TARGET_PROPERTY(LIBEXEC2 libmmg5_example2 LOCATION)
      SET( LISTEXEC ${LISTEXEC} ${LIBEXEC2} )
    ENDIF()
  ENDIF()

  # Add more tests
  INCLUDE( ${CMAKE_SOURCE_DIR}/testing.cmake )
  IF ( INSER_SINGULARITIES )
    INCLUDE( ${CMAKE_SOURCE_DIR}/testing_insert_singul.cmake )
  ENDIF ()
ENDIF(BUILD_TESTING)

###############################################################################
#####
#####         Create Package
#####
###############################################################################
INCLUDE(InstallRequiredSystemLibraries)
SET( CPACK_SOURCE_GENERATOR TGZ )
SET( CPACK_GENERATOR TGZ )

#SET(CPACK_BINARY_CYGWIN "ON") #Only on windows
#SET(CPACK_BINARY_OSX11 "ON")  #only on Mac

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MMG3D: ISOTROPIC TETRAHEDRAL REMESHER")
SET(CPACK_PACKAGE_VENDOR "Cécile Dobrzynski, Pascal Frey, Charles Dapogny, Algiane Froehly")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/Copyright.txt")
SET(CPACK_PACKAGE_VERSION "5.0.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "5")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "MMG3D ${CPACK_PACKAGE_VERSION}")
SET(CPACK_PACKAGE_EXECUTABLES "mmg3d")


SET(CPACK_IGNORE_FILES "libmmg3d5_example1 libmmg3d5_example2")
SET(CPACK_SOURCES_IGNORE_FILES "${CMAKE_SOURCE_DIR}/README_Developper.txt")
SET(CPACK_SOURCES_IGNORE_FILES "${CMAKE_SOURCE_DIR}/CMakeFiles/.* ${CPACK_SOURCES_IGNORE_FILES}")
SET(CPACK_SOURCES_IGNORE_FILES "${CMAKE_SOURCE_DIR}/libmesh5.supp ${CPACK_SOURCES_IGNORE_FILES}")
SET(CPACK_SOURCES_IGNORE_FILES "${CMAKE_SOURCE_DIR}/libSystem+libmesh5.supp ${CPACK_SOURCES_IGNORE_FILES}")
SET(CPACK_SOURCES_IGNORE_FILES "${CMAKE_SOURCE_DIR}/archives/.* ${CPACK_SOURCES_IGNORE_FILES}")
SET(CPACK_SOURCES_IGNORE_FILES "${CMAKE_SOURCE_DIR}/testing.cmake ${CPACK_SOURCES_IGNORE_FILES}")
SET( CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/.*~ ${CMAKE_SOURCE_DIR}/.*# ${CPACK_SOURCES_IGNORE_FILES}")
SET( CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/build.*/.*"  "${CPACK_SOURCES_IGNORE_FILES}"  )
SET( CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/sources/.*"  "${CPACK_SOURCES_IGNORE_FILES}"  )

SET( CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/.*"  "${CPACK_SOURCES_IGNORE_FILES}"  )

INCLUDE(CPack)
ADD_CUSTOM_TARGET(my_package_source
 COMMAND ${CMAKE_MAKE_PROGRAM} package_source
 COMMAND ${CMAKE_COMMAND} -E copy
${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz ${CMAKE_SOURCE_DIR}
 COMMAND ${CMAKE_COMMAND} -E remove ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz
)
ADD_CUSTOM_TARGET(my_package
 COMMAND ${CMAKE_MAKE_PROGRAM} package
 COMMAND ${CMAKE_COMMAND} -E copy
${CPACK_PACKAGE_FILE_NAME}.tar.gz ${CMAKE_SOURCE_DIR}
 COMMAND ${CMAKE_COMMAND} -E remove ${CPACK_PACKAGE_FILE_NAME}.tar.gz
)
