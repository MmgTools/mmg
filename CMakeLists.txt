CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (MMG3D5)
# Continuous integration
ENABLE_TESTING()

# To see falgs and options of compilation
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
  IF(APPLE)
    # Add flags to the compiler to work on old mac
    ADD_DEFINITIONS( -mmacosx-version-min=10.4 -arch x86_64)

    # To avoid pbs with binary files...
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -arch x86_64")
    #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
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

  #ELSEIF(CMAKE_COMPILER_IS_INTEL)
  # Static Link??
  #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Bstatic") # all is static
ENDIF()

# Static Link??
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")

IF(NOT CMAKE_COMPILER_IS_CLANG)
  # Compiler options for profiling... but not possible with clang
  OPTION ( PROFILING "Enable/Disable PROFILING" OFF )
  IF(PROFILING)
    ADD_DEFINITIONS(-pg)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  ENDIF(PROFILING)
ENDIF(NOT CMAKE_COMPILER_IS_CLANG)

# should we use SCOTCH
OPTION ( USE_SCOTCH "Use SCOTCH TOOL for renumbering" OFF)

# add SCOTCH library?
#
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
    #SET(CMAKE_REQUIRED_INCLUDES ${SCOTCH_INCLUDE_DIR})
    #SET(CMAKE_REQUIRED_LIBRARIES ${SCOTCH_LIBRARY})
    #SET(CMAKE_REQUIRED_LIBRARIES ${SCOTCHERR_LIBRARY})

    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_SCOTCH")

  ELSE ()
    MESSAGE("Please set the environment variable SCOTCH_DIR to your SCOTCH directory whose name should contain the \"scotch\" string" )
  ENDIF()

ENDIF ( USE_SCOTCH )


# Add openMP flags if OpenMP found
OPTION ( OPENMP "Enable/Disable OpenMP" OFF )
IF(OPENMP)
  INCLUDE(FindOpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    MESSAGE(STATUS "OpenMP added: " ${CMAKE_C_FLAGS})
  ELSE(OPENMP_FOUND)
    MESSAGE(STATUS "OpenMP not found")
  ENDIF(OPENMP_FOUND)
ENDIF(OPENMP)

# Header files
INCLUDE_DIRECTORIES(${MMG_SOURCE_DIR})
# Source files
FILE(
  GLOB
  source_files
  ${MMG_SOURCE_DIR}/*.c
  )

#compil.date : if problem with the file creation, see ADD_CUSTOM_TARGET
IF (NOT EXISTS ${MMG_SOURCE_DIR}/compil.date)

  SET_SOURCE_FILES_PROPERTIES( ${MMG_SOURCE_DIR}/compil.date PROPERTIES GENERATED TRUE )

  ADD_CUSTOM_COMMAND(
    OUTPUT ${MMG_SOURCE_DIR}/compil.date
    COMMAND ${CMAKE_COMMAND} -E echo \"\#define COMPIL \" '\"' `date` '\"' > ${MMG_SOURCE_DIR}/compil.date
    DEPENDS ${MMG_SOURCE_DIR}
    COMMENT "compil.date creation"
    )
ELSEIF(EXISTS ${MMG_SOURCE_DIR}/compil.date)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${MMG_SOURCE_DIR}/compil.date
    COMMAND ${CMAKE_COMMAND} -E echo \"\#define COMPIL \" '\"' `date` '\"' > ${MMG_SOURCE_DIR}/compil.date
    DEPENDS ${source_files}
    COMMENT "compil.date modification"
    )
ENDIF(NOT EXISTS ${MMG_SOURCE_DIR}/compil.date)

# Executable and dependencies
ADD_EXECUTABLE(mmg3d5 ${MMG_SOURCE_DIR}/compil.date ${source_files})

# Uncomment to link with mmg3d4 library
#FIND_LIBRARY(mmglib mmg3dlib4.0 "${MMG_SOURCE_DIR}/../../mmg3d4/build")
#TARGET_LINK_LIBRARIES(mmg3d5 ${mmglib})
# End uncomment

SET_TARGET_PROPERTIES(mmg3d5 PROPERTIES DEBUG_POSTFIX _debug) # in debug mode we name the executable mmg3d5_debug
SET_TARGET_PROPERTIES(mmg3d5 PROPERTIES RELEASE_POSTFIX _O3) # in Release mode we name the executable mmg3d5_O3
SET_TARGET_PROPERTIES(mmg3d5 PROPERTIES MINSIZEREL_POSTFIX _Os)
# /!\ ADD_CUSTOM_TARGET don't work with cmake-2.8.4 and 2.8.10.2 but is necessary for other versions

IF ( ${CMAKE_COMMAND} MATCHES 2.8.4  )
ELSEIF (${CMAKE_COMMAND} VERSION_GREATER 2.8.10.0 )
  ADD_CUSTOM_TARGET(compil_date DEPENDS ${MMG_SOURCE_DIR}/compil.date)
  IF (APPLE)
    EXEC_PROGRAM(uname ARGS -v  OUTPUT_VARIABLE DARWIN_VERSION)
    STRING(REGEX MATCH "[0-9]+" DARWIN_VERSION ${DARWIN_VERSION})
    MESSAGE(STATUS "DARWIN_VERSION=${DARWIN_VERSION}")
    IF (DARWIN_VERSION GREATER 9)
      MESSAGE(STATUS "Warning: no update for compil.date because cmake-2.8.10.1 and 2.8.10.2 don't work on lion/mountain lion.")
    ENDIF (DARWIN_VERSION GREATER 9)
  ENDIF(APPLE)
ELSE ()
  ADD_CUSTOM_TARGET(compil_date DEPENDS ${MMG_SOURCE_DIR}/compil.date)
ENDIF ( ${CMAKE_COMMAND} MATCHES 2.8.4 )

ADD_DEPENDENCIES(mmg3d5 compil_date)

# Libraries
FIND_LIBRARY(M_LIB m)
TARGET_LINK_LIBRARIES(mmg3d5 ${M_LIB} ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY})

# Continuous integration

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

  # Coverage
  # TARGET_LINK_LIBRARIES(mmg3d5 ${M_LIB} gcov)

  IF (CMAKE_COMPILER_IS_GNUCC)
    SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -g -Wall -fprofile-arcs -ftest-coverage --coverage")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS} -g -fprofile-arcs -ftest-coverage --coverage")
  ENDIF(CMAKE_COMPILER_IS_GNUCC)

  INCLUDE( ${CMAKE_SOURCE_DIR}/testing.cmake )
ENDIF(BUILD_TESTING)
