CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (MMG3D5)
# Continuous integration
ENABLE_TESTING()

# To see falgs and options of compilation
#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

# Executable path
SET(EXECUTABLE_OUTPUT_PATH $ENV{HOME}/bin/$ENV{ARCHI}/)
SET(MMG_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../sources/)
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
    # Add flags to the compiler
    ADD_DEFINITIONS(-static-libgcc -mmacosx-version-min=10.4 -arch x86_64)

    # To avoid pbs with binary files...
    SET(CMAKE_EXE_LINKER_FLAGS "-arch x86_64")
    # Determine if the processor supports 64bit execution
    EXECUTE_PROCESS(
      COMMAND sysctl hw.cpu64bit_capable
      ERROR_QUIET
      OUTPUT_VARIABLE 64_CMD
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    STRING(REGEX REPLACE "^hw.cpu64bit_capable: (.*)" "\\1" 64_BIT "${64_CMD}")
  ELSEIF(UNIX)# UNIX must be after APPLE becaus APPLE is UNIX too
  ENDIF(APPLE)

ELSEIF(CMAKE_COMPILER_IS_INTEL)
  ADD_DEFINITIONS(-Bstatic)
ENDIF(CMAKE_COMPILER_IS_GNUCC)

IF(NOT CMAKE_COMPILER_IS_CLANG)
  # Compiler options for profiling... but not possible with clang
  SET(PROFILING CACHE BOOL "Enable/Disable PROFILING")
  IF(PROFILING)
    ADD_DEFINITIONS(-pg)
  ENDIF(PROFILING)
ENDIF(NOT CMAKE_COMPILER_IS_CLANG)



  # Add openMP flags if OpenMP found
SET(OPENMP CACHE BOOL "Enable/Disable OpenMP")
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
# /!\ ADD_CUSTOM_TARGET don't work with cmake-2.8.4 but is necessary for other versions
IF ( ${CMAKE_COMMAND} MATCHES 2.8.4 )
ELSE ()
    ADD_CUSTOM_TARGET(compil_date DEPENDS ${MMG_SOURCE_DIR}/compil.date)
ENDIF ( ${CMAKE_COMMAND} MATCHES 2.8.4 )

ADD_DEPENDENCIES(mmg3d5 compil_date)

# Libraries
FIND_LIBRARY(M_LIB m)
TARGET_LINK_LIBRARIES(mmg3d5 ${M_LIB})

# Continuous integration

OPTION(BUILD_TESTING "Enable/Disable continuous integration" OFF)
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
    #SET(CMAKE_C_FLAGS "-Wall -fprofile-arcs -ftest-coverage")
    SET(CMAKE_C_FLAGS "-Wall")
    #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -lgcov")
    #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
  ENDIF(CMAKE_COMPILER_IS_GNUCC)

  INCLUDE( ${CMAKE_SOURCE_DIR}/testing.cmake )
ENDIF(BUILD_TESTING)
