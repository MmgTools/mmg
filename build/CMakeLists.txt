CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (MMG3D5)     

# integration continue
ENABLE_TESTING()

SET(EXECUTABLE_OUTPUT_PATH $ENV{HOME}/bin/$ENV{ARCHI}/)
SET(CMAKELISTS_PATH .) 

MESSAGE(STATUS "EXECUTABLE SAVED IN :" ${EXECUTABLE_OUTPUT_PATH})

IF(APPLE) 
   ADD_DEFINITIONS(-static-libgcc -mmacosx-version-min=10.4 -arch x86_64) 
   #to avoid pbs with binary files... 
   SET(CMAKE_EXE_LINKER_FLAGS "-arch x86_64")
   # determine if the processor supports 64bit execution
   EXECUTE_PROCESS(
      COMMAND sysctl hw.cpu64bit_capable
      ERROR_QUIET
      OUTPUT_VARIABLE 64_CMD
      OUTPUT_STRIP_TRAILING_WHITESPACE
   )
   STRING(REGEX REPLACE "^hw.cpu64bit_capable: (.*)" "\\1" 64_BIT "${64_CMD}")
   
   # display the results
   MESSAGE(STATUS "CMAKE_OSX_ARCHITECTURES: " ${CMAKE_OSX_ARCHITECTURES})
ENDIF(APPLE)  

#file sources  
FILE(
	GLOB
	source_files
	${CMAKELISTS_PATH}/sources/*.c 
)

#compil.date
IF (NOT EXISTS ${CMAKELISTS_PATH}/sources/compil.date)
  SET_SOURCE_FILES_PROPERTIES( ${CMAKELISTS_PATH}/sources/compil.date PROPERTIES GENERATED TRUE ) 
  
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKELISTS_PATH}/sources/compil.date
    COMMAND ${CMAKE_COMMAND} -E echo \"\#define COMPIL \" '\"' `date` '\"' > ${CMAKELISTS_PATH}/sources/compil.date
    DEPENDS ${sources_files}
    COMMENT "Creation de compil.date"
    )
ELSEIF(EXISTS ${CMAKELISTS_PATH}/sources/compil.date) 
 ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKELISTS_PATH}/sources/compil.date
    COMMAND ${CMAKE_COMMAND} -E echo \"\#define COMPIL \" '\"' `date` '\"' > ${CMAKELISTS_PATH}/sources/compil.date
    DEPENDS ${source_files}
    COMMENT "Modification de compil.date"
    )
ENDIF(NOT EXISTS ${CMAKELISTS_PATH}/sources/compil.date)

# Executable et dependances
ADD_EXECUTABLE(mmg3d5 ${CMAKELISTS_PATH}/sources/compil.date ${source_files})
ADD_CUSTOM_TARGET(compil_date DEPENDS ${CMAKELISTS_PATH}/sources/compil.date)
ADD_DEPENDENCIES(mmg3d5 compil_date)

# Libraries
FIND_LIBRARY(M_LIB m)
TARGET_LINK_LIBRARIES(mmg3d5 ${M_LIB})

# Tests
INCLUDE(CTest)
SET(REG_TESTS_PATH ${CMAKELISTS_PATH}/../../RegTests)
ADD_TEST(SimpleCube $ENV{HOME}/bin/$ENV{ARCHI}/mmg3d5 ${REG_TESTS_PATH}/Cube/cube.mesh)