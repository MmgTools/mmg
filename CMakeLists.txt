## =============================================================================
##  This file is part of the Mmg software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
##
##  Mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  Mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with Mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the Mmg distribution only if you accept them.
## =============================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
PROJECT (mmg)

###############################################################################
#####
#####         Release version and date
#####
###############################################################################
SET (CMAKE_RELEASE_VERSION_MAJOR "5")
SET (CMAKE_RELEASE_VERSION_MINOR "0")
SET (CMAKE_RELEASE_VERSION_PATCH "0")
SET (CMAKE_RELEASE_DATE "Jul. 6, 2012")

SET (CMAKE_RELEASE_VERSION "${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}")

# Create mmg3d.h file with the good release infos.
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/mmg3d/mmg3d.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mmg3d/mmg3d.h @ONLY)

###############################################################################
#####
#####         Main CMake compilation variables
#####
###############################################################################

# To see flags and options of compilation
#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

# Executable path
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
SET(MMG_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/mmg3d)
SET(MMGS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/mmgs)
SET(COMMON_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/common)

# Find "exotic" compilers
IF (${CMAKE_C_COMPILER} MATCHES Clang OR ${CMAKE_C_COMPILER} MATCHES clang)
    # using clang
    SET(CMAKE_COMPILER_IS_CLANG TRUE)
ELSEIF(${CMAKE_C_COMPILER} MATCHES Icc OR ${CMAKE_C_COMPILER} MATCHES icc)
    # using icc
    SET(CMAKE_COMPILER_IS_INTEL TRUE)
ENDIF()

# Hide some options
MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES
    CMAKE_OSX_DEPLOYMENT_TARGET
    CMAKE_OSX_SYSROOT)

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
#####         Choose executable target to compile
#####
###############################################################################
SET ( CMAKE_BUILD_TYPE CACHE STRING "Release" )

SET ( BUILD "MMG" CACHE STRING
    "\"MMG3D\", \"MMGS\" or \"MMG\" to compile both MMGS and MMG3D" )

IF ( BUILD MATCHES "MMGS" )
    SET(BUILD_MMG        OFF )
    SET(BUILD_MMG3D      OFF )
    SET(BUILD_MMGS       ON  )
ELSEIF ( BUILD MATCHES "MMG3D" )
    SET(BUILD_MMG3D      ON  )
    SET(BUILD_MMGS       OFF )
    SET(BUILD_MMG        OFF )
ELSEIF ( BUILD MATCHES "MMG" )
    SET(BUILD_MMG        ON  )
    SET(BUILD_MMGS       OFF )
    SET(BUILD_MMG3D      OFF )
ELSE ( )
    MESSAGE (WARNING "Nothing to build")
ENDIF()

IF ( BUILD_MMG3D OR BUILD_MMG )

    ###############################################################################
    #####
    #####         Fortran header: libmmg3df.h
    #####
    ###############################################################################
    ADD_EXECUTABLE(genheader ${CMAKE_SOURCE_DIR}/scripts/genheader.c)

    ADD_CUSTOM_COMMAND(OUTPUT ${MMG_SOURCE_DIR}/libmmg3df.h
        COMMAND genheader ${MMG_SOURCE_DIR}/libmmg3df.h ${MMG_SOURCE_DIR}/libmmg3d.h ${CMAKE_SOURCE_DIR}/scripts/genfort.pl
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        DEPENDS genheader ${MMG_SOURCE_DIR}/libmmg3d.h ${CMAKE_SOURCE_DIR}/scripts/genfort.pl
        COMMENT "Generating Fortran header"
        )
    
    ###############################################################################
    #####
    #####         Choose type of insertion: delaunay kernel or patterns
    #####
    ###############################################################################
    # Should we use patterns for insertion istead of delaunay kernel
    OPTION ( PATTERN
        "If ON, insertion by patterns, otherwise insertion by delaunay kernel"
        OFF)
    IF ( NOT PATTERN )
        MESSAGE(STATUS "Vertex insertion by delaunay kernel")
    ELSE()
        MESSAGE(STATUS "Vertex insertion by patterns")
        SET(CMAKE_C_FLAGS "-DPATTERN ${CMAKE_C_FLAGS}")
    ENDIF()
    # Remove the next line to have this option visible in basic cmake mode
    MARK_AS_ADVANCED(PATTERN)

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
    ${COMMON_SOURCE_DIR}/*.c ${COMMON_SOURCE_DIR}/*.h
    )
LIST(REMOVE_ITEM source_files
    ${MMG_SOURCE_DIR}/mmg3d.c
    ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3d.c
    ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3df.c
    ${REMOVE_FILE})
FILE(
    GLOB
    main_file
    ${MMG_SOURCE_DIR}/mmg3d.c
    )
FILE(
    GLOB
    lib_file
    ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3d.c
    ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3df.c
    )
FILE(
    GLOB
    sourcemmgs_files
    ${MMGS_SOURCE_DIR}/*.c   ${MMGS_SOURCE_DIR}/*.h
    ${COMMON_SOURCE_DIR}/*.c ${COMMON_SOURCE_DIR}/*.h
    )

INCLUDE_DIRECTORIES(${COMMON_SOURCE_DIR})

IF ((NOT WIN32) AND (NOT WIN64))
    FIND_LIBRARY(M_LIB m)
    SET( LIBRARIES ${M_LIB} ${LIBRARIES})
ENDIF()

IF ( BUILD_MMG3D OR BUILD_MMG )
    ###############################################################################
    #####
    #####         Scotch
    #####
    ###############################################################################
    # add SCOTCH library?
    SET(SCOTCH_DIR "" CACHE PATH "Installation directory for scotch")
    INCLUDE(cmake/modules/FindScotch.cmake)

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
    #####    Set the full RPATH to find libraries independently from LD_LIBRARY_PATH
    #####
    ###############################################################################
    SET(CMAKE_MACOSX_RPATH TRUE)

    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # the RPATH to be used when installing, but only if it's not a system directory
    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    IF("${isSystemDir}" STREQUAL "-1")
        SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    ENDIF("${isSystemDir}" STREQUAL "-1")

    ###############################################################################
    #####
    #####         Compile mmg3d libraries
    #####
    ###############################################################################
    SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

    OPTION ( LIBMMG3D_STATIC "Compile static library" ON)
    IF ( LIBMMG3D_STATIC )
        ADD_LIBRARY(${PROJECT_NAME}3d_a  STATIC
            ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3df.h
            ${source_files} ${lib_file} )
        SET_TARGET_PROPERTIES(${PROJECT_NAME}3d_a PROPERTIES OUTPUT_NAME ${PROJECT_NAME}3d)
        TARGET_LINK_LIBRARIES(${PROJECT_NAME}3d_a ${LIBRARIES})
        INSTALL(FILES ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3df.h DESTINATION include)
        INSTALL(FILES ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3d.h DESTINATION include)
        INSTALL(FILES ${COMMON_SOURCE_DIR}/chrono.h DESTINATION include)
        INSTALL(TARGETS ${PROJECT_NAME}3d_a
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
    ENDIF()

    # Compile shared library
    OPTION ( LIBMMG3D_SHARED "Compile dynamic library" OFF)
    IF ( LIBMMG3D_SHARED )
        ADD_LIBRARY(${PROJECT_NAME}3d_so SHARED
            ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3df.h
            ${source_files} ${lib_file})
        SET_TARGET_PROPERTIES(${PROJECT_NAME}3d_so PROPERTIES OUTPUT_NAME ${PROJECT_NAME}3d)
        TARGET_LINK_LIBRARIES(${PROJECT_NAME}3d_so ${LIBRARIES})
        INSTALL(FILES ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3df.h DESTINATION include)
        INSTALL(FILES ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3d.h DESTINATION include)
        INSTALL(FILES ${COMMON_SOURCE_DIR}/chrono.h DESTINATION include)
        INSTALL(TARGETS ${PROJECT_NAME}3d_so
            ARCHIVE DESTINATION
            LIBRARY DESTINATION lib)
    ENDIF()

    ###############################################################################
    #####
    #####         Compile program to test library
    #####
    ###############################################################################
    OPTION ( TEST_LIBMMG3D "Compile tests of library" OFF)
    IF ( TEST_LIBMMG3D )
        IF ( LIBMMG3D_STATIC )
            ADD_EXECUTABLE(libmmg_example0_a ${CMAKE_SOURCE_DIR}/libexamples/example0/example0_a/main.c)
            TARGET_LINK_LIBRARIES(libmmg_example0_a ${PROJECT_NAME}3d_a)
            INSTALL(TARGETS libmmg_example0_a RUNTIME DESTINATION bin )

            ADD_EXECUTABLE(libmmg_example0_b ${CMAKE_SOURCE_DIR}/libexamples/example0/example0_b/main.c)
            TARGET_LINK_LIBRARIES(libmmg_example0_b ${PROJECT_NAME}3d_a)
            INSTALL(TARGETS libmmg_example0_b RUNTIME DESTINATION bin )

            ADD_EXECUTABLE(libmmg_example1 ${CMAKE_SOURCE_DIR}/libexamples/example1/main.c)
            TARGET_LINK_LIBRARIES(libmmg_example1 ${PROJECT_NAME}3d_a)
            INSTALL(TARGETS libmmg_example1 RUNTIME DESTINATION bin )
            ADD_EXECUTABLE(libmmg_example2 ${CMAKE_SOURCE_DIR}/libexamples/example2/main.c)
            TARGET_LINK_LIBRARIES(libmmg_example2 ${PROJECT_NAME}3d_a)
            INSTALL(TARGETS libmmg_example2 RUNTIME DESTINATION bin )

            ENABLE_LANGUAGE (Fortran)
            ADD_EXECUTABLE(libmmg_fortran_a ${CMAKE_SOURCE_DIR}/libexamples/example0_fortran/example0_a/main.F90)
            TARGET_LINK_LIBRARIES(libmmg_fortran_a ${PROJECT_NAME}3d_a)
            INSTALL(TARGETS libmmg_fortran_a RUNTIME DESTINATION bin )
            ADD_EXECUTABLE(libmmg_fortran_b ${CMAKE_SOURCE_DIR}/libexamples/example0_fortran/example0_b/main.F90)
            TARGET_LINK_LIBRARIES(libmmg_fortran_b ${PROJECT_NAME}3d_a)
            INSTALL(TARGETS libmmg_fortran_b RUNTIME DESTINATION bin )

            ADD_EXECUTABLE(libmmg_example3  ${CMAKE_SOURCE_DIR}/libexamples/example3/mmg3d.c)
            TARGET_LINK_LIBRARIES(libmmg_example3 ${PROJECT_NAME}3d_a)
            INSTALL(TARGETS libmmg_example3 RUNTIME DESTINATION bin )
        ELSEIF ( LIBMMG3D_SHARED )
            ADD_EXECUTABLE(libmmg_example0_a ${CMAKE_SOURCE_DIR}/libexamples/example0/example0_a/main.c)
            TARGET_LINK_LIBRARIES(libmmg_example0_a ${PROJECT_NAME}3d_so)
            INSTALL(TARGETS libmmg_example0_a RUNTIME DESTINATION bin )

            ADD_EXECUTABLE(libmmg_example0_b ${CMAKE_SOURCE_DIR}/libexamples/example0/example0_b/main.c)
            TARGET_LINK_LIBRARIES(libmmg_example0_b ${PROJECT_NAME}3d_so)
            INSTALL(TARGETS libmmg_example0_b RUNTIME DESTINATION bin )

            ADD_EXECUTABLE(libmmg_example1 ${CMAKE_SOURCE_DIR}/libexamples/example1/main.c)
            TARGET_LINK_LIBRARIES(libmmg_example1 ${PROJECT_NAME}3d_so)
            INSTALL(TARGETS libmmg_example1 RUNTIME DESTINATION bin )
            ADD_EXECUTABLE(libmmg_example2  ${CMAKE_SOURCE_DIR}/libexamples/example2/main.c)
            TARGET_LINK_LIBRARIES(libmmg_example2 ${PROJECT_NAME}3d_so)
            INSTALL(TARGETS libmmg_example2 RUNTIME DESTINATION bin )

            ADD_EXECUTABLE(libmmg_example3  ${CMAKE_SOURCE_DIR}/libexamples/example3/mmg3d.c)
            TARGET_LINK_LIBRARIES(libmmg_example3 ${PROJECT_NAME}3d_so)
            INSTALL(TARGETS libmmg_example3 RUNTIME DESTINATION bin )
        ELSE ()
            MESSAGE(WARNING "You must activate the compilation of the static or shared ${PROJECT_NAME} library to compile this tests." )
        ENDIF ()
    ENDIF()

    ###############################################################################
    #####
    #####         Compile MMG3D executable
    #####
    ###############################################################################
    ADD_EXECUTABLE(${PROJECT_NAME}3d
        ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}3df.h
        ${source_files} ${main_file})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME}3d ${LIBRARIES})
    INSTALL(TARGETS ${PROJECT_NAME}3d RUNTIME DESTINATION bin)


    SET_TARGET_PROPERTIES(${PROJECT_NAME}3d PROPERTIES DEBUG_POSTFIX _debug) # in debug mode we name the executable mmg3d_debug
    SET_TARGET_PROPERTIES(${PROJECT_NAME}3d PROPERTIES RELEASE_POSTFIX _O3) # in Release mode we name the executable mmg3d_O3
    SET_TARGET_PROPERTIES(${PROJECT_NAME}3d PROPERTIES RELWITHDEBINFO_POSTFIX _O3d) # in Release mode we name the executable mmg3d_O3
    SET_TARGET_PROPERTIES(${PROJECT_NAME}3d PROPERTIES MINSIZEREL_POSTFIX _Os) # in MinSizeRel mode we name the executable mmg3d_Os

ENDIF()

###############################################################################
#####
#####         Compile MMGS executable
#####
###############################################################################
IF ( BUILD_MMGS OR BUILD_MMG )
    ADD_EXECUTABLE(${PROJECT_NAME}s    ${sourcemmgs_files} )
    TARGET_LINK_LIBRARIES(${PROJECT_NAME}s ${LIBRARIES})
    INSTALL(TARGETS ${PROJECT_NAME}s RUNTIME DESTINATION bin)


    SET_TARGET_PROPERTIES(${PROJECT_NAME}s PROPERTIES DEBUG_POSTFIX _debug) # in debug mode we name the executable mmgs_debug
    SET_TARGET_PROPERTIES(${PROJECT_NAME}s PROPERTIES RELEASE_POSTFIX _O3) # in Release mode we name the executable mmgs_O3
    SET_TARGET_PROPERTIES(${PROJECT_NAME}s PROPERTIES MINSIZEREL_POSTFIX _Os)
ENDIF()

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################
IF ( BUILD_MMG3D OR BUILD_MMG )

    OPTION ( CONTINUOUS_INTEGRATION "Enable/Disable continuous integration" OFF )

    IF(CONTINUOUS_INTEGRATION)
        ##-------------------------------------------------------------------##
        ##------- Get the continuous integration tests ----------------------##
        ##-------------------------------------------------------------------##
        # SET ( GET_MMGS_TESTS "FALSE" )
        SET ( GET_MMG3D_TESTS "FALSE" )
        # Check if the ci_tests directory exists
        IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests )
            # First download of the tests
            MESSAGE("-- Create the continuous integration directory (ci_tests)"
                " and download the test cases. May be very long...")
            FILE(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests)
            # SET ( GET_MMGS_TESTS "TRUE" )
            SET ( GET_MMG3D_TESTS "TRUE" )
        ELSE ( )
            # Check if the tests are up to date
            #IF ( EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version )
            #    FILE(READ ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version OLD_MMGS_VERSION)
            #ELSE ( )
            #    SET ( OLD_MMGS_VERSION "0" )
            #ENDIF ( )
            # File(DOWNLOAD
            #    ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version)
            # FILE(READ ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version MMGS_VERSION)
            # IF ( NOT (${OLD_MMGS_VERSION} MATCHES ${MMGS_VERSION})
            #     MESSAGE("-- Continuous integration tests cases need to be updated.")
            #     SET ( GET_MMGS_TESTS "TRUE" )
            # ENDIF ()

            IF ( EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version )
                FILE(READ ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version OLD_MMG3D_VERSION)
            ELSE ( )
                SET ( OLD_MMG3D_VERSION "0" )
            ENDIF ( )

            FILE(DOWNLOAD https://www.dropbox.com/s/s723rs0xqu2j6k8/mmg3d.version?dl=0
                ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version)
            FILE(READ ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version MMG3D_VERSION)

            IF (NOT (${OLD_MMG3D_VERSION} MATCHES ${MMG3D_VERSION}))
                MESSAGE("-- Continuous integration tests cases need to be updated.")
                SET ( GET_MMG3D_TESTS "TRUE" )
            ENDIF ()
        ENDIF()
        # Download the tests if needed
        # IF ( GET_MMGS_TESTS MATCHES "TRUE" )
        #     MESSAGE("-- Download the mmgs test cases. May be very long...")
        #     FILE(DOWNLOAD
        #         ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz
        #         SHOW_PROGRESS)
        #     IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz )
        #         MESSAGE("\n")
        #         MESSAGE(WARNING "Fail to automatically download the mmgs test cases "
        #             " Try to get it at the following link:"
        #
        #             " then untar it in the project directory (Mmg/ by default).")
        #     ELSE()
        #         EXECUTE_PROCESS(
        #             COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz
        #             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests/
        #             )
        #         IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests.tgz )
        #             MESSAGE("\n")
        #             MESSAGE(WARNING "Fail to automatically untar the mmgs "
        #                 "test cases directory (mmgs.tgz)."
        #                 "Try to untar it by hand in the project directory"
        #                 " (Mmg/ci_tests/ by default).")
        #         ENDIF()
        #         FILE(REMOVE ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz)
        #     ENDIF ()
        # ENDIF ()
        IF ( GET_MMG3D_TESTS MATCHES "TRUE" )
            MESSAGE("-- Download the mmg3d test cases. May be very long...")
            FILE(DOWNLOAD https://www.dropbox.com/s/e9yuc46617jynku/mmg3d.tgz?dl=0
                ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz
                SHOW_PROGRESS)
            IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz )
                MESSAGE("\n")
                MESSAGE(WARNING "Fail to automatically download the mmg3d test cases "
                    " Try to get it at the following link:"
                    " https://www.dropbox.com/s/e9yuc46617jynku/mmg3d.tgz?dl=0"
                    " then untar it in the project directory (Mmg/ by default).")
            ELSE()
                EXECUTE_PROCESS(
                    COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests/
                    )
                IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz )
                    MESSAGE("\n")
                    MESSAGE(WARNING "Fail to automatically untar the mmg3d "
                        "test cases directory (mmg3d.tgz)."
                        "Try to untar it by hand in the project directory"
                        " (Mmg/ci_tests/ by default).")
                ENDIF()

                FILE(REMOVE ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz)
            ENDIF()
        ENDIF()


        ##-------------------------------------------------------------------##
        ##------- Set the continuous integration options --------------------##
        ##-------------------------------------------------------------------##
        SET(MMGS_CI_TESTS ${CMAKE_SOURCE_DIR}/ci_tests/mmgs )
        SET(MMG3D_CI_TESTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d )

        # Run some tests twice with the output of the previous test as input
        OPTION ( RUN_AGAIN "Enable/Disable second run of some tests" ON )
        MARK_AS_ADVANCED(RUN_AGAIN)

        # Remove the --show-reachable=yes option from the default valgrind command
        set(VALGRIND_COMMAND_OPTIONS "-q --tool=memcheck --leak-check=yes --workaround-gcc296-bugs=yes --num-callers=100")
        # Ignore errors coming from the libmesh library
        IF(APPLE)
            # and errors from the libSystem.B library
            SET(MEMORYCHECK_SUPPRESSIONS_FILE "${CMAKE_SOURCE_DIR}/cmake/Supp/libSystem.supp")
        ENDIF(APPLE)
        # Copy config file and custom file in the good directories
        CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/config/CTestConfig.cmake ${CMAKE_SOURCE_DIR}/CTestConfig.cmake)
        CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/config/CTestCustom.cmake ./CTestCustom.cmake)

        ENABLE_TESTING()
        INCLUDE(CTest)

        # Max time for test case
        SET(DART_TESTING_TIMEOUT 2500)

        IF ( NOT APPLE)
            # Add coverage flags
            SET(CMAKE_C_FLAGS_DEBUG " -fprofile-arcs -ftest-coverage --coverage ${CMAKE_C_FLAGS_DEBUG}")
            SET(CMAKE_EXE_LINKER_FLAGS_DEBUG " -fprofile-arcs -ftest-coverage --coverage ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
        ENDIF()

        IF ( CMAKE_COMPILER_IS_GNUCC )
            SET(CMAKE_C_FLAGS_DEBUG "-g -Wall ${CMAKE_C_FLAGS_DEBUG}")
            SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "-g ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
        ELSEIF( CMAKE_COMPILER_IS_CLANG )
            SET(CMAKE_C_FLAGS_DEBUG "-g -Weverything -Wno-sign-conversion -Wno-char-subscripts -Wno-padded ${CMAKE_C_FLAGS_DEBUG}")
            SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "-g ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
        ENDIF ( )

        ##-------------------------------------------------------------------##
        ##--------------------------- Add tests and configure it ------------##
        ##-------------------------------------------------------------------##
        # Add runtime that we want to test
        IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
            SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_NAME}3d_debug)
            SET(BUILDNAME ${BUILDNAME}_debug CACHE STRING "build name variable")
        ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Release")
            SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_NAME}3d_O3)
            SET(BUILDNAME ${BUILDNAME}_O3 CACHE STRING "build name variable")
        ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Minsizerel")
            SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_NAME}3d_Os)
            SET(BUILDNAME ${BUILDNAME}_Os CACHE STRING "build name variable")
        ELSE()
            SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_NAME}3d)
            SET(BUILDNAME ${BUILDNAME} CACHE STRING "build name variable")
        ENDIF()

        SET ( LISTEXEC ${EXECUT} )
        IF ( TEST_LIBMMG3D )
            IF ( LIBMMG3D_STATIC OR LIBMMG3D_SHARED )
                SET(LIBEXEC0_a ${EXECUTABLE_OUTPUT_PATH}/libmmg_example0_a)
                SET(LIBEXEC0_b ${EXECUTABLE_OUTPUT_PATH}/libmmg_example0_b)
                SET(LIBEXEC1   ${EXECUTABLE_OUTPUT_PATH}/libmmg_example1)
                SET(LIBEXEC2   ${EXECUTABLE_OUTPUT_PATH}/libmmg_example2)

                ADD_TEST(NAME libexample0_a COMMAND ${LIBEXEC0_a})
                ADD_TEST(NAME libexample0_b COMMAND ${LIBEXEC0_b})
                ADD_TEST(NAME libexample1   COMMAND ${LIBEXEC1})
                ADD_TEST(NAME libexample2   COMMAND ${LIBEXEC2})

                SET(LIBEXEC3 ${EXECUTABLE_OUTPUT_PATH}/libmmg_example3)
                SET( LISTEXEC ${LISTEXEC} ${LIBEXEC3} )
            ENDIF()
        ENDIF()

        # Add more tests
        INCLUDE( ${CMAKE_SOURCE_DIR}/cmake/testing/mmg3d_tests.cmake )

        IF ( RUN_AGAIN )
            INCLUDE( ${CMAKE_SOURCE_DIR}/cmake/testing/mmg3d_rerun_tests.cmake )
        ENDIF()

        # Change authorisations for some files...
        IF ( EXISTS ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd7/unreadable.meshb)
            FILE(REMOVE ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd7/unwrittable.meshb)
        ENDIF()
        FILE(COPY ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd7/dout.mesh
            DESTINATION ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd7/unwrittable.meshb)

        IF ( EXISTS ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd8/unwrittable.sol)
            FILE(REMOVE ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd8/unwrittable.sol)
        ENDIF()
        FILE(COPY ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd8/dsol.sol
            DESTINATION ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd8/unwrittable.sol)

    ENDIF(CONTINUOUS_INTEGRATION)
ENDIF()

###############################################################################
#####
#####         Create API Documentation
#####
###############################################################################
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
    # MMGS Documentation
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/Doxyfile.in ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/Doxyfile @ONLY)
    ADD_CUSTOM_TARGET(mmgs_doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/Doxyfile
        # COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/latex make
        # COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/latex/refman.pdf ${CMAKE_CURRENT_SOURCE_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs
        COMMENT "Generating MMGS API documentation with Doxygen.
    Open up the doc/doxygen/mmgs/html/index.html file to see it." VERBATIM
        )
    # MMG3D Documentation
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/Doxyfile.in ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/Doxyfile @ONLY)
    ADD_CUSTOM_TARGET(mmg3d_doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/Doxyfile
        # COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/latex make
        # COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/latex/refman.pdf ${CMAKE_CURRENT_SOURCE_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d
        COMMENT "Generating MMG3D API documentation with Doxygen.
     Open up the doc/doxygen/mmg3d/html/index.html file to see it" VERBATIM
        )
    # MMGS and MMG3D Documentation
    ADD_CUSTOM_TARGET(doc
        DEPENDS mmgs_doc mmg3d_doc
        COMMENT "Generating MMGS and MMG3D API documentation with Doxygen.
    To see the MMGS documentation, open up the doc/doxygen/mmgs/html/index.html file.
    To see the MMG3D documentation, open up the doc/doxygen/mmg3d/html/index.html file."
        VERBATIM
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
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MMG3D/MMGS: ISOTROPIC TETRAHEDRAL REMESHER/SURFACE REMESHER")
SET(CPACK_PACKAGE_VENDOR "Cécile Dobrzynski, Pascal Frey, Charles Dapogny, Algiane Froehly")
SET(CPACK_PACKAGE_CONTACT "cecile.dobrzynski@inria.fr")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_VERSION ${CMAKE_RELEASE_VERSION})
SET(CPACK_PACKAGE_VERSION_MAJOR ${CMAKE_RELEASE_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${CMAKE_RELEASE_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${CMAKE_RELEASE_VERSION_PATCH})
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}")
SET(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_INSTALL_DIRECTORY}-${CMAKE_SYSTEM})
SET(CPACK_PACKAGE_EXECUTABLES "${PROJECT_NAME}")
SET(CPACK_OUTPUT_FILE_PREFIX "../archives")

SET(CPACK_IGNORE_FILES ";.\\\\\\\\.#;/\\\\\\\\.git;.~;.#;.-e;.DS_Store" ${CPACK_IGNORE_FILES})
SET(CPACK_IGNORE_FILES ";libmmg.*" ${CPACK_IGNORE_FILES})

SET(CPACK_SOURCE_IGNORE_FILES ";src/mmg3d/mmg3d.h$" ${CPACK_SOURCE_IGNORE_FILE})
SET(CPACK_SOURCE_IGNORE_FILES ";.\\\\\\\\.#;/\\\\\\\\.git;.~;.#;.-e;.DS_Store" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";lib/;bin;archives;def" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";.dir-locals.el;.tm_properties" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/example0/example0_a/mesh.o.*" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/example0/example0_b/mesh.o.*" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/example1/result0.*" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/example2/result1.*" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/example3/SphereIso.*" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";CMakeFiles" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";cmake/supp/libSystem.supp" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";CTestConfig.cmake" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";cmake/testing/testing.cmake" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";.#$"      ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";build*" ${CPACK_SOURCE_IGNORE_FILES}  )
SET(CPACK_SOURCE_IGNORE_FILES ";doc/doxygen/.*/html;doc/doxygen/.*/Doxyfile$"  ${CPACK_SOURCE_IGNORE_FILES}  )

IF ((WIN32 OR WIN64) AND NOT UNIX)
    # Warning compilation under windows not tested
    #SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/logo.pdf")
    SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\${PROJECT_NAME}.exe")
    SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} Mmg: surface and volume remeshers")
    SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.mmgtools.org‎")
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
