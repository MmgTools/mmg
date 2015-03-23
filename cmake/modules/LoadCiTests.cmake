## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
##
##  mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the mmg distribution only if you accept them.
## =============================================================================

SET ( GET_MMG_TESTS "FALSE" )
SET ( GET_MMGS_TESTS "FALSE" )
SET ( GET_MMG3D_TESTS "FALSE" )

# Check if the ci_tests directory exists
IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests )

  # First download of the tests
  MESSAGE("-- Create the continuous integration directory (ci_tests)"
    " and download the test cases. May be very long...")
  FILE(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests)

  SET ( GET_MMG_TESTS "TRUE" )
  SET ( GET_MMGS_TESTS "TRUE" )
  SET ( GET_MMG3D_TESTS "TRUE" )

ELSE ( )

  # Check if the tests are up to date
  #--------------> mmg
  IF ( EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg.version )
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmg.version OLD_MMG_MD5)
  ELSE ( )
    SET ( OLD_MMG_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://www.dropbox.com/s/yevltsxy3hxyv5u/mmg.version?dl=0
    ${CMAKE_SOURCE_DIR}/ci_tests/mmg.version
    STATUS MMG_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMG_VERSION_STATUS 0 MMG_VERSION_STATUS_0)
  LIST(GET MMG_VERSION_STATUS 1 MMG_VERSION_STATUS_1)

  IF ( MMG_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmg.version MMG_MD5)

    IF ( NOT (${OLD_MMG_MD5} MATCHES ${MMG_MD5}) )
      SET ( GET_MMG_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMG_VERSION_STATUS_1}.
 Try to get it at the following link:
        https://www.dropbox.com/s/yevltsxy3hxyv5u/mmg.version?dl=0
 then untar it in the project directory (mmg/ by default).")
  ENDIF()

    #--------------> mmgs
  IF ( EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version )
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version OLD_MMGS_MD5)
  ELSE ( )
    SET ( OLD_MMGS_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://www.dropbox.com/s/j7om7gl7r5mlnvj/mmgs.version?dl=0
    ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version
    STATUS MMGS_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMGS_VERSION_STATUS 0 MMGS_VERSION_STATUS_0)
  LIST(GET MMGS_VERSION_STATUS 1 MMGS_VERSION_STATUS_1)

  IF ( MMGS_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version MMGS_MD5)

    IF ( NOT (${OLD_MMGS_MD5} MATCHES ${MMGS_MD5}) )
      SET ( GET_MMGS_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMGS_VERSION_STATUS_1}.
 Try to get it at the following link:
        https://www.dropbox.com/s/zusur5gzz1tvmvh/mmgs.tgz?dl=0
 then untar it in the project directory (mmg/ by default).")
  ENDIF()

  #--------------> mmg3d
  IF ( EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version )
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version OLD_MMG3D_MD5)
  ELSE ( )
    SET ( OLD_MMG3D_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://www.dropbox.com/s/s723rs0xqu2j6k8/mmg3d.version?dl=0
    ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version
    STATUS MMG3D_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMG3D_VERSION_STATUS 0 MMG3D_VERSION_STATUS_0)
  LIST(GET MMG3D_VERSION_STATUS 1 MMG3D_VERSION_STATUS_1)

  IF ( MMG3D_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version MMG3D_MD5)

    IF (NOT (${OLD_MMG3D_MD5} MATCHES ${MMG3D_MD5}))
      SET ( GET_MMG3D_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMG3D_VERSION_STATUS_1}.
 Try to get it at the following link:
        https://www.dropbox.com/s/e9yuc46617jynku/mmg3d.tgz?dl=0
 then untar it in the project directory (mmg/ by default).")
  ENDIF()

ENDIF()

# Download the tests if needed
#--------------> mmg
IF ( GET_MMG_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmg test cases. May be very long...")
  FILE(DOWNLOAD https://www.dropbox.com/s/hnodvsv56mdazyx/mmg.tgz?dl=0
    ${CMAKE_SOURCE_DIR}/ci_tests/mmg.tgz
    SHOW_PROGRESS)
  IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg.tgz )
    MESSAGE("\n")
    MESSAGE(WARNING "Fail to automatically download the mmg test cases.
Try to get it at the following link:
       https://www.dropbox.com/s/hnodvsv56mdazyx/mmg.tgz?dl=0
then untar it in the project directory (mmg/ by default).")
  ELSE()
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar xzf
      ${CMAKE_SOURCE_DIR}/ci_tests/mmg.tgz
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests/
      )
    IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg.tgz )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmg "
        "test cases directory (mmg.tgz).
Try to untar it by hand in the project directory"
        " (mmg/ci_tests/ by default).")
    ENDIF()
    FILE(REMOVE ${CMAKE_SOURCE_DIR}/ci_tests/mmg.tgz)
  ENDIF ()
ENDIF ()

#--------------> mmgs
IF ( GET_MMGS_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmgs test cases. May be very long...")
  FILE(DOWNLOAD https://www.dropbox.com/s/zusur5gzz1tvmvh/mmgs.tgz?dl=0
    ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz
    SHOW_PROGRESS)
  IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz )
    MESSAGE("\n")
    MESSAGE(WARNING "Fail to automatically download the mmgs test cases.
Try to get it at the following link:
       https://www.dropbox.com/s/zusur5gzz1tvmvh/mmgs.tgz?dl=0
then untar it in the project directory (mmg/ by default).")
  ELSE()
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar xzf
      ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests/
      )
    IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmgs "
        "test cases directory (mmgs.tgz).
Try to untar it by hand in the project directory"
        " (mmg/ci_tests/ by default).")
    ENDIF()
    FILE(REMOVE ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz)
  ENDIF ()
ENDIF ()

#--------------> mmg3d
IF ( GET_MMG3D_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmg3d test cases. May be very long...")
  FILE(DOWNLOAD https://www.dropbox.com/s/e9yuc46617jynku/mmg3d.tgz?dl=0
    ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz
    SHOW_PROGRESS)
  IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz )
    MESSAGE("\n")
    MESSAGE(WARNING "Fail to automatically download the mmg3d test cases
Try to get it at the following link:
       https://www.dropbox.com/s/e9yuc46617jynku/mmg3d.tgz?dl=0
then untar it in the project directory mmg/ by default).")
  ELSE()
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar xzf
      ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests/
      )
    IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmg3d"
        "test cases directory (mmg3d.tgz).
Try to untar it by hand in the project directory"
        " (mmg/ci_tests/ by default).")
    ENDIF()

    FILE(REMOVE ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz)
  ENDIF()
ENDIF()
