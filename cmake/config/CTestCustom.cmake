## =============================================================================
##  This file is part of the MMG3D 5 software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
##
##  MMG3D 5 is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  MMG3D 5 is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with MMG3D 5 (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the MMG3D 5 distribution only if you accept them.
## =============================================================================

SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 5000)
SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 2000)

SET ( CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
  ".*: warning: array subscript has type 'char'.*")
SET ( CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
  ".*: warning: array subscript has type ‘char’.*")
SET (CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 0 )
SET (CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 0 )
