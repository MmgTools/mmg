!> @author Cecile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief Example of input output for the mmgs library for multiple solutions
!> at mesh vertices

PROGRAM main

  IMPLICIT NONE

  !> Include here the mmgs library hader file
  ! if the header file is in the "include" directory
  ! #include "libmmgsf.h"

  ! if the header file is in "include/mmg/mmgs"
#include "mmg/mmgs/libmmgsf.h"

  MMG5_DATA_PTR_T    :: mmgMesh
  MMG5_DATA_PTR_T    :: mmgSol,tmpSol
  INTEGER            :: ier,argc,i

  !! To manually recover the mesh
  INTEGER            :: nsol,np,typEntity(MMG5_NSOL_MAX),typSol(MMG5_NSOL_MAX)
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: sols

  CHARACTER(len=300) :: exec_name,filename,fileout

  PRINT*,"  -- TEST MMGSLIB"

  argc =  COMMAND_ARGUMENT_COUNT();
  CALL get_command_argument(0, exec_name)


  IF ( argc /=2 ) THEN
     PRINT*," Usage: ",TRIM(ADJUSTL(exec_name))," input_file_name output_file_name"
     CALL EXIT(1);
  ENDIF

  ! Name and path of the mesh file
  CALL get_command_argument(1, filename)
  CALL get_command_argument(2, fileout)


  !!> ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures */
  !! args of InitMesh:
  !! MMG5_ARG_start: we start to give the args of a variadic func
  !! MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
  !! &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
  !! MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
  !! &mmgSol: pointer toward your MMG5_pSol (that store your metric)

  mmgMesh = 0
  mmgSol  = 0
  tmpSol  = 0

  CALL MMGS_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end);


  !!> 2) Build initial mesh and solutions in MMG5 format
  !! Two solutions: just use the MMGS_loadMesh function that will read a .mesh(b)
  !! file formatted or manually set your mesh using the MMGS_Set* functions

  !!> Automatic loading of the mesh and multiple solutions
  CALL MMGS_loadMesh(mmgMesh,TRIM(ADJUSTL(filename)),&
       LEN(TRIM(ADJUSTL(filename))),ier)
  IF ( ier /= 1 )  CALL EXIT(102)

  CALL MMGS_loadAllSols(mmgMesh,mmgSol,TRIM(ADJUSTL(filename)),&
       LEN(TRIM(ADJUSTL(filename))),ier)
  IF ( ier /= 1 )  CALL EXIT(103)

  !!> ------------------------------ STEP II ---------------------------

  !!> 3) Transfer the solutions in a new solutions array
  !! a) Get the solutions sizes
  CALL MMGS_Get_allSolsSizes(mmgMesh,mmgSol,nsol,typEntity,np,typSol,ier)
  IF ( ier /= 1 )  CALL EXIT(104)

  !!> b) Manually set the size of the new solution: give info for the sol
  !! structure: number of solutions, type of entities on which applied the
  !! solutions, number of vertices, type of the solution */
  CALL MMGS_Set_allSolsSizes(mmgMesh,tmpSol,nsol,typEntity,np,typSol,ier)
  IF ( ier /= 1 )  CALL EXIT(105)

  !!> c) Get each solution and set it in the new structure

  !!> b) give solutions values and positions
  !! Get the entire field of a given solution
  DO i=1,nsol
    IF ( typEntity(i) /= MMG5_Vertex ) CALL EXIT(106)

    ! Get the ith solution array
    IF ( typSol(i) == MMG5_Scalar ) THEN
       ALLOCATE(sols(np))
    ELSE IF ( typSol(i) == MMG5_Vector ) THEN
       ALLOCATE(sols(3*np))
    ELSE IF ( typSol(i) == MMG5_Tensor ) THEN
       ALLOCATE(sols(6*np))
    ENDIF

    CALL MMGS_Get_ithSols_inAllSols(mmgSol,i,sols,ier)
    IF ( ier /= 1 )  CALL EXIT(107)

    ! Set the ith solution in the new structure
    CALL MMGS_Set_ithSols_inAllSols(tmpSol,i,sols,ier)
    IF ( ier /= 1 )  CALL EXIT(108)

    DEALLOCATE(sols)
  ENDDO


  !!> ------------------------------ STEP III --------------------------
  !! Save the new data
  !! Use the MMGS_saveMesh/MMGS_saveAllSols functions
  !! save the mesh
  !> 1) Automatically save the mesh
  CALL MMGS_saveMesh(mmgMesh,TRIM(ADJUSTL(fileout)),LEN(TRIM(ADJUSTL(fileout))),ier)
  IF ( ier /= 1 ) CALL EXIT(110)

  ! save the solutions array
  CALL MMGS_saveAllSols(mmgMesh,tmpSol,TRIM(ADJUSTL(fileout)), &
       LEN(TRIM(ADJUSTL(fileout))),ier)
  IF ( ier /= 1 ) CALL EXIT(111)

  !!> 3) Free the MMGS structures
  CALL MMGS_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppSols,tmpSol, &
       MMG5_ARG_ppSols,mmgSol, &
       MMG5_ARG_end)

END PROGRAM main
