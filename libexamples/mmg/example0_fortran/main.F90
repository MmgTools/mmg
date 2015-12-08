!> @author
!> Cecile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief
!>  Example for using mmglib (basic use)

!> Include the mmg library hader file
! if the "include/mmg" dir is in your include path
!#include "libmmg.h"

! if your include path do not contain the "mmg/mmg" subdirectories
! #include "mmglibmmg.h"
#include "mmg/libmmgf.h"

PROGRAM main

  IMPLICIT NONE

  MMG5_DATA_PTR_T    :: mmgMesh
  MMG5_DATA_PTR_T    :: mmgSol
  INTEGER            :: ier
  CHARACTER(len=255) :: pwd
  CHARACTER(len=300) :: tdfile,tdfileout, bdfile, bdfileout

  WRITE(*,*) "  -- TEST MMGLIB"
  ! get path of the run
  CALL getenv("PWD",pwd)

  !> ================== 2d remeshing using the mmg2d library
  !! ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures */
  !! args of InitMesh: mesh=&mmgMesh, sol=&mmgSol, input mesh name, input sol
  !! name, output mesh name */
  WRITE(bdfile,*) TRIM(pwd),"/../libexamples/mmg/example0_fortran/init"
  WRITE(bdfileout,*) TRIM(pwd),"/../libexamples/mmg/example0_fortran/result.mesh"

  mmgMesh = 0
  mmgSol  = 0
  !! Remark: %val(0) allow to pass the value 0 (i.e. NULL) instead of a pointer
  !! toward NULL.
  CALL MMG2D_Init_mesh(mmgMesh,mmgSol,%val(0))

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
  !!    file formatted or manually set your mesh using the MMG2D_Set* functions

  !> with MMG2D_loadMesh function
  !! a) (not mandatory): give the mesh name
  !!   (by default, the "mesh.mesh" file is oppened)
  !CALL MMG2D_Set_inputMeshName(mmgMesh,TRIM(ADJUSTL(bdfile)),&
  !     LEN(TRIM(ADJUSTL(bdfile))),ier)
  !IF ( ier /= 1 ) THEN
  !   CALL EXIT(101)
  !ENDIF

  !> b) function calling
  CALL MMG2D_loadMesh(mmgMesh,TRIM(ADJUSTL(bdfile)),&
       LEN(TRIM(ADJUSTL(bdfile))),ier)
  IF ( ier /= 1 )  CALL EXIT(102)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMG2D_loadSol function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMG2D_Set* functions
  CALL MMG2D_loadSol(mmgMesh,mmgSol,TRIM(ADJUSTL(bdfile)),LEN(TRIM(ADJUSTL(bdfile))),0,ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(104)
  ENDIF

  !> 4) (not mandatory): check if the number of given entities match with mesh size

  !> ------------------------------ STEP  II --------------------------
  !! library call
  CALL MMG2D_mmg2dlib(mmgMesh,mmgSol,ier)

  IF ( ier == MMG5_STRONGFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH"
     STOP 2
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG2DLIB"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results
  !! Two solutions: just use the MMG2D_saveMesh/MMG2D_saveSol functions
  !!    that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !!    using the MMG2D_getMesh/MMG2D_getSol functions

  !> 1) Automatically save the mesh
  !! a)  (not mandatory): give the ouptut mesh name using MMG2D_Set_outputMeshName
  !!   (by default, the mesh is saved in the "mesh.o.mesh" file
  !!call MMG2D_Set_outputMeshName(mmgMesh,"output.mesh",len("output.mesh"),ier)
  CALL MMG2D_Set_outputMeshName(mmgMesh,TRIM(ADJUSTL(bdfileout)),&
       LEN(TRIM(ADJUSTL(bdfile))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(105)
  ENDIF

  !! b) function calling
  CALL MMG2D_saveMesh(mmgMesh,TRIM(ADJUSTL(bdfileout)),&
       LEN(TRIM(ADJUSTL(bdfileout))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(106)
  ENDIF

  !> 2) Automatically save the solution

  !! b) save the solution in a file named bdfileout
  CALL MMG2D_saveSol(mmgMesh,mmgSol,TRIM(ADJUSTL(bdfileout)),&
       LEN(TRIM(ADJUSTL(bdfileout))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(107)
  ENDIF

  !> 3) Free the MMG2D5 structures
  CALL MMG2D_Free_all(mmgMesh,mmgSol,%val(0))

  !> ================== surface remeshing using the mmgs library

  !! ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures */
  !! args of InitMesh: mesh=&mmgMesh, sol=&mmgSol, input mesh name, input sol
  !! name, output mesh name */
  WRITE(tdfile,*) TRIM(pwd),"/../libexamples/mmg/example0_fortran/cube"
  WRITE(tdfileout,*) TRIM(pwd),"/../libexamples/mmg/example0_fortran/cube.s"

  mmgMesh = 0
  mmgSol  = 0
  !! Remark: %val(0) allow to pass the value 0 (i.e. NULL) instead of a pointer
  !! toward NULL.
  CALL MMGS_Init_mesh(mmgMesh,mmgSol,%val(0))

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMGS_loadMesh function that will read a .mesh(b)
  !!    file formatted or manually set your mesh using the MMGS_Set* functions

  !> with MMGS_loadMesh function
  !! a) (not mandatory): give the mesh name
  !!   (by default, the "mesh.mesh" file is oppened)
  CALL MMGS_Set_inputMeshName(mmgMesh,TRIM(ADJUSTL(tdfile)),&
       LEN(TRIM(ADJUSTL(tdfile))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(101)
  ENDIF

  !> b) function calling
  CALL MMGS_loadMesh(mmgMesh,ier)
  IF ( ier /= 1 )  CALL EXIT(102)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMGS_loadSol function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMGS_Set* functions

  !> With MMGS_loadSol function
  !! a) (not mandatory): give the sol name
  !!   (by default, the "mesh.sol" file is oppened)
  CALL MMGS_Set_inputSolName(mmgMesh,mmgSol,TRIM(ADJUSTL(tdfile)),&
       LEN(TRIM(ADJUSTL(tdfile))),ier)
  IF ( ier/=1 ) THEN
     CALL EXIT(103)
  ENDIF

  !> b) function calling
  CALL MMGS_loadSol(mmgMesh,mmgSol,ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(104)
  ENDIF

  !> 4) (not mandatory): check if the number of given entities match with mesh size
  CALL MMGS_Chk_meshData(mmgMesh,mmgSol,ier)
  IF ( ier /=1 ) CALL EXIT(105)

  !> ------------------------------ STEP  II --------------------------
  !! library call
  CALL MMGS_mmgslib(mmgMesh,mmgSol,ier)

  IF ( ier == MMG5_STRONGFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMGSLIB: UNABLE TO SAVE MESH"
     STOP 2
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMGSLIB"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results
  !! Two solutions: just use the MMGS_saveMesh/MMGS_saveSol functions
  !!    that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !!    using the MMGS_getMesh/MMGS_getSol functions

  !> 1) Automatically save the mesh
  !! a)  (not mandatory): give the ouptut mesh name using MMGS_Set_outputMeshName
  !!   (by default, the mesh is saved in the "mesh.o.mesh" file
  !!call MMGS_Set_outputMeshName(mmgMesh,"output.mesh",len("output.mesh"),ier)
  CALL MMGS_Set_outputMeshName(mmgMesh,TRIM(ADJUSTL(tdfileout)),&
       LEN(TRIM(ADJUSTL(tdfileout))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(105)
  ENDIF

  !! b) function calling
  CALL MMGS_saveMesh(mmgMesh,ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(108)
  ENDIF

  !> 2) Automatically save the solution
  !! a)  (not mandatory): give the ouptut sol name using MMGS_Set_outputSolName
  !!   (by default, the mesh is saved in the "mesh.o.sol" file
  CALL MMGS_Set_outputSolName(mmgMesh,mmgSol,TRIM(ADJUSTL(tdfileout)),&
       LEN(TRIM(ADJUSTL(tdfileout))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(106)
  ENDIF

  !! b) function calling
  CALL MMGS_saveSol(mmgMesh,mmgSol,ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(107)
  ENDIF

  !> 3) Free the MMGS5 structures
  CALL MMGS_Free_all(mmgMesh,mmgSol,%val(0))

  !> ================== 3d remeshing using the mmg3d library

  !! ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures */
  !! args of InitMesh: mesh=&mmgMesh, sol=&mmgSol, input mesh name, input sol
  !! name, output mesh name */
  WRITE(tdfileout,*) TRIM(pwd),"/../libexamples/mmg/example0_fortran/cube.3d"

  mmgMesh = 0
  mmgSol  = 0
  !! Remark: %val(0) allow to pass the value 0 (i.e. NULL) instead of a pointer
  !! toward NULL.
  CALL MMG3D_Init_mesh(mmgMesh,mmgSol,%val(0))

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMG3D_loadMesh function that will read a .mesh(b)
  !!    file formatted or manually set your mesh using the MMG3D_Set* functions

  !> with MMG3D_loadMesh function
  !! a) (not mandatory): give the mesh name
  !!   (by default, the "mesh.mesh" file is oppened)
  CALL MMG3D_Set_inputMeshName(mmgMesh,TRIM(ADJUSTL(tdfile)),&
       LEN(TRIM(ADJUSTL(tdfile))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(101)
  ENDIF

  !> b) function calling
  CALL MMG3D_loadMesh(mmgMesh,ier)
  IF ( ier /= 1 )  CALL EXIT(102)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMG3D_loadSol function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMG3D_Set* functions

  !> With MMG3D_loadSol function
  !! a) (not mandatory): give the sol name
  !!   (by default, the "mesh.sol" file is oppened)
  CALL MMG3D_Set_inputSolName(mmgMesh,mmgSol,TRIM(ADJUSTL(tdfile)),&
       LEN(TRIM(ADJUSTL(tdfile))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(103)
  ENDIF

  !> b) function calling
  CALL MMG3D_loadSol(mmgMesh,mmgSol,ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(104)
  ENDIF

  !> 4) (not mandatory): check if the number of given entities match with mesh size
  CALL MMG3D_Chk_meshData(mmgMesh,mmgSol,ier)
  IF ( ier /= 1 ) CALL EXIT(105)

  !> ------------------------------ STEP  II --------------------------
  !! library call
  CALL MMG3D_mmg3dlib(mmgMesh,mmgSol,ier)

  IF ( ier == MMG5_STRONGFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH"
     STOP 2
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG3DLIB"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results
  !! Two solutions: just use the MMG3D_saveMesh/MMG3D_saveSol functions
  !!    that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !!    using the MMG3D_getMesh/MMG3D_getSol functions

  !> 1) Automatically save the mesh
  !! a)  (not mandatory): give the ouptut mesh name using MMG3D_Set_outputMeshName
  !!   (by default, the mesh is saved in the "mesh.o.mesh" file
  !!call MMG3D_Set_outputMeshName(mmgMesh,"output.mesh",len("output.mesh"),ier)
  CALL MMG3D_Set_outputMeshName(mmgMesh,TRIM(ADJUSTL(tdfileout)),&
       LEN(TRIM(ADJUSTL(tdfileout))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(105)
  ENDIF

  !! b) function calling
  CALL MMG3D_saveMesh(mmgMesh,ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(108)
  ENDIF

  !> 2) Automatically save the solution
  !! a)  (not mandatory): give the ouptut sol name using MMG3D_Set_outputSolName
  !!   (by default, the mesh is saved in the "mesh.o.sol" file
  CALL MMG3D_Set_outputSolName(mmgMesh,mmgSol,TRIM(ADJUSTL(tdfileout)),&
       LEN(TRIM(ADJUSTL(tdfileout))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(106)
  ENDIF

  !! b) function calling
  CALL MMG3D_saveSol(mmgMesh,mmgSol,ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(107)
  ENDIF


  !> 3) Free the MMG3D5 structures
  CALL MMG3D_Free_all(mmgMesh,mmgSol,%val(0))

END PROGRAM main
