# Basic use of the **mmg3d** library for a **Fortran** call

  3 main steps:
    *1) Build mesh and sol at MMG5 format;
    *2) Call MMG5 library;
    *3) Get the final mesh and sol.

  Results are saved in the **_mesh.o.mesh_** and **_mesh.o.sol_** files.  

## example0_a  
  We read mesh and solution files (**_cube.mesh_** and **_cube.sol_**) using the **MMG5_loadMesh** and **MMG5_loadMet** functions.
  Results are saved using **MMG5_saveMesh** and **MMG5_saveMet** functions.

## example0_b
  The mesh and solution are hard coded.    
  They are build in MMG5 format using API functions and are recovered by the same way.  
  We show how to recover the mesh/sol by writting it in mesh.o.mesh/sol file.



