# Example of advanced use of library libmmg3d5

  We read the mesh and solution in the **_2spheres.mesh_** and **_2spheres.sol_** files.

  * First we remesh in debug mode:
    * we ask for a minimal size of 0.001, a maximal size of 40, a gradation of 2 and a global hausdorff value (applied on all the boundaries) of 0.1;
    * we save results in 2spheres_1.o.mesh/sol.

  * Second, we remesh in normal mode, with specified memory and lower verbosity:
    * in addition to previous parameters, we ask that all boundary triangles of ref 36 respect a hasdorff number of 0.01 and all boundary triangles of ref 38 a hasdorff number of 1.  
      The local hausdorff number on ref 38 has no effects because it is higher than the previous value (without local value, we apply global hausdorff (0.1)) and this value is now contained in the metric;
    * we don't save results but we reset the computed metric and reapply the initial constant metric of size 10;
    * we perform the last wave of refinment. Now we can see the effect of the local hausdorff number on ref 38;
    * we save the mesh in 2spheres_2.o.mesh/sol.

