# Progress callback example

## I/ Implementation
This example shows how to register a progress callback with the **mmg3d**
library before calling **MMG3D_mmg3dlib**.

The callback receives the current remeshing phase, the iteration count and the
number of mesh operations performed during the iteration. It returns 1 to keep
remeshing; returning 0 asks Mmg to stop.

A phase may stop before reaching its maximum iteration count when the mesh has
converged. Mmg emits a final completion notification so progress bars can close
at 100%.

The same pattern can be used with **MMG2D_Set_progressCallback** and
**MMGS_Set_progressCallback** for the 2D and surface libraries.

## II/ Compilation
  1. Build and install the **mmg3d** shared and static library. We suppose in
     the following that you have installed the **mmg3d** library in the
     **_$CMAKE_INSTALL_PREFIX_** directory;
  2. compile the main.c file specifying:
    * the **mmg3d** include directory with the **-I** option;
    * the **mmg3d** library location with the **-L** option;
    * the **mmg3d** library name with the **-l** option;
    * for the static library you must also link the executable with, if used for
      the **mmg3d** library compilation, the scotch and scotcherr libraries and
      with the math library.

> Example 1
> Command line to link the application with the **mmg3d** static library:
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include main.c -L$CMAKE_INSTALL_PREFIX/lib -L$SCOTCH_PATH -lmmg3d -lscotch -lscotcherr -lm
> ```

> Example 2
> Command line to link the application with the **mmg3d** shared library:
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include main.c -L$CMAKE_INSTALL_PREFIX/lib -lmmg3d
> export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
> ```

## III/ Run
```Shell
./a.out cube.mesh cube-progress.mesh
```

For a larger input that makes progress reporting visible, run the command-line
tool on the 2 spheres example:

```Shell
mmg3d -progress ../adaptation_example2/2spheres.mesh -out 2spheres-progress.mesh
```

The **-progress** option keeps the usual Mmg output and adds progress reporting.
With verbose output, progress is printed as regular rows to avoid overwriting
diagnostic messages.
