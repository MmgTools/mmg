# mmg - Surface and volume remeshers
mmg is an open source software for surface and volume remeshing.

It provides two applications:
  * **mmgs**: adaptation and optimization of a surface triangulation representing a piecewise linear approximation of an underlying surface geometry;
  * **mmg3d**: adaptation and optimization of a tetrahedral mesh and implicit domain meshing.

[//]: # ( comment )

## Get and compile the mmg project
  1. Get the repository:  
      `git clone https://github.com/MmgTools/mmg.git`

    The project sources are available under the **_src/_** directory, see:
      * **_src/mmgs/_**   for files related to the mmgs application;
      * **_src/mmg3d/_**  for files related to the mmg3d application;
      * **_src/common/_** for files related to the both.

  2. Fast compilation (build both **mmgs**, **mmg3d** and the mmg3d static library **libmmg3d.a**):  
      `cd mmg`  
      `mkdir build`  
      `cd build`  
      `cmake ..`  
      `make`  
      `make install`    

    The **mmgs** and **mmg3d** applications are available under the `mmgs_O3` and `mmg3d_O3` commands. 

## Documentation
### Wiki
More detailed informations about the compilation and configuration of the mmg's applications are available on the project [wiki](https://github.com/MmgTools/mmg/wiki).

[//]: ### User guide
[//]: The mmg3d user guide will be soon available on the project webpage: [http://www.mmgtools.org/doc/](http://www.mmgtools.org/doc/).

### Code documentation
Run the `make doc` command to build the Doxygen documentation.
  * To see the **mmgs** documentation, open up the **_mmg/doc/mmgs/html/index.html_** file;
  * To see the **mmg3d** documentation, open up the **_mmg/doc/mmg3d/html/index.html_** file.

## Platforms
The **mmg** applications are validated on OS X and on most of the Linux platforms. 

## Contributing
Coming soon...

## About the team
mmg's current developers and maintainers are Charles Dapogny, Cécile Dobrzynski, Pascal Frey and Algiane Froehly.

## License and copyright
Code is under the [terms of the GNU Lesser General Public License](https://raw.githubusercontent.com/MmgTools/mmg/master/LICENSE).

Copyright © Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
