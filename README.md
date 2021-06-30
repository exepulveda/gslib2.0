# GSLIB 2.0: an updated and improved version of the legendary GSLIB

This repository aims for having from the original GSLIB source code
a modern version of it.

We mean "modern" as a combination of:

* Source code migrated to support modern Fortran
* Makefiles for gfortran and ifort compilers
* Release and Debug version
* Migration of the code to avoid:
  * COMMON BLOCKS
  * GOTOs
* Migration of the code to use:
  * Modules
  * Pure/elemental functions and procedures
  * Implicit none  

Short term goals:

* Improve algorithms
  * Use KDTree in searching neigbours
  * Migrate paramter specification to use standards such as:
    * TOML
    * JSON
    * YAML
* Parallelisation by implementing Coarrays

## How to build

At this moment, we use the Fortran package manager *fpm* to build libraries and executables.
Check its documentation at https://github.com/fortran-lang/fpm.
The gfortran compiler is supported but aiming to support ifort as well in the short term and also
provide support for Windows.

To build the project just execute ```fpm build```.
To run a particular application *APP*, just execute ```fpm run APP```, where APP is one option from the list
in the Applications section.

## Applications

* kt3d

We acknowledge the original authors of the GSLIB. We base this new version of
GSLIB on the GSLIB version from 1996, which license is in the file LICENSE_original.txt
