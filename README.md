# GSLIB 2.0: an updated and improved version of the legendary GSLIB

This repository aims to create a modern version of the original GSLIB source code.

We mean "modern" as a combination of:

* Source code migrated to support modern Fortran
* Makefiles for gfortran and ifort compilers or using fpm
* Release and Debug versions
* Migration of the code to avoid as much as possible:
  * COMMON BLOCKS
  * GOTOs
* Migration of the code to use:
  * Modules
  * Pure/elemental functions and procedures
  * Implicit none

Short term goals:

* Improve algorithms
  * Use KDTree in searching neighbours
  * Migrate parameter specification to use standards such as:
    * TOML
    * JSON
    * YAML
* Parallelisation by using Coarrays

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
* sgsim

We acknowledge the original authors of the GSLIB. We base this new version on the GSLIB version from 1996,
which license is in the file LICENSE_original.txt
