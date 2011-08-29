SciD
====

SciD is a collection of numerical routines and bindings written in and for
the D programming language.  The project page can be found at

    https://github.com/cristicbz/scid

Most of the code from the original project by Lars Tandle Kyllingstad has been
rewritten during the 2011 Google Summer of Code by Cristi Cobzarenco, mentored
by David Simcha. David Simcha also contributed code to the current version of
the library.

To compile the library, BLAS and LAPACK implementations are required, the build
script defaults to deps/blaslapack.lib under Windows systems and
/usr/lib/libblas/libblas.a and usr/lib/lapack/liblapack.a under Posix systems.

These are hardcoded into the build script and if different paths are needed
the only way to change them is through modifying build.d


Authors:         Cristi Cobzarenco, Lars Tandle Kyllingstad and David Simcha
Original Author: Lars Tandle Kyllingstad
Licence:         Boost License 1.0 (see licence.txt)


How to build
============

To build the library and generate header (.di) files, run

    rdmd build

To build only the library file, run

    rdmd build lib

To only generate header files, run

    rdmd build headers

To make the documentation, run

    rdmd build html

To clean up after the build script, run

    rdmd build clean

