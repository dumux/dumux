Why CMake
=========

You can use CMake 2.6 or higher as alternative to the build system 
provided by DUNE. CMake is included in most GNU/Linux distributions 
or can be downloaded at www.cmake.org. Using CMake has several 
advantages compared to autotools:

 - Out-of-tree builds are the default way to build software: The
   directory where the source code resides won't get modified during
   compilation.
 - Much faster checking of the configuration
 - Much simpler to write your own modules
 - Better dependency handling
 - Less noisy output during compilation

But alas, it comes with a cost:

 - Parameters to the compiler are usually not exactly the same as
   those picked by DUNE's build system. This might sometimes lead to
   problems.
 - Dependencies between DUNE modules are not yet handled, i.e. all
   DUNE modules on which DuMuX depends (-> common, grid, istl, 
   localfunctions, pdelab) need to be installed already.

Preparing the installation
--------------------------

After installing cmake on your system, make sure that the 'cmake'
executable is in your PATH. To configure the project, create a 
new empty directory which will be used to store the files from
the build process and run

cd path/to/empty/build/directory
cmake path/to/DUMUX/source/directory

You might want to set a few parameters if you didn't install the
required libraries system wide. CMake cache variables can be set using
the "-D" command line switch to 'cmake'. The most important paramerters
are probably:

DUNE_DIR        Path where all DUNE modules have been compiled and can be
                found in subdirectories with the module name
DUNE_common_DIR Path where the DUNE-common module can be found. If there 
                is a 'dune-common' subdirectory in $DUNE_DIR, this is not
                required.
(All other DUNE modules are analogous)

BOOST_ROOT        Path where the BOOST libraries from www.boost.org reside. 
                  Usually this is detected automagically.
UG_DIR            Location of the UG grid manager
ALUGrid_DIR       Location of the ALUGrid grid manager
METIS_DIR         Location of the METIS graph partioning library (required
                  for the parallel version of ALUGrid)
CMAKE_BUILD_TYPE  Type of build. Either 'debug' or 'release', default is 
                  'release'.


If this was successfull, the project can be build by 

make

Finally, install everything:

make install

Examples
========

cd path/to/empty/build/directory
cmake -DCMAKE_BUILD_TYPE=debug \
      -DDUNE_DIR=/usr/local/dune/ \
      -DUG_DIR=/usr/local/ug \
      -DALUGrid_DIR=/usr/local/alugrid \
      -DMETIS_DIR=/usr/local/metis \
      path/to/DUMUX/source/directory

With gcc >= 4.5.0, compilation might fail due to an internal compiler 
error. In this case, you might want to specify the compiler explicitly 
by using the cmake options -DCMAKE_CXX_COMPILER and -DCMAKE_CC_COMPILER. 

