##############
# general stuff
cmake_minimum_required(VERSION 2.8)
set(ProjectName            "DuMuX")
set(ProjectVersion         "2.1.0")
set(ProjectMaintainer      "Bernd Flemisch")
set(ProjectMaintainerEmail "Bernd.Flemisch_at_iws dot uni-stuttgart dot de")
project(${ProjectName} CXX)
# needed for tests like pthread and BLAS
enable_language(C)

##############
# make sure our own modules will be found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

##############
# Set the policy how CMake resolves library paths to the
# policy introduced by CMake 2.6 (this does not apply for
# CMake 2.4 and below, of course). For details, see
# http://www.cmake.org/cmake/help/cmake-2.6.html#policy:CMP0003
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

##############
# Find the required packages
find_package(DUNE_grid REQUIRED)
find_package(DUNE_istl REQUIRED)
find_package(DUNE_localfunctions REQUIRED)
find_package(DUNE_common REQUIRED)

##############
# Find the optional packages
find_package(DUNE_geometry)
find_package(DUNE_pdelab)
find_package(MPI)
find_package(Boost)
find_package(Alberta)
find_package(UG)
find_package(ALUGrid)
find_package(METIS)
find_package(SuperLU)
if(SUPERLU_FOUND)
  set(SUPERLU_CPPFLAGS "-I${SUPERLU_INCLUDE_DIRS} -DENABLE_SUPERLU")
  set(SUPERLU_LIBS "-L. ${SUPERLU_LIBRARIES} ${BLAS_LIBRARIES}")
endif(SUPERLU_FOUND)

##############
# Find the required include files
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
check_include_file("malloc.h" HAVE_MALLOC_H)
check_include_file("valgrind/memcheck.h" HAVE_VALGRIND)
check_include_file_cxx("memory" HAVE_MEMORY)
check_include_file_cxx("tr1/array" HAVE_TR1_ARRAY)
check_include_file_cxx("tr1/memory" HAVE_TR1_MEMORY)

##############
# determine C++11 (former C++0x) feature support
set(CMAKE_REQUIRED_FLAGS "-std=c++0x")

include(CheckCXXSourceCompiles)

# nullptr
CHECK_CXX_SOURCE_COMPILES(
"   int main(void) {
    char* ch = nullptr;
    return 0; }
"  HAVE_NULLPTR
)

# array and fill
# CHECK_CXX_SOURCE_COMPILES("
#     #include <array>
#     
#     int main(void)
#     {
#       std::array<int,2> a;
#       a.fill(9);
#       return 0;
#     }
#     " HAVE_ARRAY
# ) # commented out due to a bug in dune-common <= 2.1.1

# __attribute__((always_inline))
CHECK_CXX_SOURCE_COMPILES(
"
   void __attribute__((always_inline)) foo(void) {}
   int main(void) { foo(); return 0; };
"  HAVE_ATTRIBUTE_ALWAYS_INLINE
)

# __attribute__((unused))
CHECK_CXX_SOURCE_COMPILES(
"
   int main(void) { int __attribute__((unused)) foo; return 0; };
"  HAVE_ATTRIBUTE_UNUSED
)

# __attribute__((deprecated))
CHECK_CXX_SOURCE_COMPILES(
"
#define DEP __attribute__((deprecated))
   class bar { bar() DEP; };
   class peng { } DEP;
   template <class T>
   class t_bar { t_bar() DEP; };
   template <class T>
   class t_peng { t_peng() {}; } DEP;
   void foo() DEP;
   void foo() {};
   int main(void) { return 0; };
"  HAVE_ATTRIBUTE_DEPRECATED
)

# __attribute__((deprecated("msg")))
CHECK_CXX_SOURCE_COMPILES("
#define DEP __attribute__((deprecated(\"message\")))
   class bar { bar() DEP; };
   class peng { } DEP;
   template <class T>
   class t_bar { t_bar() DEP; };
   template <class T>
   class t_peng { t_peng() {}; } DEP;
   void foo() DEP;
   void foo() {};
   int main(void) { return 0; };
"  HAVE_ATTRIBUTE_DEPRECATED_MSG
)

# static assert
CHECK_CXX_SOURCE_COMPILES("
   int main(void) {
      static_assert(true,\"MSG\");
    return 0; }
"  HAVE_STATIC_ASSERT
)

# variadic template support
CHECK_CXX_SOURCE_COMPILES("
   #include <cassert>

   template<typename... T>
   int addints(T... x);

   int add_ints()
   {
     return 0;
   }

   template<typename T1, typename... T>
   int add_ints(T1 t1, T... t)
   {
     return t1 + add_ints(t...);
   }

   int main(void)
   {
     assert( 5 == add_ints(9,3,-5,-2) );
     return 0;
   }
" HAVE_VARIADIC_TEMPLATES
)

# SFINAE on variadic template constructors within template classes
CHECK_CXX_SOURCE_COMPILES("
  #include <functional>

  template<typename... U>
  struct A
  {
    template<typename... T,
             typename = typename std::enable_if<(sizeof...(T) < 2)>::type
            >
    A(T... t)
    : i(1)
    {}

    template<typename... T,
             typename = typename std::enable_if<(sizeof...(T) >= 2)>::type,
             typename = void
            >
    A(T... t)
    : i(-1)
    {}

    A()
    : i(1)
    {}

    int i;
  };

  int main(void)
  {
    return (A<int>().i + A<int>(2).i + A<int>("foo",3.4).i + A<int>(8,'a',A<int>()).i == 0 ? 0 : 1);
  }
" HAVE_VARIADIC_CONSTRUCTOR_SFINAE
)

# rvalue references
CHECK_CXX_SOURCE_COMPILES("
  #include <cassert>
  #include <utility>
  int foo(int&& x) { return 1; }
  int foo(const int& x) { return -1; }

  template<typename T>
  int forward(T&& x)
  {
    return foo(std::forward<T>(x));
  }

  int main(void)
  {
    int i = 0;
    assert( forward(i) + forward(int(2)) == 0);
    return 0;
  }
" HAVE_RVALUE_REFERENCES
)

find_package(SharedPtr)

##############
# use this macros in the CMakelists of the subdirectories. 

# -> for TARGET_LINK_LIBRARIES
set(DumuxLinkLibraries 
#    "dumux" # the DUMUX library. CMake replaces this by the proper file location
    ${DUNE_common_LIBRARIES}
    ${DUNE_geometry_LIBRARIES}
    ${DUNE_grid_LIBRARIES}
    ${DUNE_mux_LIBRARIES}
    "-lm -fprofile-arcs -ftest-coverage")

# -> for INCLUDE_DIRECTORIES
set(DumuxIncludeDirectories 
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}
    ${DUNE_grid_INCLUDE_DIRS}
    ${DUNE_geometry_INCLUDE_DIRS}
    ${DUNE_common_INCLUDE_DIRS}
    ${DUNE_istl_INCLUDE_DIRS}
    ${DUNE_localfunctions_INCLUDE_DIRS}
    ${DUNE_pdelab_INCLUDE_DIRS})

if(BOOST_FOUND)
  set(DumuxLinkDirectories "${DumuxLinkDirectories} ${Boost_LIBRARY_DIRS}")
  set(DumuxLinkLibraries "${DumuxLinkLibraries} ${Boost_LIBRARIES}")
  set(DumuxIncludeDirectories "${DumuxIncludeDirectories} ${Boost_INCLUDE_DIR}")
endif(BOOST_FOUND)

##############
# set appropriate compiler flags for debug/release compilation modes
add_definitions("-std=c++0x -Wall -Wno-sign-compare --no-strict-aliasing")

if(("${CMAKE_BUILD_TYPE}" STREQUAL "debug") OR (NOT (DEFINED "${CMAKE_BUILD_TYPE}")))
  add_definitions(-DNDEBUG)
  # debug mode
  add_definitions("-g -O0 -fprofile-arcs -ftest-coverage")
#  add_definitions(-DDEBUG -DDUNE_DEVEL_MODE=1 -DDUNE_ISTL_WITH_CHECKING)
else()
  # Release mode
  add_definitions("-O3 -march=native")
endif( )

##############
# deal with the config.h include file...
macro(SetConfigHVar ConfigHName CMakeName)
  if(${CMakeName})
    set(${ConfigHName} ${${CMakeName}})
  else(${CMakeName})
    set(${ConfigHName} 0)
  endif(${CMakeName})
endmacro(SetConfigHVar)

SetConfigHVar(HAVE_BOOST               Boost_FOUND)
SetConfigHVar(HAVE_DUNE                DUNE_common_FOUND)
SetConfigHVar(HAVE_DUNE_GRID           DUNE_grid_FOUND)
SetConfigHVar(HAVE_DUNE_ISTL           DUNE_istl_FOUND)
SetConfigHVar(HAVE_DUNE_LOCALFUNCTIONS DUNE_localfunctions_FOUND)
SetConfigHVar(HAVE_DUNE_PDELAB         DUNE_pdelab_FOUND)
SetConfigHVar(HAVE_MPI                 MPI_FOUND)
SetConfigHVar(HAVE_SUPERLU             SUPERLU_FOUND)
SetConfigHVar(PROJECT_NAME             ProjectName)
SetConfigHVar(PROJECT_VERSION          ProjectVersion)
SetConfigHVar(PROJECT_MAINTAINER       ProjectMaintainer)
SetConfigHVar(PROJECT_MAINTAINER_EMAIL ProjectMaintainerEmail)
##############

##############
# add dune-common version from dune.module to config.h
if(DUNE_common_DIR)
  file(READ "${DUNE_common_DIR}/dune.module" DUNE_COMMON_MODULE)
else()
  file(READ "${DUNE_DIR}/dune-common/dune.module" DUNE_COMMON_MODULE)
endif(DUNE_common_DIR)

# find version string
string(REGEX REPLACE ".*Version:[ ]*([^ \n]+).*" "\\1" DUNE_COMMON_VERSION "${DUNE_COMMON_MODULE}")
string(REGEX REPLACE "([0-9]).*" "\\1" DUNE_COMMON_VERSION_MAJOR "${DUNE_COMMON_VERSION}")
string(REGEX REPLACE "[0-9]*\\.([0-9]).*" "\\1" DUNE_COMMON_VERSION_MINOR "${DUNE_COMMON_VERSION}")
string(REGEX REPLACE "[0-9]*\\.[0-9]*\\.([0-9]).*" "\\1" DUNE_COMMON_VERSION_REVISION "${DUNE_COMMON_VERSION}")

# remove false matches
string(REGEX MATCH "[^0-9]" NON_NUMBER_CHARACTER "${DUNE_COMMON_VERSION_MINOR}")
if(NON_NUMBER_CHARACTER)
  set(DUNE_COMMON_VERSION_MINOR "0")
endif(NON_NUMBER_CHARACTER)
string(REGEX MATCH "[^0-9]" NON_NUMBER_CHARACTER "${DUNE_COMMON_VERSION_REVISION}")
if(NON_NUMBER_CHARACTER)
  set(DUNE_COMMON_VERSION_REVISION "0")
endif(NON_NUMBER_CHARACTER)
##############

##############
# adapt build system to detected packages

# deal with UG
if(UG_FOUND)
  set(UG_CPPFLAGS "-I${UG_INCLUDE_DIRS} -DENABLE_UG")
  set(UG_LIBS "-L. ${UG_LIBRARIES}")
endif(UG_FOUND)
SetConfigHVar(HAVE_UG UG_FOUND)

# deal with ALUGrid
if(ALUGrid_FOUND)
  set(DumuxIncludeDirectories ${DumuxIncludeDirectories} ${ALUGrid_INCLUDE_DIRS})
  set(DumuxLinkLibraries ${DumuxLinkLibraries} ${ALUGrid_LIBRARIES})
endif(ALUGrid_FOUND)
SetConfigHVar(HAVE_ALUGRID ALUGrid_FOUND)

# deal with Alberta
if(Alberta_FOUND)
 set(DumuxIncludeDirectories ${DumuxIncludeDirectories} ${Alberta_INCLUDE_DIRS})
 set(DumuxLinkLibraries ${DumuxLinkLibraries} ${Alberta_LIBRARIES})
endif(Alberta_FOUND)
SetConfigHVar(HAVE_ALBERTA Alberta_FOUND)

# deal with METIS
if(METIS_FOUND)
 set(DumuxIncludeDirectories ${DumuxIncludeDirectories} ${METIS_INCLUDE_DIRS})
 set(DumuxLinkLibraries ${DumuxLinkLibraries} ${METIS_LIBRARIES})
endif(METIS_FOUND)
SetConfigHVar(HAVE_METIS METIS_FOUND)
##############

if(MPI_FOUND)
   set(TMP ${MPI_COMPILE_FLAGS})
   separate_arguments(TMP)
   add_definitions(${TMP})
   add_definitions(-DENABLE_MPI=1)
# add_definitions(-DModelP) # tell UG that the model is parallelized
   set(DumuxLinkLibraries ${DumuxLinkLibraries} ${MPI_LIBRARIES})
   set(DumuxIncludeDirectories ${DumuxIncludeDirectories} ${MPI_INCLUDE_PATH})
endif(MPI_FOUND)

# actually write the config.h file to disk 
configure_file(config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h )
#add_definitions(-DHAVE_CONFIG_H)

##############
# tell cmake that we've got a few subdirectories. (that's the
# directories where the actual programs are)
add_subdirectory("test")
add_subdirectory("tutorial")

# copy the testing script
make_directory(bin)
make_directory(references)
make_directory(parameters)
file(COPY bin/runTest.sh DESTINATION bin)
file(COPY bin/fuzzycomparevtu.py DESTINATION bin)
file(COPY test/boxmodels/1p/1ptest-reference.vtu DESTINATION references)
file(COPY test/boxmodels/1p/test_1p.input DESTINATION parameters)
file(COPY test/boxmodels/1p2c/outflow-reference.vtu DESTINATION references)
file(COPY test/boxmodels/1p2c/test_1p2c.input DESTINATION parameters)
file(COPY test/boxmodels/2p/lens-reference.vtu DESTINATION references)
file(COPY test/boxmodels/2p/test_2p.input DESTINATION parameters)
file(COPY test/boxmodels/2pni/injection2pni-reference.vtu DESTINATION references)
file(COPY test/boxmodels/2pni/test_2pni.input DESTINATION parameters)
file(COPY test/boxmodels/2p2c/injection-reference.vtu DESTINATION references)
file(COPY test/boxmodels/2p2c/test_2p2c.input DESTINATION parameters)
file(COPY test/boxmodels/2p2cni/waterair-reference.vtu DESTINATION references)
file(COPY test/boxmodels/2p2cni/test_2p2cni.input DESTINATION parameters)
file(COPY test/boxmodels/3p3c/infiltration-reference.vtu DESTINATION references)
file(COPY test/boxmodels/3p3c/test_3p3c.input DESTINATION parameters)
file(COPY test/boxmodels/3p3cni/kuevette3p3cni-reference.vtu DESTINATION references)
file(COPY test/boxmodels/3p3cni/test_3p3cni.input DESTINATION parameters)
file(COPY test/boxmodels/mpnc/obstacle-reference.vtu DESTINATION references)
file(COPY test/boxmodels/mpnc/test_mpnc.input DESTINATION parameters)
file(COPY test/boxmodels/richards/richardslens-reference.vtu DESTINATION references)
file(COPY test/boxmodels/richards/test_richards.input DESTINATION parameters)
file(COPY test/decoupled/1p/diffusion-reference.vtu DESTINATION references)
file(COPY test/decoupled/1p/test_1p-reference.vtu DESTINATION references)
file(COPY test/decoupled/1p/test_1p.input DESTINATION test/decoupled/1p/)
file(COPY test/decoupled/2p/test_transport-reference.vtu DESTINATION references)
file(COPY test/decoupled/2p/test_transport.input DESTINATION parameters)
file(COPY test/decoupled/2p/test_impes-reference.vtu DESTINATION references)
file(COPY test/decoupled/2p/test_impes.input DESTINATION parameters)
file(COPY test/decoupled/2padaptive/test_2padaptive-reference.vtu DESTINATION references)
file(COPY test/decoupled/2padaptive/test_impesadaptive.input DESTINATION parameters)
file(COPY test/decoupled/2p2c/test_dec2p2c-reference.vtu DESTINATION references)
file(COPY test/decoupled/2p2c/test_dec2p2c.input DESTINATION parameters)
file(COPY test/decoupled/2p2c/test_multiphysics2p2c-reference.vtu DESTINATION references)
file(COPY test/common/generalproblem/generallens_box-reference.vtu DESTINATION references)
file(COPY test/common/generalproblem/generallens_decoupled-reference.vtu DESTINATION references)
file(COPY tutorial/tutorial_coupled.input DESTINATION parameters)
file(COPY tutorial/tutorial_decoupled.input DESTINATION parameters)
if(SUPERLU_FOUND)
  file(COPY test/freeflow/stokes/stokes-reference.vtu DESTINATION references)
  file(COPY test/freeflow/stokes/test_stokes.input DESTINATION parameters)
  file(COPY test/freeflow/stokes2c/stokes2c-reference.vtu DESTINATION references)
  file(COPY test/freeflow/stokes2c/test_stokes2c.input DESTINATION parameters)
  file(COPY test/freeflow/stokes2cni/stokes2cni-reference.vtu DESTINATION references)
  file(COPY test/freeflow/stokes2cni/test_stokes2cni.input DESTINATION parameters)
endif(SUPERLU_FOUND)

# set up CTest 
enable_testing()
include(CTest)
add_test(test_propertysystem test/common/propertysystem/test_propertysystem)
add_test(test_general_box    bin/runTest.sh references/generallens_box-reference.vtu       generallens_box-00003.vtu       test/common/generalproblem/test_generalproblem_2p --box       1e2 2e1)
add_test(test_general_dec    bin/runTest.sh references/generallens_decoupled-reference.vtu generallens_decoupled-00003.vtu test/common/generalproblem/test_generalproblem_2p --decoupled 1e2 2e1)
add_test(test_spline         test/common/spline/test_spline)

add_test(test_fluidsystems    test/material/fluidsystems/test_fluidsystems)
add_test(test_immiscibleflash test/material/immiscibleflash/test_immiscibleflash)
add_test(test_ncpflash        test/material/ncpflash/test_ncpflash)
add_test(test_pengrobinson    test/material/pengrobinson/test_pengrobinson)
add_test(test_tabulation      test/material/tabulation/test_tabulation)

add_test(test_1p       bin/runTest.sh references/1ptest-reference.vtu        1ptest-00002.vtu        test/boxmodels/1p/test_1p             --parameter-file=parameters/test_1p.input       -gridFile test/boxmodels/1p/grids/test_1p_2d.dgf               -tEnd 1    -dtInitial 1)
add_test(test_1p2c     bin/runTest.sh references/outflow-reference.vtu       outflow-00010.vtu       test/boxmodels/1p2c/test_1p2c         --parameter-file=parameters/test_1p2c.input     -gridFile test/boxmodels/1p2c/grids/test_1p2c.dgf              -tEnd 100  -dtInitial 1)
add_test(test_2p       bin/runTest.sh references/lens-reference.vtu          lens-00009.vtu          test/boxmodels/2p/test_2p             --parameter-file=parameters/test_2p.input       -gridFile test/boxmodels/2p/grids/test_2p.dgf                  -tEnd 3000 -dtInitial 250)
add_test(test_2pni     bin/runTest.sh references/injection2pni-reference.vtu injection2pni-00009.vtu test/boxmodels/2pni/test_2pni         --parameter-file=parameters/test_2pni.input                                                                    -tEnd 1e4  -dtInitial 250)
add_test(test_2p2c     bin/runTest.sh references/injection-reference.vtu     injection-00009.vtu     test/boxmodels/2p2c/test_2p2c         --parameter-file=parameters/test_2p2c.input     -gridFile test/boxmodels/2p2c/grids/test_2p2c.dgf              -tEnd 1e4  -dtInitial 250)
add_test(test_2p2cni   bin/runTest.sh references/waterair-reference.vtu      waterair-00010.vtu      test/boxmodels/2p2cni/test_2p2cni     --parameter-file=parameters/test_2p2cni.input   -gridFile test/boxmodels/2p2cni/grids/test_2p2cni.dgf          -tEnd 1e4  -dtInitial 250)
add_test(test_3p3c      bin/runTest.sh references/infiltration-reference.vtu     infiltration-00024.vtu     test/boxmodels/3p3c/test_3p3c         --parameter-file=parameters/test_3p3c.input -gridFile test/boxmodels/3p3c/grids/test_3p3c.dgf)
set_tests_properties(test_3p3c PROPERTIES TIMEOUT 3600)
add_test(test_3p3cni   bin/runTest.sh references/kuevette3p3cni-reference.vtu      kuevette3p3cni-00039.vtu      test/boxmodels/3p3cni/test_3p3cni     --parameter-file=parameters/test_3p3cni.input -gridFile test/boxmodels/3p3cni/grids/kuev_2p2cni.dgf)
set_tests_properties(test_3p3cni PROPERTIES TIMEOUT 3600)
add_test(test_mpnc     bin/runTest.sh references/obstacle-reference.vtu      obstacle-00010.vtu      test/boxmodels/mpnc/test_mpnc         --parameter-file=parameters/test_mpnc.input     -gridFile test/boxmodels/mpnc/grids/obstacle_24x16.dgf         -tEnd 1e4  -dtInitial 250)
add_test(test_richards bin/runTest.sh references/richardslens-reference.vtu  richardslens-00008.vtu  test/boxmodels/richards/test_richards --parameter-file=parameters/test_richards.input -gridFile test/boxmodels/richards/grids/richardslens-24x16.dgf -tEnd 3000 -dtInitial 100)

add_test(test_diffusion bin/runTest.sh references/diffusion-reference.vtu      mimeticdiffusion-00001.vtu test/decoupled/1p/test_diffusion 3)
add_test(test_dec1p     bin/runTest.sh references/test_1p-reference.vtu        test_1p-00001.vtu          test/decoupled/1p/test_dec1p --parameter-file=test/decoupled/1p/test_1p.input)
add_test(test_transport bin/runTest.sh references/test_transport-reference.vtu test_transport-00006.vtu   test/decoupled/2p/test_transport --parameter-file=parameters/test_transport.input -gridFile test/decoupled/2p/grids/test_transport.dgf)
add_test(test_impes     bin/runTest.sh references/test_impes-reference.vtu     test_impes-00013.vtu       test/decoupled/2p/test_impes     --parameter-file=parameters/test_impes.input)
add_test(test_impesadaptive     bin/runTest.sh references/test_2padaptive-reference.vtu test_2padaptive-00007.vtu      test/decoupled/2padaptive/test_impesadaptive  --parameter-file=parameters/test_impesadaptive.input)
add_test(test_dec2p2c   bin/runTest.sh references/test_dec2p2c-reference.vtu   test_dec2p2c-00021.vtu     test/decoupled/2p2c/test_dec2p2c --parameter-file=parameters/test_dec2p2c.input)
add_test(test_multiphysics2p2c   bin/runTest.sh references/test_multiphysics2p2c-reference.vtu   test_multiphysics2p2c-00021.vtu     test/decoupled/2p2c/test_multiphysics2p2c)

add_test(tutorial_coupled   tutorial/tutorial_coupled --parameter-file=parameters/tutorial_coupled.input)
add_test(tutorial_decoupled tutorial/tutorial_decoupled --parameter-file=parameters/tutorial_decoupled.input)

if(SUPERLU_FOUND)
  add_test(test_stokes     bin/runTest.sh references/stokes-reference.vtu     stokes-00013.vtu     test/freeflow/stokes/test_stokes         --parameter-file=parameters/test_stokes.input     -gridFile test/freeflow/stokes/grids/test_stokes.dgf)
  add_test(test_stokes2c   bin/runTest.sh references/stokes2c-reference.vtu   stokes2c-00007.vtu   test/freeflow/stokes2c/test_stokes2c     --parameter-file=parameters/test_stokes2c.input   -gridFile test/freeflow/stokes2c/grids/test_stokes2c.dgf)
  add_test(test_stokes2cni bin/runTest.sh references/stokes2cni-reference.vtu stokes2cni-00008.vtu test/freeflow/stokes2cni/test_stokes2cni --parameter-file=parameters/test_stokes2cni.input -gridFile test/freeflow/stokes2cni/grids/test_stokes2cni.dgf)
endif(SUPERLU_FOUND)
