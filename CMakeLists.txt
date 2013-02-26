##############
# general stuff
cmake_minimum_required(VERSION 2.8)
set(ProjectName            "DuMuX")
set(ProjectVersion         "2.3-svn")
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

# tuple
CHECK_CXX_SOURCE_COMPILES(
"  #include <tuple> 
   int main(void) {
    return 0; }
"  HAVE_TUPLE
)

# constexpr
CHECK_CXX_SOURCE_COMPILES(
"  int main(void) {
    constexpr double g = 9.81;
    return 0; }
"  HAVE_CONSTEXPR
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
"  HAS_ATTRIBUTE_UNUSED
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
"  HAS_ATTRIBUTE_DEPRECATED
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
"  HAS_ATTRIBUTE_DEPRECATED_MSG
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
add_definitions("-std=c++0x -Wall -Wno-sign-compare -fno-strict-aliasing -DUSE_AMGBACKEND")

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
  set(UG_CPPFLAGS "-I${UG_INCLUDE_DIRS} -DENABLE_UG=1")
  set(UG_LIBS "-L. ${UG_LIBRARIES}")
  set(HAVE_UG "ENABLE_UG")
else(UG_FOUND)
  set(HAVE_UG 0)
endif(UG_FOUND)

# deal with ALUGrid
if(ALUGrid_FOUND)
  add_definitions(-DENABLE_ALUGRID=1)
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
file(COPY test/references/ DESTINATION references)
file(COPY bin/runTest.sh DESTINATION bin)
file(COPY bin/fuzzycomparevtu.py DESTINATION bin)
file(COPY test/implicit/1p/test_box1p.input DESTINATION test/implicit/1p)
file(COPY test/implicit/1p/test_box1pwithamg.input DESTINATION test/implicit/1p)
file(COPY test/implicit/1p2c/test_box1p2c.input DESTINATION test/implicit/1p2c)
file(COPY test/implicit/2p/test_box2p.input DESTINATION test/implicit/2p)
file(COPY test/implicit/2pni/test_box2pni.input DESTINATION test/implicit/2pni)
file(COPY test/implicit/2p2c/test_box2p2c.input DESTINATION test/implicit/2p2c)
file(COPY test/implicit/2p2cni/test_box2p2cni.input DESTINATION test/implicit/2p2cni)
file(COPY test/implicit/3p3c/test_box3p3c.input DESTINATION test/implicit/3p3c)
file(COPY test/implicit/3p3cni/test_box3p3cni.input DESTINATION test/implicit/3p3cni)
file(COPY test/implicit/co2/test_ccco2.input DESTINATION test/implicit/co2)
file(COPY test/implicit/co2/test_boxco2.input DESTINATION test/implicit/co2)
file(COPY test/implicit/co2ni/test_ccco2ni.input DESTINATION test/implicit/co2ni)
file(COPY test/implicit/co2ni/test_boxco2ni.input DESTINATION test/implicit/co2ni)
file(COPY test/implicit/mpnc/test_boxmpnc.input DESTINATION test/implicit/mpnc)
file(COPY test/implicit/mpnc/test_forchheimer1p.input DESTINATION test/implicit/mpnc)
file(COPY test/implicit/mpnc/test_forchheimer2p.input DESTINATION test/implicit/mpnc)
file(COPY test/implicit/richards/test_boxrichards.input DESTINATION test/implicit/richards)
file(COPY test/implicit/1p/test_cc1p.input DESTINATION test/implicit/1p)
file(COPY test/implicit/1p/test_cc1pwithamg.input DESTINATION test/implicit/1p)
file(COPY test/implicit/1p2c/test_cc1p2c.input DESTINATION test/implicit/1p2c)
file(COPY test/implicit/2p/test_cc2p.input DESTINATION test/implicit/2p)
file(COPY test/implicit/2pni/test_cc2pni.input DESTINATION test/implicit/2pni)
file(COPY test/implicit/2p2c/test_cc2p2c.input DESTINATION test/implicit/2p2c)
file(COPY test/implicit/2p2cni/test_cc2p2cni.input DESTINATION test/implicit/2p2cni)
file(COPY test/implicit/3p3c/test_cc3p3c.input DESTINATION test/implicit/3p3c)
file(COPY test/implicit/3p3cni/test_cc3p3cni.input DESTINATION test/implicit/3p3cni)
file(COPY test/implicit/mpnc/test_ccmpnc.input DESTINATION test/implicit/mpnc)
file(COPY test/implicit/richards/test_ccrichards.input DESTINATION test/implicit/richards)
file(COPY test/decoupled/1p/test_1p.input DESTINATION test/decoupled/1p/)
file(COPY test/decoupled/2p/test_transport.input DESTINATION test/decoupled/2p)
file(COPY test/decoupled/2p/test_impes.input DESTINATION test/decoupled/2p)
file(COPY test/decoupled/2p/test_impeswithamg.input DESTINATION test/decoupled/2p)
file(COPY test/decoupled/2p/test_mpfa2p.input DESTINATION test/decoupled/2p)
file(COPY test/decoupled/2p/test_impesadaptive.input DESTINATION test/decoupled/2p)
file(COPY test/decoupled/2p2c/test_adaptive2p2c.input DESTINATION test/decoupled/2p2c)
file(COPY test/decoupled/2p2c/test_dec2p2c.input DESTINATION test/decoupled/2p2c)
file(COPY test/common/generalproblem/test_generalproblem2p.input DESTINATION test/common/generalproblem)
file(COPY tutorial/tutorial_coupled.input DESTINATION tutorial)
file(COPY tutorial/tutorial_decoupled.input DESTINATION tutorial)
if(SUPERLU_FOUND)
  file(COPY test/freeflow/stokes/test_stokes.input DESTINATION test/freeflow/stokes)
  file(COPY test/freeflow/stokes2c/test_stokes2c.input DESTINATION test/freeflow/stokes2c)
  file(COPY test/freeflow/stokes2cni/test_stokes2cni.input DESTINATION test/freeflow/stokes2cni)
  file(COPY test/freeflow/navierstokes/test_navierstokes.input DESTINATION test/freeflow/navierstokes)
endif(SUPERLU_FOUND)

# set up CTest 
enable_testing()
include(CTest)
add_test(test_boxco2      bin/runTest.sh references/co2box-reference.vtu            heterogeneousbox-00020.vtu            test/implicit/co2/test_boxco2 -Grid.File test/implicit/co2/grids/heterogeneousSmall.dgf) 
add_test(test_propertysystem test/common/propertysystem/test_propertysystem)
add_test(test_general_box    bin/runTest.sh references/generallens_box-reference.vtu       generallens_box-00003.vtu       test/common/generalproblem/test_generalproblem2p -ModelType Box       -TimeManager.TEnd 1e2 -TimeManager.DtInitial 2e1)
add_test(test_general_dec    bin/runTest.sh references/generallens_decoupled-reference.vtu generallens_decoupled-00003.vtu test/common/generalproblem/test_generalproblem2p -ModelType Decoupled -TimeManager.TEnd 1e2 -TimeManager.DtInitial 2e1)
add_test(test_spline         test/common/spline/test_spline)

add_test(test_fluidsystems    test/material/fluidsystems/test_fluidsystems)
add_test(test_immiscibleflash test/material/immiscibleflash/test_immiscibleflash)
add_test(test_ncpflash        test/material/ncpflash/test_ncpflash)
add_test(test_pengrobinson    test/material/pengrobinson/test_pengrobinson)
add_test(test_tabulation      test/material/tabulation/test_tabulation)

add_test(test_box1p       bin/runTest.sh references/1ptestbox-reference.vtu         1ptestbox-00002.vtu         test/implicit/1p/test_box1p             -Grid.File test/implicit/1p/grids/test_1p_2d.dgf               -TimeManager.TEnd 1       -TimeManager.DtInitial 1)
add_test(test_cc1p        bin/runTest.sh references/1ptestcc-reference.vtu          1ptestcc-00002.vtu          test/implicit/1p/test_cc1p              -Grid.File test/implicit/1p/grids/test_1p_2d.dgf               -TimeManager.TEnd 1       -TimeManager.DtInitial 1)
add_test(test_box1p2c     bin/runTest.sh references/outflowbox-reference.vtu        outflowbox-00010.vtu        test/implicit/1p2c/test_box1p2c         -Grid.File test/implicit/1p2c/grids/test_1p2c.dgf              -TimeManager.TEnd 100     -TimeManager.DtInitial 1)
add_test(test_cc1p2c      bin/runTest.sh references/outflowcc-reference.vtu         outflowcc-00010.vtu         test/implicit/1p2c/test_cc1p2c          -Grid.File test/implicit/1p2c/grids/test_1p2c.dgf              -TimeManager.TEnd 100     -TimeManager.DtInitial 1)
add_test(test_box2p       bin/runTest.sh references/lensbox-reference.vtu           lensbox-00009.vtu           test/implicit/2p/test_box2p             -Grid.File test/implicit/2p/grids/test_2p.dgf                  -TimeManager.TEnd 3000    -TimeManager.DtInitial 250)
add_test(test_cc2p        bin/runTest.sh references/lenscc-reference.vtu            lenscc-00009.vtu            test/implicit/2p/test_cc2p              -Grid.File test/implicit/2p/grids/test_2p.dgf                  -TimeManager.TEnd 3000    -TimeManager.DtInitial 250)
add_test(test_box2pni     bin/runTest.sh references/injection2pnibox-reference.vtu  injection2pnibox-00009.vtu  test/implicit/2pni/test_box2pni                                                                         -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
add_test(test_cc2pni      bin/runTest.sh references/injection2pnicc-reference.vtu   injection2pnicc-00009.vtu   test/implicit/2pni/test_cc2pni                                                                         -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
add_test(test_box2p2c     bin/runTest.sh references/injectionbox-reference.vtu      injectionbox-00009.vtu      test/implicit/2p2c/test_box2p2c         -Grid.File test/implicit/2p2c/grids/test_2p2c.dgf              -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
add_test(test_cc2p2c      bin/runTest.sh references/injectioncc-reference.vtu       injectioncc-00009.vtu       test/implicit/2p2c/test_cc2p2c          -Grid.File test/implicit/2p2c/grids/test_2p2c.dgf              -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
add_test(test_box2p2cni   bin/runTest.sh references/waterairbox-reference.vtu       waterairbox-00010.vtu       test/implicit/2p2cni/test_box2p2cni     -Grid.File test/implicit/2p2cni/grids/test_2p2cni.dgf          -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
add_test(test_cc2p2cni    bin/runTest.sh references/wateraircc-reference.vtu        wateraircc-00010.vtu        test/implicit/2p2cni/test_cc2p2cni      -Grid.File test/implicit/2p2cni/grids/test_2p2cni.dgf          -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
add_test(test_box3p3c     bin/runTest.sh references/infiltrationbox-reference.vtu   infiltrationbox-00007.vtu   test/implicit/3p3c/test_box3p3c         -Grid.File test/implicit/3p3c/grids/test_3p3c_coarse.dgf       -TimeManager.TEnd 8.64e5  -TimeManager.DtInitial 8.64e4)
add_test(test_cc3p3c      bin/runTest.sh references/infiltrationcc-reference.vtu    infiltrationcc-00006.vtu    test/implicit/3p3c/test_cc3p3c          -Grid.File test/implicit/3p3c/grids/test_3p3c_coarse.dgf       -TimeManager.TEnd 8.64e5  -TimeManager.DtInitial 8.64e4)
add_test(test_box3p3cni   bin/runTest.sh references/column3p3cnibox-reference.vtu   columnxylolbox-00060.vtu    test/implicit/3p3cni/test_box3p3cni     -Grid.File test/implicit/3p3cni/grids/column.dgf)
add_test(test_cc3p3cni    bin/runTest.sh references/column3p3cnicc-reference.vtu    columnxylolcc-00054.vtu     test/implicit/3p3cni/test_cc3p3cni      -Grid.File test/implicit/3p3cni/grids/column.dgf)
add_test(test_ccco2      bin/runTest.sh references/co2cc-reference.vtu            heterogeneouscc-00019.vtu            test/implicit/co2/test_ccco2 -Grid.File test/implicit/co2/grids/heterogeneousSmall.dgf) 
add_test(test_ccco2ni    bin/runTest.sh references/co2nicc-reference.vtu          heterogeneousccni-00019.vtu          test/implicit/co2ni/test_ccco2ni -Grid.File test/implicit/co2ni/grids/heterogeneousSmall.dgf)
add_test(test_boxco2ni    bin/runTest.sh references/co2nibox-reference.vtu          heterogeneousboxni-00020.vtu          test/implicit/co2ni/test_boxco2ni -Grid.File test/implicit/co2ni/grids/heterogeneousSmall.dgf)
add_test(test_boxmpnc     bin/runTest.sh references/obstaclebox-reference.vtu       obstaclebox-00010.vtu       test/implicit/mpnc/test_boxmpnc         -Grid.File test/implicit/mpnc/grids/obstacle_24x16.dgf         -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
add_test(test_ccmpnc      bin/runTest.sh references/obstaclecc-reference.vtu        obstaclecc-00010.vtu        test/implicit/mpnc/test_ccmpnc          -Grid.File test/implicit/mpnc/grids/obstacle_24x16.dgf         -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
add_test(test_forchheimer1p    bin/runTest.sh references/forchheimer1p-reference.vtp      test_forchheimer1p-00002.vtp       test/implicit/mpnc/test_forchheimer1p  -Grid.File test/implicit/mpnc/grids/forchheimer1d.dgf)
add_test(test_forchheimer2p    bin/runTest.sh references/forchheimer2p-reference.vtu      test_forchheimer2p-00009.vtu       test/implicit/mpnc/test_forchheimer2p  -Grid.File test/implicit/mpnc/grids/obstacle_24x16.dgf)
if(MPI_FOUND)
  add_test(test_boxrichards bin/runTest.sh references/richardslensbox-reference-parallel.vtu s0002-p0000-richardslensbox-00008.vtu  ${MPIEXEC} -np 2 test/implicit/richards/test_boxrichards -Grid.File test/implicit/richards/grids/richardslens-24x16.dgf -TimeManager.TEnd 3000    -TimeManager.DtInitial 100)
else(MPI_FOUND)
  add_test(test_boxrichards bin/runTest.sh references/richardslensbox-reference.vtu   richardslensbox-00008.vtu   test/implicit/richards/test_boxrichards -Grid.File test/implicit/richards/grids/richardslens-24x16.dgf -TimeManager.TEnd 3000    -TimeManager.DtInitial 100)
  add_test(test_ccrichards  bin/runTest.sh references/richardslenscc-reference.vtu    richardslenscc-00008.vtu    test/implicit/richards/test_ccrichards  -Grid.File test/implicit/richards/grids/richardslens-24x16.dgf -TimeManager.TEnd 3000    -TimeManager.DtInitial 100)
endif(MPI_FOUND)

add_test(test_diffusion        bin/runTest.sh references/diffusion-reference.vtu             mimeticdiffusion-00001.vtu      test/decoupled/1p/test_diffusion 3)
add_test(test_dec1p            bin/runTest.sh references/test_1p-reference.vtu               test_1p-00001.vtu               test/decoupled/1p/test_dec1p -ParameterFile test/decoupled/1p/test_1p.input)
add_test(test_transport        bin/runTest.sh references/test_transport-reference.vtu        test_transport-00006.vtu        test/decoupled/2p/test_transport -Grid.File test/decoupled/2p/grids/test_transport.dgf)
add_test(test_impes            bin/runTest.sh references/test_impes-reference.vtu            test_impes-00013.vtu            test/decoupled/2p/test_impes)
add_test(test_mpfao2p          bin/runTest.sh references/test_mpfao2p-reference.vtu          test_mpfa2p-00007.vtu           test/decoupled/2p/test_mpfa2p -ModelType MPFAO)
add_test(test_mpfal2p          bin/runTest.sh references/test_mpfal2p-reference.vtu          test_mpfa2p-00007.vtu           test/decoupled/2p/test_mpfa2p -ModelType MPFAL)
add_test(test_mpfal2padaptive  bin/runTest.sh references/test_mpfal2padaptive-reference.vtu  test_mpfa2p-00007.vtu           test/decoupled/2p/test_mpfa2p -ModelType MPFALAdaptive)
add_test(test_impesadaptive    bin/runTest.sh references/test_2padaptive-reference.vtu       test_2padaptive-00007.vtu       test/decoupled/2p/test_impesadaptive)
add_test(test_adaptive2p2c     bin/runTest.sh references/test_adaptive2p2c-reference.vtu     test_adaptive2p2c-00008.vtu     test/decoupled/2p2c/test_adaptive2p2c)
add_test(test_dec2p2c          bin/runTest.sh references/test_dec2p2c-reference.vtu          test_dec2p2c-00021.vtu          test/decoupled/2p2c/test_dec2p2c)
add_test(test_multiphysics2p2c bin/runTest.sh references/test_multiphysics2p2c-reference.vtu test_multiphysics2p2c-00021.vtu test/decoupled/2p2c/test_multiphysics2p2c)

add_test(tutorial_coupled      tutorial/tutorial_coupled)
add_test(tutorial_decoupled    tutorial/tutorial_decoupled)

if(DUNE_pdelab_FOUND)
  add_test(test_box1pwithamg   bin/runTest.sh references/1ptestbox-reference.vtu             1ptestboxwithamg-00002.vtu      test/implicit/1p/test_box1pwithamg      -Grid.File test/implicit/1p/grids/test_1p_2d.dgf               -TimeManager.TEnd 1       -TimeManager.DtInitial 1)
  add_test(test_cc1pwithamg    bin/runTest.sh references/1ptestcc-reference.vtu              1ptestccwithamg-00002.vtu       test/implicit/1p/test_cc1pwithamg       -Grid.File test/implicit/1p/grids/test_1p_2d.dgf               -TimeManager.TEnd 1       -TimeManager.DtInitial 1)
  if(MPI_FOUND)
    add_test(test_impeswithamg bin/runTest.sh references/test_impes-reference-parallel.vtu   s0002-p0001-test_impeswithamg-00013.vtu ${MPIEXEC} -np 2 test/decoupled/2p/test_impeswithamg -Grid.File test/decoupled/2p/grids/test_impeswithamg.dgf)
  else(MPI_FOUND)
    add_test(test_impeswithamg bin/runTest.sh references/test_impes-reference.vtu            test_impeswithamg-00013.vtu     test/decoupled/2p/test_impeswithamg     -Grid.File test/decoupled/2p/grids/test_impeswithamg.dgf)
  endif(MPI_FOUND)
endif(DUNE_pdelab_FOUND)

if(SUPERLU_FOUND)
  add_test(test_stokes       bin/runTest.sh references/stokes-reference.vtu       stokes-00013.vtu       test/freeflow/stokes/test_stokes             -Grid.File test/freeflow/stokes/grids/test_stokes.dgf)
  add_test(test_stokes2c     bin/runTest.sh references/stokes2c-reference.vtu     stokes2c-00007.vtu     test/freeflow/stokes2c/test_stokes2c         -Grid.File test/freeflow/stokes2c/grids/test_stokes2c.dgf)
  add_test(test_stokes2cni   bin/runTest.sh references/stokes2cni-reference.vtu   stokes2cni-00008.vtu   test/freeflow/stokes2cni/test_stokes2cni     -Grid.File test/freeflow/stokes2cni/grids/test_stokes2cni.dgf)
  add_test(test_navierstokes bin/runTest.sh references/navierstokes-reference.vtu navierstokes-00008.vtu test/freeflow/navierstokes/test_navierstokes -Grid.File test/freeflow/navierstokes/grids/test_navierstokes.dgf)
endif(SUPERLU_FOUND)

add_test(test_restart  bin/runTest.sh references/co2box-reference.vtu            heterogeneousbox-00020.vtu            test/implicit/co2/test_boxco2 -Grid.File test/implicit/co2/grids/heterogeneousSmall.dgf -Restart 42008.6 -TimeManager.DtInitial 23452.1)
