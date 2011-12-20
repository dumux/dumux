##############
# general stuff
cmake_minimum_required(VERSION 2.6)
set(ProjectName            "DuMuX")
set(ProjectVersion         "2.1-svn")
set(ProjectMaintainer      "Bernd Flemisch")
set(ProjectMaintainerEmail "Bernd.Flemisch_at_iws dot uni-stuttgart dot de")
project(${ProjectName} CXX)
##############

##############
# make sure our own modules will be found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules")
##############

##############
# Set the policy how CMake resolves library paths to the
# policy introduced by CMake 2.6 (this does not apply for
# CMake 2.4 and below, of course). For details, see
# http://www.cmake.org/cmake/help/cmake-2.6.html#policy:CMP0003
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)
##############


##############
# Find the required packages
FIND_PACKAGE(DUNE_grid REQUIRED)
FIND_PACKAGE(DUNE_istl REQUIRED)
FIND_PACKAGE(DUNE_localfunctions REQUIRED)
FIND_PACKAGE(DUNE_pdelab REQUIRED)
FIND_PACKAGE(DUNE_common REQUIRED)
FIND_PACKAGE(Boost REQUIRED)
##############

##############
# Find the optional packages
FIND_PACKAGE(MPI)
#FIND_PACKAGE(UG)
FIND_PACKAGE(ALUGrid)
FIND_PACKAGE(METIS)
#FIND_PACKAGE(Alberta)
##############

##############
# Find the required include files
INCLUDE (CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("malloc.h" HAVE_MALLOC_H)
CHECK_INCLUDE_FILE_CXX("valgrind/memcheck.h" HAVE_VALGRIND)
CHECK_INCLUDE_FILE_CXX("tr1/array" HAVE_TR1_ARRAY)
##############


##############
# use this macros in the CMakelists of the subdirectories. 

# -> for TARGET_LINK_LIBRARIES
set(DumuxLinkLibraries 
#    "dumux" # the DUMUX library. CMake replaces this by the proper file location
    ${DUNE_grid_LIBRARIES}
    ${DUNE_common_LIBRARIES}
    ${DUNE_mux_LIBRARIES}
    ${Boost_LIBRARIES})

# -> for LINK_DIRECTORIES
set(DumuxLinkDirectories ${Boost_LIBRARY_DIRS})

# -> for INCLUDE_DIRECTORIES
set(DumuxIncludeDirectories 
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}
    ${DUNE_grid_INCLUDE_DIRS}
    ${DUNE_common_INCLUDE_DIRS}
    ${DUNE_istl_INCLUDE_DIRS}
    ${DUNE_localfunctions_INCLUDE_DIRS}
    ${DUNE_pdelab_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR})
##############


##############
# set appropriate compiler flags for debug/release compilation modes
add_definitions("-std=c++0x -Wall -Wno-sign-compare --no-strict-aliasing")

if(CMAKE_BUILD_TYPE STREQUAL "debug")
  # debug mode
  add_definitions("-g")
#  add_definitions(-DDEBUG -DDUNE_DEVEL_MODE=1 -DDUNE_ISTL_WITH_CHECKING)
else(CMAKE_BUILD_TYPE STREQUAL "debug")
  # Release mode
  add_definitions("-O3 -march=native")
  add_definitions(-DNDEBUG)
endif(CMAKE_BUILD_TYPE STREQUAL "debug")
##############


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
SetConfigHVar(PROJECT_NAME             ProjectName)
SetConfigHVar(PROJECT_VERSION          ProjectVersion)
SetConfigHVar(PROJECT_MAINTAINER       ProjectMaintainer)
SetConfigHVar(PROJECT_MAINTAINER_EMAIL ProjectMaintainerEmail)
##############

##############
# adapt build system to detected packages

# deal with UG
#if (UG_FOUND)
#  set(DumuxLinkLibraries ${DumuxLinkLibraries} 
#      ${UG_LIBRARIES})
#  set(DumuxIncludeDirectories ${DumuxIncludeDirectories} 
#      ${UG_INCLUDE_DIRS})
#endif(UG_FOUND)
#SetConfigHVar(HAVE_UG UG_FOUND)

# deal with ALUGrid
#if (ALUGrid_FOUND)
#  set(DumuxLinkLibraries ${DumuxLinkLibraries} 
#      ${ALUGrid_LIBRARIES})
#  set(DumuxIncludeDirectories ${DumuxIncludeDirectories} 
#      ${ALUGrid_INCLUDE_DIRS})
#endif(ALUGrid_FOUND)
#SetConfigHVar(HAVE_ALUGRID ALUGrid_FOUND)

# deal with METIS
#if (METIS_FOUND)
#  set(DumuxLinkLibraries ${DumuxLinkLibraries} 
#      ${METIS_LIBRARIES})
#  set(DumuxIncludeDirectories ${DumuxIncludeDirectories} 
#      ${METIS_INCLUDE_DIRS})
#endif(METIS_FOUND)
#SetConfigHVar(HAVE_METIS METIS_FOUND)

# deal with Alberta
#if (Alberta_FOUND)
#  set(DumuxLinkLibraries ${DumuxLinkLibraries} 
#      ${Alberta_LIBRARIES})
#  set(DumuxIncludeDirectories ${DumuxIncludeDirectories} 
#      ${Alberta_INCLUDE_DIRS})
#endif(Alberta_FOUND)
#SetConfigHVar(HAVE_ALBERTA Alberta_FOUND)
##############

if (MPI_FOUND)
   set(TMP ${MPI_COMPILE_FLAGS})
   separate_arguments(TMP)
   add_definitions(${TMP})
   add_definitions(-DENABLE_MPI=1)
#   add_definitions(-DModelP) # tell UG that the model is parallelized

   set(DumuxLinkLibraries ${DumuxLinkLibraries} ${MPI_LIBRARIES})
   set(DumuxIncludeDirectories ${DumuxIncludeDirectories} ${MPI_INCLUDE_PATH})

endif (MPI_FOUND)

# actually write the config.h file to disk 
configure_file(config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h )
#add_definitions(-DHAVE_CONFIG_H)

##############
# tell cmake that we've got a few subdirectories. (that's the
# directories where the actual programs are)
add_subdirectory("test")
add_subdirectory("tutorial")
##############

# copy the testing script
make_directory(util)
file(COPY util/runTest.sh DESTINATION util)
file(COPY test/boxmodels/1p/1ptest-reference.vtu DESTINATION test/boxmodels/1p)
file(COPY test/boxmodels/1p2c/tissue-reference.vtu DESTINATION test/boxmodels/1p2c)
file(COPY test/boxmodels/2p/lens-reference.vtu DESTINATION test/boxmodels/2p)
file(COPY test/boxmodels/2pni/injection2pni-reference.vtu DESTINATION test/boxmodels/2pni)
file(COPY test/boxmodels/2p2c/injection-reference.vtu DESTINATION test/boxmodels/2p2c)
file(COPY test/boxmodels/2p2cni/waterair-reference.vtu DESTINATION test/boxmodels/2p2cni)
file(COPY test/boxmodels/MpNc/obstacle-reference.vtu DESTINATION test/boxmodels/MpNc)
file(COPY test/boxmodels/richards/richardslens-reference.vtu DESTINATION test/boxmodels/richards)
file(COPY test/decoupled/1p/diffusion-reference.vtu DESTINATION test/decoupled/1p)
file(COPY test/decoupled/1p/test_1p-reference.vtu DESTINATION test/decoupled/1p)
file(COPY test/decoupled/2p/test_transport-reference.vtu DESTINATION test/decoupled/2p)
file(COPY test/decoupled/2p/test_impes-reference.vtu DESTINATION test/decoupled/2p)
file(COPY test/decoupled/2p2c/test_dec2p2c-reference.vtu DESTINATION test/decoupled/2p2c)

# set up CTest 
ENABLE_TESTING()
INCLUDE(CTest)
add_test(test_propertysystem test/common/propertysystem/test_propertysystem)
add_test(test_generalproblem test/common/generalproblem/test_generalproblem_2p)
add_test(test_spline test/common/spline/test_spline)

add_test(test_fluidsystems test/material/fluidsystems/test_fluidsystems)
add_test(test_ncpflash test/material/ncpflash/test_ncpflash)
add_test(test_tabulation test/material/tabulation/test_tabulation)

add_test(test_1p util/runTest.sh test/boxmodels/1p/1ptest-reference.vtu 1ptest-00002.vtu test/boxmodels/1p/test_1p test/boxmodels/1p/grids/test_1p_2d.dgf 1 1)
add_test(test_1p2c util/runTest.sh test/boxmodels/1p2c/tissue-reference.vtu tissue-00003.vtu test/boxmodels/1p2c/test_1p2c test/boxmodels/1p2c/grids/test_1p2c.dgf 2e3 1e3)
add_test(test_2p util/runTest.sh test/boxmodels/2p/lens-reference.vtu lens-00010.vtu test/boxmodels/2p/test_2p 3000 250)
add_test(test_2pni util/runTest.sh test/boxmodels/2pni/injection2pni-reference.vtu injection2pni-00009.vtu test/boxmodels/2pni/test_2pni test/boxmodels/2pni/grids/test_2pni.dgf 1e4 250)
add_test(test_2p2c util/runTest.sh test/boxmodels/2p2c/injection-reference.vtu injection-00009.vtu test/boxmodels/2p2c/test_2p2c test/boxmodels/2p2c/grids/test_2p2c.dgf 1e4 250)
add_test(test_2p2cni util/runTest.sh test/boxmodels/2p2cni/waterair-reference.vtu waterair-00011.vtu test/boxmodels/2p2cni/test_2p2cni test/boxmodels/2p2cni/grids/test_2p2cni.dgf 1e4 250)
add_test(test_MpNc util/runTest.sh test/boxmodels/MpNc/obstacle-reference.vtu obstacle-00011.vtu test/boxmodels/MpNc/test_MpNc test/boxmodels/MpNc/grids/obstacle_24x16.dgf 1e4 250)
add_test(test_richards util/runTest.sh test/boxmodels/richards/richardslens-reference.vtu richardslens-00009.vtu test/boxmodels/richards/test_richards test/boxmodels/richards/grids/richardslens.dgf 100 10)

add_test(test_diffusion util/runTest.sh test/decoupled/1p/diffusion-reference.vtu mimeticdiffusion-00001.vtu test/decoupled/1p/test_diffusion 3)
add_test(test_dec1p util/runTest.sh test/decoupled/1p/test_1p-reference.vtu test_1p-00001.vtu test/decoupled/1p/test_dec1p 3)
add_test(test_transport util/runTest.sh test/decoupled/2p/test_transport-reference.vtu test_transport-00006.vtu test/decoupled/2p/test_transport test/decoupled/2p/grids/test_transport.dgf 1e4)
add_test(test_impes util/runTest.sh test/decoupled/2p/test_impes-reference.vtu test_impes-00013.vtu test/decoupled/2p/test_impes 1e7)
add_test(test_dec2p2c util/runTest.sh test/decoupled/2p2c/test_dec2p2c-reference.vtu test_dec2p2c-00021.vtu test/decoupled/2p2c/test_dec2p2c)

add_test(tutorial_coupled tutorial/tutorial_coupled 1 1)
add_test(tutorial_decoupled tutorial/tutorial_decoupled 1)

