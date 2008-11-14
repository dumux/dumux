##############
# general stuff
cmake_minimum_required(VERSION 2.6)
set(ProjectName            "DuMuX")
set(ProjectVersion         "0.0.1")
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
FIND_PACKAGE(DUNE_disc REQUIRED)
FIND_PACKAGE(DUNE_istl REQUIRED)
#FIND_PACKAGE(DUNE_mux REQUIRED)
FIND_PACKAGE(DUNE_common REQUIRED)
FIND_PACKAGE(Boost REQUIRED)
##############

##############
# Find the optional packages
FIND_PACKAGE(MPI)
FIND_PACKAGE(UG)
FIND_PACKAGE(ALUGrid)
#FIND_PACKAGE(Alberta)
##############

##############
# Find the required include files
INCLUDE (CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("malloc.h" HAVE_MALLOC_H)
CHECK_INCLUDE_FILE_CXX("tr1/array" HAVE_TR1_ARRAY)
##############


##############
# use this macros in the CMakelists of the subdirectories. 

# -> for TARGET_LINK_LIBRARIES
set(DumuxLinkLibraries 
    "dumux" # the DUMUX library. CMake replaces this by the proper file location
    ${DUNE_grid_LIBRARIES}
    ${DUNE_disc_LIBRARIES}
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
    ${DUNE_disc_INCLUDE_DIRS}
    ${DUNE_common_INCLUDE_DIRS}
    ${DUNE_istl_INCLUDE_DIRS}
    ${DUNE_mux_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR})
##############


##############
# set appropriate compiler flags for debug/release compilation modes
set( CMAKE_CXX_FLAGS_DEBUG "-g -Wall" )
set( CMAKE_C_FLAGS_DEBUG "-g -Wall" )

set( CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -funroll-loops -finline-functions" )
set( CMAKE_C_FLAGS_RELEASE "-O3 -Wall -funroll-loops -finline-functions" )

if(CMAKE_BUILD_TYPE STREQUAL "debug")
  # debug mode
  add_definitions(-DDEBUG -DDUNE_DEVEL_MODE=1 -DDUNE_ISTL_WITH_CHECKING)
else(CMAKE_BUILD_TYPE STREQUAL "debug")
  # Release mode
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
SetConfigHVar(HAVE_DUNE_DISC           DUNE_disc_FOUND)
SetConfigHVar(HAVE_DUNE_ISTL           DUNE_istl_FOUND)
SetConfigHVar(HAVE_MPI                 MPI_FOUND)
SetConfigHVar(PROJECT_NAME             ProjectName)
SetConfigHVar(PROJECT_VERSION          ProjectVersion)
SetConfigHVar(PROJECT_MAINTAINER       ProjectMaintainer)
SetConfigHVar(PROJECT_MAINTAINER_EMAIL ProjectMaintainerEmail)
##############

##############
# adapt build system to detected packages

# deal with UG
if (UG_FOUND)
  set(DumuxLinkLibraries ${DumuxLinkLibraries} 
      ${UG_LIBRARIES})
  set(DumuxIncludeDirectories ${DumuxIncludeDirectories} 
      ${UG_INCLUDE_DIRS})
endif(UG_FOUND)
SetConfigHVar(HAVE_UG UG_FOUND)

# deal with ALUGrid
if (ALUGrid_FOUND)
  set(DumuxLinkLibraries ${DumuxLinkLibraries} 
      ${ALUGrid_LIBRARIES})
  set(DumuxIncludeDirectories ${DumuxIncludeDirectories} 
      ${ALUGrid_INCLUDE_DIRS})
endif(ALUGrid_FOUND)
SetConfigHVar(HAVE_ALUGRID ALUGrid_FOUND)

# deal with Alberta
if (Alberta_FOUND)
  set(DumuxLinkLibraries ${DumuxLinkLibraries} 
      ${Alberta_LIBRARIES})
  set(DumuxIncludeDirectories ${DumuxIncludeDirectories} 
      ${Alberta_INCLUDE_DIRS})
endif(Alberta_FOUND)
SetConfigHVar(HAVE_ALBERTA Alberta_FOUND)
##############

# actually write the config.h file to disk 
configure_file(config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h )
add_definitions(-DHAVE_CONFIG_H)

##############
# tell cmake that we've got a few subdirectories. (that's the
# directories where the actual programs are)
add_subdirectory("dumux")

add_subdirectory("test")
add_subdirectory("benchmarkproblems")
##############
