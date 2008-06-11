##############
# general stuff
set(ProjectName            "Stupid")
set(ProjectVersion         "0.0")
set(ProjectMaintainer      "Andreas Lauser")
set(ProjectMaintainerEmail "and at poware dot org")
project(${ProjectName} CXX)
##############


##############
# make sure our own modules will be found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules")
##############


##############
# Find the required packages
FIND_PACKAGE(UG REQUIRED)
FIND_PACKAGE(DUNE_grid REQUIRED)
FIND_PACKAGE(DUNE_disc REQUIRED)
FIND_PACKAGE(DUNE_istl REQUIRED)
FIND_PACKAGE(DUNE_mux REQUIRED)
FIND_PACKAGE(DUNE_common REQUIRED)
FIND_PACKAGE(Boost REQUIRED)
##############


##############
# use this macros in the CMakelists of the subdirectories. 

# -> for TARGET_LINK_LIBRARIES
set(StupidLinkLibraries 
    ${DUNE_grid_LIBRARIES}
    ${DUNE_disc_LIBRARIES}
    ${UG_LIBRARIES}
    ${DUNE_common_LIBRARIES}
    ${Boost_LIBRARIES}
    "boost_program_options")

# -> for LINK_DIRECTORIES
set(StupidLinkDirectories ${Boost_LIBRARY_DIRS})

# -> for INCLUDE_DIRECTORIES
set(StupidIncludeDirectories 
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}
    ${UG_INCLUDE_DIRS}
    ${DUNE_grid_INCLUDE_DIRS}
    ${DUNE_disc_INCLUDE_DIRS}
    ${DUNE_common_INCLUDE_DIRS}
    ${DUNE_istl_INCLUDE_DIRS}
    ${DUNE_mux_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR})
##############


##############
# set appropriate compiler flags for debug/release compilation modes
set( CMAKE_CXX_FLAGS_DEBUG "-g -Wall -DDEBUG -DDUNE_ISTL_WITH_CHECKING" )
set( CMAKE_CXX_FLAGS_RELEASE "-O1 -Wall -DNDEBUG" )

set( CMAKE_C_FLAGS_DEBUG "-g -Wall -DDEBUG -DDUNE_ISTL_WITH_CHECKING" )
set( CMAKE_C_FLAGS_RELEASE "-O1 -Wall -DNDEBUG" )
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
SetConfigHVar(HAVE_UG                  UG_FOUND)
SetConfigHVar(PROJECT_NAME             ProjectName)
SetConfigHVar(PROJECT_VERSION          ProjectVersion)
SetConfigHVar(PROJECT_MAINTAINER       ProjectMaintainer)
SetConfigHVar(PROJECT_MAINTAINER_EMAIL ProjectMaintainerEmail)
configure_file(config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/StupidConfig.hh )
##############


##############
# tell cmake that we've got a few subdirectories. (that's the
# directories where the actual programs are)
add_subdirectory("Problems")
add_subdirectory("Tests")
##############
