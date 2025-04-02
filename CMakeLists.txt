# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.22)

project(dumux LANGUAGES C CXX)

# if no path to dune-common is given, we assume
# that dune-common uses the same build directory as dumux
if(NOT (dune-common_DIR OR dune-common_ROOT OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
  string(REPLACE ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR ${PROJECT_BINARY_DIR})
endif()

find_package(dune-common REQUIRED)

# make cmake find dune-common's and dumux's cmake modules
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

# include the dune macros
include(DuneMacros)

# option for enabling Python bindings
# The default is OFF, however, note that this does not affect upstream modules
# since they are configured already before Dumux is configured. We set this option
# to have a default for DUMUX_ENABLE_PYTHONBINDINGS that is consistent
# with the Dune configuration.
option(DUNE_ENABLE_PYTHONBINDINGS "Enable Dune Python bindings" OFF)
# this option can be used to turn off Dumux Python bindings independent of Dune
option(DUMUX_ENABLE_PYTHONBINDINGS "Enable Dumux Python bindings" ${DUNE_ENABLE_PYTHONBINDINGS})

if (DUMUX_ENABLE_PYTHONBINDINGS AND NOT DUNE_ENABLE_PYTHONBINDINGS)
  message(FATAL_ERROR "Using the Python bindings requires DUNE_ENABLE_PYTHONBINDINGS=ON")
endif()

# If Python bindings are enabled we need to build shared libraries.
# Note that when building Dumux with shared libraries, we should also build all
# upstream Dune libraries as shared libraries. In the case of mixing
# upstream static and downstream shared libraries which depend on the
# upstream libraries (e.g. libdumux), there is a risk of symbols existing
# multiple times which can lead to nasty memory leaks
# (e.g. for global variables such as Dune's stdstreams).
# The option BUILD_SHARED_LIBS should be set (e.g. via the cmake.opts file)
# and used when configuring and compiling upstream modules (e.g. with dunecontrol).
if (DUMUX_ENABLE_PYTHONBINDINGS AND NOT BUILD_SHARED_LIBS)
  message(FATAL_ERROR "Using the Python bindings requires BUILD_SHARED_LIBS=ON")
endif()

# deactivate global include-directories
if (dune-common_VERSION VERSION_GREATER_EQUAL 2.11)
  dune_policy(SET DP_DEFAULT_INCLUDE_DIRS NEW)
  dune_policy(SET DP_TEST_ADD_ALL_FLAGS NEW)
endif()

# start a dune project with information from dune.module
dune_project()

# external dependencies
include(FetchContent)

option(DUMUX_ENABLE_GRIDFORMAT "Enable gridformat I/O integration" OFF)
if (DUMUX_ENABLE_GRIDFORMAT)
  find_package(gridformat QUIET)
  if (gridformat_FOUND)
    set(DUMUX_HAVE_GRIDFORMAT true)
    set(DUMUX_GRIDFORMAT_HINT_DIR ${gridformat_DIR})
  elseif (EXISTS ${CMAKE_SOURCE_DIR}/deps/gridformat/CMakeLists.txt)
    FetchContent_Declare(gridformat
      SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/gridformat"
      OVERRIDE_FIND_PACKAGE ON
    )
    FetchContent_MakeAvailable(gridformat)
    set(DUMUX_HAVE_GRIDFORMAT true)
    set(DUMUX_GRIDFORMAT_HINT_DIR ${gridformat_BINARY_DIR}/cmake)
  else ()
    FetchContent_Declare(gridformat
      GIT_REPOSITORY https://github.com/dglaeser/gridformat.git
      GIT_TAG 6e1f1963dd69443dc8da8fd14423d4f1d057654e
      GIT_PROGRESS true
      GIT_SUBMODULES_RECURSE ON
      OVERRIDE_FIND_PACKAGE ON
    )
    FetchContent_MakeAvailable(gridformat)
    set(DUMUX_GRIDFORMAT_HINT_DIR ${gridformat_BINARY_DIR}/cmake)
    set(DUMUX_HAVE_GRIDFORMAT true)
  endif()
endif()

find_package(TBB QUIET)
find_package(Kokkos QUIET)
find_package(OpenMP QUIET)

# Create the dumux library target exported as Dumux::Dumux
include(DumuxAddLibrary)
dumux_add_library(dumux
  NAMESPACE Dumux::
  EXPORT_NAME Dumux
  LINK_LIBRARIES
    Dune::Common Dune::Geometry Dune::Grid
    Dune::LocalFunctions Dune::ISTL
)

# enforce C++-20
target_compile_features(dumux PUBLIC cxx_std_20)

# check multithreading features
include(DumuxMultithreadingBackend)

if (DUMUX_HAVE_GRIDFORMAT)
  target_compile_definitions(dumux PUBLIC DUMUX_HAVE_GRIDFORMAT=1)
  target_link_libraries(dumux INTERFACE gridformat::gridformat)
endif()

if (TBB_FOUND)
  target_compile_definitions(dumux PUBLIC DUMUX_HAVE_TBB=1)
  target_link_libraries(dumux INTERFACE TBB::tbb)
endif()

if (Kokkos_FOUND)
  target_compile_definitions(dumux PUBLIC DUMUX_HAVE_KOKKOS=1)
  target_link_libraries(dumux INTERFACE Kokkos::kokkos)
endif()

if (OpenMP_FOUND)
  target_compile_definitions(dumux PUBLIC DUMUX_HAVE_OPENMP=1)
  target_link_libraries(dumux INTERFACE OpenMP::OpenMP_CXX)
endif()

if (NOT DUMUX_HAVE_STD_FORMAT)
  dumux_add_library(dumux_format NAMESPACE Dumux:: EXPORT_NAME Format)
  target_compile_features(dumux_format PUBLIC cxx_std_20)
  if (dune-common_VERSION VERSION_GREATER_EQUAL 2.11)
    dune_default_include_directories(dumux_format PUBLIC)
  endif()
  target_link_libraries(dumux INTERFACE dumux_format)
endif()

if (dune-common_VERSION VERSION_GREATER_EQUAL 2.11)
  dune_default_include_directories(dumux PUBLIC)
endif()

if(PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME)
  include(FeatureSummary)
  add_feature_info("Grid I/O via gridformat" DUMUX_HAVE_GRIDFORMAT "writing and reading grid data supporting many different VTK formats")
  set(DUMUX_MULTITHREADING_SUPPORTED true)
  if (DUMUX_MULTITHREADING_BACKEND STREQUAL "Serial")
    set(DUMUX_MULTITHREADING_SUPPORTED false)
  endif()
  add_feature_info("Shared-memory parallelism" DUMUX_MULTITHREADING_SUPPORTED "multithreaded loops, assembly and solvers if available")
  add_feature_info("Distributed-memory parallelism" MPI_FOUND "MPI parallel programs")
  add_feature_info("Python bindings for DuMux" DUMUX_ENABLE_PYTHONBINDINGS "using DuMux via its Python API (experimental)")
endif()

set(DUMUX_RUNTEST ${CMAKE_SOURCE_DIR}/bin/testing/dumux_runtest.py)

add_subdirectory(bin)
add_subdirectory(cmake/modules)
add_subdirectory(doc)
add_subdirectory(dumux)
add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(examples EXCLUDE_FROM_ALL)

# only add bindings folder if Python bindings are enabled
if(DUMUX_ENABLE_PYTHONBINDINGS)
  add_subdirectory(python)
endif()

# finalize the dune project, e.g. generating config.h etc.
# this will also generate the dumux-config.cmake file which is
# used to make Dumux available to other projects
set(DUNE_CUSTOM_PKG_CONFIG_SECTION "include(CMakeFindDependencyMacro)")
# make sure downstream projects can find gridformat too
if (DUMUX_HAVE_GRIDFORMAT)
  string(JOIN "\n" DUNE_CUSTOM_PKG_CONFIG_SECTION ${DUNE_CUSTOM_PKG_CONFIG_SECTION}
    "# Make sure gridformat build in dumux is found"
    "find_dependency(gridformat HINTS ${DUMUX_GRIDFORMAT_HINT_DIR})"
  )
endif()
# make sure downstream projects can find TBB too
if (TBB_FOUND)
  string(JOIN "\n" DUNE_CUSTOM_PKG_CONFIG_SECTION ${DUNE_CUSTOM_PKG_CONFIG_SECTION}
    "# Make sure TBB is found"
    "find_dependency(TBB)"
  )
endif()
# make sure downstream projects can find Kokkos too
if (Kokkos_FOUND)
  string(JOIN "\n" DUNE_CUSTOM_PKG_CONFIG_SECTION ${DUNE_CUSTOM_PKG_CONFIG_SECTION}
    "# Make sure Kokkos is found"
    "find_dependency(Kokkos)"
  )
endif()
# make sure downstream projects can find OpenMP too
if (OpenMP_FOUND)
  string(JOIN "\n" DUNE_CUSTOM_PKG_CONFIG_SECTION ${DUNE_CUSTOM_PKG_CONFIG_SECTION}
    "# Make sure OpenMP is found"
    "find_dependency(OpenMP)"
  )
endif()
string(JOIN "\n" DUNE_CUSTOM_PKG_CONFIG_SECTION ${DUNE_CUSTOM_PKG_CONFIG_SECTION}
  "# Propagate the multithreading backend"
  "set(DUMUX_MULTITHREADING_BACKEND ${DUMUX_MULTITHREADING_BACKEND})"
)

finalize_dune_project()

if(ENABLE_HEADERCHECK)
  get_property(dumux_header_target_list_ TARGET headercheck PROPERTY MANUALLY_ADDED_DEPENDENCIES)
  foreach(header_target_ ${dumux_header_target_list_})
    target_link_libraries(${header_target_} PRIVATE Dumux::Dumux)
  endforeach()
endif()
