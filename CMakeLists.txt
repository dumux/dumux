# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.18)

project("dumux" C CXX)

if(NOT (dune-common_DIR
        OR dune-common_ROOT
        OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME}
      dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

# find dune-common and set the module path
find_package(dune-common)

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

# start a dune project with information from dune.module
dune_project()

# external dependencies
include(FetchContent)
if (EXISTS ${CMAKE_SOURCE_DIR}/deps/gridformat/CMakeLists.txt)
  FetchContent_Declare(gridformat SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/gridformat")
  FetchContent_MakeAvailable(gridformat)
  set(DUMUX_HAVE_GRIDFORMAT true)
  message(STATUS "Using gridformat from ${gridformat_SOURCE_DIR}")
endif()

# Create the dumux library target exported as Dumux::Dumux
include(DumuxAddLibrary)
dumux_add_library(dumux NAMESPACE Dumux:: EXPORT_NAME Dumux)

# enforce C++-20
target_compile_features(dumux PUBLIC cxx_std_20)

if (DUMUX_HAVE_GRIDFORMAT)
  target_compile_definitions(dumux PUBLIC DUMUX_HAVE_GRIDFORMAT=1)
  target_link_libraries(dumux INTERFACE gridformat::gridformat)
endif()

if (NOT DUMUX_HAVE_STD_FORMAT)
  dumux_add_library(dumux_format NAMESPACE Dumux:: EXPORT_NAME Format)
  target_compile_features(dumux_format PUBLIC cxx_std_20)
  target_link_libraries(dumux INTERFACE dumux_format)
endif()

# Enable all dune registered packages on the dumux library.
dune_target_enable_all_packages(dumux)

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
set(DUNE_CUSTOM_PKG_CONFIG_SECTION "")
# make sure downstream projects can find gridformat too
if (DUMUX_HAVE_GRIDFORMAT)
  string(JOIN "\n" DUNE_CUSTOM_PKG_CONFIG_SECTION ${DUNE_CUSTOM_PKG_CONFIG_SECTION}
    "# Make sure gridformat build in dumux is found"
    "include(CMakeFindDependencyMacro)"
    "find_dependency(gridformat PATHS ${gridformat_BINARY_DIR}/cmake)"
  )
endif()

finalize_dune_project()
