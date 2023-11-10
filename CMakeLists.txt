# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.16)

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
if (DUMUX_HAVE_STD_FORMAT)
  dune_enable_all_packages(MODULE_LIBRARIES dumux)
else()
  dune_enable_all_packages(MODULE_LIBRARIES dumux dumux_fmt)
endif()

# enforce C++-17
target_compile_features(dumux PUBLIC cxx_std_17)

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
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
