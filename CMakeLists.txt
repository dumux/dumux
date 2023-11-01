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
option(DUNE_ENABLE_PYTHONBINDINGS "Enable Dune Python bindings" ON)
# this option can be used to turn off Dumux Python bindings independent of Dune
option(DUMUX_ENABLE_PYTHONBINDINGS "Enable Dumux Python bindings" ${DUNE_ENABLE_PYTHONBINDINGS})

# if Python bindings are enabled with need to build shared libraries
# we set the default depending on whether the Python bindings are enabled
option(BUILD_SHARED_LIBS "Build using shared libraries" ${DUMUX_ENABLE_PYTHONBINDINGS})
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
