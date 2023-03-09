# set up project
project("dumux" C CXX)

# general stuff
cmake_minimum_required(VERSION 3.13)

if(NOT (dune-common_DIR
        OR dune-common_ROOT
        OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME}
      dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common)

list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()
if (DUMUX_HAVE_STD_FORMAT)
  dune_enable_all_packages(MODULE_LIBRARIES dumux_parameters)
else()
  dune_enable_all_packages(MODULE_LIBRARIES dumux_fmt dumux_parameters)
endif()

# enforce C++-17
dune_require_cxx_standard(MODULE "DuMuX" VERSION 17)

# option for enabling Python bindings
option(DUMUX_ENABLE_PYTHONBINDINGS "Enable Python bindings" ${DUNE_ENABLE_PYTHONBINDINGS})

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
