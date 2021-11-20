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
dune_enable_all_packages(MODULE_LIBRARIES dumux_fmt dumux_parameters)

# enforce C++-17
dune_require_cxx_standard(MODULE "DuMuX" VERSION 17)

add_subdirectory(cmake/modules)
add_subdirectory(doc)
add_subdirectory(dumux)
add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(examples EXCLUDE_FROM_ALL)

# from Dune 2.9 on Python bindings are enabled per default
if(${dune-common_VERSION} VERSION_GREATER_EQUAL 2.9)
  add_subdirectory(python)
else()
  # with Dune 2.8, only if Python bindings are enabled
  if(DUNE_ENABLE_PYTHONBINDINGS)
    add_subdirectory(python)
    dune_python_install_package(PATH "python")
  endif()
endif()

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
