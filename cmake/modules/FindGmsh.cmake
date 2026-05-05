# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# .. cmake_module::
#
#    Find the Gmsh meshing tool
#
#    You may set the following variables to modify the
#    behaviour of this module:
#
#    :ref:`GMSH_ROOT`
#       Path list to search for gmsh in addition to the standard system paths.
#
#    Sets the following variables:
#
#    :code:`Gmsh_FOUND`
#       True if the gmsh executable was found.
#
#    :code:`GMSH_EXECUTABLE`
#       Full path to the gmsh executable.
#
#    :code:`GMSH_VERSION`
#       Version string reported by gmsh (e.g. "4.11.1"), if available.
#
# .. cmake_variable:: GMSH_ROOT
#
#   You may set this variable to have :ref:`FindGmsh` look
#   for the gmsh executable in the given path before inspecting
#   system paths.
#
include_guard(GLOBAL)

find_program(GMSH_EXECUTABLE
  NAMES gmsh
  HINTS "${GMSH_ROOT}"
  PATH_SUFFIXES "bin"
)

# Try to determine the version
if(GMSH_EXECUTABLE)
  execute_process(
    COMMAND "${GMSH_EXECUTABLE}" --version
    OUTPUT_VARIABLE GMSH_VERSION
    ERROR_VARIABLE  GMSH_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
  )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh
  REQUIRED_VARS GMSH_EXECUTABLE
  VERSION_VAR   GMSH_VERSION
)

# text for feature summary
include(FeatureSummary)
set_package_properties("Gmsh" PROPERTIES
  DESCRIPTION "Meshing tool"
  PURPOSE "Generate structured and unstructured grids")
