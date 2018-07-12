# .. cmake_module::
#
#    Find the Gmsh meshing tool
#
#    You may set the following variables to modify the
#    behaviour of this module:
#
#    :ref:`GMSH_ROOT`
#       Path list to search for gmsh.
#
#    Sets the following variables:
#
#    :code:`gmsh_FOUND`
#       True if the gmsh library was found.
#
#    :code:`GMSH_EXECUTABLE`
#       Path to gmsh executable
#
# .. cmake_variable:: GMSH_ROOT
#
#   You may set this variable to have :ref:`FindGmsh` look
#   for the gmsh library in the given path before inspecting
#   system paths.
#

# look for header files, only at positions given by the user
find_program(GMSH_EXECUTABLE
  NAMES gmsh
  PATHS "${GMSH_ROOT}"
        "${CMAKE_SOURCE_DIR}/../"
        "/usr/bin/"
  PATH_SUFFIXES "src" "external/gmsh/src" "gmsh/src" "gmsh"
  NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "gmsh"
  DEFAULT_MSG
  GMSH_EXECUTABLE
)

# set macros for config.h
set(HAVE_GMSH ${gmsh_FOUND})
set(GMSH_EXECUTABLE ${GMSH_EXECUTABLE})

# text for feature summary
set_package_properties("Gmsh" PROPERTIES
  DESCRIPTION "Meshing tool"
  PURPOSE "Generate structured and unstructured grids")
