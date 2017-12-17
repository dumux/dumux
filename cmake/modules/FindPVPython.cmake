# .. cmake_module::
#
#    Find the pvpython library
#
#    You may set the following variables to modify the
#    behaviour of this module:
#
#    :ref:`PVPYTHON_ROOT`
#       Path list to search for pvpython.
#
#    Sets the following variables:
#
#    :code:`PVPYTHON_FOUND`
#       True if the pvpython library was found.
#
#    :code:`PVPYTHON_EXECUTABLE`
#       Path to pvpython executable
#
# .. cmake_variable:: PVPYTHON_ROOT
#
#   You may set this variable to have :ref:`FindPVPython" look
#   for the pvpython library in the given path before inspecting
#   system paths.
#

# look for header files, only at positions given by the user
find_program(PVPYTHON_EXECUTABLE
  NAMES pvpython
  PATHS "${PVPYTHON_ROOT}/"
        "${CMAKE_SOURCE_DIR}/../"
        "/usr/bin/"
  PATH_SUFFIXES "pvpython"
  NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "PVPython"
  DEFAULT_MSG
  PVPYTHON_EXECUTABLE
)

# set macros for config.h
set(HAVE_PVPYTHON ${PVPYTHON_FOUND})
set(PVPYTHON_EXECUTABLE ${PVPYTHON_EXECUTABLE})

# text for feature summary
set_package_properties("PVPython" PROPERTIES
  DESCRIPTION "ParaView python client"
  PURPOSE "Extract data over line or time in post-processing")
