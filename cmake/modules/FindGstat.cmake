# .. cmake_module::
#
#    Find the gstat geostatistic library
#
#    You may set the following variables to modify the
#    behaviour of this module:
#
#    :ref:`GSTAT_ROOT`
#       Path list to search for gstat.
#
#    Sets the following variables:
#
#    :code:`GSTAT_FOUND`
#       True if the gstat library was found.
#
#    :code:`GSTAT_EXECUTABLE`
#       Path to gstat executable
#
# .. cmake_variable:: GSTAT_ROOT
#
#   You may set this variable to have :ref:`FindGstat` look
#   for the gstat library in the given path before inspecting
#   system paths.
#

# look for header files, only at positions given by the user
find_program(GSTAT_EXECUTABLE
  NAMES gstat
  PATHS "${GSTAT_ROOT}"
        "${CMAKE_SOURCE_DIR}/../"
        "/usr/bin/"
  PATH_SUFFIXES "src" "external/gstat/src" "gstat/src" "gstat"
  NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "Gstat"
  DEFAULT_MSG
  GSTAT_EXECUTABLE
)

# set macros for config.h
set(HAVE_GSTAT ${GSTAT_FOUND})
set(GSTAT_EXECUTABLE ${GSTAT_EXECUTABLE})

# text for feature summary
set_package_properties("Gstat" PROPERTIES
  DESCRIPTION "Geostatistic library"
  PURPOSE "Generate random permeability and porosity fields")
