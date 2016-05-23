#
# Module that checks whether gstat has been installed
#
# Sets the following variables:
# GSTAT_FOUND         True if gstat was found
# GSTAT_EXECUTABLE    Path to gstat executable

# look for header files, only at positions given by the user
find_program(GSTAT_EXECUTABLE
  NAMES gstat
  PATHS "${CMAKE_SOURCE_DIR}/../external/gstat/src"
        "usr/bin/"
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
