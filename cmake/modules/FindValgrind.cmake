#
# Module that checks whether Valgrind's header memcheck.h is present
#
# Variables used by this module which you may want to set:
# VALGRIND_ROOT         Path list to search for memcheck.h
#
# Sets the follwing variable:
#
# Valgrind_FOUND          True if Valgrind was found
# VALGRIND_INCLUDE_DIR    Path to Valgrind's include dirs.


# look for header files, only at positions given by the user
find_path(VALGRIND_INCLUDE_DIR
  NAMES "valgrind/memcheck.h"
  PATHS ${VALGRIND_ROOT}
  PATH_SUFFIXES "include"
  NO_DEFAULT_PATH
)
# look for header files, including default paths
find_path(VALGRIND_INCLUDE_DIR
  NAMES "valgrind/memcheck.h"
  PATH_SUFFIXES "include"
)

# handle package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "Valgrind"
  DEFAULT_MSG
  VALGRIND_INCLUDE_DIR
)

# set HAVE_VALGRIND for config.h
set(HAVE_VALGRIND ${Valgrind_FOUND})

# register all Valgrind related flags
if(Valgrind_FOUND)
  dune_register_package_flags(INCLUDE_DIRS "${VALGRIND_INCLUDE_DIR}")
endif()

# text for feature summary
set_package_properties("Valgrind" PROPERTIES
  DESCRIPTION "Memory debugging, memory leak detection, profiling"
  PURPOSE "Identify undefined variables with Memcheck")
