
#
# Module that checks whether NLOPT is available and usable.
#
# Variables used by this module which you may want to set:
# NLOPT_ROOT         Path list to search for NLOPT
#
# Sets the follwing variable:
#
# NLOPT_FOUND           True if NLOPT available and usable.
# NLOPT_INCLUDE_DIRS    Path to the NLOPT include dirs.
# NLOPT_LIBRARIES       Name to the NLOPT library.
#

# look for header files, only at positions given by the user
find_path(NLOPT_INCLUDE_DIR
  NAMES nlopt.h
  PATHS ${NLOPT_PREFIX} ${NLOPT_ROOT}
  PATH_SUFFIXES "nlopt" "include/nlopt" "include" "SRC" "src"
  NO_DEFAULT_PATH
)

# look for header files, including default paths
find_path(NLOPT_INCLUDE_DIR
  NAMES nlopt.h
  PATH_SUFFIXES "nlopt" "include/nlopt" "include" "SRC" "src"
)

# look for library, only at positions given by the user
find_library(NLOPT_LIBRARY
  NAMES "nlopt" "libnlopt"
  PATHS ${NLOPT_PREFIX} ${NLOPT_ROOT}
  PATH_SUFFIXES "lib" "lib32" "lib64" "libnlopt"
  NO_DEFAULT_PATH
)

# look for library files, including default paths
find_library(NLOPT_LIBRARY
  NAMES "nlopt" "libnlopt"
  PATH_SUFFIXES "lib" "lib32" "lib64" "libnlopt"
)

# check version specific macros
include(CheckCSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state()

# we need if clauses here because variable is set variable-NOTFOUND
# if the searches above were not successful
# Without them CMake print errors like:
# "CMake Error: The following variables are used in this project, but they are set to NOTFOUND.
# Please set them or make sure they are set and tested correctly in the CMake files:"
#
if(NLOPT_INCLUDE_DIR)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${NLOPT_INCLUDE_DIR})
endif(NLOPT_INCLUDE_DIR)
if(NLOPT_LIBRARY)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${NLOPT_LIBRARY})
endif(NLOPT_LIBRARY)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "NLOPT"
  DEFAULT_MSG
  NLOPT_INCLUDE_DIR
  NLOPT_LIBRARY
)

mark_as_advanced(NLOPT_INCLUDE_DIR NLOPT_LIBRARY)

# if both headers and library are found, store results
if(NLOPT_FOUND)
  set(NLOPT_INCLUDE_DIRS ${NLOPT_INCLUDE_DIR})
  set(NLOPT_LIBRARIES    ${NLOPT_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of NLOPT succeeded:\n"
    "Include directory: ${NLOPT_INCLUDE_DIRS}\n"
    "Library directory: ${NLOPT_LIBRARIES}\n\n")
  set(NLOPT_DUNE_COMPILE_FLAGS "-I${NLOPT_INCLUDE_DIRS}"
    CACHE STRING "Compile flags used by DUNE when compiling NLOPT programs")
  set(NLOPT_DUNE_LIBRARIES ${NLOPT_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking NLOPT programs")
else(NLOPT_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining location of NLOPT failed:\n"
    "Include directory: ${NLOPT_INCLUDE_DIRS}\n"
    "Library directory: ${NLOPT_LIBRARIES}\n\n")
endif(NLOPT_FOUND)

# set HAVE_NLOPT for config.h
set(HAVE_NLOPT ${NLOPT_FOUND})

# register all NLOPT related flags
if(NLOPT_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_NLOPT=1"
                              LIBRARIES "${NLOPT_DUNE_LIBRARIES}"
                              INCLUDE_DIRS "${NLOPT_INCLUDE_DIRS}")
endif()
