find_package(BLAS QUIET REQUIRED)
if(NOT BLAS_FOUND AND REQUIRED)
  message("BLAS not found but required for SuperLU")
  return()
endif(NOT BLAS_FOUND AND REQUIRED)

# look for header files
find_path(SUPERLU_INCLUDE_DIR
  NAMES supermatrix.h
  PATH_SUFFIXES "superlu" "include/superlu" "include" "SRC"
)

# look for library
find_library(SUPERLU_LIBRARY
  NAMES "libsuperlu.a" "libsuperlu_4.3.a" "libsuperlu_4.2.a" "libsuperlu_4.1.a" "libsuperlu_4.0.a" "libsuperlu_3.1.a" "libsuperlu_3.0.a"
  PATH_SUFFIXES "lib" "lib64"
)

# check if version is 4.3
include(CheckCSourceCompiles)
set(CMAKE_REQUIRED_INCLUDES ${SUPERLU_INCLUDE_DIR})
set(CMAKE_REQUIRED_LIBRARIES ${SUPERLU_LIBRARY})
CHECK_C_SOURCE_COMPILES("
#include <slu_ddefs.h>
int main(void)
{
  return SLU_DOUBLE;
}"
SUPERLU_MIN_VERSION_4_3)

if(SUPERLU_MIN_VERSION_4_3)
  set(SUPERLU_WITH_VERSION "SuperLU >= 4.3")
else()
  set(SUPERLU_WITH_VERSION "SuperLU <= 4.2, post 2005")
endif(SUPERLU_MIN_VERSION_4_3)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "SuperLU"
  DEFAULT_MSG
  SUPERLU_INCLUDE_DIR
  SUPERLU_LIBRARY
)

mark_as_advanced(SUPERLU_INCLUDE_DIRS SUPERLU_LIBRARIES)

# if both headers and library are found, store results
if(SUPERLU_FOUND)
  set(SUPERLU_INCLUDE_DIRS ${SUPERLU_INCLUDE_DIR})
  set(SUPERLU_LIBRARIES    ${SUPERLU_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of ${SUPERLU_WITH_VERSION} succeded:\n"
    "Include directory: ${SUPERLU_INCLUDE_DIR}\n"
    "Library directory: ${SUPERLU_LIBRARY}\n\n")
else(SUPERLU_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of SuperLU failed:\n"
    "Include directory: ${SUPERLU_INCLUDE_DIR}\n"
    "Library directory: ${SUPERLU_LIBRARY}\n\n")
endif(SUPERLU_FOUND)
