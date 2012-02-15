find_package(BLAS QUIET REQUIRED)
if(NOT BLAS_FOUND)
  message("BLAS required for SuperLU not found")
  return()
endif(NOT BLAS_FOUND)

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
  set(SUPERLU_FOUND ENABLE_SUPER_LU)
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of SuperLU succeded:\n"
    "Include directory: ${SUPERLU_INCLUDE_DIR}\n"
    "Library directory: ${SUPERLU_LIBRARY}\n\n")
else(SUPERLU_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of SuperLU failed:\n"
    "Include directory: ${SUPERLU_INCLUDE_DIR}\n"
    "Library directory: ${SUPERLU_LIBRARY}\n\n")
endif(SUPERLU_FOUND)
