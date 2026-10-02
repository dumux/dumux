# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# .. cmake_module::
#
#    Find the MPI-parallel MUMPS direct solver (double precision)
#
#    The module looks for the header dmumps_c.h and the libraries dmumps, mumps_common and,
#    if present, pord, and checks that a program calling MUMPS links with MPI. Shared MUMPS
#    libraries carry their own dependencies (ScaLAPACK, BLAS, Fortran runtime), so no Fortran
#    compiler is needed. Static MUMPS libraries are only usable if their dependencies are
#    added to MUMPS_EXTRA_LIBRARIES.
#
#    You may set the following variables to modify the
#    behaviour of this module:
#
#    :ref:`MUMPS_ROOT`
#       Installation prefix of MUMPS
#
#    :ref:`MUMPS_EXTRA_LIBRARIES`
#       Additional libraries needed to link MUMPS
#
#    Sets the following variables:
#
#    :code:`MUMPS_FOUND`
#       True if MUMPS was found and a program using it links.
#
#    :code:`MUMPS_INCLUDE_DIRS`
#       Include directories of MUMPS.
#
#    :code:`MUMPS_LIBRARIES`
#       Libraries of MUMPS.
#
#    :code:`MUMPS_VERSION`
#       Version of MUMPS.
#
#    and the imported target :code:`MUMPS::MUMPS`.
#
# .. cmake_variable:: MUMPS_ROOT
#
#   Installation prefix of MUMPS (with the header in include/ and the libraries in lib/),
#   searched before the system paths.
#
# .. cmake_variable:: MUMPS_EXTRA_LIBRARIES
#
#   Additional libraries needed to link MUMPS, e.g. ScaLAPACK, BLAS and the Fortran runtime
#   for static MUMPS libraries.
#
include_guard(GLOBAL)

find_package(MPI QUIET COMPONENTS C)

# MPI-dependent builds of MUMPS may be installed next to the MPI libraries (e.g. on Fedora)
set(_mumps_mpi_hints "")
foreach(_lib IN LISTS MPI_C_LIBRARIES)
  get_filename_component(_dir "${_lib}" DIRECTORY)
  list(APPEND _mumps_mpi_hints "${_dir}")
endforeach()

find_path(MUMPS_INCLUDE_DIR
  NAMES dmumps_c.h
  PATH_SUFFIXES MUMPS mumps)

find_library(MUMPS_DMUMPS_LIBRARY NAMES dmumps HINTS ${_mumps_mpi_hints})
find_library(MUMPS_COMMON_LIBRARY NAMES mumps_common HINTS ${_mumps_mpi_hints})
find_library(MUMPS_PORD_LIBRARY NAMES pord HINTS ${_mumps_mpi_hints})
mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_DMUMPS_LIBRARY MUMPS_COMMON_LIBRARY MUMPS_PORD_LIBRARY)

unset(MUMPS_VERSION)
if(MUMPS_INCLUDE_DIR)
  file(STRINGS "${MUMPS_INCLUDE_DIR}/dmumps_c.h" _mumps_version_line
       REGEX "^#define[ \t]+MUMPS_VERSION[ \t]+\"[0-9.]+\"")
  if(_mumps_version_line)
    string(REGEX REPLACE ".*\"([0-9.]+)\".*" "\\1" MUMPS_VERSION "${_mumps_version_line}")
  endif()
endif()

set(_mumps_libraries ${MUMPS_DMUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY})
if(MUMPS_PORD_LIBRARY)
  list(APPEND _mumps_libraries ${MUMPS_PORD_LIBRARY})
endif()
list(APPEND _mumps_libraries ${MUMPS_EXTRA_LIBRARIES})

# a program calling MUMPS has to link, which fails e.g. for static libraries without their dependencies
set(_mumps_failure "")
unset(MUMPS_LINKS CACHE)
if(MUMPS_INCLUDE_DIR AND MUMPS_DMUMPS_LIBRARY AND MUMPS_COMMON_LIBRARY AND MPI_C_FOUND)
  include(CheckCXXSourceCompiles)
  include(CMakePushCheckState)
  cmake_push_check_state(RESET)
  set(CMAKE_REQUIRED_INCLUDES "${MUMPS_INCLUDE_DIR}")
  set(CMAKE_REQUIRED_LIBRARIES ${_mumps_libraries} MPI::MPI_C)
  set(CMAKE_REQUIRED_QUIET TRUE)
  check_cxx_source_compiles("
    #include <dmumps_c.h>
    int main()
    {
      DMUMPS_STRUC_C id;
      id.job = -1;
      dmumps_c(&id);
      return 0;
    }" MUMPS_LINKS)
  cmake_pop_check_state()

  if(NOT MUMPS_LINKS)
    string(CONCAT _mumps_failure "MUMPS was found (${MUMPS_DMUMPS_LIBRARY}), but a program using it does not link. "
                                 "Static MUMPS libraries need their dependencies in MUMPS_EXTRA_LIBRARIES.")
    message(WARNING "${_mumps_failure}")
  endif()
elseif(NOT MPI_C_FOUND)
  set(_mumps_failure "MUMPS requires MPI.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS
  REQUIRED_VARS MUMPS_DMUMPS_LIBRARY MUMPS_COMMON_LIBRARY MUMPS_INCLUDE_DIR MPI_C_FOUND MUMPS_LINKS
  VERSION_VAR MUMPS_VERSION
  REASON_FAILURE_MESSAGE "${_mumps_failure}")

if(MUMPS_FOUND)
  set(MUMPS_INCLUDE_DIRS "${MUMPS_INCLUDE_DIR}")
  set(MUMPS_LIBRARIES ${_mumps_libraries})
  if(NOT TARGET MUMPS::MUMPS)
    add_library(MUMPS::MUMPS INTERFACE IMPORTED)
    set_target_properties(MUMPS::MUMPS PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${MUMPS_LIBRARIES};MPI::MPI_C")
  endif()
endif()

include(FeatureSummary)
set_package_properties("MUMPS" PROPERTIES
  DESCRIPTION "MUltifrontal Massively Parallel sparse direct Solver"
  URL "https://mumps-solver.org"
  PURPOSE "Direct linear solver backend DirectSolverMumps, sequential and parallel")
