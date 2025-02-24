# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# test if we can use parallel algorithms
include(CheckCXXSymbolExists)
check_cxx_symbol_exists(
  "std::execution::par_unseq"
  "execution"
  DUMUX_HAVE_CXX_EXECUTION_POLICY
)

if(DUMUX_HAVE_CXX_EXECUTION_POLICY)
  set(DUMUX_HAVE_CPP_PARALLEL_ALGORITHMS TRUE)
endif()

# setup multithreading backend
if(NOT DUMUX_MULTITHREADING_BACKEND)
  if(TBB_FOUND)
    set(DUMUX_MULTITHREADING_BACKEND "TBB" CACHE STRING "The multithreading backend")
  elseif(OpenMP_FOUND)
    set(DUMUX_MULTITHREADING_BACKEND "OpenMP" CACHE STRING "The multithreading backend")
  elseif(Kokkos_FOUND)
    set(DUMUX_MULTITHREADING_BACKEND "Kokkos" CACHE STRING "The multithreading backend")
  elseif(DUMUX_HAVE_CXX_EXECUTION_POLICY)
    set(DUMUX_MULTITHREADING_BACKEND "Cpp" CACHE STRING "The multithreading backend")
  else()
    set(DUMUX_MULTITHREADING_BACKEND "Serial" CACHE STRING "The multithreading backend")
  endif()

# abort if a multithreading backend has been manually selected
# but it is not available
else()
  if(DUMUX_MULTITHREADING_BACKEND STREQUAL "TBB" AND NOT TBB_FOUND)
    message(FATAL_ERROR "Selected TBB as Dumux multithreading backed but TBB has not been found")
  elseif(DUMUX_MULTITHREADING_BACKEND STREQUAL "OpenMP" AND NOT OpenMP_FOUND)
    message(FATAL_ERROR "Selected OpenMP as Dumux multithreading backed but OpenMP has not been found")
  elseif(DUMUX_MULTITHREADING_BACKEND STREQUAL "Kokkos" AND NOT Kokkos_FOUND)
    message(FATAL_ERROR "Selected Kokkos as Dumux multithreading backed but Kokkos has not been found")
  elseif(DUMUX_MULTITHREADING_BACKEND STREQUAL "Cpp" AND NOT DUMUX_HAVE_CXX_EXECUTION_POLICY)
    message(FATAL_ERROR "Selected Cpp as Dumux multithreading backed but your compiler does not implement parallel STL")
  endif()
endif()

message(STATUS "Dumux multithreading backend: ${DUMUX_MULTITHREADING_BACKEND}")
