# additional macros
include(AddGstatFileLinks)
include(AddInputFileLinks)
include(DumuxDoxygen)
include(DumuxTestMacros)

find_package(GLPK QUIET)
find_package(Gnuplot QUIET)
set(HAVE_GNUPLOT ${GNUPLOT_FOUND})
find_package(Gstat QUIET)
find_package(Gmsh QUIET)
find_package(NLOPT QUIET)
find_package(PTScotch QUIET)
include(AddPTScotchFlags)
find_package(PVPython QUIET)

find_package(Kokkos QUIET)
include(AddKokkosFlags)

# possibly link against TBB
# even if an older version is found
# otherwise we get linker errors
# because of inconsistencies with
# dune-common's TBB setup
find_package(TBB)
include(AddTBBFlags)

# in a second step make sure the
# minimum TBB version required is found
set(DUMUX_MIN_TBB_VERSION 2021)
if(TBB_FOUND)
  if(TBB_VERSION_MAJOR VERSION_LESS DUMUX_MIN_TBB_VERSION)
    find_package(TBB ${DUMUX_MIN_TBB_VERSION})
    # disable TBB manually if required version not found
    if(NOT TBB_FOUND)
      message(STATUS "Disabling TBB since version requirement not satisfied (>= ${DUMUX_MIN_TBB_VERSION}).")
      set(ENABLE_TBB FALSE)
      set(HAVE_TBB FALSE)
    endif()
  endif()
endif()

find_package(OpenMP QUIET)
include(AddOpenMPFlags)

# test if we can use parallel algorithms
check_cxx_symbol_exists(
  "std::execution::par_unseq"
  "execution"
  DUMUX_HAVE_CXX_EXECUTION_POLICY
)

if(DUMUX_HAVE_CXX_EXECUTION_POLICY)
  set(HAVE_CPP_PARALLEL_ALGORITHMS TRUE)
endif()

# setup multithreading backend
if(NOT DUMUX_MULTITHREADING_BACKEND)
  if(TBB_FOUND)
    set(DUMUX_MULTITHREADING_BACKEND "TBB" CACHE STRING "The multithreading backend")
    message(STATUS "Dumux multithreading backed: TBB")
  elseif(OpenMP_FOUND)
    set(DUMUX_MULTITHREADING_BACKEND "OpenMP" CACHE STRING "The multithreading backend")
    message(STATUS "Dumux multithreading backed: OpenMP")
  elseif(Kokkos_FOUND)
    set(DUMUX_MULTITHREADING_BACKEND "Kokkos" CACHE STRING "The multithreading backend")
    message(STATUS "Dumux multithreading backed: Kokkos")
  elseif(DUMUX_HAVE_CXX_EXECUTION_POLICY)
    set(DUMUX_MULTITHREADING_BACKEND "Cpp" CACHE STRING "The multithreading backend")
    message(STATUS "Dumux multithreading backed: Cpp")
  else()
    set(DUMUX_MULTITHREADING_BACKEND "Serial" CACHE STRING "The multithreading backend")
    message(STATUS "Dumux multithreading backed: Serial")
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
  else()
    message(STATUS "Dumux multithreading backed: ${DUMUX_MULTITHREADING_BACKEND}")
  endif()
endif()
