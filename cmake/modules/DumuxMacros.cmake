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
# beacuse of inconsistencies with
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
