# checks
include(CheckConstexpr)
# symlink commands from Dune 2.4 for use with Dune 2.3
# include only in case of Dune 2.3
if(("${DUNE_COMMON_VERSION_MAJOR}" STREQUAL "2")
    AND ("${DUNE_COMMON_VERSION_MINOR}" STREQUAL "3"))
  include(CopyOfDuneSymlinkOrCopy)
endif()
# additional macros
include(AddInputFileLinks)
include(DumuxTestMacros)

find_package(Gnuplot)
set(HAVE_GNUPLOT ${GNUPLOT_FOUND})

find_package(Valgrind)
