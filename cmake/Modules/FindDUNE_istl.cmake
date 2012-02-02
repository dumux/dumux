# -*-cmake-*-
# - Try to find tje DUNE istl library
# Once done this will define:
#  DUNE_istl_FOUND        - system has dune-istl
#  DUNE_istl_INCLUDE_DIR  - incude paths to use dune-istl
#  DUNE_istl_LIBRARIES    - Link these to use dune-istl
INCLUDE(DumuxMacros)

DumuxSetup("DUNE_istl" "dune-istl" "DUNE")

DumuxFindIncludeDir("dune/istl/io.hh")

DumuxRequiredLibsFound()
DumuxIncludeDirsFound()
DumuxCheckFound()


