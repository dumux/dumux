# -*-cmake-*-
# - Try to find tje DUNE istl library
# Once done this will define:
#  DUNE_istl_FOUND        - system has dune-istl
#  DUNE_istl_INCLUDE_DIR  - incude paths to use dune-istl
#  DUNE_istl_LIBRARIES    - Link these to use dune-istl
INCLUDE(StupidMacros)

StupidSetup("DUNE_istl" "dune-istl" "DUNE")

StupidFindIncludeDir("istl/io.hh")

StupidRequiredLibsFound()
StupidIncludeDirsFound()
StupidCheckFound()


