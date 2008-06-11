# -*-cmake-*-
# - Try to find tje DUNE common library
# Once done this will define:
#  DUNE_common_FOUND        - system has dune-common
#  DUNE_common_INCLUDE_DIR  - incude paths to use dune-common
#  DUNE_common_LIBRARIES    - Link these to use dune-common
INCLUDE(StupidMacros)

StupidSetup("DUNE_common" "dune-common" "DUNE")

StupidFindIncludeDir("dune/common/misc.hh")
StupidFindLibrary("dunecommon")

StupidRequiredLibsFound("dunecommon")
StupidIncludeDirsFound()
StupidCheckFound()


