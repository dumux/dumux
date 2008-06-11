# -*-cmake-*-
# - Try to find tje DUNE fem library
# Once done this will define:
#  Dune_fem_FOUND        - system has dune-fem
#  Dune_fem_INCLUDE_DIR  - incude paths to use dune-fem
#  Dune_fem_LIBRARIES    - Link these to use dune-fem
Include(StupidMacros)

StupidSetup("DUNE_fem" "dune-fem" "DUNE")

StupidFindIncludeDir("fem/misc/entityfunction.hh")
StupidFindLibrary("dunefem")

StupidRequiredLibsFound("dunefem")
StupidIncludeDirsFound()
StupidCheckFound()
