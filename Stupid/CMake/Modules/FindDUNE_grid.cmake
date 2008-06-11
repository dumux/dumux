# -*-cmake-*-
# - Try to find tje DUNE grid library
# Once done this will define:
#  Dune_grid_FOUND        - system has dune-grid
#  Dune_grid_INCLUDE_DIR  - incude paths to use dune-grid
#  Dune_grid_LIBRARIES    - Link these to use dune-grid
INCLUDE(StupidMacros)

StupidSetup("DUNE_grid" "dune-grid" "DUNE")

StupidFindIncludeDir("dune/grid/sgrid.hh")
StupidFindLibrary("dunegrid")

StupidRequiredLibsFound("dunegrid")
StupidIncludeDirsFound()
StupidCheckFound()
