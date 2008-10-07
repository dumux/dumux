# -*-cmake-*-
# - Try to find tje DUNE fem library
# Once done this will define:
#  Dune_fem_FOUND        - system has dune-fem
#  Dune_fem_INCLUDE_DIR  - incude paths to use dune-fem
#  Dune_fem_LIBRARIES    - Link these to use dune-fem
Include(DumuxMacros)

DumuxSetup("DUNE_fem" "dune-fem" "DUNE")

DumuxFindIncludeDir("fem/misc/entityfunction.hh")
DumuxFindLibrary("dunefem")

DumuxRequiredLibsFound("dunefem")
DumuxIncludeDirsFound()
DumuxCheckFound()
