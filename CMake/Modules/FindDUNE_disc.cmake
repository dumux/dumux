# -*-cmake-*-
# - Try to find tje DUNE disc library
# Once done this will define:
#  Dune_disc_FOUND        - system has dune-disc
#  Dune_disc_INCLUDE_DIR  - incude paths to use dune-disc
#  Dune_disc_LIBRARIES    - Link these to use dune-disc
Include(DumuxMacros)

DumuxSetup("DUNE_disc" "dune-disc" "DUNE")

DumuxFindIncludeDir("disc/functions/functions.hh")
DumuxFindLibrary("dunedisc")

DumuxRequiredLibsFound("dunedisc")
DumuxIncludeDirsFound()
DumuxCheckFound()
