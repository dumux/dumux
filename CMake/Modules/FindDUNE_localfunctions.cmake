# -*-cmake-*-
# - Try to find the DUNE localfunctions library
# Once done this will define:
#  DUNE_localfunctions_FOUND        - system has dune-localfunctions
#  DUNE_localfunctions_INCLUDE_DIR  - incude paths to use dune-localfunctions
#  DUNE_localfunctions_LIBRARIES    - Link these to use dune-localfunctions
INCLUDE(DumuxMacros)

DumuxSetup("DUNE_localfunctions" "dune-localfunctions" "DUNE")

DumuxFindIncludeDir("dune/localfunctions/lagrange.hh")

DumuxRequiredLibsFound()
DumuxIncludeDirsFound()
DumuxCheckFound()


