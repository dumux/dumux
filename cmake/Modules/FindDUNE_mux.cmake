# -*-cmake-*-
# Try to find the DUMUX library
# Once done this will define:
#  DUNE_mux_FOUND        - system has dune-mux
#  DUNE_mux_INCLUDE_DIR  - incude paths to use dune-mux
#  DUNE_mux_LIBRARIES    - Link these to use dune-mux
INCLUDE(DumuxMacros)

DumuxSetup("DUNE_mux" "dune-mux" "DUNE")

set(MyLibSuffixes 
    "dumux"
    "dumux/shapefunctions"
    "dumux/onedinndgrid")
DumuxAddPathSuffixes("" "${MyLibSuffixes}" )

#DumuxFindIncludeBaseDir("dumux/material/twophaserelations.hh" "..")
DumuxFindIncludeDir("dumux/material/twophaserelations.hh")

DumuxFindLibrary("dumux")
DumuxRequiredLibsFound("dumux")
DumuxIncludeDirsFound()
DumuxCheckFound()


