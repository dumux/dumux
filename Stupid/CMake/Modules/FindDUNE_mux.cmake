# -*-cmake-*-
# Try to find the DUMUX library
# Once done this will define:
#  DUNE_mux_FOUND        - system has dune-mux
#  DUNE_mux_INCLUDE_DIR  - incude paths to use dune-mux
#  DUNE_mux_LIBRARIES    - Link these to use dune-mux
INCLUDE(StupidMacros)

StupidSetup("DUNE_mux" "dune-mux" "DUNE")

set(MyLibSuffixes 
    "dumux"
    "dumux/hapefunctions"
    "dumux/onedinndgrid")
StupidAddPathSuffixes("" "${MyLibSuffixes}" )

#StupidFindIncludeBaseDir("dumux/material/twophaserelations.hh" "..")
StupidFindIncludeDir("dumux/material/twophaserelations.hh")

StupidFindLibrary("dumux")

StupidRequiredLibsFound("dumux")
StupidIncludeDirsFound()
StupidCheckFound()


