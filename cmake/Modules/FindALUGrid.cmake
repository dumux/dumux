# -*-cmake-*-
# - Try to find the UG grid manager
# Once done this will define:
#  ALUGrid_FOUND        - system has dune-grid
#  UG_INCLUDE_DIR  - incude paths to use dune-grid
#  UG_LIBRARIES    - Link these to use dune-grid
Include(DumuxMacros)

DumuxSetup("ALUGrid" "ALUGrid" "ALUGrid")

set(MyIncludeSuffixes 
    "include/serial"
    "include/parallel"
    "include/duneinterface")

DumuxAddPathSuffixes("${MyIncludeSuffixes}" "")

DumuxFindIncludeDir("alugrid_2d.h")
DumuxFindExtraIncludeDir("ALU_SERIAL" "serialize.h")
DumuxFindExtraIncludeDir("ALU_PARALLEL" "gitter_pll_impl.h")
DumuxFindExtraIncludeDir("ALU_DUNE" "gitter_dune_impl.h")
DumuxFindLibrary("alugrid")

DumuxRequiredLibsFound("alugrid")
DumuxIncludeDirsFound()
DumuxCheckFound()
