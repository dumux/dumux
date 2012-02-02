# -*-cmake-*-
# - Try to find the Alberta grid manager
# Once done this will define:
#  Alberta_FOUND        - system has dune-grid
#  Alberta_INCLUDE_DIR  - incude paths to use dune-grid
#  Alberta_LIBRARIES    - Link these to use dune-grid
Include(DumuxMacros)

DumuxSetup("Alberta" "Alberta" "Alberta")

#DumuxAddPathSuffixes("${MyIncludeSuffixes}" "")

DumuxFindIncludeDir("alberta.h")
DumuxFindLibrary("ALBERTA22_0")
DumuxFindLibrary("ALBERTA22_1")
DumuxFindLibrary("alberta_util")

DumuxRequiredLibsFound("ALBERTA22_0" "ALBERTA22_1" "alberta_util")
DumuxIncludeDirsFound()
DumuxCheckFound()
