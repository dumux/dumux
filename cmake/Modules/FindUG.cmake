# -*-cmake-*-
# - Try to find the UG grid manager
# Once done this will define:
#  UG_FOUND        - system has dune-grid
#  UG_INCLUDE_DIR  - incude paths to use dune-grid
#  UG_LIBRARIES    - Link these to use dune-grid
Include(DumuxMacros)


DumuxSetup("UG" "UG" "UG")

set(MyIncludeSuffixes 
    "gm"
    "np")
set(MyLibSuffixes 
    "np" 
    "np/field"
    "np/udm"
    "np/procs"
    "np/amglib"
    "np/algebra"
    "dev"
    "dev/rif"
    "dev/sif"
    "dev/ps"
    "dev/meta"
    "dev/xif"
    "dev/ppm"
    "dom/lgm/ngin2d"
    "dom/lgm/ngin"
    "dom/lgm"
    "dom/std"
    "parallel/dddif"
    "parallel/util"
    "graphics"
    "graphics/uggraph"
    "graphics/grape"
    "low"
    "gm"
    "gm/gg2"
    "gm/gg3"
    "ui")
set(MyUgLibs 
  "ugS2"
  "ugS3"
  "devS"
)

DumuxAddPathSuffixes("${MyIncludeSuffixes}" "${MyLibSuffixes}" )

DumuxFindIncludeDir("ugm.h")
DumuxFindExtraIncludeDir("GM_H" "gm.h")

foreach(tmp ${MyUgLibs})
  DumuxFindLibrary(${tmp})
endforeach(tmp)

DumuxRequiredLibsFound(${MyUgLibs})
DumuxIncludeDirsFound()
DumuxCheckFound()
