# -*-cmake-*-
# - Try to find the UG grid manager
# Once done this will define:
#  ALUGrid_FOUND        - system has dune-grid
#  UG_INCLUDE_DIR  - incude paths to use dune-grid
#  UG_LIBRARIES    - Link these to use dune-grid
Include(DumuxMacros)

DumuxSetup("ALUGrid" "ALUGrid" "ALUGrid")

# set(MyIncludeSuffixes 
#     "gm"
#     "np")
# set(MyLibSuffixes 
#     "np" 
#     "np/field"
#     "np/udm"
#     "np/procs"
#     "np/amglib"
#     "np/algebra"
#     "dev"
#     "dev/rif"
#     "dev/sif"
#     "dev/ps"
#     "dev/meta"
#     "dev/xif"
#     "dev/ppm"
#     "dom/lgm/ngin2d"
#     "dom/lgm/ngin"
#     "dom/lgm"
#     "dom/std"
#     "parallel/dddif"
#     "parallel/util"
#     "graphics"
#     "graphics/uggraph"
#     "graphics/grape"
#     "low"
#     "gm"
#     "gm/gg2"
#     "gm/gg3"
#     "ui")
# set(MyUgLibs 
# #  "algebra2"
# #  "amg2"
# #  "dddif2"
# #  "domL2"
# #  "domS2"
# #  "field2"
# #  "gg2"
# #  "graphics2"
# #  "low2"
# #  "np2"
# #  "ug_gm2"
# #  "procs2"
# #  "udm2"
# #  "uggrape2"
# #  "uggraph2"
# #  "ugL2"
# #  "ugui2"
# #  "algebra3"
# #  "amg3"
# #  "dddif3"
# #  "device_meta"
# #  "device_ppm"
# #  "device_ps"
# #  "devices"
# #  "devR"
#   "ugS2"
#   "ugS3"
#   "devS"
# #  "devX"
# #  "domL3"
# #  "field3"
# #  "gg3"
# #  "graphics3"
# #  "low"
# #  "low3"
# #  "ngin"
# #  "ngin2d"
# #  "np3"
# #  "parutil"
# #  "procs3"
# #  "udm3"
# #  "ug_gm3"
# #  "uggrape3"
# #  "uggraph3"
# #  "ugL3"
# #  "ugS3"
# #  "ugui3"
# )

#DumuxAddPathSuffixes("${MyIncludeSuffixes}" "${MyLibSuffixes}" )

DumuxFindIncludeDir("alugrid_2d.h")
DumuxFindLibrary("alugrid")

DumuxRequiredLibsFound("alugrid")
DumuxIncludeDirsFound()
DumuxCheckFound()
