dune_add_test(NAME test_embedded_1d3d_1p1p_tpfatpfa
              SOURCES test_1p_1p.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCCTypeTag LOWDIMTYPETAG=BloodFlowCCTypeTag)
dune_symlink_to_source_files(FILES "grids" "test_1p_1p.input")
set(CMAKE_BUILD_TYPE Release)
