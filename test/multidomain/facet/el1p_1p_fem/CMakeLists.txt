dune_symlink_to_source_files(FILES "grids" "params.input")
set(CMAKE_BUILD_TYPE Release)

# using FEM (1st order) and TPFA
dune_add_test(NAME test_el1p1p_fem_1storder_tpfa
              SOURCES main_fem_fv.cc
              COMPILE_DEFINITIONS BULKFLOWTYPETAG=OnePBulkTpfa
              COMPILE_DEFINITIONS FACETFLOWTYPETAG=OnePFacetTpfa
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )")

# using FEM (1st order) and MPFA
dune_add_test(NAME test_el1p1p_fem_1storder_mpfa
              SOURCES main_fem_fv.cc
              COMPILE_DEFINITIONS FEMORDER=1
              COMPILE_DEFINITIONS BULKFLOWTYPETAG=OnePBulkMpfa
              COMPILE_DEFINITIONS FACETFLOWTYPETAG=OnePFacetTpfa
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )")

# using FEM (1st order) and BOX
dune_add_test(NAME test_el1p1p_fem_1storder_box
              SOURCES main_fem_fv.cc
              COMPILE_DEFINITIONS FEMORDER=1
              COMPILE_DEFINITIONS BULKFLOWTYPETAG=OnePBulkBox
              COMPILE_DEFINITIONS FACETFLOWTYPETAG=OnePFacetBox
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )")


# using FEM (2nd order) and TPFA
dune_add_test(NAME test_el1p1p_fem_2ndorder_tpfa
              SOURCES main_fem_fv.cc
              COMPILE_DEFINITIONS FEMORDER=2
              COMPILE_DEFINITIONS BULKFLOWTYPETAG=OnePBulkTpfa
              COMPILE_DEFINITIONS FACETFLOWTYPETAG=OnePFacetTpfa
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )")

# using FEM (2nd order) and TPFA
dune_add_test(NAME test_el1p1p_fem_2ndorder_mpfa
              SOURCES main_fem_fv.cc
              COMPILE_DEFINITIONS FEMORDER=2
              COMPILE_DEFINITIONS BULKFLOWTYPETAG=OnePBulkMpfa
              COMPILE_DEFINITIONS FACETFLOWTYPETAG=OnePFacetTpfa
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )")

# using FEM (2nd order) and BOX
dune_add_test(NAME test_el1p1p_fem_2ndorder_box
             SOURCES main_fem_fv.cc
             COMPILE_DEFINITIONS FEMORDER=2
             COMPILE_DEFINITIONS BULKFLOWTYPETAG=OnePBulkBox
             COMPILE_DEFINITIONS FACETFLOWTYPETAG=OnePFacetBox
             CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )")

# using Box and TPFA
dune_add_test(NAME test_el1p1p_box_tpfa
              SOURCES main_box_fv.cc
              COMPILE_DEFINITIONS BULKFLOWTYPETAG=OnePBulkTpfa
              COMPILE_DEFINITIONS FACETFLOWTYPETAG=OnePFacetTpfa
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )")

# using Box and MPFA
dune_add_test(NAME test_el1p1p_box_mpfa
              SOURCES main_box_fv.cc
              COMPILE_DEFINITIONS BULKFLOWTYPETAG=OnePBulkMpfa
              COMPILE_DEFINITIONS FACETFLOWTYPETAG=OnePFacetTpfa
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )")
