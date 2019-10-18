dune_symlink_to_source_files(FILES "grids" "params.input" "convergencetest.py")

dumux_add_test(NAME test_md_facet_1p1p_tpfa_convergence
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND gmsh_FOUND )"
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkTpfa LOWDIMTYPETAG=OnePLowDimTpfa
              COMMAND ./convergencetest.py
              CMD_ARGS test_md_facet_1p1p_tpfa_convergence 1e-4)

dumux_add_test(NAME test_md_facet_1p1p_box_convergence
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND gmsh_FOUND )"
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkBox LOWDIMTYPETAG=OnePLowDimBox
              COMMAND ./convergencetest.py
              CMD_ARGS test_md_facet_1p1p_box_convergence 1e-4 )
