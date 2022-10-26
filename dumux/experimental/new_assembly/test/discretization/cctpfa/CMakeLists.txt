dune_symlink_to_source_files(FILES "cp_grid.grdecl")

dumux_add_test(
    NAME test_cctpfa_grid_geometry
    SOURCES test_grid_geometry.cc
    LABELS unit new_assembly
)
dumux_add_test(
    NAME test_cctpfa_grid_geometry_non_conforming
    CMAKE_GUARD dune-alugrid_FOUND
    SOURCES test_grid_geometry_non_conforming.cc
    LABELS unit new_assembly
)
dumux_add_test(
    NAME test_cctpfa_grid_geometry_network
    CMAKE_GUARD dune-foamgrid_FOUND
    SOURCES test_grid_geometry_network.cc
    LABELS unit new_assembly
)

dumux_add_test(
    NAME test_cctpfa_grid_geometry_corner_point
    CMAKE_GUARD HAVE_OPM_GRID
    COMPILE_DEFINITIONS HAVE_ECL_INPUT=1
    SOURCES test_grid_geometry_corner_point.cc
    LABELS unit new_assembly
    COMMAND ./test_cctpfa_grid_geometry_corner_point
    CMD_ARGS -Grid.File "cp_grid.grdecl"
)