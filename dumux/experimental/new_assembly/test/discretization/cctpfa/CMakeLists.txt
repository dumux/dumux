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

# todo: delete after testing phase
dumux_add_test(
    NAME test_grid_geometry_performance_new_implementation
    SOURCES test_grid_geometry_performance.cc
    LABELS unit new_assembly
    COMPILE_DEFINITIONS USE_OLD=0
)
dumux_add_test(
    NAME test_grid_geometry_performance_old_implementation
    SOURCES test_grid_geometry_performance.cc
    LABELS unit new_assembly
    COMPILE_DEFINITIONS USE_OLD=1
)
