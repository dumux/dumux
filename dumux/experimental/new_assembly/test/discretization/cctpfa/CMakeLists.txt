dumux_add_test(
    NAME test_cctpfa_grid_geometry
    SOURCES test_grid_geometry.cc
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
