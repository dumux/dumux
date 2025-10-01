add_subdirectory(analyticalsolution)
add_subdirectory(mandelsubmech)
add_subdirectory(mandelsubflow)
dune_symlink_to_source_files(FILES params.input plot_results.py)
dumux_add_test(
    NAME test_md_mandel
    SOURCES main.cc
    LABELS poromechanics
)
