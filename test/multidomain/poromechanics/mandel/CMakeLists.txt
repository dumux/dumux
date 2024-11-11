add_subdirectory(analyticalsolution)

dune_symlink_to_source_files(FILES params.input)
dumux_add_test(
    NAME test_md_mandel
    SOURCES main.cc
    LABELS poromechanics
)
