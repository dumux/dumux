# link the parameter file and the grid file to the build directory
dune_symlink_to_source_files(FILES params.input network.dgf plot.py)

# add a new executable
dumux_add_test(
    NAME example_embedded_network_1d3d
    SOURCES main.cc
)
