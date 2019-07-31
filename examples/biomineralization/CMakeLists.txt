dune_symlink_to_source_files(FILES injections_checkpoints.dat injections_type.dat params.input)

# compile MICP simplified chemistry column setup
dune_add_test(NAME example_biomineralization
              SOURCES main.cc)
