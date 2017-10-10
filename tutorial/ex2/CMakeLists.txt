# a compositional two-phase simulation program
dune_add_test(NAME exercise2
              SOURCES exercise2.cc)

# add a symlink for the input file
dune_symlink_to_source_files(FILES "exercise2.input")
