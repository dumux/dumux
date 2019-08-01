add_input_file_links()
dune_symlink_to_source_files("params.input")


#test for the shallow water model
dune_add_test(NAME shallowwaterfriction
              SOURCES main.cc)
