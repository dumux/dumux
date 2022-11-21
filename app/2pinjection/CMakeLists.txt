add_subdirectory(mesh)

dune_add_test(NAME test_gas_injection
              SOURCES main.cc)

add_input_file_links()
