dune_symlink_to_source_files(FILES convergencetest.py)
add_input_file_links()

dumux_add_test(NAME test_ff_rans_komega_analytical
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=AnalyticalKOmega
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND test_ff_rans_komega_analytical
              CMD_ARGS params.input)
