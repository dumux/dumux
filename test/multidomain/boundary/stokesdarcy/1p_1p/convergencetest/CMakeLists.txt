add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py")

dune_add_test(NAME test_md_boundary_darcy1p_stokes1p_convtest
              SOURCES main.cc
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_darcy1p_stokes1p_convtest params.input)
