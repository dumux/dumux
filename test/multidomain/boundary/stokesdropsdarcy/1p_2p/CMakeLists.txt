add_input_file_links()

dune_add_test(NAME test_md_boundary_darcy2p_darcy2p_stokes1p
              LABELS multidomain freeflow
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              )