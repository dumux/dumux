dune_symlink_to_source_files(FILES "params.input")

dune_add_test(NAME test_md_mortar_darcy1p_stokes1p
              LABELS multidomain multidomain_boundary mortar stokesdarcy 1p stokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./test_md_mortar_darcy1p_stokes1p
              CMD_ARGS params.imput)
