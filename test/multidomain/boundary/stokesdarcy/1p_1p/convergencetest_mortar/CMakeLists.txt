set(CMAKE_BUILD_TYPE Release)
dune_symlink_to_source_files(FILES "params.input" "convergencetest.py")
dune_add_test(NAME test_md_mortar_darcy1p_stokes1p_convergence
              LABELS multidomain multidomain_boundary mortar stokesdarcy 1p stokes
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ./test_md_mortar_darcy1p_stokes1p
              CMD_ARGS params.imput)
