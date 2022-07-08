dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_twolevel_pnm1p_permeabilityupscaling
               LABELS porenetwork
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK dune-foamgrid_FOUND)
