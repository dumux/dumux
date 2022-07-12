add_subdirectory(pnm)
dune_symlink_to_source_files(FILES "params.input" "sphere.raw")
dumux_add_test(NAME test_ff_stokes_permeability_upscaling
               SOURCES main.cc
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-subgrid_FOUND )")
dumux_add_test(NAME test_ff_stokes_permeability_upscaling_chorin
               SOURCES main_chorin.cc
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-subgrid_FOUND )")
