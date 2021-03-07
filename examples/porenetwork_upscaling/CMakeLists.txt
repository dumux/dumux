dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME example_pnm1p_permeabilityupscaling
               LABELS porenetwork example
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK dune-foamgrid_FOUND
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/example_pnm1p_permeabilityupscaling
               CMD_ARGS -Problem.ReferenceData "2.95910919e-13 2.91015548e-13 2.71752264e-13"
                        -Problem.TestEpsilon 1e-7)
