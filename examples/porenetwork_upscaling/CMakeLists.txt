dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME example_pnm1p_permeabilityupscaling
               LABELS porenetwork example
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/example_pnm1p_permeabilityupscaling
               CMD_ARGS -Problem.ReferenceData "3.083710e-13 2.870186e-13 2.671460e-13"
                        -Problem.TestEpsilon 1e-7)
