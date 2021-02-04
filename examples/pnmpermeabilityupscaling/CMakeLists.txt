dune_symlink_to_source_files(FILES "params.input")

dune_add_test(NAME example_pnm1p_permeabilityupscaling
              LABELS porousmediumflow 1p example
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/example_pnm1p_permeabilityupscaling
              CMD_ARGS -Problem.ReferenceData "3.310596e-13 3.258913e-13 3.057024e-13"
                       -Problem.TestEpsilon 1e-6)
