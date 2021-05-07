dune_symlink_to_source_files(FILES "params.input")

# explicit tracer tests (mass fractions)
dumux_add_test(NAME test_ff_advectiondiffusion
              LABELS freeflow advectiondiffusion
              SOURCES main.cc
              COMPILE_DEFINITIONS TRANSPORTTYPETAG=AdvectionDiffusionTpfa FLOWTYPETAG=NavierStokesTest
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_advectiondiffusion-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_ff_advectiondiffusion-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_advectiondiffusion params.input
                       -Problem.Name test_ff_advectiondiffusion")
