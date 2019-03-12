dune_add_test(NAME test_ff_navierstokes_kovasznay
              LABELS freeflow navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS UPWINDSCHEMEORDER=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_kovasznay-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay params.input
                             -Problem.Name test_ff_navierstokes_kovasznay")

dune_add_test(NAME test_ff_navierstokes_kovasznay_higherorder
              SOURCES main.cc
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS UPWINDSCHEMEORDER=2
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_kovasznay_higherorder-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay_higherorder-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay_higherorder params.input
                             -Problem.Name test_ff_navierstokes_kovasznay_higherorder
                             -Flux.TvdApproach Hou
                             -Flux.DifferencingScheme Vanleer")

dune_symlink_to_source_files(FILES "params.input")
