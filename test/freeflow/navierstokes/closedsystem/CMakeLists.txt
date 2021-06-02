add_executable(test_ff_navierstokes_closedsystem_new EXCLUDE_FROM_ALL main_new.cc)


dumux_add_test(NAME test_ff_navierstokes_closedsystem_ldc_re1_new
              TARGET test_ff_navierstokes_closedsystem_new
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_ldc_re1-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem_ldc_re1_new-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem_new params_re1.input
                             -Problem.Name test_ff_navierstokes_closedsystem_ldc_re1_new")

dumux_add_test(NAME test_ff_navierstokes_closedsystem_ldc_re1000_new
              TARGET test_ff_navierstokes_closedsystem_new
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_ldc_re1000-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem_ldc_re1000_new-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem_new params_re1000.input
                             -Problem.Name test_ff_navierstokes_closedsystem_ldc_re1000_new")

dumux_add_test(NAME test_ff_navierstokes_closedsystem_hydrostaticpressure_new
              TARGET test_ff_navierstokes_closedsystem_new
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_hydrostaticpressure-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem_hydrostaticpressure_new-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem_new params_hydrostaticpressure.input
                             -Problem.Name test_ff_navierstokes_closedsystem_hydrostaticpressure_new"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-16})

dune_symlink_to_source_files(FILES "params_re1.input" "params_re1000.input" "params_hydrostaticpressure.input")
