add_executable(test_ff_navierstokes_closedsystem EXCLUDE_FROM_ALL main.cc)

dune_add_test(NAME test_ff_navierstokes_closedsystem_ldc_re1
              TARGET test_ff_navierstokes_closedsystem
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_ldc_re1-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem_ldc_re1-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem params_re1.input
                             -Problem.Name test_ff_navierstokes_closedsystem_ldc_re1")

dune_add_test(NAME test_ff_navierstokes_closedsystem_ldc_re1000
              TARGET test_ff_navierstokes_closedsystem
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_ldc_re1000-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem_ldc_re1000-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem params_re1000.input
                             -Problem.Name test_ff_navierstokes_closedsystem_ldc_re1000")

dune_add_test(NAME test_ff_navierstokes_closedsystem_hydrostaticpressure
              TARGET test_ff_navierstokes_closedsystem
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_hydrostaticpressure-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem_hydrostaticpressure-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_closedsystem params_hydrostaticpressure.input
                             -Problem.Name test_ff_navierstokes_closedsystem_hydrostaticpressure"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-16})

add_executable(test_ff_navierstokes_densitydrivenflow_boussinesq EXCLUDE_FROM_ALL main.cc)
add_executable(test_ff_navierstokes_densitydrivenflow EXCLUDE_FROM_ALL main.cc)

target_compile_definitions(test_ff_navierstokes_densitydrivenflow PUBLIC "NONISOTHERMAL=1")
target_compile_definitions(test_ff_navierstokes_densitydrivenflow_boussinesq PUBLIC "NONISOTHERMAL=1" "BOUSSINESQ=1")

dune_symlink_to_source_files(FILES "params_re1.input" "params_re1000.input" "params_hydrostaticpressure.input" "params_densitydrivenflow.input")
