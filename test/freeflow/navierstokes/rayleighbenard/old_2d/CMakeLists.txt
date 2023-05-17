add_input_file_links()

add_executable(test_ff_stokesni_rayleighbenard_old EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_stokesni_rayleighbenard_old PUBLIC "NONISOTHERMAL=1")

dumux_add_test(NAME test_ff_stokesni_rayleighbenard_old
              TARGET test_ff_stokesni_rayleighbenard_old
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokesni_rayleighbenard-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokesni_rayleighbenard-00005.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokesni_rayleighbenard_old params.input
                             -Problem.Name test_ff_stokesni_rayleighbenard")