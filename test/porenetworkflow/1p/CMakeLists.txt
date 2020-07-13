add_subdirectory(nonisothermal)

add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dune_add_test(NAME test_pnm_1p
              SOURCES main.cc
              LABELS porenetworkflow
              COMPILE_DEFINITIONS ISOTHERMAL=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p-00000.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p")

dune_add_test(NAME test_pnm_1p_gravity
              TARGET test_pnm_1p
              LABELS porenetworkflow
              COMPILE_DEFINITIONS ISOTHERMAL=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy --zeroThreshold {"velocity_liq \(m/s\)":1e-15,"volumeFlux":1e-25}
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1p_gravity-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p_gravity-00000.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p params_hydrostatic_pressure.input")
