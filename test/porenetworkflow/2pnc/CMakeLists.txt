add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dune_add_test(NAME test_pnm2p2c_evaporation
              SOURCES test_pnm2p2c_evaporation.cc
              COMPILE_DEFINITIONS ISOTHERMAL=1
              COMMAND ${dumux_INCLUDE_DIRS}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/pnm2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm2p-00147.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm2p")
