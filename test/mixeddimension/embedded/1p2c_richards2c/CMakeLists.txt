add_dumux_test(test_rosi2c test_rosi2c test_rosi2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/rosi2c-root-reference.vtp
                         ${CMAKE_CURRENT_BINARY_DIR}/rosi2c-root-00010.vtp
                         ${CMAKE_SOURCE_DIR}/test/references/rosi2c-soil-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/rosi2c-soil-00010.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_rosi2c -TimeManager.TEnd 86400")

dune_symlink_to_source_files(FILES "grids" "test_rosi2c.input")
