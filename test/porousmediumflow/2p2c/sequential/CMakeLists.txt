add_input_file_links()

dune_add_test(NAME test_adaptive2p2c2d
              SOURCES test_adaptive2p2c2d.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_adaptive2p2c2d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c2d-00006.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c2d")

dune_add_test(NAME test_adaptive2p2c3d
              SOURCES test_adaptive2p2c3d.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_adaptive2p2c3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c3d-00013.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c3d")

dune_add_test(NAME test_dec2p2c
              SOURCES test_dec2p2c.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_dec2p2c-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_dec2p2c-00019.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_dec2p2c")

dune_add_test(NAME test_multiphysics2p2c
              SOURCES test_multiphysics2p2c.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_multiphysics2p2c-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_multiphysics2p2c-00019.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_multiphysics2p2c")
