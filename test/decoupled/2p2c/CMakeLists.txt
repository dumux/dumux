add_input_file_links()

add_dumux_test(test_adaptive2p2c2d test_adaptive2p2c2d test_adaptive2p2c2d.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
               --script fuzzy
               --files ${CMAKE_SOURCE_DIR}/test/references/test_adaptive2p2c2d-reference.vtu
                       ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c2d-00008.vtu
               --command "${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c2d")

add_dumux_test(test_adaptive2p2c3d test_adaptive2p2c3d test_adaptive2p2c3d.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
               --script fuzzy
               --files ${CMAKE_SOURCE_DIR}/test/references/test_adaptive2p2c3d-reference.vtu
                       ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c3d-00012.vtu
               --command "${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c3d")

add_dumux_test(test_dec2p2c test_dec2p2c test_dec2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
               --script fuzzy
               --files ${CMAKE_SOURCE_DIR}/test/references/test_dec2p2c-reference.vtu
                       ${CMAKE_CURRENT_BINARY_DIR}/test_dec2p2c-00021.vtu
               --command "${CMAKE_CURRENT_BINARY_DIR}/test_dec2p2c")

add_dumux_test(test_multiphysics2p2c test_multiphysics2p2c test_multiphysics2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
               --script fuzzy
               --files ${CMAKE_SOURCE_DIR}/test/references/test_multiphysics2p2c-reference.vtu
                       ${CMAKE_CURRENT_BINARY_DIR}/test_multiphysics2p2c-00021.vtu
               --command "${CMAKE_CURRENT_BINARY_DIR}/test_multiphysics2p2c")

#install sources
install(FILES
test_adaptive2p2c2d.cc
test_adaptive2p2c2dproblem.hh
test_adaptive2p2c3d.cc
test_adaptive2p2c3dproblem.hh
test_dec2p2c.cc
test_dec2p2cproblem.hh
test_dec2p2c_spatialparams.hh
test_multiphysics2p2c.cc
test_multiphysics2p2cproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/decoupled/2p2c)
