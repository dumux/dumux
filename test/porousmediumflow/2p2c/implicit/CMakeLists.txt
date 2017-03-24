add_input_file_links()

# isothermal tests
add_dumux_test(test_box2p2c test_box2p2c test_box2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/injectionbox-reference-initial.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injectionbox-00000.vtu
                         ${CMAKE_SOURCE_DIR}/test/references/injectionbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injectionbox-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2p2c
                           -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2p2c.input")

add_dumux_test(test_box2pnc test_box2pnc test_box2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/injectionbox-reference-initial.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injectionbox-00000.vtu
                         ${CMAKE_SOURCE_DIR}/test/references/injectionbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injectionbox-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pnc
                           -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2p2c.input")
target_compile_definitions(test_box2pnc PUBLIC "USE_TWOPNC=1")

add_dumux_test(test_cc2p2c test_cc2p2c test_cc2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/injectioncc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injectioncc-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2p2c")

# non-isothermal tests
add_dumux_test(test_box2p2cni test_box2p2cni test_box2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/waterairbox-reference-initial.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/waterairbox-00000.vtu
                         ${CMAKE_SOURCE_DIR}/test/references/waterairbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/waterairbox-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2p2cni
                           -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2p2cni_reference.input")

add_dumux_test(test_box2pncni test_box2pncni test_box2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/waterairbox-reference-initial.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/waterairbox-00000.vtu
                         ${CMAKE_SOURCE_DIR}/test/references/waterairbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/waterairbox-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pncni
                           -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2p2cni_reference.input")
target_compile_definitions(test_box2pncni PUBLIC "USE_TWOPNC=1")

add_dumux_test(test_cc2p2cni test_cc2p2cni test_cc2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/wateraircc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/wateraircc-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2p2cni -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc2p2cni_reference.input")

#install sources
install(FILES
injectionproblem.hh
injectionspatialparams.hh
test_box2p2c.cc
test_box2p2cni.cc
test_cc2p2c.cc
test_cc2p2cni.cc
waterairproblem.hh
waterairspatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p2c)
