add_input_file_links()

add_executable(noflowdomain_2p2cni_ed noflowdomain_2p2cni_ed.cc)
dune_symlink_to_source_files(FILES "convergencestudy.sh")

# isothermal tests
add_dumux_test(test_box2p2c test_box2p2c test_box2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/injectionbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injectionbox-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2p2c")

add_dumux_test(test_cc2p2c test_cc2p2c test_cc2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/injectioncc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injectioncc-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2p2c")

add_dumux_test(test_ccmpfa2p2c test_ccmpfa2p2c test_ccmpfa2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/injectionccmpfa-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injectionccmpfa-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccmpfa2p2c")

# non-isothermal tests
add_dumux_test(test_box2p2cni test_box2p2cni test_box2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/waterairbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/waterairbox-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2p2cni -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2p2cni_reference.input")

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
