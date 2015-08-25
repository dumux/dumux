add_input_file_links()

add_dumux_test(test_zeroeq2cni test_zeroeq2cni test_zeroeq2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq2cni"
                  --files ${CMAKE_SOURCE_DIR}/test/references/zeroeq2cni-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/zeroeq2cni-00011.vtu)

#install sources
install(FILES
test_zeroeq2cni.cc
zeroeq2cnitestproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/zeroeq2cni)
