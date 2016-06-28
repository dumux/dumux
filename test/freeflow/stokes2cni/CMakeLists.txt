add_input_file_links()

add_dumux_test(test_stokes2cni test_stokes2cni test_stokes2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni"
                  --files ${CMAKE_SOURCE_DIR}/test/references/stokes2cni-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni-00008.vtu)

add_dumux_test(test_stokes2cni_conduction test_stokes2cni test_stokes2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes2cni_conduction.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/stokes2cni_conduction-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni_conduction-00008.vtu)

#install sources
install(FILES
stokes2cnitestproblem.hh
test_stokes2cni.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/stokes2cni)
