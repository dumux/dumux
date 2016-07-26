add_input_file_links()

dune_symlink_to_source_files(FILES evaporationRates.gp)

add_dumux_test(test_2cnizeroeq2p2cni test_2cnizeroeq2p2cni test_2cnizeroeq2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_2cnizeroeq2p2cni -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnizeroeq2p2cni_reference.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/2cnizeroeq2p2cni-ff-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/zeroeq2cni-00007.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/2cnizeroeq2p2cni-pm-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/2p2cni-00007.vtu)

#install sources
install(FILES
2cnizeroeq2p2cniproblem.hh
2cnizeroeq2p2cnispatialparameters.hh
2p2cnisubproblem.hh
test_2cnizeroeq2p2cni.cc
zeroeq2cnisubproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/2cnizeroeq2p2cni)
