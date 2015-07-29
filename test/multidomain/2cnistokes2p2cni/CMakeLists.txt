add_input_file_links()

add_dumux_test(test_2cnistokes2p2cni test_2cnistokes2p2cni test_2cnistokes2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni_reference.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cni-ff-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni-00004.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cni-pm-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/2p2cni-00004.vtu)

add_dumux_test(test_2cnistokes2p2cni_boundarylayer test_2cnistokes2p2cni test_2cnistokes2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni_boundarylayer.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cniboundarylayer-ff-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni_boundarylayer-00004.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cniboundarylayer-pm-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/2p2cni_boundarylayer-00004.vtu)
