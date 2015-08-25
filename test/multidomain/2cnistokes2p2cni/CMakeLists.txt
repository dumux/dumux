add_input_file_links()

add_dumux_test(test_2cnistokes2p2cni test_2cnistokes2p2cni test_2cnistokes2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni_reference.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cni-ff-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni-00004.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cni-pm-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/2p2cni-00004.vtu
                  --zeroThreshold {"v":1e-7,"velocityN":1e-12,"velocityW":5e-10})

add_dumux_test(test_2cnistokes2p2cni_boundarylayer test_2cnistokes2p2cni test_2cnistokes2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni_boundarylayer.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cniboundarylayer-ff-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni_boundarylayer-00004.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cniboundarylayer-pm-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/2p2cni_boundarylayer-00004.vtu
                  --zeroThreshold {"v":1e-7,"velocityN":1e-12,"velocityW":5e-10})

#install sources
install(FILES
2cnistokes2p2cniproblem.hh
2cnistokes2p2cnispatialparams.hh
2p2cnisubproblem.hh
stokes2cnisubproblem.hh
test_2cnistokes2p2cni.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/2cnistokes2p2cni)
