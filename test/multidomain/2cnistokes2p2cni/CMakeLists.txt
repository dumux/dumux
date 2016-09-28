add_input_file_links()

add_dumux_test(test_2cnistokes2p2cni test_2cnistokes2p2cni test_2cnistokes2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni_reference.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cni-ff-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni-00007.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cni-pm-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/2p2cni-00007.vtu
                  --zeroThreshold {"v_1":1e-6,"velocityN_0":5e-11,"velocityW_0":5e-10,"velocityW_1":5e-9})

add_dumux_test(test_2cnistokes2p2cni_boundarylayer test_2cnistokes2p2cni test_2cnistokes2p2cni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni_boundarylayer.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cniboundarylayer-ff-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni_boundarylayer-00008.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cniboundarylayer-pm-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/2p2cni_boundarylayer-00008.vtu
                  --zeroThreshold {"v_1":1e-6,"velocityN_0":1e-8,"velocityN_1":1e-8,"velocityW_0":1e-8,"velocityW_1":1e-8,"pc":1e2,"mobN":1e-2})

#install sources
install(FILES
2cnistokes2p2cniproblem.hh
2cnistokes2p2cnispatialparams.hh
2p2cnisubproblem.hh
stokes2cnisubproblem.hh
test_2cnistokes2p2cni.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/2cnistokes2p2cni)
