add_input_file_links()

dune_add_test(NAME test_kepsilon1p2cnidarcy2p2cni
              SOURCES test_rans1p2cnidarcy2p2cni.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_kepsilon1p2cnidarcy2p2cni_kepsilon-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_kepsilon1p2cnidarcy2p2cni_kepsilon-00040.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_kepsilon1p2cnidarcy2p2cni_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_kepsilon1p2cnidarcy2p2cni_darcy-00040.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_kepsilon1p2cnidarcy2p2cni test_rans1p2cnidarcy2p2cni.input")
target_compile_definitions(test_kepsilon1p2cnidarcy2p2cni PUBLIC "KEPSILON=1")

dune_add_test(NAME test_komega1p2cnidarcy2p2cni
              SOURCES test_rans1p2cnidarcy2p2cni.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_komega1p2cnidarcy2p2cni_komega-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_komega1p2cnidarcy2p2cni_komega-00040.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_komega1p2cnidarcy2p2cni_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_komega1p2cnidarcy2p2cni_darcy-00040.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_komega1p2cnidarcy2p2cni test_rans1p2cnidarcy2p2cni.input")
target_compile_definitions(test_komega1p2cnidarcy2p2cni PUBLIC "KOMEGA=1")

dune_add_test(NAME test_lowrekepsilon1p2cnidarcy2p2cni
              SOURCES test_rans1p2cnidarcy2p2cni.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_lowrekepsilon1p2cnidarcy2p2cni_lowrekepsilon-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_lowrekepsilon1p2cnidarcy2p2cni_lowrekepsilon-00040.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_lowrekepsilon1p2cnidarcy2p2cni_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_lowrekepsilon1p2cnidarcy2p2cni_darcy-00040.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_lowrekepsilon1p2cnidarcy2p2cni test_rans1p2cnidarcy2p2cni.input")
target_compile_definitions(test_lowrekepsilon1p2cnidarcy2p2cni PUBLIC "LOWREKEPSILON=1")

dune_add_test(NAME test_oneeq1p2cnidarcy2p2cni
              SOURCES test_rans1p2cnidarcy2p2cni.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_oneeq1p2cnidarcy2p2cni_oneeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_oneeq1p2cnidarcy2p2cni_oneeq-00040.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_oneeq1p2cnidarcy2p2cni_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_oneeq1p2cnidarcy2p2cni_darcy-00040.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_oneeq1p2cnidarcy2p2cni test_rans1p2cnidarcy2p2cni.input")
target_compile_definitions(test_oneeq1p2cnidarcy2p2cni PUBLIC "ONEEQ=1")

dune_add_test(NAME test_zeroeq1p2cnidarcy2p2cni
              SOURCES test_rans1p2cnidarcy2p2cni.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_zeroeq1p2cnidarcy2p2cni_zeroeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq1p2cnidarcy2p2cni_zeroeq-00040.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_zeroeq1p2cnidarcy2p2cni_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq1p2cnidarcy2p2cni_darcy-00040.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq1p2cnidarcy2p2cni test_rans1p2cnidarcy2p2cni.input")
target_compile_definitions(test_zeroeq1p2cnidarcy2p2cni PUBLIC "ZEROEQ=1")
