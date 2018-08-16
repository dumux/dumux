#add links to input files
add_input_file_links()

dune_add_test(NAME test_2p1cni_box
              SOURCES test_2p1c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPOneCNIBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectionbox-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_boxsteaminjection-00007.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_box test_boxsteaminjection.input")

dune_add_test(NAME test_2p1cni_tpfa
              SOURCES test_2p1c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPOneCNICCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectioncc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ccsteaminjection-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_tpfa test_ccsteaminjection.input")

dune_add_test(NAME test_2p1cni_box_gaswet
              SOURCES test_2p1c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPOneCNIBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectionbox_gaswet-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_steaminjectionbox_gaswet-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_box test_boxsteaminjection.input
                                                                               -Problem.Name test_steaminjectionbox_gaswet
                                                                               -SpatialParams.GasWetting true")

dune_add_test(NAME test_2p1cni_tpfa_gaswet
              SOURCES test_2p1c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPOneCNICCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectioncc_gaswet-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_steaminjectioncc_gaswet-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_tpfa test_ccsteaminjection.input
                                                                                -Problem.Name test_steaminjectioncc_gaswet
                                                                                -SpatialParams.GasWetting true")
