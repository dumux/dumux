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
