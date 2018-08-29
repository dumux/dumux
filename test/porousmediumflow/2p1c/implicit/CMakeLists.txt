#add links to input files
add_input_file_links()

add_executable(test_2p1cni_tpfa EXCLUDE_FROM_ALL test_2p1c_fv.cc)
target_compile_definitions(test_2p1cni_tpfa PUBLIC "TYPETAG=TwoPOneCNICCTpfaTypeTag")

add_executable(test_2p1cni_box EXCLUDE_FROM_ALL test_2p1c_fv.cc)
target_compile_definitions(test_2p1cni_box PUBLIC "TYPETAG=TwoPOneCNIBoxTypeTag")


dune_add_test(NAME test_2p1cni_box_waterwet
              TARGET test_2p1cni_box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectionbox-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_boxsteaminjection-00007.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_box test_boxsteaminjection.input")

dune_add_test(NAME test_2p1cni_box_waterwet_restart
              TARGET test_2p1cni_box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectionbox-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_boxsteaminjection_restart-00002.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_box test_boxsteaminjection.input -Problem.Name test_boxsteaminjection_restart -TimeLoop.DtInitial 15.191 -Restart.Time 277.564 -Restart.File test_boxsteaminjection-00005.vtu")

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_2p1cni_box_waterwet_restart PROPERTIES DEPENDS test_2p1cni_box_waterwet)

dune_add_test(NAME test_2p1cni_tpfa_waterwet
              TARGET test_2p1cni_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectioncc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ccsteaminjection-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_tpfa test_ccsteaminjection.input")

dune_add_test(NAME test_2p1cni_box_gaswet
              TARGET test_2p1cni_box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectionbox_gaswet-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_steaminjectionbox_gaswet-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_box test_boxsteaminjection.input
                                                                               -Problem.Name test_steaminjectionbox_gaswet
                                                                               -SpatialParams.GasWetting true")

dune_add_test(NAME test_2p1cni_tpfa_gaswet
              TARGET test_2p1cni_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectioncc_gaswet-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_steaminjectioncc_gaswet-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_tpfa test_ccsteaminjection.input
                                                                                -Problem.Name test_steaminjectioncc_gaswet
                                                                                -SpatialParams.GasWetting true")
