#add links to input files
add_input_file_links()

add_executable(test_2p1cni_steaminjection_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_2p1cni_steaminjection_tpfa PUBLIC "TYPETAG=TwoPOneCNICCTpfa")

add_executable(test_2p1cni_steaminjection_box EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_2p1cni_steaminjection_box PUBLIC "TYPETAG=TwoPOneCNIBox")

dumux_add_test(NAME test_2p1cni_steaminjection_waterwet_box
              LABELS porousmediumflow 2p1c  2p1cni
              TARGET test_2p1cni_steaminjection_box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2p1cni_steaminjection_waterwet_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_waterwet_box-00007.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_box params_box.input -Problem.Name test_2p1cni_steaminjection_waterwet_box")

dumux_add_test(NAME test_2p1cni_steaminjection_waterwet_tpfa
              LABELS porousmediumflow 2p1c  2p1cni
              TARGET test_2p1cni_steaminjection_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2p1cni_steaminjection_waterwet_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_waterwet_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_tpfa params_tpfa.input -Problem.Name test_2p1cni_steaminjection_waterwet_tpfa")

dumux_add_test(NAME test_2p1cni_steaminjection_gaswet_box
              LABELS porousmediumflow 2p1c  2p1cni
              TARGET test_2p1cni_steaminjection_box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2p1cni_steaminjection_gaswet_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_gaswet_box-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_box params_box.input
                                                                               -Problem.Name test_2p1cni_steaminjection_gaswet_box
                                                                               -SpatialParams.GasWetting true")

dumux_add_test(NAME test_2p1cni_steaminjection_gaswet_tpfa
              LABELS porousmediumflow 2p1c  2p1cni
              TARGET test_2p1cni_steaminjection_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2p1cni_steaminjection_gaswet_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_gaswet_tpfa-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_tpfa params_tpfa.input
                                                                                -Problem.Name test_2p1cni_steaminjection_gaswet_tpfa
                                                                                -SpatialParams.GasWetting true")
