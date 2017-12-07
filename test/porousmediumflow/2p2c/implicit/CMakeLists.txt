add_input_file_links()

# isothermal tests
dune_add_test(NAME test_2p2c_box
              SOURCES test_2p2c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=InjectionBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/injectionbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/injection_box-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_box test_2p2c_fv.input -Problem.Name injection_box")

dune_add_test(NAME test_2p2c_tpfa
              SOURCES test_2p2c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=InjectionCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/injectioncc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/injection_tpfa-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_tpfa test_2p2c_fv.input -Problem.Name injection_tpfa")

dune_add_test(NAME test_2p2c_mpfa
              SOURCES test_2p2c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=InjectionCCMpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/injectioncc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/injection_mpfa-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_mpfa test_2p2c_fv.input -Problem.Name injection_mpfa")

# non-isothermal tests
dune_add_test(NAME test_2p2cni_box
              SOURCES test_2p2c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=WaterAirBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/waterairbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/waterair_box-00013.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_box test_2p2cni_fv.input -Problem.Name waterair_box")

dune_add_test(NAME test_2p2cni_tpfa
              SOURCES test_2p2c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=WaterAirCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/wateraircc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/waterair_tpfa-00009.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_tpfa test_2p2cni_fv.input -Problem.Name waterair_tpfa")

#install sources
install(FILES
injectionproblem.hh
injectionspatialparams.hh
waterairproblem.hh
waterairspatialparams.hh
test_2p2c_fv.cc
test_2p2c_fv.input
test_2p2cni_fv.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p2c)
