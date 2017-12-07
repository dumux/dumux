add_input_file_links()

# isothermal tests
dune_add_test(NAME test_3p3c_box
              SOURCES test_3p3c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=InfiltrationThreePThreeCBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/infiltrationbox-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/infiltration3p3cbox-00005.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3c_box test_3p3c_fv.input -Problem.Name infiltration3p3cbox")

dune_add_test(NAME test_3p3c_tpfa
              SOURCES test_3p3c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=InfiltrationThreePThreeCCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/infiltrationcc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/infiltration3p3ccc-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3c_tpfa test_3p3c_fv.input -Problem.Name infiltration3p3ccc")

# non-isothermal tests
dune_add_test(NAME test_3p3cni_columnxylol_box
              SOURCES test_3p3c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=ColumnBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/column3p3cnibox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/columnxylol_box-00062.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_columnxylol_box test_columnxylol_fv.input -Problem.Name columnxylol_box")

dune_add_test(NAME test_3p3cni_columnxylol_tpfa
              SOURCES test_3p3c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=ColumnCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/column3p3cnicc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/columnxylol_tpfa-00053.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_columnxylol_tpfa test_columnxylol_fv.input -Problem.Name columnxylol_tpfa")

dune_add_test(NAME test_3p3cni_kuevette_box
              SOURCES test_3p3c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=KuevetteBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/kuevette3p3cnibox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/kuevette_box-00004.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_kuevette_box test_kuvette_fv.input -Problem.Name kuevette_box")

dune_add_test(NAME test_3p3cni_kuevette_tpfa
              SOURCES test_3p3c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=KuevetteCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/kuevette3p3cnicc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/kuvette_tpfa-00004.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_kuevette_tpfa test_kuvette_fv.input -Problem.Name kuvette_tpfa")

#install sources
install(FILES
columnxylolproblem.hh
columnxylolspatialparams.hh
infiltrationproblem.hh
infiltrationspatialparameters.hh
kuevetteproblem.hh
kuevettespatialparams.hh
test_3p3c_fv.cc
test_3p3c_fv.input
test_columnxylol_fv.input
test_kuvette_fv.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/3p3c)
