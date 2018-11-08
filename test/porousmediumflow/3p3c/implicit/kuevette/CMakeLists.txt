add_input_file_links(FILES params.input)

dune_add_test(NAME test_3p3cni_kuevette_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=KuevetteBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3p3cni_kuevette_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_kuevette_box-00004.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_kuevette_box params.input -Problem.Name test_3p3cni_kuevette_box")

dune_add_test(NAME test_3p3cni_kuevette_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=KuevetteCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3p3cni_kuevette_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_kuevette_tpfa-00004.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_kuevette_tpfa params.input -Problem.Name test_3p3cni_kuevette_tpfa")

#install sources
install(FILES
main.cc
spatialparams.hh
problem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/3p3c/kuevette)
