add_input_file_links()

dune_add_test(SOURCES test_richardsnc_fv.cc
              NAME test_richardsnc_box
              COMPILE_DEFINITIONS TYPETAG=RichardsWellTracerBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardswelltracerbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richardsnc_box-00035.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsnc_box test_richardsnc.input -Problem.Name test_richardsnc_box")

dune_add_test(SOURCES test_richardsnc_fv.cc
              NAME test_richardsnc_tpfa
              COMPILE_DEFINITIONS TYPETAG=RichardsWellTracerCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardswelltracercc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richardsnc_tpfa-00033.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsnc_tpfa test_richardsnc.input -Problem.Name test_richardsnc_tpfa")

#install sources
install(FILES
richardswelltracerproblem.hh
richardswelltracerspatialparams.hh
test_richardsnc_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/richardsnc/implicit)
