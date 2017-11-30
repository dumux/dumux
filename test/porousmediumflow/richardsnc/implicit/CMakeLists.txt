add_input_file_links()

dune_add_test(SOURCES test_richardsnc_fv.cc
              NAME test_boxrichardsnc
              COMPILE_DEFINITIONS TYPETAG=RichardsWellTracerBoxProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardswelltracerbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsnc-00035.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsnc test_richardsnc.input -Problem.Name test_boxrichardsnc")

dune_add_test(SOURCES test_richardsnc_fv.cc
              NAME test_ccrichardsnc
              COMPILE_DEFINITIONS TYPETAG=RichardsWellTracerCCProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardswelltracercc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsnc-00033.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsnc test_richardsnc.input -Problem.Name test_ccrichardsnc")

#install sources
install(FILES
richardswelltracerproblem.hh
richardswelltracerspatialparams.hh
test_richardsnc_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/richardsnc/implicit)
