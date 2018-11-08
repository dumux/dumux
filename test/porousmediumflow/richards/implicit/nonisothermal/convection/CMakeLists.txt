add_input_file_links(FILES params.input)

# non-isothermal tests
dune_add_test(NAME test_richardsni_convection_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConvectionBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richardsni_convection_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_convection_box-00044.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_convection_box params.input -Problem.Name test_richardsni_convection_box"
                      )

dune_add_test(NAME test_richardsni_convection_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConvectionCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richardsni_convection_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_convection_tpfa-00043.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_convection_tpfa params.input -Problem.Name test_richardsni_convection_tpfa")

#install sources
install(FILES
problem.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/richards/nonisothermal/convection)
