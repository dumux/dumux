add_input_file_links(FILES params.input)

dune_add_test(NAME test_3pwateroil_sagd_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePWaterOilSagdBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMD_ARGS --script fuzzy
            --files ${CMAKE_SOURCE_DIR}/test/references/test_3pwateroil_sagd_box-reference.vtu
                   ${CMAKE_CURRENT_BINARY_DIR}/test_3pwateroil_sagd_box-00001.vtu
            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3pwateroil_sagd_box params.input -Problem.Name test_3pwateroil_sagd_box"
  )

#install sources
install(FILES
problem.hh
spatialparams.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/3pwateroil)
