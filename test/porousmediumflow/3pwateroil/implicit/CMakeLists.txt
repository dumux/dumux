add_input_file_links()

dune_add_test(NAME test_box3pwateroil
              SOURCES test_box3pwateroil.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePWaterOilSagdBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMD_ARGS --script fuzzy
            --files ${CMAKE_SOURCE_DIR}/test/references/3pwateroilbox-reference.vtu
                   ${CMAKE_CURRENT_BINARY_DIR}/sagd-00020.vtu
            --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3pwateroil"
  )

#install sources
install(FILES
3poilwatersagdproblem.hh
3poilwatersagdspatialparams.hh
test_box3pwateroil.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/3pwateroil)
