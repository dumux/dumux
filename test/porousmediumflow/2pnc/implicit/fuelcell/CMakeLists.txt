dune_symlink_to_source_files(FILES params.input)

# isothermal tests
dumux_add_test(NAME test_2pnc_fuelcell_box
              LABELS porousmediumflow 2pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FuelCellBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pnc_fuelcell_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_fuelcell_box-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_fuelcell_box -ParameterFile params.input -Problem.Name test_2pnc_fuelcell_box")

dumux_add_test(NAME test_2pnc_fuelcell_tpfa
              LABELS porousmediumflow 2pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FuelCellCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pnc_fuelcell_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_fuelcell_tpfa-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_fuelcell_tpfa -ParameterFile params.input -Problem.Name test_2pnc_fuelcell_tpfa")


dumux_add_test(NAME test_2pncni_fuelcell_box
              LABELS porousmediumflow 2pnc  2pncni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FuelCellNIBox NONISOTHERMAL=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pncni_fuelcell_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pncni_fuelcell_box-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncni_fuelcell_box -ParameterFile params.input -Problem.Name test_2pncni_fuelcell_box")
