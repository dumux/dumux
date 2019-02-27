add_input_file_links()
dune_symlink_to_source_files(FILES grids)

# build target for the CO2 test problem
# Ignore the porosity for all box models since it is defined element-wise in these test
# but the default 2p2c implementation outputs porosity per vertex.
# Depending on the order of the elements, the porosity would differ in these cases.
dune_add_test(NAME test_co2_box
              LABELS porousmediumflow co2 box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=HeterogeneousBox
              CMAKE_GUARD "( dune-alugrid_FOUND AND DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_co2_box-00019.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_box params.input -Problem.Name test_co2_box"
                       --zeroThreshold {"porosity":1})

dune_add_test(NAME test_co2_tpfa
              LABELS porousmediumflow co2 tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=HeterogeneousCCTpfa
              CMAKE_GUARD "( dune-alugrid_FOUND AND DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_co2_tpfa-00018.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_tpfa params.input -Problem.Name test_co2_tpfa")

# build target for the CO2 non-isothermal test problem
dune_add_test(NAME test_co2ni_box
              LABELS porousmediumflow co2 co2ni box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=HeterogeneousNIBox
                                  ISOTHERMAL=0
              CMAKE_GUARD "( dune-alugrid_FOUND AND DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_co2ni_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_box-00019.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_box paramsni.input -Problem.Name test_co2ni_box"
                       --zeroThreshold {"porosity":1})

dune_add_test(NAME test_co2ni_tpfa
              LABELS porousmediumflow co2 co2ni tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=HeterogeneousNICCTpfa
                                  ISOTHERMAL=0
              CMAKE_GUARD "( dune-alugrid_FOUND AND DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_co2ni_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_tpfa-00018.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_tpfa paramsni.input -Problem.Name test_co2ni_tpfa")
