add_input_file_links()
dune_symlink_to_source_files(FILES grids)

# build target for the CO2 test problem
# Ignore the porosity for all box models since it is defined element-wise in these test
# but the default 2p2c implementation outputs porosity per vertex.
# Depending on the order of the elements, the porosity would differ in these cases.
dune_add_test(NAME test_co2_box
              SOURCES test_co2_fv.cc
              COMPILE_DEFINITIONS TYPETAG=HeterogeneousBoxTypeTag
              CMAKE_GUARD "( dune-alugrid_FOUND AND DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/co2box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_box-00019.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_box test_co2_fv.input -Problem.Name co2_box"
                       --zeroThreshold {"porosity":1})

dune_add_test(NAME test_co2_tpfa
              SOURCES test_co2_fv.cc
              COMPILE_DEFINITIONS TYPETAG=HeterogeneousCCTpfaTypeTag
              CMAKE_GUARD "( dune-alugrid_FOUND AND DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/co2cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_tpfa-00018.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_tpfa test_co2_fv.input -Problem.Name co2_tpfa")

# # restart test (currently disbled because restart is not yet implemented)
# dune_add_test(NAME test_co2_box_restart
#               TARGET test_co2_box
#               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#               CMD_ARGS --script fuzzy
#                        --files ${CMAKE_SOURCE_DIR}/test/references/co2box-reference.vtu
#                                ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousbox-00019.vtu
#                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_box -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_restartco2.input"
#                        --zeroThreshold {"porosity":1})

# the restart test has to run after the test that produces the restart file
# set_tests_properties(test_restartco2 PROPERTIES DEPENDS test_co2_box)

# build target for the CO2 non-isothermal test problem
dune_add_test(NAME test_co2ni_box
              SOURCES test_co2_fv.cc
              COMPILE_DEFINITIONS TYPETAG=HeterogeneousNIBoxTypeTag
                                  ISOTHERMAL=0
              CMAKE_GUARD "( dune-alugrid_FOUND AND DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/co2nibox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2ni_box-00019.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_box test_co2ni_fv.input -Problem.Name co2ni_box"
                       --zeroThreshold {"porosity":1})

dune_add_test(NAME test_co2ni_tpfa
              SOURCES test_co2_fv.cc
              COMPILE_DEFINITIONS TYPETAG=HeterogeneousNICCTpfaTypeTag
                                  ISOTHERMAL=0
              CMAKE_GUARD "( dune-alugrid_FOUND AND DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/co2nicc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2ni_tpfa-00018.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_tpfa test_co2ni_fv.input -Problem.Name co2ni_tpfa")

#install sources
install(FILES
heterogeneousco2tables.hh
heterogeneousproblem.hh
heterogeneousproblemni.hh
heterogeneousspatialparameters.hh
test_co2_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/co2)
