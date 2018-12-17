dune_symlink_to_source_files(FILES "params.input")

# using tpfa
add_executable(test_2p_incompressible_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_2p_incompressible_tpfa PUBLIC TYPETAG=TwoPIncompressibleTpfa)

dune_add_test(NAME test_2p_incompressible_tpfa
              TARGET test_2p_incompressible_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_incompressible_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_tpfa-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_tpfa params.input -Problem.Name test_2p_incompressible_tpfa")

# using tpfa
dune_add_test(NAME test_2p_incompressible_tpfa_restart
              TARGET test_2p_incompressible_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_incompressible_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_tpfa_restart-00003.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_tpfa params.input -Problem.Name test_2p_incompressible_tpfa_restart  -TimeLoop.DtInitial 526.62 -Restart.Time 1859.95 -Restart.File test_2p_incompressible_tpfa-00005.vtu")

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_2p_incompressible_tpfa_restart PROPERTIES DEPENDS test_2p_incompressible_tpfa)

# using box
dune_add_test(NAME test_2p_incompressible_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_incompressible_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box params.input -Problem.Name test_2p_incompressible_box")

# using box with interface solver
dune_add_test(NAME test_2p_incompressible_box_ifsolver
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleBox
              COMPILE_DEFINITIONS ENABLEINTERFACESOLVER=true
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_incompressible_box_ifsolver-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box_ifsolver-00018.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box_ifsolver params.input
                                                                                                  -Problem.Name test_2p_incompressible_box_ifsolver
                                                                                                  -Problem.UseNonConformingOutput true")

# using tpfa with an oil-wet lens
dune_add_test(NAME test_2p_incompressible_tpfa_oilwet
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_incompressible_tpfa_oilwet-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_tpfa_oilwet-00009.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_tpfa_oilwet params.input
                                                                                                 -Problem.Name test_2p_incompressible_tpfa_oilwet
                                                                                                 -Problem.EnableGravity false
                                                                                                 -SpatialParams.LensIsOilWet true")

# using mpfa
dune_add_test(NAME test_2p_incompressible_mpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleMpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_incompressible_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_mpfa-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_mpfa params.input -Problem.Name test_2p_incompressible_mpfa")
