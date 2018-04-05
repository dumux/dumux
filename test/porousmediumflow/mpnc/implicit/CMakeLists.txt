add_subdirectory(2p2ccomparison)
add_input_file_links()

dune_add_test(NAME test_mpnc_box
              SOURCES test_mpnc_obstacle_fv.cc
              COMPILE_DEFINITIONS TYPETAG=ObstacleBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/obstaclebox-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_boxmpnc-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_box test_mpnc.input -Problem.Name test_boxmpnc")

dune_add_test(NAME test_mpnc_tpfa
              SOURCES test_mpnc_obstacle_fv.cc
              COMPILE_DEFINITIONS TYPETAG=ObstacleCCTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/obstaclecc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ccmpnc-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_tpfa test_mpnc.input -Problem.Name test_ccmpnc")
# build target for the full kinetic test problem
dune_add_test(NAME test_mpnckinetic_box
              SOURCES test_boxmpnckinetic.cc
              COMPILE_DEFINITIONS TYPETAG=EvaporationAtmosphereBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/evaporationatmosphere-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_boxmpnckinetic-00011.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnckinetic_box test_boxmpnckinetic.input -Problem.Name test_boxmpnckinetic")

# build target for the energy kinetic test problem, two energy balance equations
dune_add_test(COMPILE_ONLY # since it currently fails miserably with very different results on different machines
              NAME test_mpncthermalnonequil_box
              SOURCES test_boxmpncthermalnonequil.cc
              COMPILE_DEFINITIONS TYPETAG=CombustionOneComponentBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/combustion-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_boxmpncthermalnonequil-00084.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpncthermalnonequil_box test_boxmpncthermalnonequil.input -Problem.Name test_boxmpncthermalnonequil")



#install sources
install(FILES
combustionproblem1c.hh
combustionfluidsystem.hh
combustionspatialparams.hh
evaporationatmosphereproblem.hh
evaporationatmospherespatialparams.hh
obstacleproblem.hh
obstaclespatialparams.hh
test_boxmpnc.cc
test_boxmpnckinetic.cc
test_boxmpncthermalnonequil.cc
test_ccmpnc.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/mpnc)
