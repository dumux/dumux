add_input_file_links()

add_dumux_test(test_boxmpnc test_boxmpnc test_boxmpnc.cc
               ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/obstaclebox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/obstaclebox-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxmpnc")

add_dumux_test(test_ccmpnc test_ccmpnc test_ccmpnc.cc
               ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/obstaclecc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/obstaclecc-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccmpnc")

# build target for the one-phase Forchheimer test problem
add_dumux_test(test_forchheimer1p test_forchheimer1p test_forchheimer1p.cc
               ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/forchheimer1p-reference.vtp
                         ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer1p-00001.vtp
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer1p")

# build target for the two-phase Forchheimer test problem
add_dumux_test(test_forchheimer2p test_forchheimer2p test_forchheimer2p.cc
               ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/forchheimer2p-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer2p-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer2p")

# build target for the full kinetic test problem
add_dumux_test(test_boxmpnckinetic test_boxmpnckinetic test_boxmpnckinetic.cc
               ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/evaporationatmosphere-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/evaporationatmosphere-00005.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxmpnckinetic"
                 --zeroThreshold {"TnMinusTs":0.5,"TwMinusTn":0.5,"velocity_w_1":2e-12})

# build target for the energy kinetic test problem, two energy balance equations
add_dumux_test(test_boxmpncthermalnonequil test_boxmpncthermalnonequil test_boxmpncthermalnonequil.cc
               "${CMAKE_CURRENT_BINARY_DIR}/test_boxmpncthermalnonequil")

#install sources
install(FILES
combustionproblem1c.hh
combustionspatialparams.hh
evaporationatmosphereproblem.hh
evaporationatmospherespatialparams.hh
forchheimer1pproblem.hh
forchheimer2pproblem.hh
forchheimerspatialparams.hh
obstacleproblem.hh
obstaclespatialparams.hh
test_boxmpnc.cc
test_boxmpnckinetic.cc
test_boxmpncthermalnonequil.cc
test_ccmpnc.cc
test_forchheimer1p.cc
test_forchheimer2p.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/mpnc)
