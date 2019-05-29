dune_add_test(SOURCES test_thermalconductivityjohansen.cc
              LABELS unit material
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzyData --delimiter " "
                       --files ${CMAKE_SOURCE_DIR}/test/references/thermalconductivityjohansen-reference.dat
                               ${CMAKE_CURRENT_BINARY_DIR}/johansen_lambda_eff.dat
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_thermalconductivityjohansen")

dune_add_test(SOURCES test_thermalconductivitysomerton.cc
              LABELS unit material
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzyData --delimiter " "
                      --files ${CMAKE_SOURCE_DIR}/test/references/thermalconductivitysomerton-reference.dat
                              ${CMAKE_CURRENT_BINARY_DIR}/somerton_lambda_eff.dat
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_thermalconductivitysomerton")

dune_add_test(SOURCES test_splinelaw.cc
              LABELS material
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_splinelaw)
