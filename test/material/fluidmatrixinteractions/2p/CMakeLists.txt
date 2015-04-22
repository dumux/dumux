add_dumux_test(test_thermalconductivityjohansen test_thermalconductivityjohansen test_thermalconductivityjohansen.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  exact
  ${CMAKE_SOURCE_DIR}/test/references/thermalconductivityjohansen-reference.dat
  ${CMAKE_CURRENT_BINARY_DIR}/8_lambda_eff_johansen.dat
  ${CMAKE_CURRENT_BINARY_DIR}/test_thermalconductivityjohansen
   -TimeManager.DtInitial 1 # [s]
   -TimeManager.TEnd 2 # [s]
   -Grid.File ${CMAKE_SOURCE_DIR}/test/implicit/2p/grids/test_2p.dgf
   -Conductivity.File lambda_eff_johansen
)

add_dumux_test(test_thermalconductivitysomerton test_thermalconductivitysomerton test_thermalconductivitysomerton.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  exact
  ${CMAKE_SOURCE_DIR}/test/references/thermalconductivitysomerton-reference.dat
  ${CMAKE_CURRENT_BINARY_DIR}/8_lambda_eff_somerton.dat
  ${CMAKE_CURRENT_BINARY_DIR}/test_thermalconductivitysomerton
   -TimeManager.DtInitial 1 # [s]
   -TimeManager.TEnd 2 # [s]
   -Grid.File ${CMAKE_SOURCE_DIR}/test/implicit/2p/grids/test_2p.dgf
   -Conductivity.File lambda_eff_somerton
)
