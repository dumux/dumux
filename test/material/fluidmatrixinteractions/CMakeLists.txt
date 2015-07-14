add_subdirectory("2p")

add_dumux_test(test_effectivediffusivitymillingtonquirk test_effectivediffusivitymillingtonquirk test_effectivediffusivitymillingtonquirk.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  exact
  ${CMAKE_SOURCE_DIR}/test/references/effectivediffusivitymillingtonquirk-reference.dat
  ${CMAKE_CURRENT_BINARY_DIR}/millingtonquirk_d_eff.dat
  ${CMAKE_CURRENT_BINARY_DIR}/test_effectivediffusivitymillingtonquirk
   -TimeManager.DtInitial 1 # [s]
   -TimeManager.TEnd 2 # [s]
   -Grid.File ${CMAKE_SOURCE_DIR}/test/implicit/2p/grids/test_2p.dgf
   -Diffusivity.File millingtonquirk
)

add_dumux_test(test_effectivediffusivityconstanttau test_effectivediffusivityconstanttau test_effectivediffusivityconstanttau.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  exact
  ${CMAKE_SOURCE_DIR}/test/references/effectivediffusivityconstanttau-reference.dat
  ${CMAKE_CURRENT_BINARY_DIR}/constanttau_d_eff.dat
  ${CMAKE_CURRENT_BINARY_DIR}/test_effectivediffusivityconstanttau
   -TimeManager.DtInitial 1 # [s]
   -TimeManager.TEnd 2 # [s]
   -Grid.File ${CMAKE_SOURCE_DIR}/test/implicit/2p/grids/test_2p.dgf
   -Diffusivity.File constanttau
   -tau 0.3
)
