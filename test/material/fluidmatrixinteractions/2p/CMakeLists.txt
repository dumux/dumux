add_input_file_links()

add_dumux_test(test_thermalconductivityjohansen test_thermalconductivityjohansen test_thermalconductivityjohansen.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzyData --delimiter " "
                 --files ${CMAKE_SOURCE_DIR}/test/references/thermalconductivityjohansen-reference.dat
                         ${CMAKE_CURRENT_BINARY_DIR}/johansen_lambda_eff.dat
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_thermalconductivityjohansen")

add_dumux_test(test_thermalconductivitysomerton test_thermalconductivitysomerton test_thermalconductivitysomerton.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzyData --delimiter " "
                 --files ${CMAKE_SOURCE_DIR}/test/references/thermalconductivitysomerton-reference.dat
                         ${CMAKE_CURRENT_BINARY_DIR}/somerton_lambda_eff.dat
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_thermalconductivitysomerton")

#install sources
install(FILES
thermalconductivityspatialparams.hh
test_thermalconductivityjohansen.cc
test_thermalconductivitysomerton.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/fluidmatrixinteractions/2p)
