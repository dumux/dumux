add_input_file_links()

add_dumux_test(test_densitydrivenflow test_densitydrivenflow test_densitydrivenflow.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-densitydriven-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_densitydrivenflow-00012.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_densitydrivenflow")

add_dumux_test(test_purediffusion test_channel test_channel.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-purediffusion.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_purediffusion-00013.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel test_purediffusion.input")

add_dumux_test(test_advection test_channel test_channel.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-advection.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_advection-00014.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel test_advection.input")


#install sources
install(FILES
test_densitydrivenflow.cc
densityflowproblem.hh
test_channel.cc
channeltestproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/staggerednc)
