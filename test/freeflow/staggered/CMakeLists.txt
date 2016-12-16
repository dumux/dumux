add_input_file_links()

add_dumux_test(test_liddrivencavity test_closedsystem test_closedsystem.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/liddrivencavity-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_liddrivencavity-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_closedsystem test_liddrivencavity.input")

add_dumux_test(test_hydrostaticpressure test_closedsystem test_closedsystem.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/stokeshydrostaticpressure-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_hydrostaticpressure-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_closedsystem test_hydrostaticpressure.input"
                 --zeroThreshold {"v_x_neg":1e-15,"v_x_pos":1e-15,"v_y_neg":1e-15,"v_y_pos":1e-15,"velocity_0":1e-15,"velocity_1":1e-15})

add_dumux_test(test_channel_stokes test_channel_stokes test_channel.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/channel-stokes.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokes-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokes")

add_dumux_test(test_channel_navierstokes test_channel_navierstokes test_channel.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/channel-navierstokes-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_channel_navierstokes-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_navierstokes")
target_compile_definitions(test_channel_navierstokes PUBLIC "ENABLE_NAVIERSTOKES=1")


#install sources
install(FILES
channeltestproblem.hh
test_channel.cc
closedsystemtestproblem.hh
test_closedsystem.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/staggered)
