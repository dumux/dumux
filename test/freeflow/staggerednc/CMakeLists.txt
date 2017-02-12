add_input_file_links()

add_dumux_test(test_channel_stokes_nc test_channel_stokes_nc test_channel.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/channel-stokes.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokes-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokes")


#install sources
install(FILES
channeltestproblem.hh
test_channel.cc
closedsystemtestproblem.hh
test_closedsystem.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/staggered)
