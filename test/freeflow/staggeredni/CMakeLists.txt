add_input_file_links()

add_dumux_test(test_channel_ni test_channel_ni test_channel_ni.cc ./test_channel_ni)

#install sources
install(FILES
test_channel_ni.cc
channeltestproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/staggeredni)
