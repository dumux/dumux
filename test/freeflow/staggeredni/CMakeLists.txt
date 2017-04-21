add_input_file_links()

add_dumux_test(test_channel_ni test_channel_ni test_channel_ni.cc ./test_channel_ni)

add_dumux_test(test_diffusion_ni test_diffusion_ni test_diffusion_ni.cc ./test_diffusion_ni)

set(CMAKE_BUILD_TYPE Debug)

#install sources
install(FILES
test_channel_ni.cc
channeltestproblem.hh
test_diffusion_ni.cc
diffusiontestproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/staggeredni)
