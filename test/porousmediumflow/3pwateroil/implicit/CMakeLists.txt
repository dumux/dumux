add_input_file_links()

add_dumux_test(test_box3pwateroil test_box3pwateroil test_box3pwateroil.cc "${CMAKE_CURRENT_BINARY_DIR}/test_box3pwateroil" "-TimeManager.TEnd" "43200" "-TimeManager.MaxTimeStepSize" "21600")

#install sources
install(FILES
3poilwatersagdproblem.hh
3poilwatersagdspatialparams.hh
test_box3pwateroil.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/3pwateroil)
