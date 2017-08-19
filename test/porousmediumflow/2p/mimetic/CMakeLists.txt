add_subdirectory("spe10")

add_input_file_links()

add_dumux_test(test_mimetic2p test_mimetic2p test_mimetic2p.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_mimetic2p)

add_dumux_test(test_wellmimetic2p test_wellmimetic2p test_wellmimetic2p.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_wellmimetic2p)

add_dumux_test(test_cp2p_mimetic test_cp2p_mimetic test_cp2p_mimetic.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_cp2p_mimetic)


#install sources
install(FILES
lensproblem.hh
lensspatialparams.hh
test_mimetic2p.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/2p/mimetic)
