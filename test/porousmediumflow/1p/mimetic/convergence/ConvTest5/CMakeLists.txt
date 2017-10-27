add_input_file_links()

dune_symlink_to_source_files(FILES runTestMimetic.sh)
dune_symlink_to_source_files(FILES runTestMpfao.sh)

add_dumux_test(test5_conv_mimetic test5_conv_mimetic test5_conv_mimetic.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test5_conv_mimetic)

add_dumux_test(test5_conv_mpfao test5_conv_mpfao test5_conv_mpfao.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test5_conv_mpfao)
