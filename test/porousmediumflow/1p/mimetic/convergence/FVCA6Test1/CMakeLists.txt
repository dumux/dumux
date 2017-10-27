add_input_file_links()

dune_symlink_to_source_files(FILES runTestMimetic.sh)
dune_symlink_to_source_files(FILES runTestMpfao.sh)

add_dumux_test(test1_conv_mimetic test1_conv_mimetic test1_conv_mimetic.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test1_conv_mimetic)

add_dumux_test(test1_conv_mpfao test1_conv_mpfao test1_conv_mpfao.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test1_conv_mpfao)
