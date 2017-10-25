add_input_file_links()

dune_symlink_to_source_files(FILES runAllTestsAnisotropic.sh)
dune_symlink_to_source_files(FILES runAllTestsAnisotropicCase2.sh)
dune_symlink_to_source_files(FILES runAllTestsAnisotropicCase3.sh)
dune_symlink_to_source_files(FILES runTestAnisotropic.sh)
dune_symlink_to_source_files(FILES runTestAnisotropicCase2.sh)

add_dumux_test(test_1panisotropic_mimetic test_1panisotropic_mimetic test_1panisotropic_mimetic.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1panisotropic_mimetic)

add_dumux_test(test_1panisotropic_mpfao test_1panisotropic_mpfao test_1panisotropic_mpfao.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1panisotropic_mpfao)
