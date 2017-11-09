add_input_file_links()

dune_symlink_to_source_files(FILES runTestsMimetic.sh)
dune_symlink_to_source_files(FILES runTestsMPFAO.sh)
dune_symlink_to_source_files(FILES runTestsTPFA.sh)

add_dumux_test(test_spe10 test_spe10 test_spe10.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_spe10)

add_dumux_test(test_spe10_mpfao test_spe10_mpfao test_spe10_mpfao.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_spe10_mpfao)

add_dumux_test(test_spe10_tpfa test_spe10_tpfa test_spe10_tpfa.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_spe10_tpfa)