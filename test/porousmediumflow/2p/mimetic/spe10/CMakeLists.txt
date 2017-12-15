add_input_file_links()

dune_symlink_to_source_files(FILES runTestsMimetic.sh)
dune_symlink_to_source_files(FILES runTestsMPFAO.sh)
dune_symlink_to_source_files(FILES runTestsTPFA.sh)

dune_symlink_to_source_files(FILES runTestsMimeticConstDt.sh)
dune_symlink_to_source_files(FILES runTestsMPFAOConstDt.sh)
dune_symlink_to_source_files(FILES runTestsTPFAConstDt.sh)

add_dumux_test(test_spe10 test_spe10 test_spe10.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_spe10)

add_dumux_test(test_spe10_mpfao test_spe10_mpfao test_spe10_mpfao.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_spe10_mpfao)

add_dumux_test(test_spe10_tpfa test_spe10_tpfa test_spe10_tpfa.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_spe10_tpfa)

add_dumux_test(test_spe103D_tpfa test_spe103D_tpfa test_spe103D_tpfa.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_spe103D_tpfa)

add_dumux_test(test_spe103d_box test_spe103d_box test_spe103d_box.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_spe103d_box)