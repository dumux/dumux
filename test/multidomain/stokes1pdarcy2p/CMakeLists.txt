add_input_file_links()
dune_symlink_to_source_files(FILES references)

# test: test_conservation2p_2d
add_dumux_test(test_conservation2p_2d test_conservation2p_2d test_conservation.cc
  python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
    --script fuzzy
    --command "${CMAKE_CURRENT_BINARY_DIR}/test_conservation2p_2d"
    --files ${CMAKE_CURRENT_SOURCE_DIR}/references/conservation2p-stokes-reference.vtu
            ${CMAKE_CURRENT_BINARY_DIR}/conservation2p_ref_stokes-00001.vtu
            ${CMAKE_CURRENT_SOURCE_DIR}/references/conservation2p-darcy-reference.vtu
            ${CMAKE_CURRENT_BINARY_DIR}/conservation2p_ref_darcy-00001.vtu
    --zeroThreshold {"velocityN_0":1e-9}
  )

# test: test_conservation2p_1d
add_dumux_test(test_conservation2p_1d test_conservation2p_1d test_conservation.cc
  python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
    --script fuzzy
    --command "${CMAKE_CURRENT_BINARY_DIR}/test_conservation2p_1d"
    --files ${CMAKE_CURRENT_SOURCE_DIR}/references/conservation2p_1d-stokes-reference.vtp
            ${CMAKE_CURRENT_BINARY_DIR}/conservation2p_1d-ff-00030.vtp
            ${CMAKE_CURRENT_SOURCE_DIR}/references/conservation2p_1d-darcy-reference.vtp
            ${CMAKE_CURRENT_BINARY_DIR}/conservation2p_1d-pm-00030.vtp
  )
target_compile_definitions(test_conservation2p_1d PUBLIC "DUMUX_MULTIDOMAIN_DIM=1")
