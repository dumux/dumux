add_input_file_links()

add_input_file_links()

add_dumux_test(test_1pdmp_mimetic test_1pdmp_mimetic test_1pdmp_mimetic.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1pdmp_mimetic)

add_dumux_test(test_1pdmp_mpfao test_1pdmp_mpfao test_1pdmp_mpfao.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1pdmp_mpfao)

add_dumux_test(test_1pdmpwell_mimetic test_1pdmpwell_mimetic test_1pdmpwell_mimetic.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1pdmpwell_mimetic)

add_dumux_test(test_1pdmpwell_mpfao test_1pdmpwell_mpfao test_1pdmpwell_mpfao.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1pdmpwell_mpfao)
