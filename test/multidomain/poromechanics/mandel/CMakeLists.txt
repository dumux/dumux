add_subdirectory(analyticalsolution)

add_input_file_links()
dumux_add_test(NAME test_md_mandel
               SOURCES main.cc
               LABELS benchmark mandel 2d)
