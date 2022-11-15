dumux_add_test(NAME test_gas_injection
              SOURCES main.cc)

add_input_file_links()


add_compile_options(-ftemplate-backtrace-limit=0)