dumux_add_test(NAME test_uq_mlmc
               SOURCES main.cc)

set(CMAKE_BUILD_TYPE RelWithDebInfo)
target_compile_options(test_uq_mlmc PRIVATE -fsanitize=address)
target_link_options(test_uq_mlmc PRIVATE -fsanitize=address)
