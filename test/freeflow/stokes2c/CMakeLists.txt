# add build targets
add_executable("test_stokes2c" test_stokes2c.cc)
target_link_libraries("test_stokes2c" ${DumuxLinkLibraries})
add_dune_superlu_flags("test_stokes2c")

# add required libraries and includes to the build flags
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_stokes2c"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")
