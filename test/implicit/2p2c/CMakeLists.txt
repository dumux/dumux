# build target for the twophase-twocomponent test problem
add_executable("test_2p2c" test_2p2c.cc)
target_link_libraries("test_2p2c" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_2p2c"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

