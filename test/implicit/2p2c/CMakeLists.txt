# build target for the twophase-twocomponent test problem
add_executable("test_box2p2cCS2N2" test_box2p2cCS2N2.cc)
target_link_libraries("test_box2p2cCS2N2" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_box2p2cCS2N2"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

