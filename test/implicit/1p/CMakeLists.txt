# build target for the onephase-onecomponent test problem
add_executable("test_box1p" test_box1p.cc)
target_link_libraries("test_box1p" ${DumuxLinkLibraries})

add_executable("test_cc1p" test_cc1p.cc)
target_link_libraries("test_cc1p" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_box1p"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

