# build target for the CO2 test problem
add_executable("test_boxco2" test_boxco2.cc)
target_link_libraries("test_boxco2" ${DumuxLinkLibraries})

add_executable("test_ccco2" test_ccco2.cc)
target_link_libraries("test_ccco2" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_boxco2"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

