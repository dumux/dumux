# build target for the twophase-twocomponent test problem
add_executable("test_boxmpnc" test_boxmpnc.cc)
target_link_libraries("test_boxmpnc" ${DumuxLinkLibraries})

add_executable("test_ccmpnc" test_ccmpnc.cc)
target_link_libraries("test_ccmpnc" ${DumuxLinkLibraries})

# build target for the one-phase forchheimer test problem
add_executable("test_forchheimer1p" test_forchheimer1p.cc)
target_link_libraries("test_forchheimer1p" ${DumuxLinkLibraries})

# build target for the two-phase forchheimer test problem
add_executable("test_forchheimer2p" test_forchheimer2p.cc)
target_link_libraries("test_forchheimer2p" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_boxmpnc"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")
