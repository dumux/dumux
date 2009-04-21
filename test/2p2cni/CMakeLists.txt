add_definitions(-DENABLE_ALUGRID -DENABLE_UG)

#ADD_EXECUTABLE("test_2p2cni" test_2p2cni.cc)
#TARGET_LINK_LIBRARIES("test_2p2cni" ${DumuxLinkLibraries})

ADD_EXECUTABLE("new_test_2p2cni" new_test_2p2cni.cc)
TARGET_LINK_LIBRARIES("new_test_2p2cni" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_2p2cni"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

