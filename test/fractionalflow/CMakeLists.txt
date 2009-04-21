add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
#ADD_EXECUTABLE("test_fractionalflow" test_fractionalflow.cc)
#TARGET_LINK_LIBRARIES("test_fractionalflow" ${DumuxLinkLibraries})

#ADD_EXECUTABLE("test_2p2c_decoupled" test_2p2c_decoupled.cc)
#TARGET_LINK_LIBRARIES("test_2p2c_decoupled" ${DumuxLinkLibraries})

#ADD_EXECUTABLE("test_heterogeneoussoil" test_heterogeneoussoil.cc)
#TARGET_LINK_LIBRARIES("test_heterogeneoussoil" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_2p2c_decoupled"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

