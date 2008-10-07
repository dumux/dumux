add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("fivespotcase1_ff"  fivespotcase1_ff.cc)
TARGET_LINK_LIBRARIES("fivespotcase1_ff" ${DumuxLinkLibraries})

ADD_EXECUTABLE("fivespotcase2_ff"  fivespotcase1_ff.cc)
TARGET_LINK_LIBRARIES("fivespotcase2_ff" ${DumuxLinkLibraries})


# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "fivespotcase1_ff"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

