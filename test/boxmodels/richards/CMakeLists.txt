if (MPI_FOUND)
  # tell UG that the model is parallelized
  add_definitions(-DModelP)
endif (MPI_FOUND)

# build target for the twophase-twocomponent test problem
ADD_EXECUTABLE("test_richards" test_richards.cc)
TARGET_LINK_LIBRARIES("test_richards" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_richards"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

