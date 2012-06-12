# add build targets
add_executable("test_navierstokes" test_navierstokes.cc)
target_link_libraries("test_navierstokes" ${DumuxLinkLibraries})
set_target_properties("test_navierstokes"
  PROPERTIES
  COMPILE_FLAGS "${SUPERLU_CPPFLAGS} ${UG_CPPFLAGS}")
target_link_libraries("test_navierstokes" ${SUPERLU_LIBS})
target_link_libraries("test_navierstokes" ${UG_LIBS})

# set_target_properties("test_navierstokes"
#   PROPERTIES
#   COMPILE_FLAGS ${ALUGrid_CPPFLAGS})
# target_link_libraries("test_navierstokes" ${ALUGrid_LIBS})

# add required libraries and includes to the build flags
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_navierstokes"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")
