# add build targets
add_executable("test_navierstokes" test_navierstokes.cc)
if(SUPERLU_FOUND AND (UG_FOUND OR ALUGrid_FOUND))
  set_target_properties("test_navierstokes"
    PROPERTIES
    COMPILE_FLAGS "${SUPERLU_CPPFLAGS} ${UG_CPPFLAGS} ${ALUGrid_CPPFLAGS}")
  target_link_libraries("test_navierstokes" 
    ${DumuxLinkLibraries}
    ${SUPERLU_LIBS}
    ${UG_LIBS}
    ${ALUGrid_LIBS})
endif(SUPERLU_FOUND AND (UG_FOUND OR ALUGrid_FOUND))

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
