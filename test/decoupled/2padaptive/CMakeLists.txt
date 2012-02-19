# add build targets
add_executable("test_impes_adaptive" test_impes_adaptive.cc)
if(UG_FOUND)
  set_target_properties("test_impes_adaptive"
    PROPERTIES
    COMPILE_FLAGS ${UG_CPPFLAGS})
  target_link_libraries("test_impes_adaptive" ${DumuxLinkLibraries} ${UG_LIBS})
endif(UG_FOUND)

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_impes_adaptive"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

