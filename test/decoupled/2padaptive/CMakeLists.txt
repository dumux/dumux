# add build targets
add_executable("test_impesadaptive" test_impesadaptive.cc)
if(UG_FOUND)
  set_target_properties("test_impesadaptive"
    PROPERTIES
    COMPILE_FLAGS ${UG_CPPFLAGS})
  target_link_libraries("test_impesadaptive" ${DumuxLinkLibraries} ${UG_LIBS})
endif(UG_FOUND)

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_impesadaptive"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

