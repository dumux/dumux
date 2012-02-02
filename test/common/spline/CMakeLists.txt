# add build targets
add_executable("test_spline" test_spline.cc)
target_link_libraries("test_spline" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

