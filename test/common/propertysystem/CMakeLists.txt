# build the test for the property system
add_executable("test_propertysystem" test_propertysystem.cc)
target_link_libraries("test_propertysystem" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})
