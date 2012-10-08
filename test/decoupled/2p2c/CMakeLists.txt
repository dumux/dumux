# add build targets
add_executable("test_adaptive2p2c" test_adaptive2p2c.cc)
target_link_libraries("test_adaptive2p2c" ${DumuxLinkLibraries})

add_executable("test_dec2p2c" test_dec2p2c.cc)
target_link_libraries("test_dec2p2c" ${DumuxLinkLibraries})

add_executable("test_multiphysics2p2c" test_multiphysics2p2c.cc)
target_link_libraries("test_multiphysics2p2c" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})
