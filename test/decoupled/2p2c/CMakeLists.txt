add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("test_dec2p2c" test_dec2p2c.cc)
TARGET_LINK_LIBRARIES("test_dec2p2c" ${DumuxLinkLibraries})

ADD_EXECUTABLE("test_multiphysics2p2c" test_multiphysics2p2c.cc)
TARGET_LINK_LIBRARIES("test_multiphysics2p2c" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})
