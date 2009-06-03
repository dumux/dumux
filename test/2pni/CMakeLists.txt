ADD_EXECUTABLE("test_2pni" test_2pni.cc)
TARGET_LINK_LIBRARIES("test_2pni" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})
