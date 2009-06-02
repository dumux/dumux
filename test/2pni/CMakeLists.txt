ADD_EXECUTABLE("new_test_2pni" new_test_2pni.cc)
TARGET_LINK_LIBRARIES("new_test_2pni" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})
