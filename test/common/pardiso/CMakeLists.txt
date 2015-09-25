# add build targets
ADD_EXECUTABLE("test_pardiso" test_pardiso.cc)
TARGET_LINK_LIBRARIES("test_pardiso" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

