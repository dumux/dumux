# add build targets
ADD_EXECUTABLE("test_spline" test_spline.cc)
TARGET_LINK_LIBRARIES("test_spline" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

