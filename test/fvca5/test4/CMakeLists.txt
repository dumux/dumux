add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("fvca5_test4" fvca5_test4.cc)
TARGET_LINK_LIBRARIES("fvca5_test4" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})
