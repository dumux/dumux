# build the test for the property system
ADD_EXECUTABLE("test_propertysystem" test_propertysystem.cc)

# add required libraries and includes to the build flags 
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})
