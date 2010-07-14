# build the test for the property system
ADD_EXECUTABLE("test_properties" test_properties.cc)

# add required libraries and includes to the build flags 
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})
