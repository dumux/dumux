# build the test for the property system
add_executable("test_propertysystem" test_propertysystem.cc)

# add required libraries and includes to the build flags 
include_directories(${DumuxIncludeDirectories})
