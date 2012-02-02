# build target for the simple twophase lens problem in a general definition for both box and decoupled model
add_executable("test_generalproblem_2p" test_generalproblem_2p.cc)
target_link_libraries("test_generalproblem_2p" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})
