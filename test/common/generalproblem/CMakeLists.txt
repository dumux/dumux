add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# build target for the simple twophase lens problem in a general definition for both box and decoupled model
ADD_EXECUTABLE("test_generalproblem_2p" test_generalproblem_2p.cc)
TARGET_LINK_LIBRARIES("test_generalproblem_2p" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})
