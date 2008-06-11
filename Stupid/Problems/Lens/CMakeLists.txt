ADD_EXECUTABLE("StupidLens" LensMain.cc)

LINK_DIRECTORIES(${StupidLinkDirectories})
INCLUDE_DIRECTORIES(${StupidIncludeDirectories})
TARGET_LINK_LIBRARIES("StupidLens" ${StupidLinkLibraries})

