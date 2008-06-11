ADD_EXECUTABLE("LenhardExp1" LenhardExp1Main.cc)
ADD_EXECUTABLE("LenhardExp2" LenhardExp2Main.cc)

LINK_DIRECTORIES(${StupidLinkDirectories})
INCLUDE_DIRECTORIES(${StupidIncludeDirectories})
TARGET_LINK_LIBRARIES("LenhardExp1" ${StupidLinkLibraries})
TARGET_LINK_LIBRARIES("LenhardExp2" ${StupidLinkLibraries})

