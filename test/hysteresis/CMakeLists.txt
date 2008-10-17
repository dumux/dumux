ADD_EXECUTABLE("test_pc_hysteresis_lenhard1" lenhardexp1main.cc)
ADD_EXECUTABLE("test_pc_hysteresis_lenhard2" lenhardexp2main.cc)
ADD_EXECUTABLE("test_pc_hysteresis_lens" lensmain.cc)

LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})
TARGET_LINK_LIBRARIES("test_pc_hysteresis_lenhard1" ${DumuxLinkLibraries})
TARGET_LINK_LIBRARIES("test_pc_hysteresis_lenhard2" ${DumuxLinkLibraries})
TARGET_LINK_LIBRARIES("test_pc_hysteresis_lens" ${DumuxLinkLibraries})

