add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("test_diffusion_multiscale" test_diffusion_multiscale.cc)
TARGET_LINK_LIBRARIES("test_diffusion_multiscale" ${DumuxLinkLibraries})

ADD_EXECUTABLE("test_upscaledvelocity" test_upscaledvelocity.cc)
TARGET_LINK_LIBRARIES("test_upscaledvelocity" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

