add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

add_library("dumux" STATIC ./io/readstarformat.cc ./shapefunctions/CRshapefunctions.cc)

# add includes to the build flags 
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

