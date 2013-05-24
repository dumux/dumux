# set up project
project("dumux" C CXX)

# general stuff
cmake_minimum_required(VERSION 2.8.6)

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH})
message("dune-common_MODULE_PATH=${dune-common_MODULE_PATH}")
#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

#add_subdirectory("m4")
#add_subdirectory("dumux")
#add_subdirectory("doc")
add_subdirectory("test")
add_subdirectory("tutorial")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)



# # constexpr
# CHECK_CXX_SOURCE_COMPILES(
# "  int main(void) {
#     constexpr double g = 9.81;
#     return 0; }
# "  HAVE_CONSTEXPR
# )
# 
# # __attribute__((always_inline))
# CHECK_CXX_SOURCE_COMPILES(
# "
#    void __attribute__((always_inline)) foo(void) {}
#    int main(void) { foo(); return 0; };
# "  HAVE_ATTRIBUTE_ALWAYS_INLINE
# )



# ##############
# # Check for patched DUNE-PDELab
# include(CMakePushCheckState)
# cmake_push_check_state()
# set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -DHAVE_NULLPTR=${HAVE_NULLPTR}")
# set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${DumuxIncludeDirectories})
# set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${DumuxLinkLibraries})
# 
# CHECK_CXX_SOURCE_COMPILES("
#   #include <dune/pdelab/backend/istlvectorbackend.hh>
# 
#   int main(void)
#   {
#     Dune::PDELab::ISTLBlockVectorContainer
#       <std::vector<double>, double, 3> blockVectorContainer;
#     return blockVectorContainer.size();
#   }
# " DUNE_PDELAB_IS_PATCHED_FOR_DUMUX
# )
# cmake_pop_check_state()
