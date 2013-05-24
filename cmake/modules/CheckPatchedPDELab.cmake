# Check for patched DUNE-PDELab
include(CMakePushCheckState)
cmake_push_check_state()
set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS} ${DUNE_FLAGS})
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${DUNE_INCLUDES})
set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${DUNE_LIBRARIES})

include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES("
  #include <dune/pdelab/backend/istlvectorbackend.hh>

  int main(void)
  {
    Dune::PDELab::ISTLBlockVectorContainer
      <std::vector<double>, double, 3> blockVectorContainer;
    return blockVectorContainer.size();
  }"
  DUNE_PDELAB_IS_PATCHED_FOR_DUMUX)
cmake_pop_check_state()
