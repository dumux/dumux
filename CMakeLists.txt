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

message("<<<< DUNE_GRID_VERSION: ${DUNE_GRID_VERSION}")


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


# ##############
# # tell cmake that we've got a few subdirectories. (that's the
# # directories where the actual programs are)
# add_subdirectory("test")
# add_subdirectory("tutorial")
# 
# # copy the testing script
# make_directory(bin)
# make_directory(references)
# file(COPY test/references/ DESTINATION references)
# file(COPY bin/runTest.sh DESTINATION bin)
# file(COPY bin/fuzzycomparevtu.py DESTINATION bin)
# file(COPY test/implicit/1p/test_box1p.input DESTINATION test/implicit/1p)
# file(COPY test/implicit/1p/test_box1pwithamg.input DESTINATION test/implicit/1p)
# file(COPY test/implicit/1p2c/test_box1p2c.input DESTINATION test/implicit/1p2c)
# file(COPY test/implicit/2p/test_box2p.input DESTINATION test/implicit/2p)
# file(COPY test/implicit/2pni/test_box2pni.input DESTINATION test/implicit/2pni)
# file(COPY test/implicit/2p2c/test_box2p2c.input DESTINATION test/implicit/2p2c)
# file(COPY test/implicit/2p2cni/test_box2p2cni.input DESTINATION test/implicit/2p2cni)
# file(COPY test/implicit/2pdfm/test_2pdfm.input DESTINATION test/implicit/2pdfm)
# file(COPY test/implicit/3p3c/test_box3p3c.input DESTINATION test/implicit/3p3c)
# file(COPY test/implicit/3p3cni/test_box3p3cni.input DESTINATION test/implicit/3p3cni)
# file(COPY test/implicit/co2/test_ccco2.input DESTINATION test/implicit/co2)
# file(COPY test/implicit/co2/test_boxco2.input DESTINATION test/implicit/co2)
# file(COPY test/implicit/co2ni/test_ccco2ni.input DESTINATION test/implicit/co2ni)
# file(COPY test/implicit/co2ni/test_boxco2ni.input DESTINATION test/implicit/co2ni)
# file(COPY test/implicit/mpnc/test_boxmpnc.input DESTINATION test/implicit/mpnc)
# file(COPY test/implicit/mpnc/test_forchheimer1p.input DESTINATION test/implicit/mpnc)
# file(COPY test/implicit/mpnc/test_forchheimer2p.input DESTINATION test/implicit/mpnc)
# file(COPY test/implicit/richards/test_boxrichards.input DESTINATION test/implicit/richards)
# file(COPY test/implicit/1p/test_cc1p.input DESTINATION test/implicit/1p)
# file(COPY test/implicit/1p/test_cc1pwithamg.input DESTINATION test/implicit/1p)
# file(COPY test/implicit/1p2c/test_cc1p2c.input DESTINATION test/implicit/1p2c)
# file(COPY test/implicit/2p/test_cc2p.input DESTINATION test/implicit/2p)
# file(COPY test/implicit/2pni/test_cc2pni.input DESTINATION test/implicit/2pni)
# file(COPY test/implicit/2p2c/test_cc2p2c.input DESTINATION test/implicit/2p2c)
# file(COPY test/implicit/2p2cni/test_cc2p2cni.input DESTINATION test/implicit/2p2cni)
# file(COPY test/implicit/3p3c/test_cc3p3c.input DESTINATION test/implicit/3p3c)
# file(COPY test/implicit/3p3cni/test_cc3p3cni.input DESTINATION test/implicit/3p3cni)
# file(COPY test/implicit/mpnc/test_ccmpnc.input DESTINATION test/implicit/mpnc)
# file(COPY test/implicit/richards/test_ccrichards.input DESTINATION test/implicit/richards)
# file(COPY test/decoupled/1p/test_1p.input DESTINATION test/decoupled/1p/)
# file(COPY test/decoupled/2p/test_transport.input DESTINATION test/decoupled/2p)
# file(COPY test/decoupled/2p/test_impes.input DESTINATION test/decoupled/2p)
# file(COPY test/decoupled/2p/test_impeswithamg.input DESTINATION test/decoupled/2p)
# file(COPY test/decoupled/2p/test_mpfa2p.input DESTINATION test/decoupled/2p)
# file(COPY test/decoupled/2p/test_impesadaptive.input DESTINATION test/decoupled/2p)
# file(COPY test/decoupled/2p2c/test_adaptive2p2c.input DESTINATION test/decoupled/2p2c)
# file(COPY test/decoupled/2p2c/test_dec2p2c.input DESTINATION test/decoupled/2p2c)
# file(COPY test/common/generalproblem/test_generalproblem2p.input DESTINATION test/common/generalproblem)
# file(COPY tutorial/tutorial_coupled.input DESTINATION tutorial)
# file(COPY tutorial/tutorial_decoupled.input DESTINATION tutorial)
# if(SUPERLU_FOUND)
#   file(COPY test/freeflow/stokes/test_stokes.input DESTINATION test/freeflow/stokes)
#   file(COPY test/freeflow/stokes2c/test_stokes2c.input DESTINATION test/freeflow/stokes2c)
#   file(COPY test/freeflow/stokes2cni/test_stokes2cni.input DESTINATION test/freeflow/stokes2cni)
#   file(COPY test/freeflow/navierstokes/test_navierstokes.input DESTINATION test/freeflow/navierstokes)
# endif(SUPERLU_FOUND)
# 
# # set up CTest 
# enable_testing()
# include(CTest)
# add_test(test_boxco2      bin/runTest.sh references/co2box-reference.vtu            heterogeneousbox-00020.vtu            test/implicit/co2/test_boxco2 -Grid.File test/implicit/co2/grids/heterogeneousSmall.dgf) 
# add_test(test_propertysystem test/common/propertysystem/test_propertysystem)
# add_test(test_general_box    bin/runTest.sh references/generallens_box-reference.vtu       generallens_box-00003.vtu       test/common/generalproblem/test_generalproblem2p -ModelType Box       -TimeManager.TEnd 1e2 -TimeManager.DtInitial 2e1)
# add_test(test_general_dec    bin/runTest.sh references/generallens_decoupled-reference.vtu generallens_decoupled-00003.vtu test/common/generalproblem/test_generalproblem2p -ModelType Decoupled -TimeManager.TEnd 1e2 -TimeManager.DtInitial 2e1)
# add_test(test_spline         test/common/spline/test_spline)
# 
# add_test(test_fluidsystems    test/material/fluidsystems/test_fluidsystems)
# add_test(test_immiscibleflash test/material/immiscibleflash/test_immiscibleflash)
# add_test(test_ncpflash        test/material/ncpflash/test_ncpflash)
# add_test(test_pengrobinson    test/material/pengrobinson/test_pengrobinson)
# add_test(test_tabulation      test/material/tabulation/test_tabulation)
# 
# add_test(test_box1p       bin/runTest.sh references/1ptestbox-reference.vtu         1ptestbox-00002.vtu         test/implicit/1p/test_box1p             -Grid.File test/implicit/1p/grids/test_1p_2d.dgf               -TimeManager.TEnd 1       -TimeManager.DtInitial 1)
# add_test(test_cc1p        bin/runTest.sh references/1ptestcc-reference.vtu          1ptestcc-00002.vtu          test/implicit/1p/test_cc1p              -Grid.File test/implicit/1p/grids/test_1p_2d.dgf               -TimeManager.TEnd 1       -TimeManager.DtInitial 1)
# add_test(test_box1p2c     bin/runTest.sh references/outflowbox-reference.vtu        outflowbox-00010.vtu        test/implicit/1p2c/test_box1p2c         -Grid.File test/implicit/1p2c/grids/test_1p2c.dgf              -TimeManager.TEnd 100     -TimeManager.DtInitial 1)
# add_test(test_cc1p2c      bin/runTest.sh references/outflowcc-reference.vtu         outflowcc-00010.vtu         test/implicit/1p2c/test_cc1p2c          -Grid.File test/implicit/1p2c/grids/test_1p2c.dgf              -TimeManager.TEnd 100     -TimeManager.DtInitial 1)
# add_test(test_box2p       bin/runTest.sh references/lensbox-reference.vtu           lensbox-00009.vtu           test/implicit/2p/test_box2p             -Grid.File test/implicit/2p/grids/test_2p.dgf                  -TimeManager.TEnd 3000    -TimeManager.DtInitial 250)
# add_test(test_cc2p        bin/runTest.sh references/lenscc-reference.vtu            lenscc-00009.vtu            test/implicit/2p/test_cc2p              -Grid.File test/implicit/2p/grids/test_2p.dgf                  -TimeManager.TEnd 3000    -TimeManager.DtInitial 250)
# add_test(test_box2pni     bin/runTest.sh references/injection2pnibox-reference.vtu  injection2pnibox-00009.vtu  test/implicit/2pni/test_box2pni                                                                         -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
# add_test(test_cc2pni      bin/runTest.sh references/injection2pnicc-reference.vtu   injection2pnicc-00009.vtu   test/implicit/2pni/test_cc2pni                                                                         -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
# add_test(test_box2p2c     bin/runTest.sh references/injectionbox-reference.vtu      injectionbox-00009.vtu      test/implicit/2p2c/test_box2p2c         -Grid.File test/implicit/2p2c/grids/test_2p2c.dgf              -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
# add_test(test_cc2p2c      bin/runTest.sh references/injectioncc-reference.vtu       injectioncc-00009.vtu       test/implicit/2p2c/test_cc2p2c          -Grid.File test/implicit/2p2c/grids/test_2p2c.dgf              -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
# add_test(test_box2p2cni   bin/runTest.sh references/waterairbox-reference.vtu       waterairbox-00010.vtu       test/implicit/2p2cni/test_box2p2cni     -Grid.File test/implicit/2p2cni/grids/test_2p2cni.dgf          -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
# add_test(test_cc2p2cni    bin/runTest.sh references/wateraircc-reference.vtu        wateraircc-00010.vtu        test/implicit/2p2cni/test_cc2p2cni      -Grid.File test/implicit/2p2cni/grids/test_2p2cni.dgf          -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
# add_test(test_2pdfm       bin/runTest.sh references/2pdfm-reference.vtu             2pdfm-00011.vtu             test/implicit/2pdfm/test_2pdfm          -Grid.File test/implicit/2pdfm/grids/2pdfmartmesh.net)
# add_test(test_box3p3c     bin/runTest.sh references/infiltrationbox-reference.vtu   infiltrationbox-00007.vtu   test/implicit/3p3c/test_box3p3c         -Grid.File test/implicit/3p3c/grids/test_3p3c_coarse.dgf       -TimeManager.TEnd 8.64e5  -TimeManager.DtInitial 8.64e4)
# add_test(test_cc3p3c      bin/runTest.sh references/infiltrationcc-reference.vtu    infiltrationcc-00006.vtu    test/implicit/3p3c/test_cc3p3c          -Grid.File test/implicit/3p3c/grids/test_3p3c_coarse.dgf       -TimeManager.TEnd 8.64e5  -TimeManager.DtInitial 8.64e4)
# add_test(test_box3p3cni   bin/runTest.sh references/column3p3cnibox-reference.vtu   columnxylolbox-00060.vtu    test/implicit/3p3cni/test_box3p3cni     -Grid.File test/implicit/3p3cni/grids/column.dgf)
# add_test(test_cc3p3cni    bin/runTest.sh references/column3p3cnicc-reference.vtu    columnxylolcc-00054.vtu     test/implicit/3p3cni/test_cc3p3cni      -Grid.File test/implicit/3p3cni/grids/column.dgf)
# add_test(test_ccco2      bin/runTest.sh references/co2cc-reference.vtu            heterogeneouscc-00019.vtu            test/implicit/co2/test_ccco2 -Grid.File test/implicit/co2/grids/heterogeneousSmall.dgf) 
# add_test(test_ccco2ni    bin/runTest.sh references/co2nicc-reference.vtu          heterogeneousccni-00019.vtu          test/implicit/co2ni/test_ccco2ni -Grid.File test/implicit/co2ni/grids/heterogeneousSmall.dgf)
# add_test(test_boxco2ni    bin/runTest.sh references/co2nibox-reference.vtu          heterogeneousboxni-00020.vtu          test/implicit/co2ni/test_boxco2ni -Grid.File test/implicit/co2ni/grids/heterogeneousSmall.dgf)
# add_test(test_boxmpnc     bin/runTest.sh references/obstaclebox-reference.vtu       obstaclebox-00010.vtu       test/implicit/mpnc/test_boxmpnc         -Grid.File test/implicit/mpnc/grids/obstacle_24x16.dgf         -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
# add_test(test_ccmpnc      bin/runTest.sh references/obstaclecc-reference.vtu        obstaclecc-00010.vtu        test/implicit/mpnc/test_ccmpnc          -Grid.File test/implicit/mpnc/grids/obstacle_24x16.dgf         -TimeManager.TEnd 1e4     -TimeManager.DtInitial 250)
# add_test(test_forchheimer1p    bin/runTest.sh references/forchheimer1p-reference.vtp      test_forchheimer1p-00002.vtp       test/implicit/mpnc/test_forchheimer1p  -Grid.File test/implicit/mpnc/grids/forchheimer1d.dgf)
# add_test(test_forchheimer2p    bin/runTest.sh references/forchheimer2p-reference.vtu      test_forchheimer2p-00009.vtu       test/implicit/mpnc/test_forchheimer2p  -Grid.File test/implicit/mpnc/grids/obstacle_24x16.dgf)
# if(MPI_FOUND)
#   add_test(test_boxrichards bin/runTest.sh references/richardslensbox-reference-parallel.vtu s0002-p0000-richardslensbox-00008.vtu  ${MPIEXEC} -np 2 test/implicit/richards/test_boxrichards -Grid.File test/implicit/richards/grids/richardslens-24x16.dgf -TimeManager.TEnd 3000    -TimeManager.DtInitial 100)
# else(MPI_FOUND)
#   add_test(test_boxrichards bin/runTest.sh references/richardslensbox-reference.vtu   richardslensbox-00008.vtu   test/implicit/richards/test_boxrichards -Grid.File test/implicit/richards/grids/richardslens-24x16.dgf -TimeManager.TEnd 3000    -TimeManager.DtInitial 100)
#   add_test(test_ccrichards  bin/runTest.sh references/richardslenscc-reference.vtu    richardslenscc-00008.vtu    test/implicit/richards/test_ccrichards  -Grid.File test/implicit/richards/grids/richardslens-24x16.dgf -TimeManager.TEnd 3000    -TimeManager.DtInitial 100)
# endif(MPI_FOUND)
# 
# add_test(test_diffusion        bin/runTest.sh references/diffusion-reference.vtu             mimeticdiffusion-00001.vtu      test/decoupled/1p/test_diffusion 3)
# add_test(test_dec1p            bin/runTest.sh references/test_1p-reference.vtu               test_1p-00001.vtu               test/decoupled/1p/test_dec1p -ParameterFile test/decoupled/1p/test_1p.input)
# add_test(test_transport        bin/runTest.sh references/test_transport-reference.vtu        test_transport-00006.vtu        test/decoupled/2p/test_transport -Grid.File test/decoupled/2p/grids/test_transport.dgf)
# add_test(test_impes            bin/runTest.sh references/test_impes-reference.vtu            test_impes-00013.vtu            test/decoupled/2p/test_impes)
# add_test(test_mpfao2p          bin/runTest.sh references/test_mpfao2p-reference.vtu          test_mpfa2p-00007.vtu           test/decoupled/2p/test_mpfa2p -ModelType MPFAO)
# add_test(test_mpfal2p          bin/runTest.sh references/test_mpfal2p-reference.vtu          test_mpfa2p-00007.vtu           test/decoupled/2p/test_mpfa2p -ModelType MPFAL)
# add_test(test_mpfal2padaptive  bin/runTest.sh references/test_mpfal2padaptive-reference.vtu  test_mpfa2p-00007.vtu           test/decoupled/2p/test_mpfa2p -ModelType MPFALAdaptive)
# add_test(test_impesadaptive    bin/runTest.sh references/test_2padaptive-reference.vtu       test_2padaptive-00007.vtu       test/decoupled/2p/test_impesadaptive)
# add_test(test_adaptive2p2c     bin/runTest.sh references/test_adaptive2p2c-reference.vtu     test_adaptive2p2c-00008.vtu     test/decoupled/2p2c/test_adaptive2p2c)
# add_test(test_dec2p2c          bin/runTest.sh references/test_dec2p2c-reference.vtu          test_dec2p2c-00021.vtu          test/decoupled/2p2c/test_dec2p2c)
# add_test(test_multiphysics2p2c bin/runTest.sh references/test_multiphysics2p2c-reference.vtu test_multiphysics2p2c-00021.vtu test/decoupled/2p2c/test_multiphysics2p2c)
# 
# add_test(tutorial_coupled      tutorial/tutorial_coupled)
# add_test(tutorial_decoupled    tutorial/tutorial_decoupled)
# 
# if(DUNE_pdelab_FOUND)
#   add_test(test_box1pwithamg   bin/runTest.sh references/1ptestbox-reference.vtu             1ptestboxwithamg-00002.vtu      test/implicit/1p/test_box1pwithamg      -Grid.File test/implicit/1p/grids/test_1p_2d.dgf               -TimeManager.TEnd 1       -TimeManager.DtInitial 1)
#   add_test(test_cc1pwithamg    bin/runTest.sh references/1ptestcc-reference.vtu              1ptestccwithamg-00002.vtu       test/implicit/1p/test_cc1pwithamg       -Grid.File test/implicit/1p/grids/test_1p_2d.dgf               -TimeManager.TEnd 1       -TimeManager.DtInitial 1)
#   if(MPI_FOUND)
#     add_test(test_impeswithamg bin/runTest.sh references/test_impes-reference-parallel.vtu   s0002-p0001-test_impeswithamg-00013.vtu ${MPIEXEC} -np 2 test/decoupled/2p/test_impeswithamg -Grid.File test/decoupled/2p/grids/test_impeswithamg.dgf)
#   else(MPI_FOUND)
#     add_test(test_impeswithamg bin/runTest.sh references/test_impes-reference.vtu            test_impeswithamg-00013.vtu     test/decoupled/2p/test_impeswithamg     -Grid.File test/decoupled/2p/grids/test_impeswithamg.dgf)
#   endif(MPI_FOUND)
# endif(DUNE_pdelab_FOUND)
# 
# if(SUPERLU_FOUND)
#   add_test(test_stokes       bin/runTest.sh references/stokes-reference.vtu       stokes-00013.vtu       test/freeflow/stokes/test_stokes             -Grid.File test/freeflow/stokes/grids/test_stokes.dgf)
#   add_test(test_stokes2c     bin/runTest.sh references/stokes2c-reference.vtu     stokes2c-00007.vtu     test/freeflow/stokes2c/test_stokes2c         -Grid.File test/freeflow/stokes2c/grids/test_stokes2c.dgf)
#   add_test(test_stokes2cni   bin/runTest.sh references/stokes2cni-reference.vtu   stokes2cni-00008.vtu   test/freeflow/stokes2cni/test_stokes2cni     -Grid.File test/freeflow/stokes2cni/grids/test_stokes2cni.dgf)
#   add_test(test_navierstokes bin/runTest.sh references/navierstokes-reference.vtu navierstokes-00008.vtu test/freeflow/navierstokes/test_navierstokes -Grid.File test/freeflow/navierstokes/grids/test_navierstokes.dgf)
# endif(SUPERLU_FOUND)
# 
# add_test(test_restart  bin/runTest.sh references/co2box-reference.vtu            heterogeneousbox-00020.vtu            test/implicit/co2/test_boxco2 -Grid.File test/implicit/co2/grids/heterogeneousSmall.dgf -Restart 42008.6 -TimeManager.DtInitial 23452.1)
