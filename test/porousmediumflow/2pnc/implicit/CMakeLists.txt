dune_symlink_to_source_files(FILES test_2pnc_fv.input test_2pnc_diffusion.input)

# isothermal tests
dune_add_test(NAME test_2pnc_box
              SOURCES test_2pnc_fv.cc
              COMPILE_DEFINITIONS TYPETAG=FuelCellBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pncbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fuelcell_box-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_box -ParameterFile test_2pnc_fv.input -Problem.Name fuelcell_box")

dune_add_test(NAME test_2pnc_box_restart
              TARGET test_2pnc_box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pncbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fuelcell_box_restart-00003.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_box test_2pnc_fv.input -Problem.Name fuelcell_box_restart -TimeLoop.DtInitial 218.212 -Restart.Time 326.63 -Restart.File fuelcell_box-00012.vtu")

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_2pnc_box_restart PROPERTIES DEPENDS test_2pnc_box)

dune_add_test(NAME test_2pnc_tpfa
              SOURCES test_2pnc_fv.cc
              COMPILE_DEFINITIONS TYPETAG=FuelCellCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pnccc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fuelcell_tpfa-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_tpfa -ParameterFile test_2pnc_fv.input -Problem.Name fuelcell_tpfa")

dune_add_test(NAME test_cc2pnc_maxwellstefan
              SOURCES test_cc2pnc_diffusion.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPNCDiffusionTypeTag DIFFUSIONTYPE=MaxwellStefansLaw<TypeTag>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/2pncdiffusioncc-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_cc2pnc_maxwellstefan-00026.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pnc_maxwellstefan test_2pnc_diffusion.input -Problem.Name test_cc2pnc_maxwellstefan")

dune_add_test(NAME test_cc2pnc_fickslaw
              SOURCES test_cc2pnc_diffusion.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPNCDiffusionTypeTag DIFFUSIONTYPE=FicksLaw<TypeTag>
              COMMAND  ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/2pncdiffusioncc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_cc2pnc_fickslaw-00026.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pnc_fickslaw test_2pnc_diffusion.input -Problem.Name test_cc2pnc_fickslaw")

dune_add_test(NAME test_2pncni_box
              SOURCES test_2pnc_fv.cc
              COMPILE_DEFINITIONS TYPETAG=FuelCellNIBoxTypeTag NONISOTHERMAL=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pncboxni-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fuelcell_ni_box-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncni_box -ParameterFile test_2pnc_fv.input -Problem.Name fuelcell_ni_box")

#install sources
install(FILES
fuelcellproblem.hh
fuelcellniproblem.hh
fuelcellspatialparams.hh
test_2pnc_fv.cc
test_2pnc_fvni.cc
test_2pnc_fv.input
test_cc2pnc_diffusion.cc
maxwellstefandiffusion.hh
maxwellstefandiffusionspatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pnc)
set(CMAKE_BUILD_TYPE Release)
