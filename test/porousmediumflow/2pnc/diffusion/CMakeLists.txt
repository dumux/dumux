dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_2pnc_maxwellstefan_tpfa
              LABELS porousmediumflow 2pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPNCDiffusion DIFFUSIONTYPE=MaxwellStefansLaw<TypeTag>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2pnc_diffusion_tpfa-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_maxwellstefan_tpfa-00027.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_maxwellstefan_tpfa params.input -Problem.Name test_2pnc_maxwellstefan_tpfa")

dumux_add_test(NAME test_2pnc_fickslaw_tpfa
              LABELS porousmediumflow 2pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPNCDiffusion DIFFUSIONTYPE=FicksLaw<TypeTag>
              COMMAND  ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pnc_diffusion_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_fickslaw_tpfa-00027.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_fickslaw_tpfa params.input -Problem.Name test_2pnc_fickslaw_tpfa")
