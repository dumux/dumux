dune_add_test(NAME test_embedded_1d3d_1p1p_tpfatpfa_average
              SOURCES test_1p_1p.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCCTypeTag LOWDIMTYPETAG=BloodFlowCCTypeTag COUPLINGMODE=EmbeddedCouplingMode::average
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p1p_tpfatpfa_average_1d.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_tpfatpfa_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p1p_tpfatpfa_average_3d.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_tpfatpfa_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_tpfatpfa_average test_1p_1p.input \
                                   -Problem.Name test_embedded_1d3d_1p1p_tpfatpfa_average")

dune_add_test(NAME test_embedded_1d3d_1p1p_boxtpfa_average
              SOURCES test_1p_1p.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueBoxTypeTag LOWDIMTYPETAG=BloodFlowCCTypeTag COUPLINGMODE=EmbeddedCouplingMode::average
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p1p_tpfatpfa_average_1d.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_boxtpfa_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p1p_boxbox_average_3d.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_boxtpfa_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_boxtpfa_average test_1p_1p.input \
                                   -Problem.Name test_embedded_1d3d_1p1p_boxtpfa_average -Tissue.Grid.Cells \"19 19 19\"")

dune_add_test(NAME test_embedded_1d3d_1p1p_tpfabox_average
              SOURCES test_1p_1p.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCCTypeTag LOWDIMTYPETAG=BloodFlowBoxTypeTag COUPLINGMODE=EmbeddedCouplingMode::average
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p1p_boxbox_average_1d.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_tpfabox_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p1p_tpfatpfa_average_3d.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_tpfabox_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_tpfabox_average test_1p_1p.input \
                                   -Problem.Name test_embedded_1d3d_1p1p_tpfabox_average")

dune_add_test(NAME test_embedded_1d3d_1p1p_boxbox_average
              SOURCES test_1p_1p.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueBoxTypeTag LOWDIMTYPETAG=BloodFlowBoxTypeTag COUPLINGMODE=EmbeddedCouplingMode::average
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p1p_boxbox_average_1d.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_boxbox_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p1p_boxbox_average_3d.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_boxbox_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p1p_boxbox_average test_1p_1p.input \
                                   -Problem.Name test_embedded_1d3d_1p1p_boxbox_average -Tissue.Grid.Cells \"19 19 19\"")

# make sure these configurations compile and run too
dune_add_test(NAME test_embedded_1d3d_1p1p_tpfatpfa_cylindersources
              SOURCES test_1p_1p.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCCTypeTag LOWDIMTYPETAG=BloodFlowCCTypeTag COUPLINGMODE=EmbeddedCouplingMode::cylindersources
              COMMAND ./test_embedded_1d3d_1p1p_tpfatpfa_cylindersources
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS test_1p_1p.input)

dune_add_test(NAME test_embedded_1d3d_1p1p_tpfatpfa_line
              SOURCES test_1p_1p.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCCTypeTag LOWDIMTYPETAG=BloodFlowCCTypeTag COUPLINGMODE=EmbeddedCouplingMode::line
              COMMAND ./test_embedded_1d3d_1p1p_tpfatpfa_line
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS test_1p_1p.input)

dune_add_test(NAME test_embedded_1d3d_1p1p_tpfatpfa_kernel
              SOURCES test_1p_1p.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCCTypeTag LOWDIMTYPETAG=BloodFlowCCTypeTag COUPLINGMODE=EmbeddedCouplingMode::kernel
              COMMAND ./test_embedded_1d3d_1p1p_tpfatpfa_kernel
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS test_1p_1p.input)

dune_symlink_to_source_files(FILES "grids" "test_1p_1p.input")
set(CMAKE_BUILD_TYPE Release)
