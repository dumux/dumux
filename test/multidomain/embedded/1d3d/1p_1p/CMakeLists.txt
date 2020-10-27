dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_average
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Average
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_tpfatpfa_average_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfatpfa_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_tpfatpfa_average_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfatpfa_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfatpfa_average params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfatpfa_average")

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_boxtpfa_average
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Average
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_tpfatpfa_average_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxtpfa_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_boxbox_average_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxtpfa_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxtpfa_average params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p1p_boxtpfa_average -Tissue.Grid.Cells \"19 19 19\"")

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfabox_average
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowBox COUPLINGMODE=Embedded1d3dCouplingMode::Average
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_boxbox_average_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfabox_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_tpfatpfa_average_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfabox_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfabox_average params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfabox_average")

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_boxbox_average
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=BloodFlowBox COUPLINGMODE=Embedded1d3dCouplingMode::Average
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_boxbox_average_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxbox_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_boxbox_average_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxbox_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxbox_average params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p1p_boxbox_average -Tissue.Grid.Cells \"19 19 19\"")

# make sure these configurations compile and run too
dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_surface
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Surface
              COMMAND ./test_md_embedded_1d3d_1p1p_tpfatpfa_surface
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfatpfa_surface)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_line
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Line
              COMMAND ./test_md_embedded_1d3d_1p1p_tpfatpfa_line
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfatpfa_line)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_kernel
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
              COMMAND ./test_md_embedded_1d3d_1p1p_tpfatpfa_kernel
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfatpfa_kernel)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_convergence
              LABELS multidomain multidomain_embedded 1p
              TARGET test_md_embedded_1d3d_1p1p_tpfatpfa_kernel
              COMMAND ./convergence.py
              CMAKE_GUARD dune-foamgrid_FOUND)

set_tests_properties(test_md_embedded_1d3d_1p1p_tpfatpfa_convergence PROPERTIES DEPENDS "test_md_embedded_1d3d_1p1p_tpfatpfa_surface;test_md_embedded_1d3d_1p1p_tpfatpfa_line")

dune_symlink_to_source_files(FILES "params.input" "convergence.py")
