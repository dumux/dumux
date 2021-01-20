dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES grids)

# using tpfa and analytical Jacobian
add_executable(test_1p_incompressible_tpfa_anadiff EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1p_incompressible_tpfa_anadiff PUBLIC "TYPETAG=OnePIncompressibleTpfa" "NUMDIFFMETHOD=DiffMethod::analytic")

# using mpfa and analytical Jacobian
add_executable(test_1p_incompressible_mpfa_anadiff EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1p_incompressible_mpfa_anadiff PUBLIC "TYPETAG=OnePIncompressibleMpfa" "NUMDIFFMETHOD=DiffMethod::analytic")

# using box and analytical Jacobian
add_executable(test_1p_incompressible_box_anadiff EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1p_incompressible_box_anadiff PUBLIC "TYPETAG=OnePIncompressibleBox" "NUMDIFFMETHOD=DiffMethod::analytic")

# using tpfa and analytical Jacobian
dumux_add_test(NAME test_1p_incompressible_tpfa
              TARGET test_1p_incompressible_tpfa_anadiff
              LABELS porousmediumflow 1p
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_anadiff params.input -Problem.Name test_1p_incompressible_tpfa")

# using tpfa, analytical Jacobian and quad precision
dumux_add_test(NAME test_1p_incompressible_tpfa_quad
              LABELS porousmediumflow 1p
              SOURCES main.cc
              CMAKE_GUARD HAVE_QUADMATH
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfaQuad NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_quad-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_quad params.input -Problem.Name test_1p_incompressible_tpfa_quad")

# using mpfa and analytical Jacobian
dumux_add_test(NAME test_1p_incompressible_mpfa
              TARGET test_1p_incompressible_mpfa_anadiff
              LABELS porousmediumflow 1p
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_mpfa-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_1p_mpfa_scvvelocities-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/mpfa_scv_velocities.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_mpfa_anadiff params.input -Problem.Name test_1p_incompressible_mpfa -IO.WriteMpfaVelocities true")

# using box and analytical Jacobian
dumux_add_test(NAME test_1p_incompressible_box
              TARGET test_1p_incompressible_box_anadiff
              LABELS porousmediumflow 1p
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_box_anadiff params.input -Problem.Name test_1p_incompressible_box")

# using tpfa and numeric differentiation
dumux_add_test(NAME test_1p_incompressible_tpfa_numdiff
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa NUMDIFFMETHOD=DiffMethod::numeric
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_numdiff-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_numdiff params.input
                                                                                  -Problem.Name test_1p_incompressible_tpfa_numdiff
                                                                                  -Assembly.NumericDifference.PriVarMagnitude 1e5")

# using box and numeric differentiation
dumux_add_test(NAME test_1p_incompressible_box_numdiff
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox NUMDIFFMETHOD=DiffMethod::numeric
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_box_numdiff-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_box_numdiff params.input
                                                                                  -Problem.Name test_1p_incompressible_box_numdiff
                                                                                  -Assembly.NumericDifference.PriVarMagnitude 1e5")

# using tpfa and analytical Jacobian with extrusion factor
dumux_add_test(NAME test_1p_incompressible_tpfa_extrude
              TARGET test_1p_incompressible_tpfa_anadiff
              LABELS porousmediumflow 1p
              COMMAND ./test_1p_incompressible_tpfa_anadiff
              CMD_ARGS params.input -Problem.Name test_1p_incompressible_tpfa_extrude
                                    -Problem.ExtrusionFactor 10
                                    -Vtk.AddVelocity 1
                                    -Problem.CheckIsConstantVelocity true
                                    -Problem.EnableGravity false)

# using mpfa and analytical Jacobian with extrusion factor
dumux_add_test(NAME test_1p_incompressible_mpfa_extrude
              TARGET test_1p_incompressible_mpfa_anadiff
              LABELS porousmediumflow 1p
              COMMAND ./test_1p_incompressible_mpfa_anadiff
              CMD_ARGS params.input -Problem.Name test_1p_incompressible_mpfa_extrude
                                    -Problem.ExtrusionFactor 10
                                    -Vtk.AddVelocity 1
                                    -Problem.CheckIsConstantVelocity true
                                    -Problem.EnableGravity false)

# using box and analytical Jacobian with extrusion factor
dumux_add_test(NAME test_1p_incompressible_box_extrude
              TARGET test_1p_incompressible_box_anadiff
              LABELS porousmediumflow 1p
              COMMAND ./test_1p_incompressible_box_anadiff
              CMD_ARGS params.input -Problem.Name test_1p_incompressible_box_extrude
                                    -Problem.ExtrusionFactor 10
                                    -Vtk.AddVelocity 1
                                    -Problem.CheckIsConstantVelocity true
                                    -Problem.EnableGravity false)

# using box and analytical Jacobian with extrusion factor on distorted grid
dumux_add_test(NAME test_1p_incompressible_box_extrude_distorted
              LABELS porousmediumflow 1p
              SOURCES main.cc
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ./test_1p_incompressible_box_extrude_distorted
              CMD_ARGS params.input -Problem.Name test_1p_incompressible_box_extrude_distorted
                                    -Problem.ExtrusionFactor 10
                                    -Vtk.AddVelocity 1
                                    -Problem.CheckIsConstantVelocity true
                                    -Problem.EnableGravity false
                                    -Grid.File ./grids/randomlydistorted.dgf)

# using mpfa and analytical Jacobian with extrusion factor on distorted grid
dumux_add_test(NAME test_1p_incompressible_mpfa_extrude_distorted
              LABELS porousmediumflow 1p
              SOURCES main.cc
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleMpfa NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ./test_1p_incompressible_mpfa_extrude_distorted
              CMD_ARGS params.input -Problem.Name test_1p_incompressible_mpfa_extrude_distorted
                                    -Problem.ExtrusionFactor 10
                                    -Vtk.AddVelocity 1
                                    -Problem.CheckIsConstantVelocity true
                                    -Problem.EnableGravity false
                                    -Grid.File ./grids/randomlydistorted.dgf)
