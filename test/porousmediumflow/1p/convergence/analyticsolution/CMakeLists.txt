add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py")

# if UG is found, we use UG for all tests. Otherwise, we use
# YaspGrid for the tests on unstructured grids and enforce
# the tests on unstructured grids to be skipped
add_executable(test_1p_convergence_analytic_tpfa EXCLUDE_FROM_ALL main.cc)
add_executable(test_1p_convergence_analytic_mpfa EXCLUDE_FROM_ALL main.cc)
add_executable(test_1p_convergence_analytic_box EXCLUDE_FROM_ALL main.cc)

if (dune-uggrid_FOUND)
    target_compile_definitions(test_1p_convergence_analytic_tpfa PUBLIC "TYPETAG=OnePConvergenceTpfa" PUBLIC "GRIDTYPE=Dune::UGGrid<2>")
    target_compile_definitions(test_1p_convergence_analytic_mpfa PUBLIC "TYPETAG=OnePConvergenceMpfa" PUBLIC "GRIDTYPE=Dune::UGGrid<2>")
    target_compile_definitions(test_1p_convergence_analytic_box PUBLIC "TYPETAG=OnePConvergenceBox" PUBLIC "GRIDTYPE=Dune::UGGrid<2>")
else()
    target_compile_definitions(test_1p_convergence_analytic_tpfa PUBLIC "TYPETAG=OnePConvergenceTpfa" PUBLIC "GRIDTYPE=Dune::YaspGrid<2>")
    target_compile_definitions(test_1p_convergence_analytic_mpfa PUBLIC "TYPETAG=OnePConvergenceMpfa" PUBLIC "GRIDTYPE=Dune::YaspGrid<2>")
    target_compile_definitions(test_1p_convergence_analytic_box PUBLIC "TYPETAG=OnePConvergenceBox" PUBLIC "GRIDTYPE=Dune::YaspGrid<2>")
endif()

# using tpfa and structured grid with zero off-diagonal permeability entries
dumux_add_test(NAME test_1p_convergence_analytic_tpfa_structured
               TARGET test_1p_convergence_analytic_tpfa
               LABELS porousmediumflow 1p
               COMMAND ./convergencetest.py
               CMD_ARGS test_1p_convergence_analytic_tpfa
                        tpfa_structured params.input -Problem.C 0.0)

# using mpfa and structured grid
dumux_add_test(NAME test_1p_convergence_analytic_mpfa_structured
               TARGET test_1p_convergence_analytic_mpfa
               LABELS porousmediumflow 1p
               COMMAND ./convergencetest.py
               CMD_ARGS test_1p_convergence_analytic_mpfa mpfa_structured)

# using box and structured grid
dumux_add_test(NAME test_1p_convergence_analytic_box_structured
               TARGET test_1p_convergence_analytic_box
               LABELS porousmediumflow 1p
               COMMAND ./convergencetest.py
               CMD_ARGS test_1p_convergence_analytic_box box_structured)

# tests on ustructured grids in case UG is found
if (dune-uggrid_FOUND)
    # using mpfa and unstructured grid
    dumux_add_test(NAME test_1p_convergence_analytic_mpfa_unstructured
                   TARGET test_1p_convergence_analytic_mpfa
                   LABELS porousmediumflow 1p
                   COMMAND ./convergencetest.py
                   CMD_ARGS test_1p_convergence_analytic_mpfa mpfa_unstructured
                            params.input -Grid.File ../../incompressible/grids/randomlydistorted.dgf)

    # using box and unstructured grid
    dumux_add_test(NAME test_1p_convergence_analytic_box_unstructured
                   TARGET test_1p_convergence_analytic_box
                   LABELS porousmediumflow 1p
                   COMMAND ./convergencetest.py
                   CMD_ARGS test_1p_convergence_analytic_box box_unstructured
                            params.input -Grid.File ../../incompressible/grids/randomlydistorted.dgf)

# otherwise define test with the same name but enforce skip
else()
    dumux_add_test(NAME test_1p_convergence_analytic_mpfa_unstructured
                   SOURCES main.cc
                   CMAKE_GUARD dune-uggrid_FOUND
                   LABELS porousmediumflow 1p)

    dumux_add_test(NAME test_1p_convergence_analytic_box_unstructured
                   SOURCES main.cc
                   CMAKE_GUARD dune-uggrid_FOUND
                   LABELS porousmediumflow 1p)
endif()
