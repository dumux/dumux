dune_symlink_to_source_files(FILES "grids" "test_facetcoupling_1p1p.input")

dune_add_test(NAME test_facetcoupling_tpfa_1p1p
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              SOURCES test_facetcoupling_tpfa_1p1p.cc
              COMMAND ./test_facetcoupling_tpfa_1p1p
              CMD_ARGS test_facetcoupling_1p1p.input)

dune_add_test(NAME test_facetcoupling_tpfa_1p1p_dirichlet
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              TARGET test_facetcoupling_tpfa_1p1p
              COMMAND ./test_facetcoupling_tpfa_1p1p
              CMD_ARGS test_facetcoupling_1p1p.input
                       -Problem.UseInteriorDirichletBCs true
                       -Bulk.Problem.Name test_1p_1p_bulk_dirichlet
                       -LowDim.Problem.Name test_1p_1p_lowdim_dirichlet)

dune_add_test(NAME test_facetcoupling_tpfa_1p1p_barrier
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              TARGET test_facetcoupling_tpfa_1p1p
              COMMAND ./test_facetcoupling_tpfa_1p1p
              CMD_ARGS test_facetcoupling_1p1p.input
                       -Problem.UseAnalyticProblem false
                       -Bulk.Problem.Name test_1p_1p_bulk_barrier
                       -LowDim.Problem.Name test_1p_1p_lowdim_barrier
                       -LowDim.Problem.UseDirichletBCs false)

dune_add_test(NAME test_facetcoupling_tpfa_1p1p_barrier_dirichlet
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              TARGET test_facetcoupling_tpfa_1p1p
              COMMAND ./test_facetcoupling_tpfa_1p1p
              CMD_ARGS test_facetcoupling_1p1p.input
                       -Problem.UseAnalyticProblem false
                       -Problem.UseInteriorDirichletBCs true
                       -Bulk.Problem.Name test_1p_1p_bulk_barrier_dirichlet
                       -LowDim.Problem.Name test_1p_1p_lowdim_barrier_dirichlet
                       -LowDim.Problem.UseDirichletBCs false)

dune_add_test(NAME test_facetcoupling_tpfa_1p1p_boundary
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              TARGET test_facetcoupling_tpfa_1p1p
              COMMAND ./test_facetcoupling_tpfa_1p1p
              CMD_ARGS test_facetcoupling_1p1p.input
                       -Problem.UseAnalyticProblem false
                       -Bulk.Problem.Name test_1p_1p_bulk_boundary
                       -LowDim.Problem.Name test_1p_1p_lowdim_boundary
                       -Grid.File grids/hybridgrid_boundary.msh)

dune_add_test(NAME test_facetcoupling_tpfa_1p1p_boundary_dirichlet
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              TARGET test_facetcoupling_tpfa_1p1p
              COMMAND ./test_facetcoupling_tpfa_1p1p
              CMD_ARGS test_facetcoupling_1p1p.input
                       -Problem.UseAnalyticProblem false
                       -Problem.UseInteriorDirichletBCs true
                       -Bulk.Problem.Name test_1p_1p_bulk_boundary_dirichlet
                       -LowDim.Problem.Name test_1p_1p_lowdim_boundary_dirichlet
                       -Grid.File grids/hybridgrid_boundary.msh)

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
test_facetcoupling_tpfa_1p1p.cc
bulkProblem.hh
lowdimproblem.hh
spatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/facet/1p_1p/analytical)
