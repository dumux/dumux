// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Parallel == sequential for a 1p problem with a mixed Neumann/Dirichlet boundary on a
 *        distributed network. Uses a connected path (single component) so the mixed BC is
 *        well-posed; checks that the overlap fringe does not perturb the owned Neumann/Dirichlet
 *        assembly.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/type.hh>

#include "properties.hh"

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/parallel/distribute.hh>

#ifndef TYPETAG
#define TYPETAG OnePNetworkCCTpfa
#endif

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;

    auto solve = [](const Grid& grid) -> SolutionVector
    {
        const auto gridView = grid.leafGridView();
        auto gridGeometry = std::make_shared<GridGeometry>(gridView);
        auto problem = std::make_shared<Problem>(gridGeometry);
        SolutionVector x(gridGeometry->numDofs());
        problem->applyInitialSolution(x);
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x);
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
        auto linearSolver = std::make_shared<LinearSolver>(gridView, gridGeometry->dofMapper());
        NewtonSolver<Assembler, LinearSolver> newton(assembler, linearSolver);
        newton.solve(x);
        return x;
    };

    // a connected path of N line elements along x (single component, well-posed mixed BC:
    // Neumann at the x=0 end, Dirichlet at the x=N end via Problem.NeumannXMax)
    const int N = 16;
    auto makePath = []
    {
        Dune::GridFactory<Grid> factory;
        for (int i = 0; i <= N; ++i)
            factory.insertVertex(Dune::FieldVector<double, 3>{double(i), 0.0, 0.0});
        for (unsigned int i = 0; i < unsigned(N); ++i)
            factory.insertElement(Dune::GeometryTypes::line, {i, i + 1});
        return std::shared_ptr<Grid>(factory.createGrid().release());
    };

    // reference: full grid, sequential (self communicator)
    auto fullGrid = makePath();
    const auto xFull = solve(*fullGrid);

    // distribute by element centroid and solve in parallel
    Dune::Communication<Dune::MPIHelper::MPICommunicator> comm(mpiHelper.getCommunicator());
    std::vector<int> owners(fullGrid->leafGridView().indexSet().size(0));
    for (const auto& e : elements(fullGrid->leafGridView()))
        owners[fullGrid->leafGridView().indexSet().index(e)] =
            std::min(int(e.geometry().center()[0] * comm.size() / N), comm.size() - 1);

    auto local = Dune::FoamGridParallel::distributeFromRoot(fullGrid.get(), comm, owners, /*overlapWidth=*/1);
    const auto xLocal = solve(*local);

    // compare owned dofs against the reference by global id
    const auto localGV = local->leafGridView();
    const auto par = local->parallelData();
    static constexpr int dofCodim =
        (GridGeometry::discMethod == DiscretizationMethods::box) ? int(Grid::dimension) : 0;
    double maxDiff = 0.0, maxRef = 0.0;
    for (const auto& entity : entities(localGV, Dune::Codim<dofCodim>{}))
    {
        const auto idx = localGV.indexSet().index(entity);
        const auto gid = (dofCodim == 0) ? par->elementGlobalId(idx) : par->vertexGlobalId(idx);
        maxDiff = std::max(maxDiff, std::abs(xLocal[idx][0] - xFull[gid][0]));
        maxRef = std::max(maxRef, std::abs(xFull[gid][0]));
    }
    maxDiff = comm.max(maxDiff); maxRef = comm.max(maxRef);
    const double relErr = maxDiff / (maxRef > 0.0 ? maxRef : 1.0);

    const bool ok = relErr < 1e-7;
    if (comm.rank() == 0)
        std::cout << "1p Neumann/Dirichlet network: rel diff parallel vs sequential = " << relErr
                  << (ok ? "  PASS\n" : "  FAIL\n");
    return ok ? 0 : 1;
}
