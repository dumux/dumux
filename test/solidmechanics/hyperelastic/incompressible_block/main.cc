// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::HyperelasticIncompressibleBlock;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    vtkWriter.addField(x, "d");
    vtkWriter.write(0.0);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // Pseudo-static load stepping: increase load from 0 to q in nLoadSteps steps.
    // Each step uses the previous converged solution as initial guess.
    const int nLoadSteps = getParam<int>("LoadStepping.NSteps", 10);
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    std::cout << "Starting pseudo-static load stepping with " << nLoadSteps << " steps\n";

    for (int step = 1; step <= nLoadSteps; ++step)
    {
        const Scalar loadFactor = static_cast<Scalar>(step) / nLoadSteps;
        problem->setLoadFactor(loadFactor);

        std::cout << "Load step " << step << "/" << nLoadSteps
                  << "  (factor = " << loadFactor << ")\n";

        nonLinearSolver.solve(x);

        // Re-initialise grid variables with the new solution so the next step
        // starts from the correct state.
        gridVariables->update(x);
        vtkWriter.write(static_cast<double>(step));
    }

    // Point P = centre of the load patch = symmetry corner of the quarter model
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    const GlobalPosition P{ 50.0, 50.0, 50.0 };
    const auto tipElements = intersectingEntities(P, gridGeometry->boundingBoxTree());
    if (!tipElements.empty())
    {
        const auto element = gridGeometry->element(tipElements[0]);
        const auto elemSol = elementSolution(element, x, *gridGeometry);
        const auto u = evalSolution(element, element.geometry(), *gridGeometry, elemSol, P);
        std::cout << "Displacement at P (50,50,50):  u_z = " << u[2]
                  << " mm  |u_z| = " << std::abs(u[2]) << " mm\n"
                  << "  SPP-1748 reference (standard Q1 FEM, H1): |u_z| ≈ 7.78 mm\n"
                  << "  Note: DuMux box FVM differs from Galerkin Q1 FEM.\n";
    }

    return 0;
}
