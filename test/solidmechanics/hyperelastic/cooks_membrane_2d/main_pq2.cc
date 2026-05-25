// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// SPP-1748 Cook's membrane — 2D T₂ (PQ2 single-field, pure quadratic displacement).
// Locking-free without a pressure split: P2 has enough displacement modes.
// Compare to MINI (PQ1Bubble+P1 mixed) and TH (PQ2+P1 mixed) for reference.
//
#include <config.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/assembly/assembler.hh>

#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include "properties_pq2.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::CooksMembranePQ2;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

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
    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits,
                                           LinearAlgebraTraitsFromAssembler<Assembler>>;
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    const auto lambdaValues = getParam<std::vector<Scalar>>("SpatialParams.LambdaValues");
    const auto mu = getParam<Scalar>("SpatialParams.Mu");

    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    const GlobalPosition tipPoint{ 48.0, 60.0 };

    struct Result { Scalar lambda, nu, uy; bool converged; };
    std::vector<Result> results;
    results.reserve(lambdaValues.size());

    for (const Scalar lambda : lambdaValues)
    {
        const Scalar nu = lambda / (2.0*(lambda + mu));
        problem->spatialParams().setLambda(lambda);

        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x);

        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
        auto linearSolver = std::make_shared<LinearSolver>();
        NewtonSolver nonLinearSolver(assembler, linearSolver);

        bool converged = true;
        try { nonLinearSolver.solve(x); }
        catch (const Dumux::NumericalProblem&)
        {
            converged = false;
            x = 0.0;
        }

        if (converged)
        {
            VtkOutputModule<GridVariables, SolutionVector>
                vtkWriter(*gridVariables, x,
                          "test_hyperelastic_cooks_membrane_pq2_lambda"
                          + std::to_string(static_cast<int>(lambda)));
            vtkWriter.addVolumeVariable([](const auto& v){
                return Dune::FieldVector<double,2>{ v.displacement(0), v.displacement(1) };
            }, "d");
            vtkWriter.write(1.0);
        }

        Scalar uy = 0.0;
        const auto tipElements = intersectingEntities(tipPoint, gridGeometry->boundingBoxTree());
        if (!tipElements.empty())
        {
            const auto element = gridGeometry->element(tipElements[0]);
            const auto elemSol = elementSolution(element, x, *gridGeometry);
            uy = evalSolution(element, element.geometry(), *gridGeometry, elemSol, tipPoint)[1];
        }

        results.push_back({ lambda, nu, uy, converged });
    }

    // 3D p-FEM references from SPP-1748 Table 2.9
    const std::vector<Scalar> refs3D = { 11.381, 11.156, 10.783, 10.746, 10.741 };
    std::cout << "\n"
              << "SPP-1748 Cook's membrane — 2D T₂/PQ2 single-field (locking-free without pressure split), ψ₁\n"
              << std::string(80, '-') << "\n"
              << std::setw(14) << "λ [MPa]"
              << std::setw(10) << "ν"
              << std::setw(18) << "u_y PQ2 [mm]"
              << std::setw(14) << "u/u_ref 3D"
              << "\n"
              << std::string(80, '-') << "\n";
    for (std::size_t i = 0; i < results.size(); ++i)
    {
        const auto& r = results[i];
        const Scalar ratio = (i < refs3D.size() && refs3D[i] > 0.0) ? r.uy/refs3D[i] : 0.0;
        std::cout << std::setw(14) << std::fixed << std::setprecision(3) << r.lambda
                  << std::setw(10) << std::setprecision(5) << r.nu
                  << std::setw(18) << std::setprecision(5) << r.uy
                  << std::setw(12) << std::setprecision(1) << 100.0*ratio << "%"
                  << (r.converged ? "" : "  [diverged]") << "\n";
    }
    std::cout << std::string(80, '-') << "\n"
              << "  PQ2 single-field is locking-free (pure P2 has enough displacement modes).\n";

    return 0;
}
