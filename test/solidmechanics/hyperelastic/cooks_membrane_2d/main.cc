// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
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

#include <dumux/assembly/fvassembler.hh>

#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::HyperelasticCooksMembrane;
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
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    // read λ sweep values — SPP-1748 Table 2.7
    const auto lambdaValues = getParam<std::vector<Scalar>>("SpatialParams.LambdaValues");
    const auto mu = getParam<Scalar>("SpatialParams.Mu");

    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    const GlobalPosition tipPoint{ 48.0, 60.0 };

    // store sweep results and print table at the end (avoids interleaving with DuMux output)
    struct Result { Scalar lambda, nu, uy; bool converged; };
    std::vector<Result> results;
    results.reserve(lambdaValues.size());

    // solution vector carried across the sweep: using the previous λ solution as
    // warm start for the next helps convergence for nearly-incompressible cases
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

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
            x = 0.0;  // reset to zero if diverged so the next step starts fresh
        }

        if (converged)
        {
            VtkOutputModule<GridVariables, SolutionVector> vtkWriter(
                *gridVariables, x,
                problem->name() + "_lambda" + std::to_string(static_cast<int>(lambda)));
            vtkWriter.addField(x, "d");
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

    // print results table — after all DuMux output is done
    // Reference values from SPP-1748 Table 2.9 (3D p-FEM, converged).
    // No 2D equivalent of Table 2.9 is provided in the paper; only Table 2.3 gives the
    // 2D reference for λ=432.099 (~10.596 mm). Figure 2.7(c) normalises 2D T2 results
    // against the 3D references — our 2D box results converge to ~93% of the 3D value,
    // consistent with the paper (active σ_zz in 3D plane strain increases deformation).
    // 3D converged references from SPP-1748 Table 2.9.
    // The paper provides no 2D sweep refs; 2D locking-free limit ≈ 93% of 3D ref.
    // Showing u_y / u_y_ref_3D makes locking visible: values below ~93% indicate locking.
    const std::vector<Scalar> refs3D = { 11.381, 11.156, 10.783, 10.746, 10.741 };

    std::cout << "\n"
              << "SPP-1748 Cook's membrane hyperelastic — 2D plane strain, ψ₁ (Table 2.7)\n"
              << "  u_y/u_y,ref normalised against 3D p-FEM (Table 2.9); locking-free 2D ≈ 93%\n"
              << std::string(62, '-') << "\n"
              << std::setw(14) << "λ [MPa]"
              << std::setw(10) << "ν"
              << std::setw(18) << "u_y A [mm]"
              << std::setw(12) << "u/u_ref"
              << "\n"
              << std::string(62, '-') << "\n";

    for (std::size_t i = 0; i < results.size(); ++i)
    {
        const auto& r = results[i];
        const Scalar ratio = (i < refs3D.size() && refs3D[i] > 0.0)
            ? r.uy / refs3D[i] : 0.0;
        std::cout << std::setw(14) << std::fixed << std::setprecision(3) << r.lambda
                  << std::setw(10) << std::setprecision(5) << r.nu
                  << std::setw(18) << std::setprecision(5) << r.uy
                  << std::setw(11) << std::setprecision(1) << 100.0*ratio << "%"
                  << (r.converged ? "" : "  [diverged]") << "\n";
    }
    std::cout << std::string(62, '-') << "\n";

    return 0;
}
