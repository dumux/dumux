// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// SPP-1748 Cook's membrane — 3D MINI (PQ1Bubble displacement + P1 pressure).
// Direct comparison to SPP-1748 Table 2.9 (3D p-FEM reference values).
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

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/assembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include "mixed/properties_mixed3d.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using MomTypeTag = Properties::TTag::CooksMembraneMixedMomentum3D;
    using PresTypeTag = Properties::TTag::CooksMembraneMixedPressure3D;
    using MDTraits = Dumux::Properties::CooksMembraneMixedMDTraits3D;
    using Scalar = typename MDTraits::Scalar;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<MomTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using MomGG  = GetPropType<MomTypeTag,  Properties::GridGeometry>;
    using PresGG = GetPropType<PresTypeTag, Properties::GridGeometry>;
    auto momGG  = std::make_shared<MomGG>(leafGridView);
    auto presGG = std::make_shared<PresGG>(leafGridView);

    using CouplingManager = HyperelasticVolIsoCouplingManager<MDTraits>;
    auto couplingManager = std::make_shared<CouplingManager>(momGG, presGG);

    using MomProblem  = GetPropType<MomTypeTag,  Properties::Problem>;
    using PresProblem = GetPropType<PresTypeTag, Properties::Problem>;
    auto momProblem  = std::make_shared<MomProblem>(momGG,   couplingManager);
    auto presProblem = std::make_shared<PresProblem>(presGG, couplingManager);

    using SolVec = typename MDTraits::SolutionVector;
    SolVec x;
    x[CouplingManager::momentumIdx].resize(momGG->numDofs());
    x[CouplingManager::pressureIdx].resize(presGG->numDofs());
    x = 0.0;

    using MomGV  = GetPropType<MomTypeTag,  Properties::GridVariables>;
    using PresGV = GetPropType<PresTypeTag, Properties::GridVariables>;
    using Assembler = Dumux::Experimental::MultiDomainAssembler<MDTraits, CouplingManager, DiffMethod::numeric>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits,
                                           LinearAlgebraTraitsFromAssembler<Assembler>>;
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;

    const auto lambdaValues = getParam<std::vector<Scalar>>("SpatialParams.LambdaValues");
    const auto mu = getParam<Scalar>("SpatialParams.Mu");

    // Tip point: (48, 60, 0) — top-right corner on z=0 symmetry face
    using GlobalPosition = typename MomGG::GlobalCoordinate;
    const GlobalPosition tipPoint{ 48.0, 60.0, 0.0 };

    struct Result { Scalar lambda, nu, uy; bool converged; };
    std::vector<Result> results;
    results.reserve(lambdaValues.size());

    for (const Scalar lambda : lambdaValues)
    {
        const Scalar nu = lambda / (2.0*(lambda + mu));

        momProblem->spatialParams().setLambda(lambda);
        presProblem->spatialParams().setLambda(lambda);

        auto momGV  = std::make_shared<MomGV>(momProblem,   momGG);
        auto presGV = std::make_shared<PresGV>(presProblem, presGG);
        momGV->init(x[CouplingManager::momentumIdx]);
        presGV->init(x[CouplingManager::pressureIdx]);

        couplingManager->init(momProblem, presProblem, x);
        couplingManager->computeColorsForAssembly();

        auto assembler = std::make_shared<Assembler>(
            std::make_tuple(momProblem, presProblem),
            std::make_tuple(momGG, presGG),
            std::make_tuple(momGV, presGV),
            couplingManager);

        auto linearSolver = std::make_shared<LinearSolver>();
        NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

        bool converged = true;
        try { nonLinearSolver.solve(x); }
        catch (const Dumux::NumericalProblem&)
        {
            converged = false;
            x = 0.0;
        }

        if (converged)
        {
            VtkOutputModule<MomGV, typename MDTraits::template SubDomain<0>::SolutionVector>
                vtkWriter(*momGV, x[CouplingManager::momentumIdx],
                          "test_hyperelastic_cooks_membrane_mini3d_lambda"
                          + std::to_string(static_cast<int>(lambda)));
            vtkWriter.addVolumeVariable([](const auto& v){
                return Dune::FieldVector<double,3>{ v.displacement(0), v.displacement(1), v.displacement(2) };
            }, "d");
            vtkWriter.write(1.0);
        }

        Scalar uy = 0.0;
        const auto& momSol = x[CouplingManager::momentumIdx];
        const auto tipElements = intersectingEntities(tipPoint, momGG->boundingBoxTree());
        if (!tipElements.empty())
        {
            const auto element = momGG->element(tipElements[0]);
            const auto elemSol = elementSolution(element, momSol, *momGG);
            uy = evalSolution(element, element.geometry(), *momGG, elemSol, tipPoint)[1];
        }

        results.push_back({ lambda, nu, uy, converged });
    }

    // SPP-1748 Table 2.9 — 3D p-FEM reference values
    const std::vector<Scalar> refs3D = { 11.38113537, 11.15639736, 10.78267246, 10.74641606, 10.74061099 };
    // T2 / T2P0 from Table 2.5 (finest mesh 125856 DOFs, λ=432 only)
    const std::vector<Scalar> refsT2   = { 11.37441, -1, -1, -1, -1 };
    const std::vector<Scalar> refsT2P0 = { 11.38386, -1, -1, -1, -1 };

    std::cout << "\n"
              << "SPP-1748 Cook's membrane — 3D MINI (PQ1Bubble+P1), ψ₁\n"
              << "  T2/T2P0: Galerkin FEM from Table 2.5 (125856 DOFs, λ=432 only).\n"
              << std::string(90, '-') << "\n"
              << std::setw(14) << "λ [MPa]"
              << std::setw(10) << "ν"
              << std::setw(14) << "MINI3D FVM"
              << std::setw(10) << "/%ref"
              << std::setw(12) << "T2 FEM"
              << std::setw(12) << "T2P0 FEM"
              << std::setw(14) << "p-FEM ref"
              << "\n"
              << std::string(90, '-') << "\n";
    for (std::size_t i = 0; i < results.size(); ++i)
    {
        const auto& r = results[i];
        const Scalar ratio = (i < refs3D.size() && refs3D[i] > 0.0) ? r.uy/refs3D[i] : 0.0;
        std::cout << std::setw(14) << std::fixed << std::setprecision(3) << r.lambda
                  << std::setw(10) << std::setprecision(5) << r.nu
                  << std::setw(14) << std::setprecision(3) << r.uy
                  << std::setw(9)  << std::setprecision(1) << 100.0*ratio << "%"
                  << std::setw(12) << (refsT2[i]   > 0 ? std::to_string(refsT2[i]).substr(0,7)   : "  -    ")
                  << std::setw(12) << (refsT2P0[i] > 0 ? std::to_string(refsT2P0[i]).substr(0,7) : "  -    ")
                  << std::setw(14) << std::setprecision(5) << refs3D[i]
                  << (r.converged ? "" : "  [diverged]") << "\n";
    }
    std::cout << std::string(90, '-') << "\n";

    return 0;
}
