// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// SPP-1748 incompressible block mixed formulations in one executable source:
// - P1BP1 (mini/mixed): PQ1Bubble displacement + P1 pressure
// - P2P1 (Taylor-Hood): PQ2 displacement + P1 pressure
//
#include <config.h>
#include <iomanip>
#include <iostream>
#include <string>

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
#include <dumux/discretization/concepts.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dune/common/fmatrix.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtk/fieldtype.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include "properties.hh"

#ifndef MOMENTUM_DISCRETIZATION
#define MOMENTUM_DISCRETIZATION IncompressibleBlockMomentumP1B
#endif

#ifndef PRESSURE_DISCRETIZATION
#define PRESSURE_DISCRETIZATION IncompressibleBlockPressure
#endif

int main(int argc, char** argv)
{
    using namespace Dumux;

    using MomTypeTag = Properties::TTag::MOMENTUM_DISCRETIZATION;
    using PresTypeTag = Properties::TTag::PRESSURE_DISCRETIZATION;
    using MDTraits = MultiDomainTraits<MomTypeTag, PresTypeTag>;

    using Scalar = typename MDTraits::Scalar;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<MomTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using MomGG = GetPropType<MomTypeTag, Properties::GridGeometry>;
    using PresGG = GetPropType<PresTypeTag, Properties::GridGeometry>;
    using MomDiscMethod = typename MomGG::DiscretizationMethod;
    using MomElementDisc = typename MomGG::LocalView;
    const std::string thisRunName = MomDiscMethod::name() + []()
    {
        using namespace Dumux::Experimental::Concepts;
        if constexpr (HybridElementDiscretization<MomElementDisc>)
            return std::string{" FVM-Hybrid (this run)"};
        else if constexpr (FVElementDiscretization<MomElementDisc>)
            return std::string{" FVM (this run)"};
        else if constexpr (FEElementDiscretization<MomElementDisc>)
            return std::string{" FE (this run)"};

        return std::string{" (this run)"};
    }();
    auto momGG = std::make_shared<MomGG>(leafGridView);
    auto presGG = std::make_shared<PresGG>(leafGridView);

    using CouplingManager = HyperelasticVolIsoCouplingManager<MDTraits>;
    auto couplingManager = std::make_shared<CouplingManager>(momGG, presGG);

    using MomProblem = GetPropType<MomTypeTag, Properties::Problem>;
    using PresProblem = GetPropType<PresTypeTag, Properties::Problem>;
    auto momProblem = std::make_shared<MomProblem>(momGG, couplingManager);
    auto presProblem = std::make_shared<PresProblem>(presGG, couplingManager);

    using SolVec = typename MDTraits::SolutionVector;
    SolVec x;
    x[CouplingManager::momentumIdx].resize(momGG->numDofs());
    x[CouplingManager::pressureIdx].resize(presGG->numDofs());
    x = 0.0;

    using MomGV = GetPropType<MomTypeTag, Properties::GridVariables>;
    using PresGV = GetPropType<PresTypeTag, Properties::GridVariables>;
    using Assembler = Dumux::Experimental::MultiDomainAssembler<MDTraits, CouplingManager, DiffMethod::numeric>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits,
                                           LinearAlgebraTraitsFromAssembler<Assembler>>;
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;

    const int nLoadSteps = getParam<int>("LoadStepping.NSteps", 1);
    std::cout << "Starting pseudo-static load stepping with " << nLoadSteps << " steps\n";

    auto momGV = std::make_shared<MomGV>(momProblem, momGG);
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

    VtkOutputModule<MomGV, typename MDTraits::template SubDomain<0>::SolutionVector>
        vtkWriter(*momGV, x[CouplingManager::momentumIdx],
                  getParam<std::string>("Problem.Name"));
    vtkWriter.addVolumeVariable([](const auto& v){
        return Dune::FieldVector<double,3>{ v.displacement(0), v.displacement(1), v.displacement(2) };
    }, "d");

    vtkWriter.write(0.0);

    for (int step = 1; step <= nLoadSteps; ++step)
    {
        const Scalar loadFactor = static_cast<Scalar>(step) / nLoadSteps;
        momProblem->setLoadFactor(loadFactor);

        std::cout << "Load step " << step << "/" << nLoadSteps
                  << " (factor = " << loadFactor << ")\n";

        nonLinearSolver.solve(x);

        momGV->update(x[CouplingManager::momentumIdx]);
        presGV->update(x[CouplingManager::pressureIdx]);
        vtkWriter.write(static_cast<double>(step));
    }

    // Point P = (w/2, l/2, h) = (50,50,50) mm: center of the full load patch (Fig. 3.1).
    // Maximum |u_z| occurs here; SPP-1748 Tables 3.3 & 3.4 reference ≈ 20 mm.
    using GlobalPosition = typename MomGG::GlobalCoordinate;
    const GlobalPosition P{ 50.0 - 1e-6, 50.0 - 1e-6, 50.0 - 1e-6 };
    const auto& momSol = x[CouplingManager::momentumIdx];
    const auto tipElements = intersectingEntities(P, momGG->boundingBoxTree());
    if (!tipElements.empty())
    {
        const auto element = momGG->element(tipElements[0]);
        const auto elemSol = elementSolution(element, momSol, *momGG);
        const auto u = evalSolution(element, element.geometry(), *momGG, elemSol, P);
        const Scalar uz = std::abs(u[2]);
        const int nDofs = momGG->numDofs() + presGG->numDofs();

        // ---------------------------------------------------------------
        // SPP-1748 Benchmark 3 (Kollmannsberger et al. 2020, Tables 3.3 & 3.4)
        // |u_z(P)| [mm] at P=(w/2, l/2, h) = (50,50,50) mm (center of full load patch).
        // Full model: w=100mm, l=100mm, h=50mm; quarter [0,50]³ simulated.
        // Load: q=3 MPa, a=b=25mm (Table 3.1); λ=499.92568 MPa, µ=1.61148 MPa (Table 3.2).
        // ---------------------------------------------------------------
        // Table 3.3: h-refinement (n_dofs for H1; same node count for other elements)
        // n_dofs |   H1  | H1P0  | O2P1  | H1EI9 | TSCG
        //    260 |  7.78 | 19.87 | 19.80 | 19.09 | 20.14
        //   1800 | 13.17 | 20.02 | 19.93 | 19.98 | 20.10
        //  13328 | 17.54 | 20.01 | 19.97 | 20.01 | 20.03
        // 102432 | 19.52 | 20.00 | 19.98 | 20.00 | 20.01
        // ---------------------------------------------------------------
        // Table 3.4: p-extension on 4×4×4 mesh (isotropic trunk space)
        // p=1 (375 dofs): 7.656 mm, p=2 (1275): 19.637, p=3 (2175): 19.998
        // Converged reference: |u_z(P)| ≈ 20.00 mm
        // ---------------------------------------------------------------

        constexpr int nRef = 4;
        constexpr int refDofs[nRef] = { 260, 1800, 13328, 102432 };
        constexpr double refH1[nRef] = { 7.78, 13.17, 17.54, 19.52 };
        constexpr double refH1P0[nRef] = { 19.87, 20.02, 20.01, 20.00 };
        constexpr double refO2P1[nRef] = { 19.80, 19.93, 19.97, 19.98 };
        constexpr double refH1EI9[nRef] = { 19.09, 19.98, 20.01, 20.00 };
        constexpr double refTSCG[nRef] = { 20.14, 20.10, 20.03, 20.01 };

        int row = 0;
        for (int i = 1; i < nRef; ++i)
            if (std::abs(nDofs - refDofs[i]) < std::abs(nDofs - refDofs[row]))
                row = i;

        std::cout << "\n";
        std::cout << " SPP-1748 Benchmark 3 — incompressible block |u_z(P)| comparison\n";
        std::cout << " Reference: Schröder et al. (2020), Tables 3.3 & 3.4\n";
        std::cout << " Converged reference: |u_z(P)| ≈ 20.00 mm\n\n";
        std::cout << std::left
              << " " << std::setw(28) << "Formulation"
                  << std::right << std::setw(10) << "n_dofs"
                  << std::setw(14) << "|u_z(P)| mm" << "\n";
        std::cout << " " << std::string(52, '-') << "\n";

        auto row2 = [&](const std::string& name, int nd, double val) {
            std::cout << " " << std::left << std::setw(28) << name
                      << std::right << std::setw(10) << nd
                      << std::setw(14) << std::fixed << std::setprecision(2) << val << "\n";
        };

        row2("H1 (locked)", refDofs[row], refH1[row]);
        row2("H1P0 (mixed hex)", refDofs[row], refH1P0[row]);
        row2("O2P1 (mixed tet)", refDofs[row], refO2P1[row]);
        row2("H1EI9 (enh. strain)", refDofs[row], refH1EI9[row]);
        row2("TSCG (stab. enh.)", refDofs[row], refTSCG[row]);
        std::cout << " " << std::string(52, '-') << "\n";
        row2(thisRunName, nDofs, uz);
        std::cout << "\n";
    }

    return 0;
}