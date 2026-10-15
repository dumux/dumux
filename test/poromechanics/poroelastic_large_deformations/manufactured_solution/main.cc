// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Plane-strain manufactured solution for the large-deformation poroelastic model.
 *
 * The exact solution compresses the square locally to J = 0.55 and stretches it to J = 1.29.
 * Its storage term is linear in time, so implicit Euler adds no time discretization error
 * and the L2 errors at the end time measure the spatial discretization error alone.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>

#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/multistagemultidomainassembler.hh>

#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include "properties.hh"

namespace Dumux {

//! L2 norm of the difference between a discrete solution and an exact field
template<class GridGeometry, class Solution, class ExactSolution>
double l2Error(const GridGeometry& gridGeometry, const Solution& sol, ExactSolution&& exact)
{
    double error = 0.0;
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto geometry = element.geometry();
        const auto elemSol = elementSolution(element, sol, gridGeometry);
        for (const auto& qp : Dune::QuadratureRules<double, 2>::rule(geometry.type(), 8))
        {
            const auto numeric = evalSolutionAtLocalPos(element, geometry, gridGeometry, elemSol, qp.position());
            const auto analytic = exact(geometry.global(qp.position()));
            double diff2 = 0.0;
            for (int i = 0; i < numeric.size(); ++i)
                diff2 += (numeric[i] - analytic[i])*(numeric[i] - analytic[i]);
            error += diff2*qp.weight()*geometry.integrationElement(qp.position());
        }
    }
    return std::sqrt(error);
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    using MomentumTypeTag = Properties::TTag::PoroElasticLargeDefMMSMomentum;
    using SolidPressureTypeTag = Properties::TTag::PoroElasticLargeDefMMSSolidPressure;
    using FluidPressureTypeTag = Properties::TTag::PoroElasticLargeDefMMSFluidPressure;
    using Traits = Properties::PoroElasticLargeDefMMSTraits;
    using Scalar = typename Traits::Scalar;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<MomentumTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    auto momentumGG = std::make_shared<GetPropType<MomentumTypeTag, Properties::GridGeometry>>(leafGridView);
    auto solidPressureGG = std::make_shared<GetPropType<SolidPressureTypeTag, Properties::GridGeometry>>(leafGridView);
    auto fluidPressureGG = std::make_shared<GetPropType<FluidPressureTypeTag, Properties::GridGeometry>>(leafGridView);

    using CouplingManager = PoroElasticLargeDefCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>(momentumGG, solidPressureGG, fluidPressureGG);

    auto momentumProblem = std::make_shared<GetPropType<MomentumTypeTag, Properties::Problem>>(momentumGG, couplingManager);
    auto solidPressureProblem = std::make_shared<GetPropType<SolidPressureTypeTag, Properties::Problem>>(solidPressureGG, couplingManager);
    auto fluidPressureProblem = std::make_shared<GetPropType<FluidPressureTypeTag, Properties::Problem>>(fluidPressureGG, couplingManager);

    constexpr auto momentumIdx = CouplingManager::momentumIdx;
    constexpr auto solidPressureIdx = CouplingManager::solidPressureIdx;
    constexpr auto fluidPressureIdx = CouplingManager::fluidPressureIdx;

    typename Traits::SolutionVector x;
    x[momentumIdx].resize(momentumGG->numDofs());
    x[solidPressureIdx].resize(solidPressureGG->numDofs());
    x[fluidPressureIdx].resize(fluidPressureGG->numDofs());
    x = 0.0;
    auto xOld = x;

    auto momentumGV = std::make_shared<GetPropType<MomentumTypeTag, Properties::GridVariables>>(momentumProblem, momentumGG);
    auto solidPressureGV = std::make_shared<GetPropType<SolidPressureTypeTag, Properties::GridVariables>>(solidPressureProblem, solidPressureGG);
    auto fluidPressureGV = std::make_shared<GetPropType<FluidPressureTypeTag, Properties::GridVariables>>(fluidPressureProblem, fluidPressureGG);
    momentumGV->init(x[momentumIdx]);
    solidPressureGV->init(x[solidPressureIdx]);
    fluidPressureGV->init(x[fluidPressureIdx]);

    couplingManager->init(momentumProblem, solidPressureProblem, fluidPressureProblem, x);

    auto method = std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();
    using Assembler = Experimental::MultiStageMultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(momentumProblem, solidPressureProblem, fluidPressureProblem),
        std::make_tuple(momentumGG, solidPressureGG, fluidPressureGG),
        std::make_tuple(momentumGV, solidPressureGV, fluidPressureGV),
        couplingManager, method, xOld
    );
    assembler->setLinearSystem();

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, std::make_shared<LinearSolver>(), couplingManager);
    Experimental::MultiStageTimeStepper<NewtonSolver> timeStepper(newtonSolver, method);

    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.Dt");
    TimeLoop<Scalar> timeLoop(0.0, dt, tEnd);
    timeLoop.start();
    do
    {
        xOld = x;
        assembler->setPreviousSolution(xOld);
        timeStepper.step(x, timeLoop.time(), timeLoop.timeStepSize());
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();
    } while (!timeLoop.finished());

    const Scalar t = timeLoop.time();
    const Scalar K = momentumProblem->spatialParams().bulkModulus();
    const auto displacementError = l2Error(*momentumGG, x[momentumIdx], [&](const auto& pos) {
        return PoroElasticLargeDefMMS::displacement(pos[0], pos[1], t);
    });
    const auto solidPressureError = l2Error(*solidPressureGG, x[solidPressureIdx], [&](const auto& pos) {
        return std::array<Scalar, 1>{PoroElasticLargeDefMMS::solidPressure(pos[0], pos[1], t, K)};
    });
    const auto fluidPressureError = l2Error(*fluidPressureGG, x[fluidPressureIdx], [&](const auto& pos) {
        return std::array<Scalar, 1>{PoroElasticLargeDefMMS::fluidPressure(pos[0], pos[1], t)};
    });

    Scalar hMax = 0.0;
    for (const auto& element : elements(leafGridView))
        for (int i = 0; i < element.geometry().corners(); ++i)
            for (int j = i+1; j < element.geometry().corners(); ++j)
                hMax = std::max(hMax, (element.geometry().corner(i) - element.geometry().corner(j)).two_norm());

    std::cout << "Max element diameter: " << hMax << "\n"
              << "L2 error displacement: " << displacementError << "\n"
              << "L2 error solid pressure: " << solidPressureError << "\n"
              << "L2 error fluid pressure: " << fluidPressureError << std::endl;

    return 0;
}
