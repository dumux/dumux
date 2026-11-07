// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup AssemblyTests
 * \brief Test of the multi-stage assembler on a grid that is periodic in both directions
 *
 * The heat equation on \f$ [0, 2\pi]^2 \f$ has the solution
 * \f$ u = 1 + \mathrm{e}^{-5t} \sin(x + 0.3) \sin(2y + 0.7) \f$, which satisfies no boundary
 * condition on the domain boundary other than periodicity. The test checks that
 * - implicit Euler with the multi-stage assembler reproduces the stationary-scheme assembler,
 * - the residual assembled without the Jacobian equals the one assembled with it,
 * - the degrees of freedom on both sides of the periodic boundaries keep the same value,
 * - the heat in the control volumes is conserved to round-off.
 * It writes the L2 error of the third-order DIRK solution for the convergence test driver.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <numbers>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/multistageassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/discretization/fem/interpolationpointdata.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"

namespace Dumux {

using TypeTag = Properties::TTag::TYPETAG;
using Scalar = GetPropType<TypeTag, Properties::Scalar>;
using HostGrid = GetPropType<TypeTag, Properties::Grid>;
using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
using Problem = GetPropType<TypeTag, Properties::Problem>;
using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

Scalar exactSolution(const GlobalPosition& p, Scalar t)
{ return 1.0 + std::exp(-5.0*t)*std::sin(p[0] + 0.3)*std::sin(2.0*p[1] + 0.7); }

SolutionVector initialSolution(const GridGeometry& gridGeometry)
{
    SolutionVector x(gridGeometry.numDofs());
    auto elemDisc = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        elemDisc.bind(element);
        for (const auto& localDof : localDofs(elemDisc))
            x[localDof.dofIndex()] = exactSolution(ipData(elemDisc, localDof).global(), 0.0);
    }
    return x;
}

//! The heat in the control volumes, the quantity balanced by the control-volume equations
Scalar controlVolumeHeat(const GridGeometry& gridGeometry, const SolutionVector& x)
{
    Scalar heat = 0.0;
    auto elemDisc = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        elemDisc.bind(element);
        for (const auto& scv : scvs(elemDisc))
            heat += scv.volume()*x[scv.dofIndex()][0];
    }
    return heat;
}

Scalar l2Error(const GridGeometry& gridGeometry, const SolutionVector& x, Scalar t)
{
    Scalar error = 0.0;
    auto elemDisc = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        elemDisc.bind(element);
        const auto geometry = element.geometry();
        const auto& localBasis = elemDisc.feLocalBasis();
        for (const auto& qp : Dune::QuadratureRules<Scalar, 2>::rule(geometry.type(), 8))
        {
            const auto global = geometry.global(qp.position());
            const FEInterpolationPointData<GlobalPosition, std::decay_t<decltype(localBasis)>>
                shapeData(geometry, qp.position(), global, localBasis);
            Scalar u = 0.0;
            for (const auto& localDof : localDofs(elemDisc))
                u += shapeData.shapeValues()[localDof.index()][0]*x[localDof.dofIndex()][0];
            error += (u - exactSolution(global, t))*(u - exactSolution(global, t))
                     *qp.weight()*geometry.integrationElement(qp.position());
        }
    }
    return std::sqrt(error);
}

Scalar maxPeriodicMismatch(const GridGeometry& gridGeometry, const SolutionVector& x)
{
    Scalar mismatch = 0.0;
    for (const auto& [dofIdx, mappedDofs] : gridGeometry.periodicDofMap())
        for (const auto mappedDofIdx : mappedDofs)
            mismatch = std::max(mismatch, std::abs(x[dofIdx][0] - x[mappedDofIdx][0]));
    return mismatch;
}

std::shared_ptr<GridGeometry> makeGridGeometry(GridManager<HostGrid>& gridManager)
{
    gridManager.init();
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    if (!gridGeometry->isPeriodic())
        DUNE_THROW(Dune::InvalidStateException, "The grid geometry is not periodic");
    return gridGeometry;
}

using MultiStageAssembler = Experimental::MultiStageAssembler<TypeTag, DiffMethod::numeric>;
using MultiStageLinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<MultiStageAssembler>>;

//! Time integration with the multi-stage assembler and a fixed time step size
SolutionVector solveMultiStage(std::shared_ptr<const GridGeometry> gridGeometry,
                               std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> method,
                               Scalar dt, int numSteps)
{
    auto problem = std::make_shared<Problem>(gridGeometry);
    auto x = initialSolution(*gridGeometry);
    auto xOld = x;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    auto assembler = std::make_shared<MultiStageAssembler>(problem, gridGeometry, gridVariables, method, xOld);
    using Solver = NewtonSolver<MultiStageAssembler, MultiStageLinearSolver>;
    auto solver = std::make_shared<Solver>(assembler, std::make_shared<MultiStageLinearSolver>());
    Experimental::MultiStageTimeStepper<Solver> timeStepper(solver, method);

    for (int stepIdx = 0; stepIdx < numSteps; ++stepIdx)
    {
        timeStepper.step(x, stepIdx*dt, dt);
        xOld = x;
        gridVariables->advanceTimeStep();
    }
    return x;
}

//! Implicit Euler with the assembler of single-stage schemes
SolutionVector solveImplicitEuler(std::shared_ptr<const GridGeometry> gridGeometry, Scalar dt, int numSteps)
{
    auto problem = std::make_shared<Problem>(gridGeometry);
    auto x = initialSolution(*gridGeometry);
    auto xOld = x;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, numSteps*dt, false);
    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);
    using Solver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    NewtonSolver<Assembler, Solver> nonLinearSolver(assembler, std::make_shared<Solver>());

    for (int stepIdx = 0; stepIdx < numSteps; ++stepIdx)
    {
        nonLinearSolver.solve(x);
        xOld = x;
        gridVariables->advanceTimeStep();
    }
    return x;
}

//! Largest relative difference between the residual assembled without and with the Jacobian over the stages of a step
Scalar residualAssemblyDifference(std::shared_ptr<const GridGeometry> gridGeometry,
                                  std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> method, Scalar dt)
{
    auto problem = std::make_shared<Problem>(gridGeometry);
    auto x = initialSolution(*gridGeometry);
    auto xOld = x;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    auto assembler = std::make_shared<MultiStageAssembler>(problem, gridGeometry, gridVariables, method, xOld);
    NewtonSolver<MultiStageAssembler, MultiStageLinearSolver> nonLinearSolver(assembler, std::make_shared<MultiStageLinearSolver>());

    Scalar maxRelDiff = 0.0;
    for (std::size_t stageIdx = 1; stageIdx <= method->numStages(); ++stageIdx)
    {
        assembler->prepareStage(x, std::make_shared<Experimental::MultiStageParams<Scalar>>(*method, stageIdx, 0.0, dt));
        assembler->assembleJacobianAndResidual(x);
        auto difference = assembler->residual();
        assembler->assembleResidual(x);
        difference -= assembler->residual();
        maxRelDiff = std::max(maxRelDiff, difference.two_norm()/assembler->residual().two_norm());
        nonLinearSolver.solve(x);
    }
    return maxRelDiff;
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);
    Parameters::init(argc, argv);

    const auto dirk3 = std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<Scalar>>();
    const auto implicitEuler = std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();
    const Scalar tEnd = 0.2;
    const int numSteps = 20;
    const Scalar dt = tEnd/numSteps;

    bool passed = true;
    const auto check = [&](bool condition, const auto& message)
    {
        std::cout << (condition ? "passed: " : "FAILED: ") << message << std::endl;
        passed = passed && condition;
    };

    GridManager<HostGrid> gridManager;
    const auto gridGeometry = makeGridGeometry(gridManager);

    if (getParam<bool>("Test.CompareAssemblers"))
    {
        const auto multiStage = solveMultiStage(gridGeometry, implicitEuler, dt, numSteps);
        auto difference = solveImplicitEuler(gridGeometry, dt, numSteps);
        difference -= multiStage;
        const auto relDiff = difference.infinity_norm()/multiStage.infinity_norm();
        check(relDiff < 1e-12, Fmt::format("implicit Euler with both assemblers, relative difference {:.3e}", relDiff));

        const auto residualDiff = residualAssemblyDifference(gridGeometry, dirk3, dt);
        check(residualDiff < 1e-12, Fmt::format("residual with and without Jacobian, relative difference {:.3e}", residualDiff));
    }

    const auto x = solveMultiStage(gridGeometry, dirk3, 2.0*dt, numSteps/2);
    const auto heat0 = controlVolumeHeat(*gridGeometry, initialSolution(*gridGeometry));
    const auto heatChange = std::abs(controlVolumeHeat(*gridGeometry, x) - heat0)/heat0;
    check(heatChange < 1e-12, Fmt::format("relative change of the heat {:.3e}", heatChange));
    const auto mismatch = maxPeriodicMismatch(*gridGeometry, x);
    check(mismatch < 1e-12, Fmt::format("largest difference of periodic partners {:.3e}", mismatch));

    const auto& gridView = gridGeometry->gridView();
    const auto h = (gridGeometry->bBoxMax()[0] - gridGeometry->bBoxMin()[0])/std::sqrt(static_cast<Scalar>(gridView.size(0)));
    const auto error = l2Error(*gridGeometry, x, tEnd);
    std::cout << Fmt::format("h = {:.6e}, L2 error = {:.6e}", h, error) << std::endl;
    std::ofstream(getParam<std::string>("Problem.Name") + "-errors.csv") << "h,L2(u)\n" << h << "," << error << "\n";

    return passed ? 0 : 1;
}
