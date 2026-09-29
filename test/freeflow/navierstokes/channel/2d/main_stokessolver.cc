// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Transient channel flow on the staggered grid solved with the Stokes solver and its options for
 *        transient problems: a pressure operator for small time steps and a kept matrix.
 */

#include <config.h>

#ifndef NONISOTHERMAL
#define NONISOTHERMAL 0
#endif

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/linear/stokes_solver.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/momentum/velocityoutput.hh>

#include "properties.hh"
#include "transientpressureoperator.hh"

template<class Vector, class MomGG, class MassGG, class MomP, class MomIdx, class MassIdx>
auto dirichletDofs(std::shared_ptr<MomGG> momentumGridGeometry,
                   std::shared_ptr<MassGG> massGridGeometry,
                   std::shared_ptr<MomP> momentumProblem,
                   MomIdx momentumIdx, MassIdx massIdx)
{
    Vector dirichletDofs;
    dirichletDofs[momentumIdx].resize(momentumGridGeometry->numDofs());
    dirichletDofs[massIdx].resize(massGridGeometry->numDofs());
    dirichletDofs = 0.0;

    auto fvGeometry = localView(*momentumGridGeometry);
    for (const auto& element : elements(momentumGridGeometry->gridView()))
    {
        fvGeometry.bind(element);
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary() || !scvf.isFrontal())
                continue;
            const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
            if (scv.boundary())
            {
                const auto bcTypes = momentumProblem->boundaryTypes(element, scvf);
                if (bcTypes.isDirichlet(scv.dofAxis()))
                    dirichletDofs[momentumIdx][scv.dofIndex()][0] = 1.0;
            }
        }
    }

    return dirichletDofs;
}

/*!
 * \brief Solve with a kept matrix and compare with a solve that is given the matrix
 *
 * The right-hand side carries values in the Dirichlet rows, which the symmetrization of the constraints moves
 * into the other rows; the kept matrix has to do the same. The second right-hand side is not a multiple of
 * the first and is solved with the preconditioner built for the first.
 */
template<class LinearSolver, class Matrix, class Vector>
void checkKeptMatrix(LinearSolver& linearSolver, const Matrix& matrix, const Vector& residual, const Vector& dirichletDofs)
{
    using namespace Dune::Indices;
    linearSolver.setMatrix(matrix);

    auto rhs = residual;
    for (std::size_t i = 0; i < rhs[_0].size(); ++i)
        if (dirichletDofs[_0][i][0] > 0.5)
            rhs[_0][i] = 1.0;

    for (int k = 0; k < 2; ++k)
    {
        auto givenMatrix = rhs; givenMatrix = 0.0;
        auto keptMatrix = rhs; keptMatrix = 0.0;
        linearSolver.solve(matrix, givenMatrix, rhs);
        linearSolver.solve(keptMatrix, rhs);

        auto difference = givenMatrix;
        difference -= keptMatrix;
        const auto relativeDifference = linearSolver.norm(difference)/linearSolver.norm(givenMatrix);
        if (Dune::MPIHelper::instance().rank() == 0)
            std::cout << "Kept matrix, right-hand side " << k << ": relative difference " << relativeDifference << std::endl;
        if (!(relativeDifference < 1e-6))
            DUNE_THROW(Dune::Exception, "Solving with the kept matrix gives a different solution");

        rhs[_0] *= 2.0;
        rhs[_1] *= -0.5;
    }
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    using MomentumTypeTag = Properties::TTag::ChannelTestMomentum;
    using MassTypeTag = Properties::TTag::ChannelTestMass;

    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    GridManager<GetPropType<MassTypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    using Scalar = GetPropType<MassTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);
    auto xOld = x;

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    massProblem->setTime(timeLoop->time());
    momentumProblem->setTime(timeLoop->time());

    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);
    momentumGridVariables->init(x[momentumIdx]);
    massGridVariables->init(x[massIdx]);

    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.write(0.0);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager, timeLoop, xOld);

    using Matrix = typename Assembler::JacobianMatrix;
    using Vector = typename Assembler::ResidualType;
    using LinearSolver = StokesSolver<Matrix, Vector, MomentumGridGeometry, MassGridGeometry>;
    const auto dDofs = dirichletDofs<Vector>(momentumGridGeometry, massGridGeometry, momentumProblem, momentumIdx, massIdx);
    auto linearSolver = std::make_shared<LinearSolver>(momentumGridGeometry, massGridGeometry, dDofs);

    // for small time steps the pressure Schur complement is dominated by the storage term of the velocity block
    const bool useTransientPressureOperator = getParam<std::string>("LinearSolver.PressureOperator", "MassMatrix") == "Transient";
    TransientPressureOperator<MassGridGeometry> pressureOperator(massGridGeometry,
                                                                 getParam<Scalar>("Component.LiquidDensity"),
                                                                 getParam<Scalar>("Component.LiquidDynamicViscosity"));
    if (useTransientPressureOperator)
    {
        linearSolver->setPressureMatrix(pressureOperator.matrix());
        linearSolver->setPressureDiagonal(pressureOperator.viscousDiagonal());
    }

    if (getParam<bool>("LinearSolver.CheckKeptMatrix", false))
    {
        pressureOperator.update(dt);
        assembler->assembleJacobianAndResidual(x);
        checkKeptMatrix(*linearSolver, assembler->jacobian(), assembler->residual(), dDofs);
    }

    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    timeLoop->start(); do
    {
        pressureOperator.update(timeLoop->timeStepSize());
        nonLinearSolver.solve(x, *timeLoop);

        xOld = x;
        momentumGridVariables->advanceTimeStep();
        massGridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();
        vtkWriter.write(timeLoop->time());
        timeLoop->reportTimeStep();

        massProblem->setTime(timeLoop->time());
        momentumProblem->setTime(timeLoop->time());
    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
