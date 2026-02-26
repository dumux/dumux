// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid Navier-Stokes model (Donea 2003, \cite Donea2003).
 */

#include <config.h>

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/stokes_solver.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/momentum/velocityoutput.hh>
#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors.hh>
#include <test/freeflow/navierstokes/errors_cvfe.hh>

#include "properties.hh"

#if HAVE_UMFPACK
#include <umfpack.h>
#endif

#ifndef USE_STOKES_SOLVER
#define USE_STOKES_SOLVER 0
#endif

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

        for (const auto& intersection : intersections(momentumGridGeometry->gridView(), element))
        {
            if(intersection.boundary() == false)
                continue;

            const auto bcTypes = momentumProblem->boundaryTypes(fvGeometry, intersection);
            if (bcTypes.hasDirichlet())
            {
                for (auto&& localDof : localDofs(fvGeometry, intersection))
                {
                        for (int i = 0; i < bcTypes.size(); ++i)
                            if (bcTypes.isDirichlet(i))
                                dirichletDofs[momentumIdx][localDof.index()][i] = 1.0;
                }
            }
        }
    }

    return dirichletDofs;
}

namespace Dumux {

template<class Error>
void writeError_(std::ofstream& logFile, const Error& error)
{
    for (const auto& e : error)
        logFile << ", " << e;
}

template<class MomentumProblem, class MassProblem,
         class MomentumGridVariables, class MassGridVariables,
         class SolutionVector, class MomentumIdx, class MassIdx>
void printErrors(std::shared_ptr<MomentumProblem> momentumProblem,
                 std::shared_ptr<MassProblem> massProblem,
                 const MomentumGridVariables& momentumGridVariables,
                 const MassGridVariables& massGridVariables,
                 const SolutionVector& x,
                 const MomentumIdx momentumIdx,
                 const MassIdx massIdx)
{
    using MomentumGridGeometry = std::decay_t<decltype(std::declval<MomentumProblem>().gridGeometry())>;
    using MassGridGeometry = std::decay_t<decltype(std::declval<MassProblem>().gridGeometry())>;
    static constexpr int dim = MomentumGridGeometry::GridView::dimension;

    // print discrete L2 and Linfity errors
    const bool printErrors = getParam<bool>("Problem.PrintErrors", false);
    const bool printConvergenceTestFile = getParam<bool>("Problem.PrintConvergenceTestFile", false);

    if (!printErrors && !printConvergenceTestFile)
        return;

    if constexpr (DiscretizationMethods::isCVFE<typename MomentumGridGeometry::DiscretizationMethod>
                  && DiscretizationMethods::isCVFE<typename MassGridGeometry::DiscretizationMethod>)
    {
        // first print momentum
        {
        const auto [totalVolume, errors] = calculateL2AndH1Errors(*momentumProblem, momentumGridVariables, x[momentumIdx]);

        if (printErrors)
            std::cout << "[Errors] velocity L2 = " << errors[0] << " H1 = " << errors[1] << std::endl;

        if (printConvergenceTestFile)
        {
            std::ofstream logFile(momentumProblem->name() + "_errors_velocity.csv", std::ios::app);
            auto numDofs = momentumProblem->gridGeometry().numDofs();
            logFile << numDofs << ", ";
            logFile << std::pow(totalVolume / numDofs, 1.0/dim);
            writeError_(logFile, errors);
            logFile << "\n";
        }
        }

        // now print mass
        {
        const auto [totalVolume, errors] = calculateL2AndH1Errors(*massProblem, massGridVariables, x[massIdx]);

        if (printErrors)
            std::cout << "[Errors] pressure L2 = " << errors[0] << " H1 = " << errors[1] << std::endl;

        if (printConvergenceTestFile)
        {
            std::ofstream logFile(massProblem->name() + "_errors_pressure.csv", std::ios::app);
            auto numDofs = massProblem->gridGeometry().numDofs();
            logFile << numDofs << ", ";
            logFile << std::pow(totalVolume / numDofs, 1.0/dim);
            writeError_(logFile, errors);
            logFile << "\n";
        }
        }
    }
    else
    {
        NavierStokesTest::Errors errors(momentumProblem, massProblem, x);
        NavierStokesTest::ErrorCSVWriter(
            momentumProblem, massProblem, std::to_string(x[massIdx].size())
        ).printErrors(errors);

        if (printConvergenceTestFile)
            NavierStokesTest::convergenceTestAppendErrors(momentumProblem, massProblem, errors);
    }
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tags for this problem
    using MomentumTypeTag = Properties::TTag::DoneaTestMomentum;
    using MassTypeTag = Properties::TTag::TYPETAG_MASS;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    using GridManager = Dumux::GridManager<GetPropType<MassTypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    // start timer
    Dune::Timer timer;

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometries
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;

    // share the basic grid geometry
    using BasicGridGeometry = typename MassGridGeometry::BasicGridGeometry;
    auto basicGridGeometry = std::make_shared<BasicGridGeometry>(leafGridView);

    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(basicGridGeometry);
    auto massGridGeometry = std::make_shared<MassGridGeometry>(basicGridGeometry);

    // mass-momentum coupling manager
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);

    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);

    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // initialize the coupling stencils
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    // initializing the gridvariables requires the coupling manager to be set up
    momentumGridVariables->init(x[momentumIdx]);
    massGridVariables->init(x[massIdx]);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());

    NavierStokesTest::AnalyticalSolutionVectors analyticalSolVectors(momentumProblem, massProblem);
    vtkWriter.addField(analyticalSolVectors.analyticalPressureSolution(), "pressureExact");
    vtkWriter.addField(analyticalSolVectors.analyticalVelocitySolution(), "velocityExact");
    //vtkWriter.addFaceField(analyticalSolVectors.analyticalVelocitySolutionOnFace(), "faceVelocityExact");
    vtkWriter.write(0.0);

    // use the multidomain FV assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);

    // the linear solver
#if USE_STOKES_SOLVER
    using Matrix = typename Assembler::JacobianMatrix;
    using Vector = typename Assembler::ResidualType;
    using LinearSolver = StokesSolver<Matrix, Vector, MomentumGridGeometry, MassGridGeometry>;
    auto dDofs = dirichletDofs<Vector>(momentumGridGeometry, massGridGeometry, momentumProblem, momentumIdx, massIdx);
    auto linearSolver = std::make_shared<LinearSolver>(momentumGridGeometry, massGridGeometry, dDofs);
#else
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
#endif

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // linearize & solve
    nonLinearSolver.solve(x);

    Dumux::printErrors(momentumProblem, massProblem, *momentumGridVariables, *massGridVariables, x, momentumIdx, massIdx);

    // write vtk output
    vtkWriter.write(1.0);

    timer.stop();

    const auto& comm = Dune::MPIHelper::getCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
