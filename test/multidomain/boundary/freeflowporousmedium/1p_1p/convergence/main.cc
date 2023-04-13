// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief A convergence test for the coupled FreeFlow/Darcy problem (1p).
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/partial.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/format.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/multidomain/staggeredtraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>

#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors.hh>

#include "testcase.hh"
#include "properties.hh"

/*!
* \brief Creates analytical solution vector for output
* Returns a tuple of the analytical solution for the pressure and velocity at cell centers
* \param problem the problem for which to evaluate the analytical solution
*/
template<class Scalar, class Problem>
auto createDarcyAnalyticalSolution(const Problem& problem)
{
    const auto& gridGeometry = problem.gridGeometry();
    using GridView = typename std::decay_t<decltype(gridGeometry)>::GridView;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    std::vector<Scalar> analyticalPressure;
    std::vector<VelocityVector> analyticalVelocity;

    analyticalPressure.resize(gridGeometry.numDofs());
    analyticalVelocity.resize(gridGeometry.numDofs());

    for (const auto& element : elements(gridGeometry.gridView()))
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto ccDofIdx = scv.dofIndex();
            const auto ccDofPosition = scv.dofPosition();
            const auto analyticalSolutionAtCc = problem.fullAnalyticalSolution(ccDofPosition);
            analyticalPressure[ccDofIdx] = analyticalSolutionAtCc[dim];

            for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
                analyticalVelocity[ccDofIdx][dirIdx] = analyticalSolutionAtCc[dirIdx];
        }
    }

    return std::make_tuple(analyticalPressure, analyticalVelocity);
}

template<class MassP, class MomP, class SolutionVector>
void printFreeFlowL2Error(std::shared_ptr<MomP> momentumProblem,
                          std::shared_ptr<MassP> massProblem,
                          const SolutionVector& sol)
{
    using namespace Dumux;

    NavierStokesTest::Errors errors(momentumProblem, massProblem, sol);
    const int numCellCenterDofs = massProblem->gridGeometry().numDofs();
    const int numFaceDofs = momentumProblem->gridGeometry().numDofs();
    const auto absL2 = errors.l2Absolute();
    const auto relL2 = errors.l2Relative();

    std::cout << Fmt::format("** L2 error (abs/rel) for {} cc dofs and {} face dofs (total: {}): ",
                             numCellCenterDofs, numFaceDofs, numCellCenterDofs + numFaceDofs)
              << Fmt::format("L2(p) = {:.8e} / {:.8e}", absL2[0], relL2[0])
              << Fmt::format(", L2(vx) = {:.8e} / {:.8e}", absL2[1], relL2[1])
              << Fmt::format(", L2(vy) = {:.8e} / {:.8e}", absL2[2], relL2[2])
              << std::endl;

    // write the norm into a log file
    std::ofstream logFile(massProblem->name() + ".log", std::ios::app);
    logFile << "[ConvergenceTest] L2(p) = " << absL2[0]
            << " L2(vx) = " << absL2[1]
            << " L2(vy) = " << absL2[2]
            << std::endl;
}

template<class Problem, class SolutionVector>
void printDarcyL2Error(std::shared_ptr<Problem> problem, const SolutionVector& x)
{
    using namespace Dumux;
    using Scalar = double;

    Scalar l2error = 0.0;
    for (const auto& element : elements(problem->gridGeometry().gridView()))
    {
        auto fvGeometry = localView(problem->gridGeometry());
        fvGeometry.bindElement(element);

        for (auto&& scv : scvs(fvGeometry))
        {
            const auto dofIdx = scv.dofIndex();
            const Scalar delta = x[dofIdx] - problem->fullAnalyticalSolution(scv.center())[2/*pressureIdx*/];
            l2error += scv.volume()*(delta*delta);
        }
    }
    using std::sqrt;
    l2error = sqrt(l2error);

    const auto numDofs = problem->gridGeometry().numDofs();
    std::cout << Fmt::format("** L2 error (abs) for {} cc dofs", numDofs)
              << Fmt::format("L2 error = {:.8e}", l2error)
              << std::endl;

    // write the norm into a log file
    std::ofstream logFile(problem->name() + ".log", std::ios::app);
    logFile << "[ConvergenceTest] L2(p) = " << l2error << std::endl;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using FreeFlowMomentumTypeTag = Properties::TTag::FreeFlowOnePMomentum;
    using FreeFlowMassTypeTag = Properties::TTag::FreeFlowOnePMass;
    using DarcyTypeTag = Properties::TTag::DarcyOneP;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using DarcyGridManager = Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy"); // pass parameter group

    using FreeFlowGridManager = Dumux::GridManager<GetPropType<FreeFlowMomentumTypeTag, Properties::Grid>>;
    FreeFlowGridManager freeFlowGridManager;
    freeFlowGridManager.init("FreeFlow"); // pass parameter group

    // we compute on the leaf grid view
    const auto& darcyGridView = darcyGridManager.grid().leafGridView();
    const auto& freeFlowGridView = freeFlowGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FreeFlowMomentumGridGeometry = GetPropType<FreeFlowMomentumTypeTag, Properties::GridGeometry>;
    auto freeFlowMomentumGridGeometry = std::make_shared<FreeFlowMomentumGridGeometry>(freeFlowGridView);
    using FreeFlowMassGridGeometry = GetPropType<FreeFlowMassTypeTag, Properties::GridGeometry>;
    auto freeFlowMassGridGeometry = std::make_shared<FreeFlowMassGridGeometry>(freeFlowGridView);
    using DarcyFVGridGeometry = GetPropType<DarcyTypeTag, Properties::GridGeometry>;
    auto darcyGridGeometry = std::make_shared<DarcyFVGridGeometry>(darcyGridView);

    using Traits = MultiDomainTraits<FreeFlowMomentumTypeTag, FreeFlowMassTypeTag, DarcyTypeTag>;
    using CouplingManager = FreeFlowPorousMediumCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the indices
    constexpr auto freeFlowMomentumIndex = CouplingManager::freeFlowMomentumIndex;
    constexpr auto freeFlowMassIndex = CouplingManager::freeFlowMassIndex;
    constexpr auto porousMediumIndex = CouplingManager::porousMediumIndex;

    // the problem (initial and boundary conditions)
    const auto testCaseName = getParam<std::string>("Problem.TestCase");
    const auto testCase = [&]()
    {
        if (testCaseName == "ShiueExampleOne")
            return DarcyStokesTestCase::ShiueExampleOne;
        else if (testCaseName == "ShiueExampleTwo")
            return DarcyStokesTestCase::ShiueExampleTwo;
        else if (testCaseName == "Rybak")
            return DarcyStokesTestCase::Rybak;
        else if (testCaseName == "Schneider")
            return DarcyStokesTestCase::Schneider;
        else
            DUNE_THROW(Dune::InvalidStateException, testCaseName + " is not a valid test case");
    }();

    // the problem (initial and boundary conditions)
    using FreeFlowMomentumProblem = GetPropType<FreeFlowMomentumTypeTag, Properties::Problem>;
    auto freeFlowMomentumProblem = std::make_shared<FreeFlowMomentumProblem>(freeFlowMomentumGridGeometry, couplingManager, testCase, testCaseName);
    using FreeFlowMassProblem = GetPropType<FreeFlowMassTypeTag, Properties::Problem>;
    auto freeFlowMassProblem = std::make_shared<FreeFlowMassProblem>(freeFlowMassGridGeometry, couplingManager, testCase, testCaseName);
    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto spatialParams = std::make_shared<typename DarcyProblem::SpatialParams>(darcyGridGeometry, testCase);
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyGridGeometry, couplingManager, spatialParams, testCase, testCaseName);

    // the solution vector
    Traits::SolutionVector sol;
    sol[freeFlowMomentumIndex].resize(freeFlowMomentumGridGeometry->numDofs());
    sol[freeFlowMassIndex].resize(freeFlowMassGridGeometry->numDofs());
    sol[porousMediumIndex].resize(darcyGridGeometry->numDofs());

    // the grid variables
    using FreeFlowMomentumGridVariables = GetPropType<FreeFlowMomentumTypeTag, Properties::GridVariables>;
    auto freeFlowMomentumGridVariables = std::make_shared<FreeFlowMomentumGridVariables>(freeFlowMomentumProblem, freeFlowMomentumGridGeometry);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    using FreeFlowMassGridVariables = GetPropType<FreeFlowMassTypeTag, Properties::GridVariables>;
    auto freeFlowMassGridVariables = std::make_shared<FreeFlowMassGridVariables>(freeFlowMassProblem, freeFlowMassGridGeometry);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables = std::make_shared<DarcyGridVariables>(darcyProblem, darcyGridGeometry);

    couplingManager->init(freeFlowMomentumProblem, freeFlowMassProblem, darcyProblem,
                          std::make_tuple(freeFlowMomentumGridVariables, freeFlowMassGridVariables, darcyGridVariables),
                          sol);

    freeFlowMomentumGridVariables->init(sol[freeFlowMomentumIndex]);
    freeFlowMassGridVariables->init(sol[freeFlowMassIndex]);
    darcyGridVariables->init(sol[porousMediumIndex]);

    // initialize the vtk output module
    VtkOutputModule freeflowVtkWriter(*freeFlowMassGridVariables, sol[freeFlowMassIndex], freeFlowMassProblem->name());
    GetPropType<FreeFlowMassTypeTag, Properties::IOFields>::initOutputModule(freeflowVtkWriter); // Add model specific output fields
    freeflowVtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<FreeFlowMassGridVariables>>());
    NavierStokesTest::AnalyticalSolutionVectors freeFlowAnalyticSol(freeFlowMomentumProblem, freeFlowMassProblem);
    freeflowVtkWriter.addField(freeFlowAnalyticSol.analyticalPressureSolution(), "pressureExact");
    freeflowVtkWriter.addField(freeFlowAnalyticSol.analyticalVelocitySolution(), "velocityExact");
    freeflowVtkWriter.write(0.0);

    VtkOutputModule darcyVtkWriter(*darcyGridVariables, sol[porousMediumIndex],  darcyProblem->name());
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(darcyVtkWriter);
    darcyVtkWriter.addVelocityOutput(std::make_shared<GetPropType<DarcyTypeTag, Properties::VelocityOutput>>(*darcyGridVariables));
    const auto darcyAnalyticalSolution = createDarcyAnalyticalSolution<double>(*darcyProblem);
    darcyVtkWriter.addField(std::get<0>(darcyAnalyticalSolution), "pressureExact");
    darcyVtkWriter.addField(std::get<1>(darcyAnalyticalSolution), "velocityExact");
    darcyVtkWriter.write(0.0);

    // the assembler for a stationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(freeFlowMomentumProblem, freeFlowMassProblem, darcyProblem),
                                                 std::make_tuple(freeFlowMomentumGridGeometry,
                                                                 freeFlowMassGridGeometry,
                                                                 darcyGridGeometry),
                                                 std::make_tuple(freeFlowMomentumGridVariables,
                                                                 freeFlowMassGridVariables,
                                                                 darcyGridVariables),
                                                 couplingManager);

    // the linear solver
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    const auto dofsVS = freeFlowMomentumGridGeometry->numDofs();
    const auto dofsPS = freeFlowMassGridGeometry->numDofs();
    const auto dofsPD = darcyGridGeometry->numDofs();
    std::cout << "Stokes velocity dofs: " << dofsVS << std::endl;
    std::cout << "Stokes pressure dofs: " << dofsPS << std::endl;
    std::cout << "Darcy pressure dofs: " << dofsPD << std::endl;
    std::cout << "Total number of dofs: " << dofsVS + dofsPS + dofsPD << std::endl;

    // solve the non-linear system
    nonLinearSolver.solve(sol);

    // write vtk output
    freeflowVtkWriter.write(1.0);
    darcyVtkWriter.write(1.0);

    // print L2 errors
    printFreeFlowL2Error(freeFlowMomentumProblem, freeFlowMassProblem, partial(sol, freeFlowMomentumIndex, freeFlowMassIndex));
    printDarcyL2Error(darcyProblem, sol[porousMediumIndex]);

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
}
