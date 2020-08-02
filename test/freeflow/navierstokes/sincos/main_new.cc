// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the instationary staggered grid Navier-Stokes model
 *        with analytical solution.
 */

#include <config.h>

#include <ctime>
#include <iostream>
#include <type_traits>
#include <tuple>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/istlsolverfactorybackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "../l2error.hh"
#include "../analyticalsolution.hh"
#include "problem_new.hh"

template<class MomentumProblem, class MassProblem, class SolutionVector>
void printL2Error(const MomentumProblem& momentumProblem,
                  const MassProblem& massProblem,
                  const SolutionVector& x)
{
    const auto velocityL2error = calculateL2Error(momentumProblem, x[Dune::index_constant<0>()]);
    const auto pressureL2error = calculateL2Error(massProblem, x[Dune::index_constant<1>()]);
    static constexpr auto pressureIdx = MassProblem::Indices::pressureIdx;
    static constexpr auto velocityXIdx = MomentumProblem::Indices::velocityXIdx;
    static constexpr auto velocityYIdx = MomentumProblem::Indices::velocityYIdx;
    const auto numCellCenterDofs = massProblem.gridGeometry().numDofs();
    const auto numFaceDofs = momentumProblem.gridGeometry().numDofs();

    std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                << std::scientific
                << "L2(p) = " << pressureL2error.absolute[pressureIdx] << " / " << pressureL2error.relative[pressureIdx]
                << " , L2(vx) = " << velocityL2error.absolute[velocityXIdx] << " / " << velocityL2error.relative[velocityXIdx]
                << " , L2(vy) = " << velocityL2error.absolute[velocityYIdx] << " / " << velocityL2error.relative[velocityYIdx]
                << std::endl;

    // write the norm into a log file
    std::ofstream logFile;
    logFile.open(momentumProblem.name() + ".log", std::ios::app);
    logFile << "[ConvergenceTest] L2(p) = " << pressureL2error.absolute[pressureIdx] << " L2(vx) = "
            << velocityL2error.absolute[velocityXIdx] << " L2(vy) = " << velocityL2error.absolute[velocityYIdx] << std::endl;
    logFile.close();
}

namespace Dumux::Properties{

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::SincosTest>
{
private:
    using Traits = MultiDomainTraits<TTag::SincosTestMomentum, TTag::SincosTestMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::SincosTestMomentum;
    using MassTypeTag = Properties::TTag::SincosTestMass;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    momentumGridGeometry->update();

    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);
    massGridGeometry->update();

    // the coupling manager
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using CouplingManager = StaggeredFreeFlowCouplingManager<Traits>;

    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);

    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // get some time loop parameters
    using Scalar = typename Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the problem (initial and boundary conditions)
    momentumProblem->updateTimeStepSize(timeLoop->timeStepSize());
    massProblem->updateTimeStepSize(timeLoop->timeStepSize());

    // the solution vector
    constexpr auto momentumIdx = Dune::index_constant<0>();
    constexpr auto massIdx = Dune::index_constant<1>();
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);
    auto xOld = x;

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);

    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    const bool isStationary = getParam<bool>("Problem.IsStationary");

    if (isStationary)
        couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    else
        couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);

    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    // intialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    auto exactPressure = getScalarAnalyticalSolution(*massProblem)[GetPropType<MassTypeTag, Properties::ModelTraits>::Indices::pressureIdx];
    auto exactVelocity = getVelocityAnalyticalSolution(*momentumProblem);
    vtkWriter.addField(exactPressure, "pressureExact");
    vtkWriter.addField(exactVelocity, "velocityExact");

    vtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = isStationary ? std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                                std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                                std::make_tuple(momentumGridVariables, massGridVariables),
                                                                couplingManager)
                                  : std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                                                std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                                                std::make_tuple(momentumGridVariables, massGridVariables),
                                                                                couplingManager, timeLoop, xOld);

    // the linear solver
    using LinearSolver = LINEARSOLVER;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    const bool shouldPrintL2Error = getParam<bool>("Problem.PrintL2Error");

    if (isStationary)
    {
        // linearize & solve
        Dune::Timer timer;
        nonLinearSolver.solve(x);

        if (shouldPrintL2Error)
            printL2Error(*momentumProblem, *massProblem, x);

        // write vtk output
        vtkWriter.write(1.0);

        timer.stop();

        const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
        std::cout << "Simulation took " << timer.elapsed() << " seconds on "
                << comm.size() << " processes.\n"
                << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";
    }
    else
    {
        // time loop
        timeLoop->start(); do
        {
            // solve the non-linear system with time step control
            nonLinearSolver.solve(x, *timeLoop);

            // make the new solution the old solution
            xOld = x;
            momentumGridVariables->advanceTimeStep();
            massGridVariables->advanceTimeStep();

            if (shouldPrintL2Error)
                printL2Error(*momentumProblem, *massProblem, x);

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();
            momentumProblem->updateTime(timeLoop->time());
            massProblem->updateTime(timeLoop->time());

            // update the analytical solution for the new time step
            exactPressure = getScalarAnalyticalSolution(*massProblem)[GetPropType<MassTypeTag, Properties::ModelTraits>::Indices::pressureIdx];
            exactVelocity = getVelocityAnalyticalSolution(*momentumProblem);

            // write vtk output
            vtkWriter.write(timeLoop->time());

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // set new dt as suggested by newton solver
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
            momentumProblem->updateTimeStepSize(timeLoop->timeStepSize());
            massProblem->updateTimeStepSize(timeLoop->timeStepSize());

        } while (!timeLoop->finished());

        timeLoop->finalize(leafGridView.comm());
    }

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
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
