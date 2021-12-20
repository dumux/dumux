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
 * \brief Test for the staggered grid Navier-Stokes model (Donea 2003, \cite Donea2003).
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/integrate.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/incompressiblestokessolver.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors.hh>

#include "properties.hh"

namespace Dumux {

template<class MassGridProblem, class MomentumProblem, class CouplingManager, class SolutionVector, class AnalyticSolVecs>
void printErrors(std::shared_ptr<MassGridProblem> massProblem,
                 std::shared_ptr<MomentumProblem> momentumProblem,
                 std::shared_ptr<CouplingManager> couplingManager,
                 const SolutionVector& x,
                 const AnalyticSolVecs& analyticalSolVectors)
{
    const bool printErrors = getParam<bool>("Problem.PrintErrors", false);
    const bool printConvergenceTestFile = getParam<bool>("Problem.PrintConvergenceTestFile", false);
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;

    if (printErrors || printConvergenceTestFile)
    {
        using MomentumGridGeometry = std::decay_t<decltype(momentumProblem->gridGeometry())>;
        if constexpr (MomentumGridGeometry::discMethod == DiscretizationMethods::fcstaggered)
        {
            // print discrete L2 and Linfity errors
            NavierStokesTest::Errors errors(momentumProblem, massProblem, x);
            NavierStokesTest::ErrorCSVWriter csvWriter(
                momentumProblem, massProblem, std::to_string(x[massIdx].size())
            );
            csvWriter.printErrors(errors);

            if (printConvergenceTestFile)
                NavierStokesTest::convergenceTestAppendErrors(momentumProblem, massProblem, errors);
        }
        else if (MomentumGridGeometry::discMethod == DiscretizationMethods::fcdiamond)
        {
            // convert format to solution vector type
            const auto& vOnFace = analyticalSolVectors.analyticalVelocitySolutionOnFace();
            const auto& pOnDof = analyticalSolVectors.analyticalPressureSolution();
            auto analyticVel = x[momentumIdx];
            auto analyticPressure = x[massIdx];
            auto momFVGeometry = localView(momentumProblem->gridGeometry());
            auto massFVGeometry = localView(massProblem->gridGeometry());
            for (const auto& element : elements(massProblem->gridGeometry().gridView()))
            {
                momFVGeometry.bindElement(element);
                for (const auto& scv : scvs(momFVGeometry))
                    analyticVel[scv.dofIndex()] = vOnFace[scv.dofIndex()];

                massFVGeometry.bindElement(element);
                for (const auto& scv : scvs(massFVGeometry))
                    analyticPressure[scv.dofIndex()] = pOnDof[scv.dofIndex()];
            }

            const auto l2ErrorVel = integrateL2Error(
                momentumProblem->gridGeometry(), x[momentumIdx], analyticVel, 3
            );

            auto zeroVel = analyticVel; zeroVel = 0.0;
            const auto l2ErrorNormalizerVel = integrateL2Error(
                momentumProblem->gridGeometry(), zeroVel, analyticVel, 3
            );

            std::cout << "Velocity L2 error: " << l2ErrorVel/l2ErrorNormalizerVel << std::endl;

            const auto l2ErrorPressure = integrateL2Error(
                massProblem->gridGeometry(), x[massIdx], analyticPressure, 3
            );

            auto zeroP = analyticPressure; zeroP = 0.0;
            const auto l2ErrorNormalizerP = integrateL2Error(
                massProblem->gridGeometry(), zeroP, analyticPressure, 3
            );

            std::cout << "Pressure L2 error: " << l2ErrorPressure/l2ErrorNormalizerP << std::endl;
        }
    }
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tags for this problem
    using MomentumTypeTag = Properties::TTag::DoneaTestMomentum;
    using MassTypeTag = Properties::TTag::DoneaTestMass;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

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
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);

    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // mass-momentum coupling manager
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using CouplingManager = GetPropType<MassTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);

    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
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
    // intializing the gridvariables requires the coupling manager to be set up
    momentumGridVariables->init(x[momentumIdx]);
    massGridVariables->init(x[massIdx]);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    std::string discSuffix = std::string("_") + MomentumGridGeometry::DiscretizationMethod::name();
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name() + discSuffix);
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
    using LinearSolver = IncompressibleStokesSolver<typename Assembler::JacobianMatrix, typename Assembler::ResidualType>;
    auto linearSolver = std::make_shared<LinearSolver>(*couplingManager);

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // linearize & solve
    nonLinearSolver.solve(x);

    // maybe print errors
    printErrors(massProblem, momentumProblem, couplingManager, x, analyticalSolVectors);

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
