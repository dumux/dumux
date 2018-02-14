// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief Test for the 1d-3d embedded mixed-dimension model coupling two
 *        one-phase porous medium flow problems
 */
#include <config.h>

#include <ctime>
#include <iostream>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/geometry/diameter.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtoncontroller.hh>
#include <dumux/mixeddimension/embedded/integrationpointsource.hh>

#include "matrixproblem.hh"
#include "fractureproblem.hh"
#include "couplingmanager.hh"

namespace Dumux {
namespace Properties {

SET_PROP(MatrixTypeTag, CouplingManager)
{
    using Traits = MultiDomainTraits<TypeTag, TTAG(FractureTypeTag)>;
    using type = Dumux::EmbeddedFractureCouplingManager<Traits>;
};

SET_PROP(FractureTypeTag, CouplingManager)
{
    using Traits = MultiDomainTraits<TTAG(MatrixTypeTag), TypeTag>;
    using type = Dumux::EmbeddedFractureCouplingManager<Traits>;
};

SET_TYPE_PROP(MatrixTypeTag, PointSource, IntegrationPointSource<TypeTag>);
SET_TYPE_PROP(FractureTypeTag, PointSource, IntegrationPointSource<TypeTag>);

SET_TYPE_PROP(MatrixTypeTag, PointSourceHelper, IntegrationPointSourceHelper);
SET_TYPE_PROP(FractureTypeTag, PointSourceHelper, IntegrationPointSourceHelper);

SET_STRING_PROP(MatrixTypeTag, ModelParameterGroup, "Matrix");
SET_STRING_PROP(FractureTypeTag, ModelParameterGroup, "Fracture");

} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using BulkTypeTag = TTAG(MatrixTypeTag);
    using LowDimTypeTag = TTAG(FractureTypeTag);

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using BulkGridCreator = typename GET_PROP_TYPE(BulkTypeTag, GridCreator);
    BulkGridCreator::makeGrid("Matrix"); // pass parameter group

    using LowDimGridCreator = typename GET_PROP_TYPE(LowDimTypeTag, GridCreator);
    LowDimGridCreator::makeGrid("Fracture"); // pass parameter group

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& bulkGridView = BulkGridCreator::grid().leafGridView();
    const auto& lowDimGridView = LowDimGridCreator::grid().leafGridView();

    // create the finite volume grid geometry
    using BulkFVGridGeometry = typename GET_PROP_TYPE(BulkTypeTag, FVGridGeometry);
    auto bulkFvGridGeometry = std::make_shared<BulkFVGridGeometry>(bulkGridView);
    bulkFvGridGeometry->update();
    using LowDimFVGridGeometry = typename GET_PROP_TYPE(LowDimTypeTag, FVGridGeometry);
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);
    lowDimFvGridGeometry->update();

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template DomainIdx<0>();
    constexpr auto lowDimIdx = Traits::template DomainIdx<1>();

    // the coupling manager
    using CouplingManager = typename GET_PROP_TYPE(BulkTypeTag, CouplingManager);
    auto couplingManager = std::make_shared<CouplingManager>(bulkFvGridGeometry, lowDimFvGridGeometry);

    // the problem (initial and boundary conditions)
    using BulkProblem = typename GET_PROP_TYPE(BulkTypeTag, Problem);
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, couplingManager);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimTypeTag, Problem);
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, couplingManager);

    // the solution vector
    Traits::SolutionVector sol;
    sol[bulkIdx].resize(bulkFvGridGeometry->numDofs());
    sol[lowDimIdx].resize(lowDimFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(sol[bulkIdx]);
    lowDimProblem->applyInitialSolution(sol[lowDimIdx]);
    auto oldSol = sol;

    couplingManager->init(bulkProblem, lowDimProblem, sol);
    bulkProblem->computePointSourceMap();
    lowDimProblem->computePointSourceMap();

    // the grid variables
    using BulkGridVariables = typename GET_PROP_TYPE(BulkTypeTag, GridVariables);
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    bulkGridVariables->init(sol[bulkIdx], oldSol[bulkIdx]);
    using LowDimGridVariables = typename GET_PROP_TYPE(LowDimTypeTag, GridVariables);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    lowDimGridVariables->init(sol[lowDimIdx], oldSol[lowDimIdx]);

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDivisions = getParam<int>("TimeLoop.MaxTimeStepDivisions");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // intialize the vtk output module
    VtkOutputModule<BulkTypeTag> bulkVtkWriter(*bulkProblem, *bulkFvGridGeometry, *bulkGridVariables, sol[bulkIdx], bulkProblem->name());
    GET_PROP_TYPE(BulkTypeTag, VtkOutputFields)::init(bulkVtkWriter);
    bulkVtkWriter.write(0.0);

    VtkOutputModule<LowDimTypeTag> lowDimVtkWriter(*lowDimProblem, *lowDimFvGridGeometry, *lowDimGridVariables, sol[lowDimIdx], lowDimProblem->name());
    GET_PROP_TYPE(LowDimTypeTag, VtkOutputFields)::init(lowDimVtkWriter);
    lowDimVtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager, timeLoop);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonController = MultiDomainNewtonController<double, CouplingManager>;
    auto newtonController = std::make_shared<NewtonController>(timeLoop, couplingManager);
    using NewtonMethod = NewtonMethod<NewtonController, Assembler, LinearSolver>;
    NewtonMethod nonLinearSolver(newtonController, assembler, linearSolver);

    // time loop
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(oldSol);

        // try solving the non-linear system
        for (int i = 0; i < maxDivisions; ++i)
        {
            // assemble
            // assembler->assembleJacobianAndResidual(sol);
            // auto& b = assembler->residual(); b *= -1.0;

            // solve
            // linearSolver->template solve<2>(assembler->jacobian(), sol, b);
            auto converged = nonLinearSolver.solve(sol);
            // auto converged = true;

            // update
            // couplingManager->updateSolution(sol);
            // bulkGridVariables->update(sol[bulkIdx]);
            // lowDimGridVariables->update(sol[lowDimIdx]);

            const auto& jac = assembler->jacobian();
            const auto& res = assembler->residual();
            using namespace Dune::Indices;
            if (jac[_0][_0].N() < 28)
            {
                Dune::printmatrix(std::cout, jac[_0][_0], "", "", 15, 15);
                Dune::printmatrix(std::cout, jac[_0][_1], "", "", 15, 15);
                Dune::printmatrix(std::cout, jac[_1][_0], "", "", 15, 15);
                Dune::printmatrix(std::cout, jac[_1][_1], "", "", 15, 15);
                Dune::printvector(std::cout, res[_0], "", "", 15, 15);
                Dune::printvector(std::cout, res[_1], "", "", 15, 15);
                Dune::printvector(std::cout, sol[_0], "", "", 15, 15);
                Dune::printvector(std::cout, sol[_1], "", "", 15, 15);
            }

            if (converged)
                break;

            if (!converged && i == maxDivisions-1)
                DUNE_THROW(Dune::MathError,
                           "Newton solver didn't converge after "
                           << maxDivisions
                           << " time-step divisions. dt="
                           << timeLoop->timeStepSize()
                           << ".\nThe solutions of the current and the previous time steps "
                           << "have been saved to restart files.");
        }

        // make the new solution the old solution
        oldSol = sol;
        bulkGridVariables->advanceTimeStep();
        lowDimGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // output the source terms
        bulkProblem->computeSourceIntegral(sol[bulkIdx], *bulkGridVariables);
        lowDimProblem->computeSourceIntegral(sol[lowDimIdx], *lowDimGridVariables);

        // write vtk output
        bulkVtkWriter.write(timeLoop->time());
        lowDimVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        // timeLoop->setTimeStepSize(newtonController->suggestTimeStepSize(timeLoop->timeStepSize()));
    }

    timeLoop->finalize(mpiHelper.getCollectiveCommunication());

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
