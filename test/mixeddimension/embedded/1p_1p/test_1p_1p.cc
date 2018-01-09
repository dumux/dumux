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

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/nonlinear/mixeddimensionnewtoncontroller.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/mixeddimension/embedded/cellcentered/bboxtreecouplingmanagersimple.hh>
#include <dumux/mixeddimension/embedded/integrationpointsource.hh>

#include "bloodflowproblem.hh"
#include "tissueproblem.hh"

namespace Dumux {
namespace Properties {

SET_PROP(TissueCCTypeTag, CouplingManager)
{
    using Traits = MultiDomainTraits<TypeTag, TTAG(BloodFlowCCTypeTag)>;
    using type = Dumux::CCBBoxTreeEmbeddedCouplingManagerSimple<Traits>;
};

SET_PROP(BloodFlowCCTypeTag, CouplingManager)
{
    using Traits = MultiDomainTraits<TTAG(TissueCCTypeTag), TypeTag>;
    using type = Dumux::CCBBoxTreeEmbeddedCouplingManagerSimple<Traits>;
};

SET_TYPE_PROP(TissueCCTypeTag, PointSource, IntegrationPointSource<TypeTag>);
SET_TYPE_PROP(BloodFlowCCTypeTag, PointSource, IntegrationPointSource<TypeTag>);

SET_TYPE_PROP(TissueCCTypeTag, PointSourceHelper, IntegrationPointSourceHelper);
SET_TYPE_PROP(BloodFlowCCTypeTag, PointSourceHelper, IntegrationPointSourceHelper);

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
    using BulkTypeTag = TTAG(TissueCCTypeTag);
    using LowDimTypeTag = TTAG(BloodFlowCCTypeTag);

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using BulkGridCreator = typename GET_PROP_TYPE(BulkTypeTag, GridCreator);
    BulkGridCreator::makeGrid("Tissue"); // pass parameter group

    using LowDimGridCreator = typename GET_PROP_TYPE(LowDimTypeTag, GridCreator);
    LowDimGridCreator::makeGrid("Vessel"); // pass parameter group

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
    using CouplingManager = CCBBoxTreeEmbeddedCouplingManagerSimple<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>(bulkFvGridGeometry, lowDimFvGridGeometry);

    // the problem (initial and boundary conditions)
    using BulkProblem = typename GET_PROP_TYPE(BulkTypeTag, Problem);
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, couplingManager);
    bulkProblem->computePointSourceMap();
    using LowDimProblem = typename GET_PROP_TYPE(LowDimTypeTag, Problem);
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, couplingManager);
    lowDimProblem->computePointSourceMap();

    // the solution vector
    auto sol = std::make_shared<Traits::SolutionVector>();
    (*sol)[bulkIdx].resize(bulkFvGridGeometry->numDofs());
    (*sol)[lowDimIdx].resize(lowDimFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution((*sol)[bulkIdx]);
    lowDimProblem->applyInitialSolution((*sol)[lowDimIdx]);
    auto oldSol = std::make_shared<Traits::SolutionVector>(*sol);

    couplingManager->init(bulkProblem, lowDimProblem, sol);

    // the grid variables
    using BulkGridVariables = typename GET_PROP_TYPE(BulkTypeTag, GridVariables);
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    bulkGridVariables->init((*sol)[bulkIdx], (*oldSol)[bulkIdx]);
    using LowDimGridVariables = typename GET_PROP_TYPE(LowDimTypeTag, GridVariables);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    lowDimGridVariables->init((*sol)[lowDimIdx], (*oldSol)[lowDimIdx]);

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDivisions = getParam<int>("TimeLoop.MaxTimeStepDivisions");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // intialize the vtk output module
    VtkOutputModule<BulkTypeTag> bulkVtkWriter(*bulkProblem, *bulkFvGridGeometry, *bulkGridVariables, (*sol)[bulkIdx], bulkProblem->name());
    GET_PROP_TYPE(BulkTypeTag, VtkOutputFields)::init(bulkVtkWriter);
    bulkProblem->addVtkOutputFields(bulkVtkWriter);
    bulkVtkWriter.write(0.0);

    VtkOutputModule<LowDimTypeTag> lowDimVtkWriter(*lowDimProblem, *lowDimFvGridGeometry, *lowDimGridVariables, (*sol)[lowDimIdx], lowDimProblem->name());
    GET_PROP_TYPE(LowDimTypeTag, VtkOutputFields)::init(lowDimVtkWriter);
    lowDimVtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    // using Assembler = FVAssembler<Traits, DiffMethod::numeric>;
    // auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonController = MixedDimensionNewtonController<double>;
    auto newtonController = std::make_shared<NewtonController>(timeLoop);
    // using NewtonMethod = NewtonMethod<NewtonController, Assembler, LinearSolver>;
    // NewtonMethod nonLinearSolver(newtonController, assembler, linearSolver);

    typename Traits::JacobianMatrix M;

    auto& A11 = M[bulkIdx][bulkIdx];
    auto& A12 = M[bulkIdx][lowDimIdx];
    auto& A22 = M[lowDimIdx][bulkIdx];
    auto& A21 = M[lowDimIdx][lowDimIdx];

    A11.setSize(3, 3);
    A12.setSize(3, 8);
    A22.setSize(8, 8);
    A21.setSize(8, 3);

    {
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(3, 3);
        for (int i = 0; i < occupationPattern.rows(); ++i)
            occupationPattern.add(i,i);
        occupationPattern.exportIdx(A11);
    }

    {
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(3, 8);
        for (int i = 0; i < occupationPattern.rows(); ++i)
            occupationPattern.add(i,i);
        occupationPattern.exportIdx(A12);
    }

    {
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(8, 3);
        for (int i = 0; i < 3; ++i)
            occupationPattern.add(i,i);
        occupationPattern.exportIdx(A21);
    }

    {
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(8, 8);
        for (int i = 0; i < occupationPattern.rows(); ++i)
            occupationPattern.add(i,i);
        occupationPattern.exportIdx(A22);
    }

    M = 1.0;

    Dune::printmatrix(std::cout, A11, "", "");
    Dune::printmatrix(std::cout, A12, "", "");
    Dune::printmatrix(std::cout, A21, "", "");
    Dune::printmatrix(std::cout, A22, "", "");


    // time loop
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // set previous solution for storage evaluations
        // assembler->setPreviousSolution(solOld);

        // try solving the non-linear system
        for (int i = 0; i < maxDivisions; ++i)
        {
            // linearize & solve
            // auto converged = nonLinearSolver.solve(sol);

            // if (converged)
            //     break;
            //
            // if (!converged && i == maxDivisions-1)
            //     DUNE_THROW(Dune::MathError,
            //                "Newton solver didn't converge after "
            //                << maxDivisions
            //                << " time-step divisions. dt="
            //                << timeLoop->timeStepSize()
            //                << ".\nThe solutions of the current and the previous time steps "
            //                << "have been saved to restart files.");
        }

        // make the new solution the old solution
        (*oldSol) = (*sol);
        bulkGridVariables->advanceTimeStep();
        lowDimGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // output the source terms
        bulkProblem->computeSourceIntegral((*sol)[bulkIdx], *bulkGridVariables);
        lowDimProblem->computeSourceIntegral((*sol)[lowDimIdx], *lowDimGridVariables);

        // write vtk output
        bulkVtkWriter.write(timeLoop->time());
        lowDimVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(newtonController->suggestTimeStepSize(timeLoop->timeStepSize()));
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
