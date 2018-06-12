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
 * \brief test for the 2p tracer CC model
 */
#include <config.h>

#include "2ptestproblem.hh"
#include "2ptracertestproblem.hh"

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/valgrind.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/methods.hh>

#include <dumux/io/vtkoutputmodule.hh>

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */



void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe list of mandatory arguments for this program is:\n"
                                        "\t-TimeManager.TEnd               End of the simulation [s] \n"
                                        "\t-TimeManager.DtInitial          Initial timestep size [s] \n"
                                        "\t-Grid.LowerLeft                 Lower left corner coordinates\n"
                                        "\t-Grid.UpperRight                Upper right corner coordinates\n"
                                        "\t-Grid.Cells                     Number of cells in respective coordinate directions\n"
                                        "\t                                definition in DGF format\n"
                                        "\t-SpatialParams.LensLowerLeft   coordinates of the lower left corner of the lens [m] \n"
                                        "\t-SpatialParams.LensUpperRight  coordinates of the upper right corner of the lens [m] \n"
                                        "\t-SpatialParams.Permeability     Permeability of the domain [m^2] \n"
                                        "\t-SpatialParams.PermeabilityLens Permeability of the lens [m^2] \n";

        std::cout << errorMessageOut
                  << "\n";
    }
}


int main(int argc, char** argv) try
{

    std::cout << "Halloooooooo:";

    using namespace Dumux;

    //! define the type tags for this problem
    using TwoPTypeTag = TTAG(TwoPIncompressibleTpfa);
    using TracerTypeTag = TTAG(TwoPTracerTestCCTypeTag);

    //! initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    //! print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    //! parse command line arguments
    Parameters::init(argc, argv, usage);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    // only create the grid once using the 1p type tag
    // try to create a grid (from the given grid file or the input file)
    using GridCreator = typename GET_PROP_TYPE(TwoPTypeTag, GridCreator);
    try { GridCreator::makeGrid(); }
    catch (...) {
        std::cout << "\n\t -> Creation of the grid failed! <- \n\n";
        throw;
    }
    GridCreator::loadBalance();


    // get some time loop parameters
    using Scalar = typename GET_PROP_TYPE(TwoPTypeTag, Scalar);
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))
        restartTime = getParam<Scalar>("TimeLoop.Restart");

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);  // restartTime anstelle 0.0
    // auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);


    ////////////////////////////////////////////////////////////
    // set 2p Problem
    ////////////////////////////////////////////////////////////


    // we compute on the leaf grid view
    const auto& leafGridView = GridCreator::grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = typename GET_PROP_TYPE(TwoPTypeTag, FVGridGeometry);
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (initial and boundary conditions)
    using TwoPProblem = typename GET_PROP_TYPE(TwoPTypeTag, Problem);
    auto twoPProblem = std::make_shared<TwoPProblem>(fvGridGeometry);

    // the solution vector
    using TwoPSolutionVector = typename GET_PROP_TYPE(TwoPTypeTag, SolutionVector);
    TwoPSolutionVector p(fvGridGeometry->numDofs());
    twoPProblem->applyInitialSolution(p);
    auto pOld = p;

    // maybe update the interface parameters
    if (ENABLEINTERFACESOLVER)
        twoPProblem->spatialParams().updateMaterialInterfaceParams(p);

    // the grid variables
    using TwoPGridVariables = typename GET_PROP_TYPE(TwoPTypeTag, GridVariables);
    auto twoPGridVariables = std::make_shared<TwoPGridVariables>(twoPProblem, fvGridGeometry);
    twoPGridVariables->init(p, pOld);

    // intialize the vtk output module
    using TwoPVtkOutputFields = typename GET_PROP_TYPE(TwoPTypeTag, VtkOutputFields);

    // use non-conforming output for the test with interface solver
    const auto ncOutput = getParam<bool>("Problem.UseNonConformingOutput", false);
    VtkOutputModule<TwoPTypeTag> twoPVtkWriter(*twoPProblem, *fvGridGeometry, *twoPGridVariables, p, twoPProblem->name(), "",
                                       (ncOutput ? Dune::VTK::nonconforming : Dune::VTK::conforming));

    TwoPVtkOutputFields::init(twoPVtkWriter); //!< Add model specific output fields
    twoPVtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using TwoPAssembler = FVAssembler<TwoPTypeTag, DiffMethod::numeric>;
    auto twoPAssembler = std::make_shared<TwoPAssembler>(twoPProblem, fvGridGeometry, twoPGridVariables, timeLoop);

    // the linear solver
    using LinearSolver = AMGBackend<TwoPTypeTag>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, fvGridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<TwoPAssembler, LinearSolver>;
    NewtonSolver nonLinearSolver(twoPAssembler, linearSolver);


    ////////////////////////////////////////////////////////////
    // set tracer Problem
    ////////////////////////////////////////////////////////////

     //! the problem (initial and boundary conditions)
    using TracerProblem = typename GET_PROP_TYPE(TracerTypeTag, Problem);
    auto tracerProblem = std::make_shared<TracerProblem>(fvGridGeometry);

    //! the solution vector
    using TracerSolutionVector = typename GET_PROP_TYPE(TracerTypeTag, SolutionVector);
    TracerSolutionVector x(leafGridView.size(0));
    tracerProblem->applyInitialSolution(x);
    auto xOld = x;

    //! the grid variables
    using TracerGridVariables = typename GET_PROP_TYPE(TracerTypeTag, GridVariables);
    auto tracerGridVariables = std::make_shared<TracerGridVariables>(tracerProblem, fvGridGeometry);
    tracerGridVariables->init(x, xOld);

     //! the linear system
    using JacobianMatrix = typename GET_PROP_TYPE(TracerTypeTag, JacobianMatrix);
    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<TracerSolutionVector>();

    //! the assembler with time loop for instationary problem
    using TracerAssembler = FVAssembler<TracerTypeTag, DiffMethod::analytic, /*implicit=*/false>;
    auto tracerAssembler = std::make_shared<TracerAssembler>(tracerProblem, fvGridGeometry, tracerGridVariables, timeLoop);
    tracerAssembler->setLinearSystem(A, r);

    //! initialize the vtk output module
    VtkOutputModule<TracerTypeTag> vtkWriter(*tracerProblem, *fvGridGeometry, *tracerGridVariables, x, tracerProblem->name());
    using TracerVtkOutputFields = typename GET_PROP_TYPE(TracerTypeTag, VtkOutputFields);
    TracerVtkOutputFields::init(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);

    //! initialize the flux vector
    using ModelTraits = typename GET_PROP_TYPE(TwoPTypeTag, ModelTraits);
    static const int numPhases = ModelTraits::numPhases();
    using ScvfVector  = std::vector<Scalar>;
    using FieldVector = std::vector<Dune::FieldVector<ScvfVector, numPhases>>;
    ScvfVector phaseFlux;
    FieldVector volumeFlux;

    for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        volumeFlux[phaseIdx] = ScvfVector(fvGridGeometry->numScvf(), 0.0);
    }


    //! set some check points for the time loop
    timeLoop->setPeriodicCheckPoint(tEnd/10.0);

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem
    ////////////////////////////////////////////////////////////


    // time loop
    timeLoop->start(); do
    {
        // set previous solution for storage evaluations
        twoPAssembler->setPreviousSolution(pOld);

        // solve the non-linear system with time step control
        nonLinearSolver.solve(p, *timeLoop);

        // make the new solution the old solution
        pOld = p;
        twoPGridVariables->advanceTimeStep();

        // // advance to the time loop to the next step
        // timeLoop->advanceTimeStep();

        // write vtk output
        twoPVtkWriter.write(timeLoop->time());

        // report statistics of this time step              ??
        timeLoop->reportTimeStep();

        // set new dt as suggested by the Newton solver     ??
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        // } while (!timeLoop->finished());

        // timeLoop->finalize(leafGridView.comm());

        ////////////////////////////////////////////////////////////
        // compute volume fluxes for the tracer model
        ///////////////////////////////////////////////////////////

        using FluxVariables =  typename GET_PROP_TYPE(TwoPTypeTag, FluxVariables);
        auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };
        for (const auto& element : elements(leafGridView))
        {
            auto fvGeometry = localView(*fvGridGeometry);
            fvGeometry.bind(element);

            auto elemVolVars = localView(twoPGridVariables->curGridVolVars());
            elemVolVars.bind(element, fvGeometry, p);

            auto elemFluxVars = localView(twoPGridVariables->gridFluxVarsCache());
            elemFluxVars.bind(element, fvGeometry, elemVolVars);

            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    const auto idx = scvf.index();

                    if (!scvf.boundary())
                    {
                        FluxVariables fluxVars;
                        fluxVars.init(*twoPProblem, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                        phaseFlux[idx] = fluxVars.advectiveFlux(phaseIdx, upwindTerm);
                    }
                    else
                    {
                        const auto bcTypes = twoPProblem->boundaryTypes(element, scvf);
                        if (bcTypes.hasOnlyDirichlet())
                        {
                            FluxVariables fluxVars;
                            fluxVars.init(*twoPProblem, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                            phaseFlux[idx] = fluxVars.advectiveFlux(phaseIdx, upwindTerm);
                        }
                    }
                }
                volumeFlux[phaseIdx] = phaseFlux;

            }
        }

        ////////////////////////////////////////////////////////////
        // solve tracer problem on the same grid / run instationary non-linear simulation
        ////////////////////////////////////////////////////////////

        // set the flux from the 1p problem
        tracerProblem->spatialParams().setVolumeFlux(volumeFlux);

        //! get some time loop parameters
        //const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
        //auto dt = getParam<Scalar>("TimeLoop.DtInitial");
        //const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");

        // //! instantiate time loop
        // auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
        // timeLoop->setMaxTimeStepSize(maxDt);

        // //! set some check points for the time loop
        // timeLoop->setPeriodicCheckPoint(tEnd/10.0);

        // //! start the time loop
        // timeLoop->start(); do
        // {

        // set previous solution for storage evaluations
        tracerAssembler->setPreviousSolution(xOld);

        Dune::Timer tracerAssembleTimer;
        tracerAssembler->assembleJacobianAndResidual(x);
        tracerAssembleTimer.stop();

        // solve the linear system A(xOld-xNew) = r
        Dune::Timer solveTimer;
        TracerSolutionVector xDelta(x);
        linearSolver->solve(*A, xDelta, *r);
        solveTimer.stop();

        // update solution and grid variables
        Dune::Timer updateTimer;
        updateTimer.reset();
        x -= xDelta;
        tracerGridVariables->update(x);
        updateTimer.stop();

        // statistics
        Dune::Timer assembleTimer;
        const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
        std::cout << "Assemble/solve/update time: "
                  <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                  <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                  <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                  <<  std::endl;

        // make the new solution the old solution
        xOld = x;
        tracerGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output on check points
        if (timeLoop->isCheckPoint())
            twoPVtkWriter.write(timeLoop->time());
            vtkWriter.write(timeLoop->time());        // CHECKPOINTS MIT BEIDEN .WRITER ... ?!

        //  // report statistics of this time step
        //  timeLoop->reportTimeStep();

        //  // set new dt
        //  timeLoop->setTimeStepSize(dt);

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    //! print dumux end message
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
