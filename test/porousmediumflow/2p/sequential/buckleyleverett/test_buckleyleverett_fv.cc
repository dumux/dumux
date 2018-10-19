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
 * \brief test for the two-phase porousmedium flow model
 */
#define PROBLEM 3

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include "pressureproblem.hh"
#include "transportproblem.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/methods.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

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
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////
    using TwoPImpesTT = TTAG(TwoPImpes);
    using TwoPTransportTT = TTAG(TwoPTransport);

    // we simply extract the grid creator from one of the type tags
    using GridManager = Dumux::GridManager<typename GET_PROP_TYPE(TwoPImpesTT, Grid)>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometries
    using TwoPImpesFVGridGeometry = typename GET_PROP_TYPE(TwoPImpesTT, FVGridGeometry);
    using TwoPTransportFVGridGeometry = typename GET_PROP_TYPE(TwoPTransportTT, FVGridGeometry);
    auto twoPImpesFvGridGeometry = std::make_shared<TwoPImpesFVGridGeometry>(leafGridView);
    auto twoPTransportFvGridGeometry = std::make_shared<TwoPTransportFVGridGeometry>(leafGridView);
    twoPImpesFvGridGeometry->update();
    twoPTransportFvGridGeometry->update();

    // the problems (boundary conditions)
    using TwoPImpes = typename GET_PROP_TYPE(TwoPImpesTT, Problem);
    using TwoPTransportProblem = typename GET_PROP_TYPE(TwoPTransportTT, Problem);
    auto twoPImpesProblem = std::make_shared<TwoPImpes>(twoPImpesFvGridGeometry);
    auto twoPTransportProblem = std::make_shared<TwoPTransportProblem>(twoPTransportFvGridGeometry);

    // the solution vector
    using SVTwoPImpes = typename GET_PROP_TYPE(TwoPImpesTT, SolutionVector);
    SVTwoPImpes x_p(twoPImpesFvGridGeometry->numDofs());
    twoPImpesProblem->applyInitialSolution(x_p);
    auto x_p_old = x_p;

    using SVTwoPTransport = typename GET_PROP_TYPE(TwoPTransportTT, SolutionVector);
    SVTwoPTransport x_s(twoPTransportFvGridGeometry->numDofs());
    twoPTransportProblem->applyInitialSolution(x_s);
    auto x_s_old = x_s;

    // the grid variables
    using GVTwoPImpes = typename GET_PROP_TYPE(TwoPImpesTT, GridVariables);
    auto gridVarTwoPImpes = std::make_shared<GVTwoPImpes>(twoPImpesProblem, twoPImpesFvGridGeometry);
    gridVarTwoPImpes->init(x_p, x_p_old);

    using GVTwoPTransport = typename GET_PROP_TYPE(TwoPTransportTT, GridVariables);
    auto gridVarTwoPTransport = std::make_shared<GVTwoPTransport>(twoPTransportProblem, twoPTransportFvGridGeometry);
    gridVarTwoPTransport->init(x_s, x_s_old);

    // get some time loop parameters
    using Scalar = typename GET_PROP_TYPE(TwoPImpesTT, Scalar);
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))
        restartTime = getParam<Scalar>("TimeLoop.Restart");

    // intialize the vtk output module
    using VtkOutputFields = typename GET_PROP_TYPE(TwoPImpesTT, VtkOutputFields);

    VtkOutputModule<GVTwoPImpes, SVTwoPImpes> vtkWriter(*gridVarTwoPImpes, x_p, twoPImpesProblem->name(), "",Dune::VTK::conforming);
    VtkOutputFields::init(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);

    using VtkOutputFieldsTransport = typename GET_PROP_TYPE(TwoPTransportTT, VtkOutputFields);

    std::string nameTransport = twoPImpesProblem->name() + "_transport";
    VtkOutputModule<GVTwoPTransport, SVTwoPTransport> vtkWriterTransport(*gridVarTwoPTransport, x_s, nameTransport, "",Dune::VTK::conforming);
    VtkOutputFieldsTransport::init(vtkWriterTransport); //!< Add model specific output fields
    vtkWriterTransport.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    //! the assembler with time loop for instationary problem
    using TwoPImpesAssembler = FVAssembler<TwoPImpesTT, DiffMethod::analytic, true>;
    auto twoPImpesAssembler = std::make_shared<TwoPImpesAssembler>(twoPImpesProblem,
                                                                   twoPImpesFvGridGeometry,
                                                                   gridVarTwoPImpes,
                                                                   timeLoop);
    using JMTwoPImpes = typename GET_PROP_TYPE(TwoPImpesTT, JacobianMatrix);
    auto Ap = std::make_shared<JMTwoPImpes>();
    auto rp = std::make_shared<SVTwoPImpes>();
    twoPImpesAssembler->setLinearSystem(Ap, rp);


    using TwoPTransportAssembler = FVAssembler<TwoPTransportTT, DiffMethod::analytic, false>;
    auto twoPTransportAssembler = std::make_shared<TwoPTransportAssembler>(twoPTransportProblem,
                                                                           twoPTransportFvGridGeometry,
                                                                           gridVarTwoPTransport,
                                                                           timeLoop);
    using JMTwoPTransport = typename GET_PROP_TYPE(TwoPTransportTT, JacobianMatrix);
    auto As = std::make_shared<JMTwoPTransport>();
    auto rs = std::make_shared<SVTwoPTransport>();
    twoPTransportAssembler->setLinearSystem(As, rs);

    // the linear solver
    using LinearSolver = Dumux::ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    //! set some check points for the time loop
    timeLoop->setPeriodicCheckPoint(getParam<Scalar>("TimeLoop.OutputTimeInterval"));

    //! start the time loop
    timeLoop->start();
    while (!timeLoop->finished())
    {
        std::vector<double> satVals(x_s.size() + twoPTransportFvGridGeometry->numBoundaryScvf(),0.0);
        for (const auto& element : elements(leafGridView))
        {
            auto fvGeometry = localView(*twoPTransportFvGridGeometry);
            fvGeometry.bind(element);

            auto elemVolVars = localView(gridVarTwoPTransport->curGridVolVars());
            elemVolVars.bind(element, fvGeometry, x_s);

            const auto globalI = twoPTransportFvGridGeometry->elementMapper().index(element);
            auto&& scvI = fvGeometry.scv(globalI);
            const auto index = scvI.dofIndex();
            satVals[index] = elemVolVars[scvI].saturation(0);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                {
                    const auto bcTypes = twoPTransportProblem->boundaryTypes(element, scvf);
                    if (bcTypes.hasOnlyDirichlet())
                    {
                        const auto idx = scvf.outsideScvIdx();
                        satVals[idx] = elemVolVars[idx].saturation(0);
                    }
                }
            }
        }

        twoPImpesProblem->spatialParams().setWettingSaturation(satVals);

        // This additional update is required to guarantee that the correct
        // saturation values are stored in the pressure volVars.
        // This has only an effect if chaching is enabled.
        gridVarTwoPImpes->update(x_p);

        // set previous solution for storage evaluations
        twoPImpesAssembler->setPreviousSolution(x_p_old);

        twoPImpesAssembler->assembleJacobianAndResidual(x_p);

        SVTwoPImpes xpDelta(x_p);
        linearSolver->solve(*Ap, xpDelta, *rp);

        x_p -= xpDelta;
        gridVarTwoPImpes->update(x_p);

        // make the new solution the old solution
        x_p_old = x_p;
        gridVarTwoPImpes->advanceTimeStep();

        using Scalar =  typename GET_PROP_TYPE(TwoPImpesTT, Scalar);
        std::vector<Scalar> volumeFlux(twoPImpesFvGridGeometry->numScvf(), 0.0);

        using FluxVariables =  typename GET_PROP_TYPE(TwoPImpesTT, FluxVariables);
        auto upwindTerm = [](const auto& volVars) { return 1.0; };
        for (const auto& element : elements(leafGridView))
        {
            auto fvGeometry = localView(*twoPImpesFvGridGeometry);
            fvGeometry.bind(element);

            auto elemVolVars = localView(gridVarTwoPImpes->curGridVolVars());
            elemVolVars.bind(element, fvGeometry, x_p);

            auto elemFluxVars = localView(gridVarTwoPImpes->gridFluxVarsCache());
            elemFluxVars.bind(element, fvGeometry, elemVolVars);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto idx = scvf.index();

                if (!scvf.boundary())
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*twoPImpesProblem, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                    volumeFlux[idx] = fluxVars.advectiveFlux(0, upwindTerm);
                }
                else
                {
                    const auto bcTypes = twoPImpesProblem->boundaryTypes(element, scvf);
                    if (bcTypes.hasOnlyDirichlet())
                    {
                        FluxVariables fluxVars;
                        fluxVars.init(*twoPImpesProblem, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                        volumeFlux[idx] = fluxVars.advectiveFlux(0, upwindTerm);
                    }
                    else
                    {
                        //volumeFlux[idx] = -0.1;
                    }
                }
            }
        }

        twoPTransportProblem->spatialParams().setVolumeFlux(volumeFlux);

        //Now update saturation
        twoPTransportAssembler->setPreviousSolution(x_s_old);

        twoPTransportAssembler->assembleJacobianAndResidual(x_s);

        SVTwoPTransport xsDelta(x_s);
        linearSolver->solve(*As, xsDelta, *rs);

        x_s -= xsDelta;
        gridVarTwoPTransport->update(x_s);

        // make the new solution the old solution
        x_s_old = x_s;
        gridVarTwoPTransport->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // write vtk output on check points
        if (timeLoop->isCheckPoint() || timeLoop->finished() || timeLoop->timeStepIndex() == 1)
        {
            vtkWriter.write(timeLoop->time());
            vtkWriterTransport.write(timeLoop->time());
        }

        // set new dt
        timeLoop->setTimeStepSize(dt);

    }

    timeLoop->finalize(leafGridView.comm());

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
