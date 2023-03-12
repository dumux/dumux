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
 * \brief test for the pore network model
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/common/defaultusagemessage.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>
#include <dumux/porenetwork/2p/newtonsolver.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include "problem_dynamic.hh"


int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::DrainageProblem;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    using GridManager = Dumux::PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);

    // the spatial parameters
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometry);

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, spatialParams);

    // the solution vector
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(leafGridView.size(GridView::dimension));
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    const auto poresToPlot = getParam<std::vector<std::size_t>>("PlotMaterialLaw.PoresToPlot", std::vector<std::size_t>{});
    const auto sRange = getParam<std::array<double, 2>>("PlotMaterialLaw.SaturationRange", std::array<double, 2>{0, 1});
    spatialParams->plotPcSw<typename GridVariables::VolumeVariables>(poresToPlot, *problem, sRange[0], sRange[1]);

    // const auto throatsToPlot = getParam<std::vector<std::size_t>>("PlotMaterialLaw.ThroatsToPlot", std::vector<std::size_t>{});
    // spatialParams->plotTransmissibilities<typename GridVariables::VolumeVariables>(throatsToPlot, *problem, gridVariables->gridFluxVarsCache());

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // specify whether to drainage shall occur with step-wise applied pressures or using constant BCs
    const bool stepWiseDrainage = getParam<int>("Problem.NumSteps") > 1;

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))
        restartTime = getParam<Scalar>("TimeLoop.Restart");

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    PoreNetwork::VtkOutputModule<GridVariables, GetPropType<TypeTag, Properties::FluxVariables>, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); //! Add model specific output fields

    vtkWriter.write(0.0);

    Dumux::PoreNetwork::AveragedValues<GridVariables, SolutionVector> avgValues(*gridVariables, x);
    using FS = typename GridVariables::VolumeVariables::FluidSystem;
    avgValues.addAveragedQuantity([](const auto& v){ return v.saturation(FS::phase0Idx); }, [](const auto& v){ return v.poreVolume(); }, "avgSat");
    avgValues.addAveragedQuantity([](const auto& v){ return v.pressure(FS::phase0Idx); }, [](const auto& v){ return v.saturation(FS::phase0Idx)*v.poreVolume(); }, "avgPw");
    avgValues.addAveragedQuantity([](const auto& v){ return v.pressure(FS::phase1Idx); }, [](const auto& v){ return v.saturation(FS::phase1Idx)*v.poreVolume(); }, "avgPn");
    std::vector<std::size_t> dofsToNeglect;

    for (const auto& vertex : vertices(leafGridView))
    {
        using Labels = GetPropType<TypeTag, Properties::Labels>;
        const auto vIdx = gridGeometry->vertexMapper().index(vertex);
        if (gridGeometry->poreLabel(vIdx) == Labels::inlet || gridGeometry->poreLabel(vIdx) == Labels::outlet)
            dofsToNeglect.push_back(vIdx);
    }

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = PoreNetwork::TwoPNewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        // try solving the non-linear system
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // calculate the averaged values
        avgValues.eval(dofsToNeglect);
        problem->postTimeStep(timeLoop->time(), avgValues, gridVariables->gridFluxVarsCache().invasionState().numThroatsInvaded(), timeLoop->timeStepSize());

        // write vtk output
        if(problem->shouldWriteOutput(timeLoop->timeStepIndex(), *gridVariables))
            vtkWriter.write(timeLoop->time());

        // check if all drainge steps have been performed
        if(problem->simulationFinished())
            timeLoop->setFinished();

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        static const Scalar increasedTimeStepShiftThreshold = getParam<Scalar>("Problem.IncreasedTimeStepShiftThreshold", 1e-8);
        static const Scalar maxDtNonEquilibrium = getParam<Scalar>("TimeLoop.MaxDtNonEquilibrium", 1e-5);
        const Scalar newTimeStepSize = problem->dSwDt() > increasedTimeStepShiftThreshold ? maxDtNonEquilibrium : nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize());
        if (problem->dSwDt() > increasedTimeStepShiftThreshold)
        {
            std::cout << "Temporal change in saturation " << problem->dSwDt() << " is greater than threshold " << increasedTimeStepShiftThreshold << ". ";
            std::cout <<"Setting new time step to " << newTimeStepSize << std::endl;
        }
        else
        {
            std::cout << "Temporal change in saturation " << problem->dSwDt() << " is smaller than threshold " << increasedTimeStepShiftThreshold << ". ";
            std::cout <<"Setting new time step to " << newTimeStepSize << std::endl;
        }
        timeLoop->setTimeStepSize(newTimeStepSize);
        // timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    //plot the pc-S curve, if desired
#ifdef HAVE_GNUPLOT
    if(stepWiseDrainage)
        problem->plotPcS();
#endif


    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    timeLoop->finalize(leafGridView.comm());

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

}
