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
 * \brief Test for the shallow water model.
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dune/grid/uggrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/amgbackend.hh>
//#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/discretization/methods.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/hdf5reader.hh>
#include <dumux/io/xdmfwriter.hh>

#include "sweproblem.hh"

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
                    errorMessageOut += "\n\nThe list of mandatory options for this program is:\n"
                                        "\t-TimeManager.TEnd      End of the simulation [s] \n"
                                        "\t-TimeManager.DtInitial Initial timestep size [s] \n"
                                        "\t-Grid.File             Name of the file containing the grid \n"
                                        "\t                       definition in DGF format\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

/*!
 * \brief Converts a list of doubles into an vector
 *
 * \param string  string that should be splited
 * \return vector vector containing the time steps where to plot
 *
 */

std::vector<double> getDoublesFromString(std::string const & numberstring){

    std::istringstream myISS(numberstring);

    return std::vector<double>{
        std::istream_iterator<double>(myISS),
        std::istream_iterator<double>()
        };
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv) try
{

    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = TTAG(TYPETAG);

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);

    //Read a triangular grid from xdmf/HDF5
    typedef Dune::UGGrid<2> GridType;
    //typedef Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming> GridType;
    Dune::Dumux::HDF5Reader<GridType> reader("init.hdf5");
    reader.readGrid();
    reader.readAllCellData("/timesteps/start");
    reader.loadBalance();
    auto& elementdata = reader.getCellData();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    //xdmf-way:
    const auto& leafGridView  = reader.grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    auto spatialParams = std::make_shared<typename Problem::SpatialParams>(fvGridGeometry);
    auto problem = std::make_shared<Problem>(fvGridGeometry, spatialParams);

    // set the inital values
    problem->setInputData(elementdata);
    problem->spatialParams().setElementdata(elementdata);

    // the solution vector
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    SolutionVector x(fvGridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x, xOld);

    // get some time loop parameters
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDivisions = getParam<int>("TimeLoop.MaxTimeStepDivisions");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    std::string plottimesString =  getParam<std::string>("TimeLoop.Plottimes");
    std::vector<double> plottimes = getDoublesFromString(plottimesString);
    plottimes.push_back(tEnd + 1.0);


    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))
        restartTime = getParam<Scalar>("TimeLoop.Restart");

    /*// intialize the vtk output module
    using VtkOutputFields = typename GET_PROP_TYPE(TypeTag, VtkOutputFields);
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    VtkOutputFields::init(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);
    */

    //intialize xdmf output module
    Dune::Swf::XDMFWriter<typename GridType::LeafGridView>
        writer("result", leafGridView,"Dune-SWF implicit", 0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop);

    // the linear solver AMG works parallel and seriell (parallel AMG is slow for triangular grids)
    using LinearSolver = Dumux::AMGBackend<TypeTag>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, fvGridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // read in the boundary values
    problem->setBoundaryValues();


    // write a first output file with no data
    auto& plotMap = problem->xdmfGetVariable(x, *gridVariables, timeLoop->time());
    writer.beginTimeStep(0.0);
    writer.writeCellData(plotMap["h"],"h","m");
    writer.writeCellData(plotMap["u"],"u","m/s");
    writer.writeCellData(plotMap["v"],"v","m/s");
    writer.writeCellData(plotMap["z"],"z","m");
    writer.writeCellData(plotMap["theta"],"theta","m");
    writer.endTimeStep();


    bool doPlot = false;
    int timePos = 0;
    double timestepSuggestedBefore = timeLoop->timeStepSize();

    // time loop
    timeLoop->start(); do
    {
        //preprocessing
        problem->preTimeStep(x, *gridVariables, timeLoop->time(),timeLoop->timeStepSize());

        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        problem->postTimeStep(x, *gridVariables, timeLoop->timeStepSize());
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write output
        if (doPlot){
            plotMap = problem->xdmfGetVariable(x, *gridVariables, timeLoop->time());
            writer.beginTimeStep(timeLoop->time());
            writer.writeCellData(plotMap["h"],"h","m");
            writer.writeCellData(plotMap["u"],"u","m/s");
            writer.writeCellData(plotMap["v"],"v","m/s");
            writer.writeCellData(plotMap["z"],"z","m");
            writer.writeCellData(plotMap["theta"],"theta","m");
            writer.endTimeStep();
            doPlot = false;
        }

        // set new dt as suggested by newton controller and limited by printout definitions
        auto optTime = nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize());
        optTime = std::max(timestepSuggestedBefore,optTime);

        if ((std::max(plottimes[timePos]-timeLoop->time(),1.0E-3) < optTime)&&(plottimes[timePos]-timeLoop->time() > 1.0E-9))
        {
            timestepSuggestedBefore = optTime;
            optTime = std::max(plottimes[timePos]-timeLoop->time(),1.0E-3);
        }

        timeLoop->setTimeStepSize(optTime);

        //check if time is in plottimes
        if (timeLoop->time() >= plottimes[timePos]){
          timePos += 1;
          doPlot = true;
        }
        if (mpiHelper.rank() == 0){
            //simple output
            std::cout << "\n=====================================" << std::endl;
            std::cout << "time " << timeLoop->time() << " dt " << timeLoop->timeStepSize() << "\n" << std::endl;

            //show the fluxes over all boundaries
            problem->printBoundaryFluxes();
        }


    } while (!timeLoop->finished());

    // do the last plot after the timeloop
    plotMap = problem->xdmfGetVariable(x, *gridVariables, timeLoop->time());
    writer.beginTimeStep(timeLoop->time());
    writer.writeCellData(plotMap["h"],"h","m");
    writer.writeCellData(plotMap["u"],"u","m/s");
    writer.writeCellData(plotMap["v"],"v","m/s");
    writer.writeCellData(plotMap["z"],"z","m");
    writer.writeCellData(plotMap["theta"],"theta","m");
    writer.endTimeStep();

    // output some Newton statistics
    nonLinearSolver.report();
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
}

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
