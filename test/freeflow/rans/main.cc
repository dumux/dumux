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
 * \ingroup RANSTests
 * \brief Pipe flow test for the staggered grid RANS model,
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "problem.hh"

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
                    errorMessageOut += "\nPlease use the provided input files.\n";
        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::PipeLauferProblem;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    problem->setTimeLoop(timeLoop);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x[FVGridGeometry::cellCenterIdx()].resize(fvGridGeometry->numCellCenterDofs());
    x[FVGridGeometry::faceIdx()].resize(fvGridGeometry->numFaceDofs());
    problem->applyInitialSolution(x);
    problem->updateStaticWallProperties();
    problem->updateDynamicWallProperties(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop);

    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // update wall properties
        problem->updateDynamicWallProperties(x);

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

#if HAVE_PVPYTHON
    static const bool plotLawOfTheWall = getParam<bool>("Output.PlotLawOfTheWall", false);
    static const bool plotVelocityProfile = getParam<bool>("Output.PlotVelocityProfile", false);
    if (plotLawOfTheWall || plotVelocityProfile)
    {
        char fileName[255];
        std::string fileNameFormat = "%s-%05d";
        sprintf(fileName, fileNameFormat.c_str(), problem->name().c_str(), timeLoop->timeStepIndex());
        std::cout << fileName << std::endl;
        std::string vtuFileName = std::string(fileName) + ".vtu";
        std::string script = std::string(DUMUX_SOURCE_DIR) + "/bin/postprocessing/extractlinedata.py";
        std::string syscom;

        // execute the pvpython script
        std::string command = std::string(PVPYTHON_EXECUTABLE) + " " + script
                              + " -f " + vtuFileName
                              + " -v 0"
                              + " -r 10000";
        syscom =  command + " -p1 8.0 0.0 0.0"
                          + " -p2 8.0 0.2469 0.0"
                          + " -of " + std::string(fileName) + "\n";

        if (!system(syscom.c_str()))
        {
            Dumux::GnuplotInterface<Scalar> gnuplotLawOfTheWall;
            Dumux::GnuplotInterface<Scalar> gnuplotVelocityProfile;
            char gnuplotFileName[255];
            sprintf(gnuplotFileName, fileNameFormat.c_str(), "lawOfTheWall", timeLoop->timeStepIndex());
            gnuplotLawOfTheWall.setOpenPlotWindow(plotLawOfTheWall);
            gnuplotLawOfTheWall.setDatafileSeparator(',');
            gnuplotLawOfTheWall.resetPlot();
            gnuplotLawOfTheWall.setXlabel("y^+ [-]");
            gnuplotLawOfTheWall.setYlabel("u_+ [-]");
            gnuplotLawOfTheWall.setYRange(0.0, 30.0);
            gnuplotLawOfTheWall.setOption("set log x");
            gnuplotLawOfTheWall.setOption("set xrange [1:3000]");
            gnuplotLawOfTheWall.addFileToPlot("references/laufer_re50000_u+y+.csv", "u 1:2 w p t 'Laufer 1954, Re=50000'");
#if LOWREKEPSILON
            gnuplotLawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 11:12 w l lc 7");
#elif KEPSILON
            gnuplotLawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 11:12 w l lc 7 t 'with u_{tau}'");
            gnuplotLawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 15:16 w l lc 8 t 'with u_{tau,nom}'");
            gnuplotLawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 11:12 w l lc 7");
#elif KOMEGA
            gnuplotLawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 11:12 w l lc 7");
#elif ONEEQ
            gnuplotLawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 11:12 w l lc 7");
#else
            gnuplotLawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 11:12 w l lc 7");
#endif
            gnuplotLawOfTheWall.plot(std::string(gnuplotFileName));

            sprintf(gnuplotFileName, fileNameFormat.c_str(), "velProfile", timeLoop->timeStepIndex());
            gnuplotVelocityProfile.setOpenPlotWindow(plotVelocityProfile);
            gnuplotVelocityProfile.setDatafileSeparator(',');
            gnuplotVelocityProfile.resetPlot();
            gnuplotVelocityProfile.setXlabel("v_x/v_{x,max} [-]");
            gnuplotVelocityProfile.setYRange(0.0, 1.0);
            gnuplotVelocityProfile.setYlabel("y [-]");
            gnuplotVelocityProfile.addFileToPlot("references/laufer_re50000.csv", "u 2:1 w p t 'Laufer 1954, Re=50000'");
#if LOWREKEPSILON
            gnuplotVelocityProfile.addFileToPlot(std::string(fileName) + ".csv", "u 6:($25/0.2456) w l lc 7");
#elif KEPSILON
            gnuplotVelocityProfile.addFileToPlot(std::string(fileName) + ".csv", "u 6:($29/0.2456) w l lc 7");
#elif KOMEGA
            gnuplotVelocityProfile.addFileToPlot(std::string(fileName) + ".csv", "u 6:($25/0.2456) w l lc 7");
#elif ONEEQ
            gnuplotVelocityProfile.addFileToPlot(std::string(fileName) + ".csv", "u 6:($24/0.2456) w l lc 7");
#else
            gnuplotVelocityProfile.addFileToPlot(std::string(fileName) + ".csv", "u 6:($23/0.2456) w l lc 7");
#endif
            gnuplotVelocityProfile.plot(std::string(gnuplotFileName));
        }
        else
        {
            std::cerr << "An error occurred when calling pvpython.";
        }
    }
#endif

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
