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
 * \brief Pipe flow test for the staggered grid RANS model
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */
#include <config.h>

#define IS_TURBULENT 1

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include "pipelauferproblem.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/methods.hh>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>

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
    using TypeTag = TTAG(PipeLauferProblem);

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);

    // try to create a grid (from the given grid file or the input file)
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    GridCreator::makeGrid();
    GridCreator::loadBalance();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = GridCreator::grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // get some time loop parameters
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))
        restartTime = getParam<Scalar>("TimeLoop.Restart");

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    problem->setTimeLoop(timeLoop);

    // the solution vector
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    const auto numDofsCellCenter = leafGridView.size(0);
    const auto numDofsFace = leafGridView.size(1);
    SolutionVector x;
    x[FVGridGeometry::cellCenterIdx()].resize(numDofsCellCenter);
    x[FVGridGeometry::faceIdx()].resize(numDofsFace);
    problem->applyInitialSolution(x);
    problem->updateStaticWallProperties();
    problem->updateDynamicWallProperties(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x, xOld);

    // intialize the vtk output module
    using VtkOutputFields = typename GET_PROP_TYPE(TypeTag, VtkOutputFields);
    StaggeredVtkOutputModule<TypeTag, GET_PROP_VALUE(TypeTag, PhaseIdx)> vtkWriter(*problem, *fvGridGeometry, *gridVariables, x, problem->name());
    VtkOutputFields::init(vtkWriter); //!< Add model specific output fields
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

    Dumux::GnuplotInterface<Scalar> gnuplot_lawOfTheWall;
    Dumux::GnuplotInterface<Scalar> gnuplot_velocityProfile;

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
            Dumux::GnuplotInterface<Scalar> gnuplot_lawOfTheWall;
            Dumux::GnuplotInterface<Scalar> gnuplot_velocityProfile;
            char gnuplotFileName[255];
            sprintf(gnuplotFileName, fileNameFormat.c_str(), "lawOfTheWall", timeLoop->timeStepIndex());
            gnuplot_lawOfTheWall.setOpenPlotWindow(plotLawOfTheWall);
            gnuplot_lawOfTheWall.setDatafileSeparator(',');
            gnuplot_lawOfTheWall.resetPlot();
            gnuplot_lawOfTheWall.setXlabel("y^+ [-]");
            gnuplot_lawOfTheWall.setYlabel("u_+ [-]");
            gnuplot_lawOfTheWall.setOption("set log x");
            gnuplot_lawOfTheWall.setOption("set xrange [1:3000]");
            gnuplot_lawOfTheWall.addFileToPlot("laufer_re50000_u+y+.csv", "u 1:2 w p t 'Laufer 1954, Re=50000'");
#if LOWREKEPSILON
            gnuplot_lawOfTheWall.addFileToPlot("pdelab-lowrekepsilon.csv", "u 23:22 w l lw 2 t 'PDELab Low-Re k-epsilon'");
            gnuplot_lawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 12:13 w l lc 7");
#elif KEPSILON
            gnuplot_lawOfTheWall.addFileToPlot("pdelab-kepsilon_twolayers.csv", "u 34:33 w l lw 2 t 'PDELab k-epsilon (two layers)'");
            gnuplot_lawOfTheWall.addFileToPlot("pdelab-kepsilon_wallfunction.csv", "u 36:35 w l lw 2 t 'PDELab k-epsilon (wall function)'");
            gnuplot_lawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 12:13 w l lc 7");
#elif KOMEGA
            gnuplot_lawOfTheWall.addFileToPlot("pdelab-komega.csv", "u 24:23 w l lw 2 t 'PDELab k-omega'");
            gnuplot_lawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 12:13 w l lc 7");
#else
            gnuplot_lawOfTheWall.addFileToPlot("pdelab-zeroeq.csv", "u 22:21 w l lw 2 t 'PDELab 0-Eq.'");
            gnuplot_lawOfTheWall.addFileToPlot(std::string(fileName) + ".csv", "u 12:13 w l lc 7");
#endif
            gnuplot_lawOfTheWall.plot(std::string(gnuplotFileName));

            sprintf(gnuplotFileName, fileNameFormat.c_str(), "velProfile", timeLoop->timeStepIndex());
            gnuplot_velocityProfile.setOpenPlotWindow(plotVelocityProfile);
            gnuplot_velocityProfile.setDatafileSeparator(',');
            gnuplot_velocityProfile.resetPlot();
            gnuplot_velocityProfile.setXlabel("v_x/v_{x,max} [-]");
            gnuplot_velocityProfile.setYRange(0.0, 1.0);
            gnuplot_velocityProfile.setYlabel("y [-]");
            gnuplot_velocityProfile.addFileToPlot("laufer_re50000.csv", "u 2:1 w p t 'Laufer 1954, Re=50000'");
#if LOWREKEPSILON
            gnuplot_velocityProfile.addFileToPlot("pdelab-lowrekepsilon.csv", "u 5:($27/0.2456) w l lw 2 t 'PDELab Low-Re k-epsilon'");
            gnuplot_velocityProfile.addFileToPlot(std::string(fileName) + ".csv", "u 7:($26/0.2456) w l lc 7");
#elif KEPSILON
            gnuplot_velocityProfile.addFileToPlot("pdelab-kepsilon_twolayers.csv", "u ($5):($40/0.2469) w l lw 2 t 'PDELab k-epsilon (two layers)'");
            gnuplot_velocityProfile.addFileToPlot("pdelab-kepsilon_wallfunction.csv", "u ($5):($40/0.2469) w l lw 2 t 'PDELab k-epsilon (wall function)'");
            gnuplot_velocityProfile.addFileToPlot(std::string(fileName) + ".csv", "u 7:($28/0.2456) w l lc 7");
#elif KOMEGA
            gnuplot_velocityProfile.addFileToPlot("pdelab-komega.csv", "u 5:($29/0.2456) w l lw 2 t 'PDELab k-omega'");
            gnuplot_velocityProfile.addFileToPlot(std::string(fileName) + ".csv", "u 7:($30/0.2456) w l lc 7");
#else
            gnuplot_velocityProfile.addFileToPlot("pdelab-zeroeq.csv", "u 5:($26/0.2456) w l lw 2 t 'PDELab 0-Eq.'");
            gnuplot_velocityProfile.addFileToPlot(std::string(fileName) + ".csv", "u 7:($24/0.2456) w l lc 7");
#endif
            gnuplot_velocityProfile.plot(std::string(gnuplotFileName));
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
