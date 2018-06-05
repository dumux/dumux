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
 *
 * \brief Backwards Facing Step Benchmark Problem for the staggered grid RANS model
 *
 *  This test simulates the experiments performed by
 *  Driver and Seegmiller in 1985 \cite DriverSeegmiller1985.
 *  More information available at https://turbmodels.larc.nasa.gov/backstep_val.html
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
#include <dumux/io/grid/subgridgridcreator.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "problem.hh"


/*!
 * \brief A method providing an () operator in order to select elements for a subgrid.
 */
template<class GlobalPosition>
class StepSelector
{
public:
    StepSelector(const GlobalPosition& step) : steptip_(step) {}

    //! Select all elements within a circle around a center point.
    int operator() (const auto& element) const
    {
        const auto x = element.geometry().center()[0];
        const auto y = element.geometry().center()[1];
        return ((x > steptip_[0]) || (y > steptip_[1]));
    }
private:
    const GlobalPosition steptip_;
};

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

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

   // The hostgrid
    constexpr int dim = 2;
    using GlobalPosition = Dune::FieldVector<double, dim>;
    using HostGrid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, dim> >;
    using HostGridManager = GridManager<HostGrid>;
    HostGridManager hostGridManager;
    hostGridManager.init();
    auto& hostGrid = hostGridManager.grid();


    const GlobalPosition step{getParam<double>("Grid.StepX");), getParam<double>("Problem.StepY")};
    StepSelector<GlobalPosition> subgridSelector(step);

    auto subGridPtr = SubgridGridCreator<HostGrid>::makeGrid(hostGrid, selector);

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = subGridPtr->leafGridView();

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
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
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

    // initialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVolumeVariable([](const auto& v){ return v.pressure() - 1.1e+5; }, "DeltaP");
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
        assembler->updateGridVariables(x);

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


 #if HAVE_PVPYTHON
// Plot Pre-Step Velocity Profile at 4H Before the Step
     static const bool plotPreStepVelocityProfile = getParam<bool>("Output.PlotPreStepVelocityProfile", false);
     if (plotPreStepVelocityProfile)
    {
        char fileName[255];
        std::string fileNameFormat = "%s-%05d";
        sprintf(fileName, fileNameFormat.c_str(), problem->name().c_str(), timeLoop->timeStepIndex());
        std::cout << fileName << std::endl;
        std::string vtuFileName = std::string(fileName) + ".vtu";
        std::string script = std::string(DUMUX_SOURCE_DIR) + "/bin/postprocessing/extractlinedata.py";
        std::string refPath = "../../../../test/references/BackwardsFacingStep_VelocityDistribution_Experiment.csv";
        std::string syscom;

        // execute the pvpython script
        std::string command = std::string(PVPYTHON_EXECUTABLE) + " " + script
                              + " -f " + vtuFileName
                              + " -v 2"
                              + " -r 10000";
        syscom =  command + " -p1 106 1.0 0.0"
                          + " -p2 106 9.0 0.0"
                          + " -of " + std::string(fileName)+ "_pre04 \n";
        if (!system(syscom.c_str()))
        {
            Dumux::GnuplotInterface<Scalar> gnuplot_preStepVelocityProfile;
            char gnuplotFileName[255];
            gnuplot_preStepVelocityProfile.resetAll();
            gnuplot_preStepVelocityProfile.setTerminalType("pngcairo size 1300,1000 solid");
            gnuplot_preStepVelocityProfile.setOption("set output 'VelDist_PreStep.png'");
            gnuplot_preStepVelocityProfile.setDatafileSeparator(',');
            gnuplot_preStepVelocityProfile.setXlabel("v_x/v_{max}");
            gnuplot_preStepVelocityProfile.setYlabel("y");
            gnuplot_preStepVelocityProfile.setOption("set xlabel font ',20'");
            gnuplot_preStepVelocityProfile.setOption("set ylabel font ',20'");
            gnuplot_preStepVelocityProfile.setOption("set xrange [0:1.01]");
            gnuplot_preStepVelocityProfile.setOption("set yrange [1:5]");
            gnuplot_preStepVelocityProfile.setOption("set grid ytics lc 9 lw 1 lt 0");
            gnuplot_preStepVelocityProfile.setOption("set grid xtics lc 9 lw 1 lt 0");
            gnuplot_preStepVelocityProfile.setOption("set key font ',20'");
            gnuplot_preStepVelocityProfile.setOption("set key left Left reverse center samplen 1");
            gnuplot_preStepVelocityProfile.addFileToPlot(std::string(fileName) + "_pre04.csv", "using 7:27 w l lc 7t 'Simulation x/H=-4'");
            gnuplot_preStepVelocityProfile.addFileToPlot(refPath, "using 3:2 w p ps 3 pt 3 lc 8 t 'Experiment x/H=-4'");
            gnuplot_preStepVelocityProfile.plot(std::string(gnuplotFileName));

        }
        else
        {
            std::cerr << "An error occurred when calling pvpython.";
        }
    }

// Plot the Post-Step Velocity Profile in 4 Locations
    static const bool plotPostStepVelocityProfile = getParam<bool>("Output.PlotPostStepVelocityProfile", false);
if (plotPostStepVelocityProfile)
    {
        char fileName[255];
        std::string fileNameFormat = "%s-%05d";
        sprintf(fileName, fileNameFormat.c_str(), problem->name().c_str(), timeLoop->timeStepIndex());
        std::string vtuFileName = std::string(fileName) + ".vtu";
        std::string script = std::string(DUMUX_SOURCE_DIR) + "/bin/postprocessing/extractlinedata.py";
        std::string refPath = "../../../../test/references/BackwardsFacingStep_VelocityDistribution_Experiment.csv";

        // execute the pvpython script
        std::string command = std::string(PVPYTHON_EXECUTABLE) + " " + script
                              + " -f " + vtuFileName
                              + " -v 2"
                              + " -r 10000";
        std::string syscom01 = command + " -p1 111 0.0 0.0"
                                       + " -p2 111 9.0 0.0"
                                       + " -of " + std::string(fileName) + "_post01 \n";
        std::string syscom04 = command + " -p1 114 0.0 0.0"
                                       + " -p2 114 9.0 0.0"
                                       + " -of " + std::string(fileName) + "_post04 \n";
        std::string syscom06 = command + " -p1 116 0.0 0.0"
                                       + " -p2 116 9.0 0.0"
                                       + " -of " + std::string(fileName) + "_post06 \n";
        std::string syscom10 = command + " -p1 120 0.0 0.0"
                                       + " -p2 120 9.0 0.0"
                                       + " -of " + std::string(fileName) + "_post10 \n";

        if (!system(syscom01.c_str()) && !system(syscom04.c_str()) && !system(syscom06.c_str()) && !system(syscom10.c_str()) )
        {
            Dumux::GnuplotInterface<Scalar> gnuplot_postStepVelocityProfile;
            char gnuplotFileName[255];
            gnuplot_postStepVelocityProfile.resetAll();
            gnuplot_postStepVelocityProfile.setTerminalType("pngcairo size 1300,1000 solid");
            gnuplot_postStepVelocityProfile.setOption("set output 'VelDist_PostStep.png'");
            gnuplot_postStepVelocityProfile.setDatafileSeparator(',');
            gnuplot_postStepVelocityProfile.setXlabel("v_x/v_{max}");
            gnuplot_postStepVelocityProfile.setYlabel("y");
            gnuplot_postStepVelocityProfile.setOption("set xlabel font ',20'");
            gnuplot_postStepVelocityProfile.setOption("set ylabel font ',20'");
            gnuplot_postStepVelocityProfile.setOption("set xrange [1:9]");
            gnuplot_postStepVelocityProfile.setOption("set yrange [0:2.5]");
            gnuplot_postStepVelocityProfile.setOption("set grid ytics lc 9 lw 1 lt 0");
            gnuplot_postStepVelocityProfile.setOption("set grid xtics lc 9 lw 1 lt 0");
            gnuplot_postStepVelocityProfile.setOption("set key font ',20'");
            gnuplot_postStepVelocityProfile.setOption("set arrow from 1,0 to 1,2.5 nohead");
            gnuplot_postStepVelocityProfile.setOption("set arrow from 3,-0.12 to 3,2.5 nohead");
            gnuplot_postStepVelocityProfile.setOption("set arrow from 5,-0.12 to 5,2.5 nohead");
            gnuplot_postStepVelocityProfile.setOption("set arrow from 7,-0.12 to 7,2.5 nohead");
            gnuplot_postStepVelocityProfile.setOption("set arrow from 9,0 to 9,2.5 nohead");
            gnuplot_postStepVelocityProfile.setOption("set xtics ('-1' 1, '0' 2, '1   -1' 3, '0' 4, '1   -1' 5, '0' 6, '1   -1' 7, '0' 8, '1' 9)");
            gnuplot_postStepVelocityProfile.setOption("set label 1 at  1.75,1 'x/H = 1' font ',20' center");
            gnuplot_postStepVelocityProfile.setOption("set label 2 at  3.75,1 'x/H = 4' font ',20' center");
            gnuplot_postStepVelocityProfile.setOption("set label 3 at  5.75,1 'x/H = 6' font ',20' center");
            gnuplot_postStepVelocityProfile.setOption("set label 4 at  7.75,1 'x/H = 10' font ',20' center");
            gnuplot_postStepVelocityProfile.setOption("set label 5 at 0.7,1.025 '(Step' font ',12' center");
            gnuplot_postStepVelocityProfile.setOption("set label 6 at 0.7,0.975 'Height)' font ',12' center");
            gnuplot_postStepVelocityProfile.setOption("set key at 2,2 box reverse opaque center");
            gnuplot_postStepVelocityProfile.setOption("set border back");
            gnuplot_postStepVelocityProfile.addFileToPlot(std::string(fileName) + "_post01.csv", "using 1:1 w l lc 8 t 'Simulation Results'");
            gnuplot_postStepVelocityProfile.addFileToPlot(std::string(fileName) + "_post01.csv", "using (2+$7):26 w l lc rgb '#0072bd' t ''");
            gnuplot_postStepVelocityProfile.addFileToPlot(std::string(fileName) + "_post04.csv", "using (4+$7):26 w l lc rgb '#d95319' t ''");
            gnuplot_postStepVelocityProfile.addFileToPlot(std::string(fileName) + "_post06.csv", "using (6+$7):26 w l lc rgb '#7e2f8e' t ''");
            gnuplot_postStepVelocityProfile.addFileToPlot(std::string(fileName) + "_post10.csv", "using (8+$7):26 w l lc rgb '#77ac30' t ''");
            gnuplot_postStepVelocityProfile.addFileToPlot(refPath, "using 1:1        w p ps 3 pt 3 lc 8 t 'Experimental Data'");
            gnuplot_postStepVelocityProfile.addFileToPlot(refPath, "using (2+$7):6   w p ps 3 pt 3 lc rgb '#0072bd' t ''");
            gnuplot_postStepVelocityProfile.addFileToPlot(refPath, "using (4+$11):10 w p ps 3 pt 3 lc rgb '#d95319' t ''");
            gnuplot_postStepVelocityProfile.addFileToPlot(refPath, "using (6+$15):14 w p ps 3 pt 3 lc rgb '#7e2f8e' t ''");
            gnuplot_postStepVelocityProfile.addFileToPlot(refPath, "using (8+$19):18 w p ps 3 pt 3 lc rgb '#77ac30' t ''");
            gnuplot_postStepVelocityProfile.plot(std::string(gnuplotFileName));
        }
        else
        {
            std::cerr << "An error occurred when calling pvpython.";
        }
    }

    // Plot Coefficient of Friction Along the Post-Step Base Wall
     static const bool plotBaseWallFrictionCoeff = getParam<bool>("Output.PlotBaseWallFrictionCoeff", false);
     if (plotBaseWallFrictionCoeff)
    {
        char fileName[255];
        std::string fileNameFormat = "%s-%05d";
        sprintf(fileName, fileNameFormat.c_str(), problem->name().c_str(), timeLoop->timeStepIndex());
        std::cout << fileName << std::endl;
        std::string vtuFileName = std::string(fileName) + ".vtu";
        std::string script = std::string(DUMUX_SOURCE_DIR) + "/bin/postprocessing/extractlinedata.py";
        std::string refPath = "../../../../test/references/BackwardsFacingStep_VelocityDistribution_Experiment.csv";
        std::string syscom;

        // execute the pvpython script
        std::string command = std::string(PVPYTHON_EXECUTABLE) + " " + script
                              + " -f " + vtuFileName
                              + " -v 2"
                              + " -r 10000";
        syscom =  command + " -p1 110 0.0 0.0"
                          + " -p2 150 0.0 0.0"
                          + " -of " + std::string(fileName) + "_BaseWall \n";

        if (!system(syscom.c_str()))
        {
            Dumux::GnuplotInterface<Scalar> gnuplot_baseWallCoeffFriction;
            char gnuplotFileName[255];
            gnuplot_baseWallCoeffFriction.resetAll();
            gnuplot_baseWallCoeffFriction.setTerminalType("pngcairo size 1300,1000 solid");
            gnuplot_baseWallCoeffFriction.setOption("set output 'FrictionCoeff_BaseWall.png'");
            gnuplot_baseWallCoeffFriction.setDatafileSeparator(',');
            gnuplot_baseWallCoeffFriction.setXlabel("x");
            gnuplot_baseWallCoeffFriction.setYlabel("Coefficient of Friction C_{f}");
            gnuplot_baseWallCoeffFriction.setOption("set xlabel font ',20'");
            gnuplot_baseWallCoeffFriction.setOption("set ylabel font ',20'");
            gnuplot_baseWallCoeffFriction.setOption("set xrange [0:20]");
            gnuplot_baseWallCoeffFriction.setOption("set yrange [-0.003:0.003]");
            gnuplot_baseWallCoeffFriction.setOption("set xtics 0,1,20");
            gnuplot_baseWallCoeffFriction.setOption("set grid ytics lc 9 lw 1 lt 0");
            gnuplot_baseWallCoeffFriction.setOption("set grid xtics lc 9 lw 1 lt 0");
            gnuplot_baseWallCoeffFriction.setOption("set arrow from 0,0 to 20,0 nohead lc 8 lw 2 lt 0");
            gnuplot_baseWallCoeffFriction.setOption("set arrow from 6.26,-0.003 to 6.26,0.003 nohead lc 8 lw 2 lt 0");
            gnuplot_baseWallCoeffFriction.setOption("set label 1 '(Reattachment Length)' at  6.5,-0.0015 font ',20' rotate by 90 center");
            gnuplot_baseWallCoeffFriction.setOption("set key font ',20'");
            gnuplot_baseWallCoeffFriction.setOption("set key at 3,0.002 box reverse opaque center");
            gnuplot_baseWallCoeffFriction.addFileToPlot(std::string(fileName) + "_BaseWall.csv", "using ($25-110):(($9)*($3+$16)/(0.5*($19/$7)*($19/$7))) w l lw 3 lt 1 lc 7 t 'Simulation C_{f} '");
            gnuplot_baseWallCoeffFriction.addFileToPlot(refPath, "using 21:22 w p ps 3 pt 3 lc 8 t 'Experimental C_{f} '");
            gnuplot_baseWallCoeffFriction.plot(std::string(gnuplotFileName));
        }
        else
        {
            std::cerr << "An error occurred when calling pvpython.";
        }
#endif
    }

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
