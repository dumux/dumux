// $Id: test_impes.cc 4144 2010-08-24 10:10:47Z bernd $
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff, Klaus Mosthaf                       *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#include "config.h"

#include "external_interface.hh"
#include "buckleyleverettproblem.hh"

#include <dune/grid/common/gridinfo.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

#include <iomanip>

//#include "dumux/decoupled/2p/old_fractionalflow/variableclass2p.hh"
//#include "dumux/decoupled/2p/old_fractionalflow/define2pmodel.hh"
//#include "buckleyleverettproblem.hh"
//#include "impes_buckleyleverett_analytic.hh"
//#include "dumux/timedisc/expliciteulerstep.hh"
//#include "dumux/timedisc/timeloop.hh"

void usage(const char *progname)
{
    std::cout << boost::format("usage: %s tEnd\n")%progname;
    exit(1);
}


////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    try
    {
        typedef TTAG(BuckleyLeverettProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;

        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);


        //load interface-file
        Dumux::InterfaceProblemProperties interfaceProbProps("interface_BL.xml");
        double discretizationLength = interfaceProbProps.IPP_DiscretizationLength;

        // define the problem dimensions
        Dune::FieldVector<Scalar, dim> lowerLeft(0);
        Dune::FieldVector<Scalar, dim> upperRight(300);
        upperRight[1] = 75;

        int cellNumberX = static_cast<int>(300/discretizationLength);
        int cellNumberY = static_cast<int>(75/discretizationLength);

        Dune::FieldVector<int, dim> cellNumbers(cellNumberX);
        cellNumbers[1] = cellNumberY;
        if (argc != 2)
            usage(argv[0]);

        std::string arg1(argv[1]);
        std::istringstream is1(arg1);
        double tEnd;
        is1 >> tEnd;

        double dt = tEnd;

        // create a grid object
        Grid grid(cellNumbers, lowerLeft, upperRight);
        Dune::gridinfo(grid);

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////

        Dumux::InterfaceFluidProperties interfaceFluidProps("interface_BL.xml");
        typedef GET_PROP_TYPE(TypeTag, PTAG(WettingPhase)) WettingPhase;
        typedef GET_PROP_TYPE(TypeTag, PTAG(NonwettingPhase)) NonwettingPhase;

        WettingPhase::Component::setViscosity(interfaceFluidProps.IFP_ViscosityWettingFluid);
        NonwettingPhase::Component::setViscosity(interfaceFluidProps.IFP_ViscosityNonWettingFluid);

        WettingPhase::Component::setDensity(interfaceFluidProps.IFP_DensityWettingFluid);
        NonwettingPhase::Component::setDensity(interfaceFluidProps.IFP_DensityNonWettingFluid);

        Problem problem(grid.leafView(), lowerLeft, upperRight);

        problem.timeManager().init(problem, 0, dt, tEnd);
        problem.timeManager().run();
        return 0;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }

// OLD STUFF -> PROBLEM
//        // IMPES parameters
//        int iterFlag = 0;
//        int nIter = 30;
//        double maxDefect = 1e-5;
//
//        // plotting parameters
//        const char* fileName = "buckleyleverett";
//        int modulo = 1;
//
//        typedef Dumux::FVVelocity2P<GridView, Scalar, VariableType,
//                ProblemType> DiffusionType;
//        DiffusionType diffusion(gridView, problem, modelDef);
//
//        typedef Dumux::FVSaturation2P<GridView, Scalar, VariableType,
//                ProblemType> TransportType;
//        TransportType transport(gridView, problem, modelDef);
//
//        typedef Dune::IMPESBLAnalytic<GridView, DiffusionType, TransportType,
//                VariableType> IMPESType;
//        IMPESType impes(diffusion, transport, iterFlag, nIter, maxDefect);
//        Dumux::ExplicitEulerStep<GridType, IMPESType> timestep;
//        Dumux::TimeLoop<GridView, IMPESType> timeloop(gridView, tStart, tEnd, fileName,
//                modulo, cFLFactor, maxDt, firstDt, timestep);
//        Dune::Timer timer;
//        timer.reset();
//        timeloop.execute(impes, false);
//        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds"
//                << std::endl;
//        //printvector(std::cout, *fractionalflow, "saturation", "row", 200, 1);
}
