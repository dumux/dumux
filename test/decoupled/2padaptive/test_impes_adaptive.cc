/*****************************************************************************
 *   Copyright (C) 20010 by Markus Wolff                                     *
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
/*!
 * \file
 *
 * \ingroup IMPETtests
 * \brief test for the sequential 2p model
 */
#include "config.h"

#include "test_impes_adaptive_problem.hh"

#include <dune/grid/common/gridinfo.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>

#include <iostream>
#include <boost/format.hpp>


////////////////////////
// the main function
////////////////////////
void usage(const char *progname)
{
    std::cout << boost::format("usage: %s [--restart restartTime] InputFileName\n")%progname;
    exit(1);
}

int main(int argc, char** argv)
{
    try {
        typedef TTAG(TestIMPESAdaptiveProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;
        typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

        //todo: diese zwei Zeile nach dem testen entfernen
        typedef GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
        typedef typename GridView::Codim<0>::Iterator ElementLeafIterator;

//        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc < 2)
            usage(argv[0]);

        // deal with the restart stuff
        int argPos = 1;
        bool restart = false;
        double restartTime = 0;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> restartTime;
        }

        if (argc - argPos != 1) {
            usage(argv[0]);
        }


        std::string inputFileName;
        inputFileName = argv[argPos++];


        ////////////////////////////////////////////////////////////
        // Read Input file and create grid
        ////////////////////////////////////////////////////////////

        Dune::ParameterTreeParser::readINITree(inputFileName, Params::tree());

        std::string fileName = Params::tree().get<std::string>("gridFile");
        Dune::GridPtr<Grid> gridPtr(fileName);


        gridPtr->globalRefine(Params::tree().get<int>("levelInit"));

        // read the initial time step and the end time
        double tEnd, dt;
        tEnd = Params::tree().get<double>("tEnd");
        dt = tEnd;

//		std::cout << "size before adapt: " << gridPtr->leafView().size(0) << std::endl;
////		Einzelne Zellen verfeinern: links
//		for (ElementLeafIterator it=gridPtr->leafView().begin<0>();
//				it!=gridPtr->leafView().end<0>(); ++it)
//		{
//			if ((it->geometry().corner(0)[0]>=100))
//			{
//			gridPtr->mark(1,*it);
//			}
//		}
//
//		gridPtr->preAdapt();
//		gridPtr->adapt();
//		gridPtr->postAdapt();
//
//		std::cout << "size after first adapt: " << gridPtr->leafView().size(0) << std::endl;

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        TimeManager timeManager;
        Problem problem(timeManager, gridPtr->leafView());
        problem.setGrid(*gridPtr);
        problem.gridAdapt().setLevels(Params::tree().get<int>("levelMin"), Params::tree().get<int>("levelMax"));

        // load restart file if necessarry
        if (restart)
            problem.deserialize(restartTime);

        timeManager.init(problem, 0, dt, tEnd, !restart);
        timeManager.run();
        return 0;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }

    return 3;
}
