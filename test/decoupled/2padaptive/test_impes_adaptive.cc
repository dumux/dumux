// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 20010 by Markus Wolff                                     *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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

#if !HAVE_UG
#warning "You need to have an UGGrid installed to run this test"

#include <iostream>

int main()
{
    std::cerr << "You need to have an UGGrid installed to run this test\n";
    return 1;
}
#else

#include "test_impes_adaptive_problem.hh"

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>

#include <iostream>


////////////////////////
// the main function
////////////////////////
void usage(const char *progname)
{
    std::cout << "usage: " << progname << " [--restart restartTime] InputFileName\n";
    exit(1);
}

int main(int argc, char** argv)
{
    try {
        typedef TTAG(TestIMPESAdaptiveProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;
        typedef typename GET_PROP(TypeTag, ParameterTree) Params;

        static const int dim = Grid::dimension;

        //todo: diese zwei Zeile nach dem testen entfernen
        typedef GET_PROP_TYPE(TypeTag, GridView) GridView;
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
        double startTime = 0;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> startTime;
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

        Dune::array< unsigned int, dim > numberOfCells;
        numberOfCells[0] = 2;
        numberOfCells[1] = 1;
        Dune::FieldVector<double, dim> lowerLeftCorner(0);
        Dune::FieldVector<double, dim> domainSize(300);
        domainSize[1] = 100;

        Dune::shared_ptr<Grid> grid(Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeftCorner, domainSize, numberOfCells));
        grid->setClosureType(Grid::ClosureType::NONE);

        grid->globalRefine(Params::tree().get<int>("MaxLevel"));

        // read the initial time step and the end time
        double tEnd, dt;
        tEnd = Params::tree().get<double>("tEnd");
        dt = tEnd;

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        TimeManager timeManager;
        Problem problem(timeManager, grid->leafView());
        problem.setGrid(*grid);

        timeManager.init(problem, startTime, dt, tEnd, restart);
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
#endif
