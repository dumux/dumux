// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \brief test for the decoupled one-phase model.
 */
#include "config.h"
#include <iostream>
#include <boost/format.hpp>
#include <dune/common/parametertreeparser.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/gridinfo.hh>

#include "groundwater_problem.hh"

////////////////////////
// the main function
////////////////////////
void usage(const char *progname)
{
    std::cout << boost::format("usage: %s #inputFileName\n")%progname;
    exit(1);
}

int main(int argc, char** argv)
{
    try {
        typedef TTAG(GroundwaterProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        static const int dim = Grid::dimension;
        typedef Dune::FieldVector<Scalar, dim> GlobalPosition;
        typedef GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc != 2)
            usage(argv[0]);

        std::string inputFileName;
        inputFileName = argv[1];

        ////////////////////////////////////////////////////////////
        // Read Input file
        ////////////////////////////////////////////////////////////

        Dune::ParameterTreeParser::readINITree(inputFileName, Params::tree());

        Dune::FieldVector<int,2> N = Params::tree().get<Dune::FieldVector<int,2> >("Geometry.numberOfCells");
        GlobalPosition H = Params::tree().get<GlobalPosition>("Geometry.domainSize");

        if (N[0]*N[1]>10000)
        {
        	std::cout << "Number of cells exceeds limit. Choose less than 10000 cells!\n";
        	exit(2);
        }

        if (N[0]<2 || N[1]<2)
        {
			std::cout << "Number of cells too low. Choose at least 2 cell in each direction!\n";
			exit(2);

        }

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        GlobalPosition L(0.0);
        Grid grid(N,L,H);

        ////////////////////////////////////////////////////////////
        // adjust fluid properties
        ////////////////////////////////////////////////////////////
        typedef GET_PROP_TYPE(TypeTag, PTAG(Fluid)) Fluid;
		Fluid::Component::setViscosity(0.001);
		Fluid::Component::setDensity(1000);

//      Fluid::Component::setViscosity(inputParameters.get<double>("Fluid.viscosity"));
//		Fluid::Component::setDensity(inputParameters.get<double>("Fluid.density"));

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////

        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        Problem problem(grid.leafView());
        problem.init();
        problem.writeOutput();

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
