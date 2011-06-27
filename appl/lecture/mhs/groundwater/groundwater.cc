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

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/gridinfo.hh>

#include "external_interface.hh"
#include "groundwater_problem.hh"

////////////////////////
// the main function
////////////////////////
void usage(const char *progname)
{
    std::cout << boost::format("usage: %s #refine [delta]\n")%progname;
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

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
//        if (argc != 2 && argc != 3)
//            usage(argv[0]);
//
//        int numRefine;
//        std::istringstream(argv[1]) >> numRefine;



        ////////////////////////////////////////////////////////////
        // Read Input file
        ////////////////////////////////////////////////////////////

        Dumux::InterfaceProblemProperties interfaceProbProps("interface_groundwater.xml");
        Dumux::InterfaceFluidProperties interfaceFluidProps("interface_groundwater.xml");

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        GlobalPosition L(0.0);
        Grid grid(interfaceProbProps.resolution,L,interfaceProbProps.size);

        ////////////////////////////////////////////////////////////
        // adjust fluid properties
        ////////////////////////////////////////////////////////////
        typedef GET_PROP_TYPE(TypeTag, PTAG(Fluid)) Fluid;
        Fluid::Component::setViscosity(interfaceFluidProps.viscosity);
        Fluid::Component::setDensity(interfaceFluidProps.density);

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        Problem problem(grid.leafView(), interfaceProbProps);
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
