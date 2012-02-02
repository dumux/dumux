// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
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
 * \brief test for the two-phase box model
 */
#include "config.h"
#include "lensproblem.hh"
#include <dumux/common/start.hh>

#include <dune/grid/common/gridinfo.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/uggrid.hh>

#include <iostream>

void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::cout << errorMsg << "\n"
                  << "\n";
    }
    std::cout
        << "Usage: " << progName << " [options]\n"
        << "Mandatory options are:\n"
        << "\t--t-end=ENDTIME                  The time of the end of the simlation [s]\n"
        << "\t--dt-initial=STEPSIZE            The initial time step size [s]\n"
        << "\n"
        << "Alternativ supported syntax:\n"
        << "\t-tEnd ENDTIME                    The time of the end of the simlation [s]\n"
        << "\t-dtInitial STEPSIZE              The initial time step size [s]\n"
        << "\n"
        << "If --parameter-file is specified parameters can also be defined there. In this case,\n"
        << "camel case is used for the parameters (e.g.: --grid-file becomes gridFile). Parameters\n"
        << "specified on the command line have priority over those in the parameter file.\n"
        << "Important optional options include:\n"
        << "\t--help,-h                        Print this usage message and exit\n"
        << "\t--print-parameters[=true|false]  Print the run-time modifiable parameters _after_ \n"
        << "\t                                 the simulation [default: true]\n"
        << "\t--print-properties[=true|false]  Print the compile-time parameters _before_ \n"
        << "\t                                 the simulation [default: true]\n"
        << "\t--parameter-file=FILENAME        File with parameter definitions\n"
        << "\t--restart=RESTARTTIME            Restart simulation from a restart file\n"
        << "\n"
        << "For the case of no arguments given, the input parameter file is expected to be named './parameter.input' \n"
        << "\n";
}


//! \cond INTERNAL
////////////////////////
// helper class for grid instantiation
////////////////////////
template <class Grid, class Scalar>
class CreateGrid
{
};

#if HAVE_UG
template <class Scalar>
class CreateGrid<Dune::UGGrid<2>, Scalar>
{
public:
    static Dune::UGGrid<2> *create(const Dune::FieldVector<Scalar, 2> &upperRight,
                                   const Dune::FieldVector<int, 2> &cellRes)
    {
        Dune::GridFactory<Dune::UGGrid<2> > factory;
        for (int i=0; i<=cellRes[0]; i++) {
            for (int j=0; j<=cellRes[1]; j++) {
                Dune::FieldVector<double,2> pos;
                pos[0] = upperRight[0]*double(i)/cellRes[0];
                pos[1] = upperRight[1]*double(j)/cellRes[1];
                factory.insertVertex(pos);
            }
        }

        for (int i=0; i<cellRes[0]; i++) {
            for (int j=0; j<cellRes[1]; j++) {
#if CUBES
                std::vector<unsigned int> v(4);
#else
                std::vector<unsigned int> v(3);
#endif

                int i0 = i*(cellRes[1]+1) + j;
                int i1 = i*(cellRes[1]+1) + j+1;
                int i2 = (i+1)*(cellRes[1]+1) + j;
                int i3 = (i+1)*(cellRes[1]+1) + j+1;

#if CUBES
                v[0] = i0;
                v[1] = i1;
                v[2] = i2;
                v[3] = i3;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), v);
#else
                v[0] = i0;
                v[1] = i1;
                v[2] = i2;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), v);

                v[0] = i1;
                v[1] = i2;
                v[2] = i3;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), v);
#endif
            }
        }

        Dune::UGGrid<2> *grid = factory.createGrid();
        grid->loadBalance();
        return grid;
    }
};
#endif

template <class Scalar>
class CreateGrid<Dune::YaspGrid<2>, Scalar>
{
public:
    static Dune::YaspGrid<2> *create(const Dune::FieldVector<Scalar, 2> &upperRight,
                                     const Dune::FieldVector<int, 2> &cellRes)
    {
        return new Dune::YaspGrid<2>(
#ifdef HAVE_MPI
            Dune::MPIHelper::getCommunicator(),
#endif
            upperRight, // upper right
            cellRes, // number of cells
            Dune::FieldVector<bool,2>(false), // periodic
            0); // overlap
    };
};

template <class Scalar>
class CreateGrid<Dune::SGrid<2, 2>, Scalar>
{
public:
    static Dune::SGrid<2, 2> *create(const Dune::FieldVector<Scalar, 2> &upperRight,
                                     const Dune::FieldVector<int, 2> &cellRes)
    {
        return new Dune::SGrid<2,2>(cellRes, // number of cells
                                    Dune::FieldVector<Scalar, 2>(0.0), // lower left
                                    upperRight); // upper right

    };
};
//! \endcond






////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
#ifdef NDEBUG
    try {
#endif
        typedef TTAG(LensProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;

        static const int dim = Grid::dimension;

        // print all properties
        Dumux::Properties::print<TypeTag>();

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        Scalar tEnd, dt;
        Scalar restartTime=0;
        bool restart = false;
        Dumux::startWithParametersProvideMyOwnGrid<TypeTag, Scalar>(argc,
                                                                    argv,
                                                                    usage,
                                                                    tEnd,
                                                                    dt,
                                                                    restart,
                                                                    restartTime);

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        GlobalPosition lowerLeft(0.0);
        GlobalPosition upperRight;
        Dune::FieldVector<int,dim> res; // cell resolution
        upperRight[0] = 6.0;
        upperRight[1] = 4.0;

        /*
        res[0] = 48*4;
        res[1] = 32*4;
        */

        /*
        res[0] = 96;
        res[1] = 64;
        */

        res[0] = 48;
        res[1] = 32;

        /*
        res[0] = 24;
        res[1] = 16;
        */

        /*
        res[0] = 6;
        res[1] = 4;
        */

        /*
        res[0] = 1;
        res[1] = 2;
        */

        std::auto_ptr<Grid> grid(CreateGrid<Grid, Scalar>::create(upperRight, res));
        //grid->globalRefine(2);

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////

        // specify dimensions of the low-permeable lens
        GlobalPosition lowerLeftLens, upperRightLens;
        lowerLeftLens[0] = 1.0;
        lowerLeftLens[1] = 2.0;
        upperRightLens[0] = 4.0;
        upperRightLens[1] = 3.0;

        // instantiate and run the concrete problem
        TimeManager timeManager;
        Problem problem(timeManager, grid->leafView(), lowerLeftLens, upperRightLens);
        timeManager.init(problem, restartTime, dt, tEnd, restart);
        timeManager.run();

        return 0;

#ifdef NDEBUG
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }
#endif // NDEBUG

    return 3;
}
