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

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progname  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
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
        << "Important optional options include:\n"
        << "\t--help,-h                        Print this usage message and exit\n"
        << "\t--print-parameters[=true|false]  Print the run-time modifiable parameters _after_ \n"
        << "\t                                 the simulation [default: true]\n"
        << "\t--print-properties[=true|false]  Print the compile-time parameters _before_ \n"
        << "\t                                 the simulation [default: true]\n"
        << "\t--parameter-file=FILENAME        File with parameter definitions\n"
        << "\t--restart=RESTARTTIME            Restart simulation from a restart file\n"
        << "\t--cells-x=NUM                    Number of cells in horizontal direction\n"
        << "\t--cells-y=NUM                    Number of cells in vertical direction\n"
        << "\t--lens-lower-left-x=VALUE        X-Coordinate of the lower-left corner\n"
        << "                                   of the low permeability lens\n"
        << "\t--lens-lower-left-y=VALUE        Y-Coordinate of the lower-left corner\n"
        << "                                   of the low permeability lens\n"
        << "\t--lens-upper-right-x=VALUE       X-Coordinate of the upper-right corner\n"
        << "                                   of the low permeability lens\n"
        << "\t--lens-upper-right-y=VALUE       Y-Coordinate of the upper-right corner\n"
        << "                                   of the low permeability lens\n"
        << "\n"
        << "All parameters can also be specified using the alternative syntax:\n"
        << "\t-tEnd ENDTIME                    The time of the end of the simlation [s]\n"
        << "\t-dtInitial STEPSIZE              The initial time step size [s]\n"
        << "\n"
        << "If --parameter-file is specified, parameters can also be defined there. In this case,\n"
        << "camel case is used for the parameters (e.g.: --grid-file becomes GridFile). Parameters\n"
        << "specified on the command line have priority over those in the parameter file.\n"
        << "\n"
        << "For the case of no arguments given, the input parameter file is expected to be named './parameter.input' \n"
        << "\n";
}


//! \cond INTERNAL
////////////////////////
// helper class for grid instantiation
////////////////////////
template <class TypeTag, class Grid = typename GET_PROP_TYPE(TypeTag, Grid)>
class LensGridCreator;

#if HAVE_UG
template <class TypeTag>
class LensGridCreator<TypeTag, Dune::UGGrid<2> >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::UGGrid<2> Grid;

public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        Dune::FieldVector<int, 2> cellRes;
        Dune::FieldVector<Scalar, 2> upperRight;
        Dune::FieldVector<Scalar, 2> lowerLeft;

        lowerLeft[0] = 0.0;
        lowerLeft[1] = 0.0;
        upperRight[0] = 6.0;
        upperRight[1] = 4.0;

        cellRes[0] = GET_RUNTIME_PARAM(TypeTag, int, CellsX);
        cellRes[1] = GET_RUNTIME_PARAM(TypeTag, int, CellsY);
        
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

        grid_ = factory.createGrid();
        grid_->loadBalance();
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *grid_;
    };

private:
    static Grid *grid_;
};

template <class TypeTag>
Dune::UGGrid<2> *LensGridCreator<TypeTag, Dune::UGGrid<2> >::grid_;
#endif

template <class TypeTag>
class LensGridCreator<TypeTag, Dune::YaspGrid<2> >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::YaspGrid<2> Grid;

public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        Dune::FieldVector<int, 2> cellRes;
        Dune::FieldVector<Scalar, 2> upperRight;
        Dune::FieldVector<Scalar, 2> lowerLeft;

        lowerLeft[0] = 0.0;
        lowerLeft[1] = 0.0;
        upperRight[0] = 6.0;
        upperRight[1] = 4.0;

        cellRes[0] = GET_RUNTIME_PARAM(TypeTag, int, CellsX);
        cellRes[1] = GET_RUNTIME_PARAM(TypeTag, int, CellsY);

        grid_ = new Dune::YaspGrid<2>(
#ifdef HAVE_MPI
            Dune::MPIHelper::getCommunicator(),
#endif
            upperRight, // upper right
            cellRes, // number of cells
            Dune::FieldVector<bool,2>(false), // periodic
            0); // overlap
    };

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *grid_;
    };

private:
    static Grid *grid_;
};

template <class TypeTag>
Dune::YaspGrid<2> *LensGridCreator<TypeTag, Dune::YaspGrid<2> >::grid_;

// set the GridCreator property
namespace Dumux {
namespace Properties {
SET_TYPE_PROP(LensProblem, GridCreator, LensGridCreator<TypeTag>);
}}

//! \endcond


////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    typedef TTAG(LensProblem) TypeTag;

    SET_RUNTIME_DEFAULT(TypeTag, int, CellsX, "48");
    SET_RUNTIME_DEFAULT(TypeTag, int, CellsY, "32");
    
    return Dumux::startWithParameters<TypeTag>(argc, argv, usage);
}
