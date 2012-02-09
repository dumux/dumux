// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Melanie Darcis                               *
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
 * \brief test for the 2pni box model
 */
#include "config.h"
#include "injectionproblem2pni.hh"
#include <dumux/common/start.hh>


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
                    errorMessageOut += "\n\nThe List of Mandatory arguments for this program is:\n"
                                        "\t-tEnd                          The end of the simulation. [s] \n"
                                        "\t-dtInitial                     The initial timestep size. [s] \n"
                                        "\t-Grid.numberOfCellsX           Resolution in x-direction [-]\n"
                                        "\t-Grid.numberOfCellsY           Resolution in y-direction [-]\n"
                                        "\t-Grid.upperRightX              Dimension of the grid [m]\n"
                                        "\t-Grid.upperRightY              Dimension of the grid [m]\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    typedef TTAG(InjectionProblem2PNI) ProblemTypeTag;
    return Dumux::start<ProblemTypeTag>(argc, argv, usage);
}

//! \cond INTERNAL
////////////////////////
// helper class for grid instantiation
////////////////////////
template <class TypeTag, class Grid = typename GET_PROP_TYPE(TypeTag, Grid)>
class TwoPNiGridCreator;

#if HAVE_UG
template <class TypeTag>
class TwoPNiGridCreator<TypeTag, Dune::UGGrid<2> >
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
        upperRight[0] = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.upperRightX);
        upperRight[1] = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.upperRightY);

        cellRes[0] = GET_RUNTIME_PARAM(TypeTag, int, Grid.numberOfCellsX);
        cellRes[1] = GET_RUNTIME_PARAM(TypeTag, int, Grid.numberOfCellsY);

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
Dune::UGGrid<2> *TwoPNiGridCreator<TypeTag, Dune::UGGrid<2> >::grid_;
#endif

template <class TypeTag>
class TwoPNiGridCreator<TypeTag, Dune::YaspGrid<2> >
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
        upperRight[0] = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.upperRightX);
        upperRight[1] = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.upperRightY);

        cellRes[0] = GET_RUNTIME_PARAM(TypeTag, int, Grid.numberOfCellsX);
        cellRes[1] = GET_RUNTIME_PARAM(TypeTag, int, Grid.numberOfCellsY);

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
Dune::YaspGrid<2> *TwoPNiGridCreator<TypeTag, Dune::YaspGrid<2> >::grid_;

// set the GridCreator property
namespace Dumux {
namespace Properties {
SET_TYPE_PROP(InjectionProblem2PNI, GridCreator, TwoPNiGridCreator<TypeTag>);
}}

//! \endcond



