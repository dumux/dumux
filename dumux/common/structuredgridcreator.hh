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
#include <dune/grid/utility/structuredgridfactory.hh>
//#include <dumux/../util/structuredgridfactorysg.hh> //under discussion!


namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);
}

//! \cond INTERNAL
////////////////////////
// helper class for grid instantiation
////////////////////////


template <class TypeTag>
class CubeGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid)  Grid;
    typedef Dune::shared_ptr<Grid> GridPointer;

    enum
    {
        dim = Grid::dimension
    };

public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        Dune::array< unsigned int, dim > cellRes;
        Dune::FieldVector<Scalar, 2> upperRight;
        Dune::FieldVector<Scalar, 2> lowerLeft;

        lowerLeft[0] = 0.0;
        upperRight[0] = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.upperRightX);
        cellRes[0] = GET_RUNTIME_PARAM(TypeTag, int, Grid.numberOfCellsX);
        if (dim > 1)
        {
            lowerLeft[1] = 0.0;
            upperRight[1] = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.upperRightY);
            cellRes[1] = GET_RUNTIME_PARAM(TypeTag, int, Grid.numberOfCellsY);
        }
        if (dim > 2)
        {
            lowerLeft[2] = 0.0;
            upperRight[2] = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.upperRightZ);
            cellRes[2] = GET_RUNTIME_PARAM(TypeTag, int, Grid.numberOfCellsZ);
        }

        grid_ = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, cellRes);
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *grid_;
    };

    /*!
     * \brief Call loadBalance() function of Grid.
     */
    static void loadBalance()
    {
        grid_->loadBalance();
    };

private:
    static GridPointer grid_;
};

template <class TypeTag>
typename CubeGridCreator<TypeTag>::GridPointer CubeGridCreator<TypeTag>::grid_;

template <class TypeTag>
class SimplexGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid)  Grid;
    typedef Dune::shared_ptr<Grid> GridPointer;

    enum
    {
        dim = Grid::dimension
    };

public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        Dune::array< unsigned int, dim > cellRes;
        Dune::FieldVector<Scalar, 2> upperRight;
        Dune::FieldVector<Scalar, 2> lowerLeft;

        lowerLeft[0] = 0.0;
        upperRight[0] = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.upperRightX);
        cellRes[0] = GET_RUNTIME_PARAM(TypeTag, int, Grid.numberOfCellsX);
        if (dim > 1)
        {
            lowerLeft[1] = 0.0;
            upperRight[1] = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.upperRightY);
            cellRes[1] = GET_RUNTIME_PARAM(TypeTag, int, Grid.numberOfCellsY);
        }
        if (dim > 2)
        {
            lowerLeft[2] = 0.0;
            upperRight[2] = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.upperRightZ);
            cellRes[2] = GET_RUNTIME_PARAM(TypeTag, int, Grid.numberOfCellsZ);
        }

        grid_ = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lowerLeft, upperRight, cellRes);
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *grid_;
    };

    /*!
     * \brief Call loadBalance() function of Grid.
     */
    static void loadBalance()
    {
        grid_->loadBalance();
    };

private:
    static GridPointer grid_;
};

template <class TypeTag>
typename SimplexGridCreator<TypeTag>::GridPointer SimplexGridCreator<TypeTag>::grid_;

}

//! \endcond



