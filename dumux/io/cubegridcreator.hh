// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Provides a grid creator which a regular grid made of
 *        quadrilaterals.
 */
#ifndef DUMUX_CUBE_GRID_CREATOR_HH
#define DUMUX_CUBE_GRID_CREATOR_HH

#include <dune/grid/utility/structuredgridfactory.hh>

#include <dumux/common/basicproperties.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);
}

/*!
 * \brief Provides a grid creator which a regular grid made of
 *        quadrilaterals.
 *
 * A quadirlateral is a line segment in 1D, a rectangle in 2D and a
 * cube in 3D.
 */
template <class TypeTag>
class CubeGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid)  Grid;
    typedef Dune::shared_ptr<Grid> GridPointer;

    enum { dim = Grid::dimension };

public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        Dune::array< unsigned int, dim > cellRes;
        Dune::FieldVector<Scalar, dim> upperRight;
        Dune::FieldVector<Scalar, dim> lowerLeft;

        lowerLeft[0] = 0.0;
        upperRight[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, UpperRightX);
        cellRes[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, NumberOfCellsX);
        if (dim > 1)
        {
            lowerLeft[1] = 0.0;
            upperRight[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, UpperRightY);
            cellRes[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, NumberOfCellsY);
        }
        if (dim > 2)
        {
            lowerLeft[2] = 0.0;
            upperRight[2] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, UpperRightZ);
            cellRes[2] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, NumberOfCellsZ);
        }

        gridPtr() = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, cellRes);
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *gridPtr();
    }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     */
    static void loadBalance()
    {
        gridPtr()->loadBalance();
    }

    /*!
     * \brief Returns a reference to the shared pointer to the grid.
     */
    static GridPointer &gridPtr()
    {
        static GridPointer cubeGrid;
        return cubeGrid;
    }
};

}

#endif
