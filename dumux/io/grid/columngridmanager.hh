// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup InputOutput
 * \brief Provides a grid manager for a piece of a column grid
 */
#ifndef DUMUX_COLUMN_GRID_MANAGER_HH
#define DUMUX_COLUMN_GRID_MANAGER_HH

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <dune/common/dynvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/math.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for creating a two-dimensional
 *        grid with a horizontal and a vertical column.
 */
template <class Grid>
class ColumnGridManager
{
    using Scalar = typename Grid::ctype;

    using GridFactory = Dune::GridFactory<Grid>;
    using GridPointer = std::shared_ptr<Grid>;

    enum { dim = Grid::dimension,
           dimWorld = Grid::dimensionworld };

public:
    /*!
     * \brief Make the grid.
     */
    void init(const std::string& modelParamGroup = "")
    {
        static_assert(dim == 2, "The ColumnGridManager is only implemented for 2D.");

        const bool verbose = getParamFromGroup<bool>(modelParamGroup, "Grid.Verbosity", false);

        gridPtr() = createColumnGrid(modelParamGroup, verbose);

        loadBalance();
    }

    /*!
     * \brief Creates a grid consisting of a horizontal and a vertical column.
     *
     * \param modelParamGroup name of the model parameter group
     * \param verbose if the output should be verbose
     */
    std::unique_ptr<Grid> createColumnGrid(const std::string& modelParamGroup,
                                           bool verbose = false)
    {
        const auto columnX = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.HorizontalColumn");
        const auto columnY = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.VerticalColumn");
        const auto cells = getParamFromGroup<std::vector<unsigned>>(modelParamGroup, "Grid.Cells");

        // grid size in each coordinate direction
        const auto hX = (columnX[1] - columnX[0])/cells[0];
        const auto hY = (columnY[1] - columnY[0])/cells[1];

        // get the x-coordinate of the column intersection, by default the start of the horizontal column
        auto intersectX = getParamFromGroup<Scalar>(modelParamGroup, "Grid.IntersectionX", columnX[0]);

        // number of cells before the vertical column
        const unsigned intersectIdx = (intersectX - columnX[0])/hX;

        // evtl. adapt the intersection coordinate such that it coincides with a vertex location
        intersectX = columnX[0] + intersectIdx*hX;

        GridFactory gridFactory;
        constexpr auto type = Dune::GeometryTypes::cube(dim);

        // insert vertices for the horizontal column
        for (auto i = 0u; i <= cells[0]; ++i)
        {
            Dune::FieldVector <double, dim> v(0.0);

            // determine x-coordinate
            v[0] = columnX[0] + i*hX;

            // lower vertex in horizontal column
            v[1] = columnY[0] - hY;
            if(verbose)
                printCoordinate(v);
            gridFactory.insertVertex(v);

            // upper vertex in horizontal column
            v[1] = columnY[0];
            if(verbose)
                printCoordinate(v);
            gridFactory.insertVertex(v);
        }

        // insert vertices for the vertical column;
        // the bottom two vertices have been inserted
        // already above as part of the horizontal column
        for (auto i = 1u; i <= cells[1]; ++i)
        {
            Dune::FieldVector <double, dim> v(0.0);

            // determine y-coordinate
            v[1] = columnY[0] + i*hY;

            // left vertex in vertical column
            v[0] = intersectX;
            if(verbose)
                printCoordinate(v);
            gridFactory.insertVertex(v);

            // right vertex in vertical column
            v[0] = intersectX + hX;
            if(verbose)
                printCoordinate(v);
            gridFactory.insertVertex(v);
        }

        // insert elements in horizontal column
        for (auto i = 0u; i < cells[0]; ++i)
        {
            const std::vector<unsigned> vid({2*i, 2*i+2, 2*i+1, 2*i+3});
            if (verbose)
                printIndices(vid);

            gridFactory.insertElement(type, vid);
        }

        // the indices of the lowest element of the vertical column
        const auto offset = 2*cells[0];
        const std::vector<unsigned> vid({2*intersectIdx+1, 2*intersectIdx+3, offset+2, offset+3});
        if (verbose)
            printIndices(vid);
        gridFactory.insertElement(type, vid);

        // insert remaining elements in vertical column
        for (auto i = 1u; i < cells[1]; ++i)
        {
            const std::vector<unsigned> vid({offset+2*i, offset+2*i+1, offset+2*i+2, offset+2*i+3});
            if (verbose)
                printIndices(vid);

            gridFactory.insertElement(type, vid);
        }

        // return the grid pointer
        return std::unique_ptr<Grid>(gridFactory.createGrid());
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    {
        return *gridPtr();
    }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     */
    void loadBalance()
    {
        gridPtr()->loadBalance();
    }

protected:
    static void printCoordinate(const Dune::FieldVector <double, dim>& v)
    {
        std::cout << "Coordinates of : ";
        for (int k = 0; k < v.size(); ++k)
            std::cout << v[k] << " ";
        std::cout << std::endl;
    }

    static void printIndices(const std::vector<unsigned>& vid)
    {
        std::cout << "element vertex indices: ";
        for (int k = 0; k < vid.size(); ++k)
            std::cout << vid[k] << " ";
        std::cout << std::endl;
    }

    /*!
     * \brief Returns a reference to the shared pointer to the grid.
     */
    GridPointer& gridPtr()
    {
        return columnGrid_;
    }

private:
    GridPointer columnGrid_;
};

} // end namespace Dumux

#endif
