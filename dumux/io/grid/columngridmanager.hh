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
 * \brief Provides a grid manager with a method for creating vectors
 *        with coordinates for columns.
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

        std::array<std::vector<Scalar>, dim> coordinates;

        createVectors(coordinates, modelParamGroup, verbose);
        gridPtr() = createColumnGrid(coordinates, modelParamGroup, verbose);

        loadBalance();
    }

    /*!
     * \brief Create vectors containing coordinates of all points.
     *
     * All keys are expected to be in group GridParameterGroup.
     * The following keys are recognized:
     * - Cells : number of cells array for x-coordinate (Again, an added 0, 1 or 2 specifies x, y or z)
     * - Grading : grading factor array for x-coordinate (Same here)
     * - Verbosity : whether the grid construction should output to standard out
     *
     * The grading factor \f$ g \f$ specifies the ratio between the next and the current cell size:
     * \f$ g = \frac{h_{i+1}}{h_i} \f$.
     * Negative grading factors are converted to
     * \f$ g = -\frac{1}{g_\textrm{negative}} \f$
     * to avoid issues with imprecise fraction numbers.
     */
    static void createVectors(std::array<std::vector<Scalar>, dim> &coordinates,
                              const std::string& modelParamGroup,
                              bool verbose = false)
    {
        // The positions
        std::array<std::vector<Scalar>, dim> positions;
        for (int i = 0; i < dim; ++i)
        {
            positions[i] = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.Pos" + std::to_string(i));

            if (!std::is_sorted(positions[i].begin(), positions[i].end()))
                DUNE_THROW(Dune::GridError, "Make sure to specify a monotone increasing \"Positions\" array");

            if(positions[i].size() < 2)
                DUNE_THROW(Dune::GridError, "Make sure to specify position arrays with at least two entries (min and max value).");
        }

        // the number of cells (has a default)
        std::array<std::vector<int>, dim> cells;
        for (int i = 0; i < dim; ++i)
        {
            cells[i].resize(positions[i].size()-1, 1.0);
            cells[i] = getParamFromGroup<std::vector<int>>(modelParamGroup, "Grid.Cells" +  std::to_string(i), cells[i]);
            if (cells[i].size() + 1 != positions[i].size())
                DUNE_THROW(Dune::RangeError, "Make sure to specify the correct length of \"Cells\" and \"Positions\" arrays");

        }
        std::cout << " positions" << std::endl;

        // grading factor (has a default)
        std::array<std::vector<Scalar>, dim> grading;
        for (int i = 0; i < dim; ++i)
        {
            grading[i].resize(positions[i].size()-1, 1.0);
            grading[i] = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.Grading" +  std::to_string(i), grading[i]);
            if (grading[i].size() + 1 != positions[i].size())
                DUNE_THROW(Dune::RangeError, "Make sure to specify the correct length of \"Grading\" and \"Positions\" arrays");
        }

        // make the grid
        std::array<std::vector<Scalar>, dim> globalPositions;
        using std::pow;
        using Dune::power;
        for (int dimIdx = 0; dimIdx < dim; dimIdx++)
        {
            // Each grid direction is subdivided into (numCells + 1) points
            std::size_t numGlobalPositions = 1;
            for (int zoneIdx = 0; zoneIdx < cells[dimIdx].size(); ++zoneIdx)
                numGlobalPositions += cells[dimIdx][zoneIdx];

            globalPositions[dimIdx].resize(numGlobalPositions);
            std::size_t posIdx = 0;
            for (int zoneIdx = 0; zoneIdx < cells[dimIdx].size(); ++zoneIdx)
            {
                const Scalar lower = positions[dimIdx][zoneIdx];
                const Scalar upper = positions[dimIdx][zoneIdx+1];
                const int numCells = cells[dimIdx][zoneIdx];
                Scalar gradingFactor = grading[dimIdx][zoneIdx];
                const Scalar length = upper - lower;
                Scalar height = 1.0;
                bool increasingCellSize = false;

                if (verbose)
                {
                    std::cout << "dim " << dimIdx
                              << " lower "  << lower
                              << " upper "  << upper
                              << " numCells "  << numCells
                              << " grading "  << gradingFactor;
                }

                if (gradingFactor > 1.0)
                    increasingCellSize = true;

                // take absolute values and reverse cell size increment to achieve
                // reverse behavior for negative values
                if (gradingFactor < 0.0)
                {
                    using std::abs;
                    gradingFactor = abs(gradingFactor);

                    if (gradingFactor < 1.0)
                        increasingCellSize = true;
                }

                const bool useGrading = Dune::FloatCmp::eq(gradingFactor, 1.0) ? false : true;

                // if the grading factor is exactly 1.0 do equal spacing
                if (!useGrading)
                {
                    height = 1.0 / numCells;
                    if (verbose)
                        std::cout << " -> h "  << height * length << std::endl;
                }
                // if grading factor is not 1.0, do power law spacing
                else
                {
                    height = (1.0 - gradingFactor) / (1.0 - power(gradingFactor, numCells));

                    if (verbose)
                    {
                        std::cout << " -> grading_eff "  << gradingFactor
                                  << " h_min "  << height * power(gradingFactor, 0) * length
                                  << " h_max "  << height * power(gradingFactor, numCells-1) * length
                                  << std::endl;
                    }
                }

                // the positions inside a global segment
                Dune::DynamicVector<Scalar> localPositions(numCells, 0.0);
                for (int i = 1; i < numCells; i++)
                {
                    Scalar hI = height;
                    if (useGrading)
                    {
                        if (increasingCellSize)
                            hI *= power(gradingFactor, i-1);

                        else
                            hI *= power(gradingFactor, numCells-i);
                    }
                    localPositions[i] = localPositions[i-1] + hI;
                }

                localPositions *= length;
                localPositions += lower;

                for (int i = 0; i < numCells; ++i)
                    globalPositions[dimIdx][posIdx++] = localPositions[i];
            }

            globalPositions[dimIdx][posIdx] = positions[dimIdx].back();
        }
        for (int dimIdx = 0; dimIdx < dim; dimIdx++)
        {
            coordinates[dimIdx] = globalPositions[dimIdx];
        }

    }

    /*!
     * \brief Creates Cartesian grid from coordinates.
     *
     * \param coordinates Vector containing coordinates
     * \param modelParamGroup name of the model parameter group
     * \param verbose if the output should be verbose
     */
    std::unique_ptr<Grid> createColumnGrid(std::array<std::vector<Scalar>, dim> &coordinates,
                                         const std::string& modelParamGroup,
                                         bool verbose = false)
    {
        std::vector<Scalar> x = coordinates[0];
        std::vector<Scalar> y = coordinates[1];

        int maxy = y.size() - 1;

        GridFactory gridFactory;
        constexpr auto type = Dune::GeometryTypes::cube(dim);

        for (int j = 0; j <= maxy; ++j)
        {
            for (int i = 0; i <= x.size()- 1; ++i)
            {
                // transform into Cartesian coordinates
                Dune::FieldVector <double, dim> v(0.0);

                v[0] = x[i];
                v[1] = y[j];
                if(verbose)
                    printCoordinate(v);
                gridFactory.insertVertex(v);
            }
        }

        std::cout << "Filled node vector" << std::endl;

        // assign nodes
        unsigned int z = 0;
        unsigned int t = 0;
        unsigned int xSize = x.size();
        for (int j = 0; j < y.size() - 1; ++j)
        {
            if (j < maxy)
            {
                for (int i = 0; i < x.size() - 1; ++i)
                {
                    std::vector<unsigned int> vid({z, z+1, z+xSize, z+xSize+1});
                    if (verbose)
                        printIndices(vid);

                    gridFactory.insertElement(type, vid);
                    z++;
                }
                z++;
            }

            else
            {
                for (int i = 0; i < x.size() - 1; ++i)
                {

                    std::vector<unsigned int> vid({z, z+1, t, t+1});

                    if (verbose)
                        printIndices(vid);

                    gridFactory.insertElement(type, vid);
                    t++;
                    z++;
                }
                t++;
                z++;

                if (verbose)
                    std::cout << "assign nodes ends..." << std::endl;
            }
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

    static void printIndices(const std::vector<unsigned int>& vid)
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
