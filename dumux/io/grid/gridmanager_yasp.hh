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
 * \brief Grid manager specialization for YaspGrid
 */
#ifndef DUMUX_IO_GRID_MANAGER_YASP_HH
#define DUMUX_IO_GRID_MANAGER_YASP_HH

#include <dune/common/math.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#ifndef DUMUX_IO_GRID_MANAGER_BASE_HH
#include <dumux/io/grid/gridmanager_base.hh>
#endif

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for YaspGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - File : a DGF file to load the coarse grid from
 * - LowerLeft : lower left corner of the domain (only for EquidistantOffsetCoordinates, otherwise 0-origin)
 * - UpperRight : upper right corner of the domain
 * - Cells : the number of cells in each direction
 * - Periodic : true or false for each direction
 * - Overlap : overlap size in cells
 * - Partitioning : a non-standard load-balancing, number of processors per direction
 * - KeepPyhsicalOverlap : whether to keep the physical overlap
 *     in physical size or in number of cells upon refinement
 * - Refinement : the number of global refines to apply initially.
 *
 */
template<class Coordinates, int dim>
class GridManager<Dune::YaspGrid<dim, Coordinates>>
: public GridManagerBase<Dune::YaspGrid<dim, Coordinates>>
{
    using ct = typename Dune::YaspGrid<dim, Coordinates>::ctype;
    using GlobalPosition = Dune::FieldVector<ct, dim>;
public:
    using Grid = typename Dune::YaspGrid<dim, Coordinates>;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {
        // First try to create it from a DGF file in GridParameterGroup.File
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromDgfFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"));
            postProcessing_(modelParamGroup);
            return;
        }

        // Then look for the necessary keys to construct from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.UpperRight"))
        {
            // get the upper right corner coordinates
            const auto upperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight");

            // number of cells in each direction
            std::array<int, dim> cells; cells.fill(1);
            cells = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Cells", cells);

            // \todo TODO periodic boundaries with yasp (the periodicity concept of yasp grid is currently not supported, use dune-spgrid)
            // const auto periodic = getParamFromGroup<std::bitset<dim>>(modelParamGroup, "Grid.Periodic", std::bitset<dim>());
            const std::bitset<dim> periodic;

            // get the overlap
            const int overlap =  getParamFromGroup<int>(modelParamGroup, "Grid.Overlap", 1);

            if constexpr (std::is_same_v<Dune::EquidistantCoordinates<ct, dim>, Coordinates>)
                init(upperRight, cells, modelParamGroup, overlap, periodic);
            else
            {
                const auto lowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.LowerLeft", GlobalPosition(0.0));
                init(lowerLeft, upperRight, cells, modelParamGroup, overlap, periodic);
            }
        }

        // Didn't find a way to construct the grid
        else
        {
            const auto prefix = modelParamGroup.empty() ? modelParamGroup : modelParamGroup + ".";
            DUNE_THROW(ParameterException, "Please supply one of the parameters "
                                           << prefix + "Grid.UpperRight"
                                           << ", or a grid file in " << prefix + "Grid.File");

        }
    }

    /*!
     * \brief Make the grid using input data not read from the input file.
     * \note Use this function for EquidistantCoordinates where lowerLeft is always zero.
     */
    void init(const GlobalPosition& upperRight,
              const std::array<int, dim>& cells,
              const std::string& modelParamGroup = "",
              const int overlap = 1,
              const std::bitset<dim> periodic = std::bitset<dim>{})
    {
        static_assert(std::is_same_v<Dune::EquidistantCoordinates<ct, dim>, Coordinates>,
                      "Use init function taking lowerLeft as argument when working with EquidistantOffsetCoordinates");

        if (!hasParamInGroup(modelParamGroup, "Grid.Partitioning"))
        {
            // construct using default load balancing
            ParentType::gridPtr() = std::make_unique<Grid>(upperRight, cells, periodic, overlap);
        }
        else
        {
            // construct using user defined partitioning
            const auto partitioning = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Partitioning");
            Dune::YaspFixedSizePartitioner<dim> lb(partitioning);
            ParentType::gridPtr() = std::make_unique<Grid>(upperRight, cells, periodic, overlap, typename Grid::CollectiveCommunicationType(), &lb);
        }

        postProcessing_(modelParamGroup);
    }

    /*!
     * \brief Make the grid using input data not read from the input file.
     * \note Use this function for EquidistantOffsetCoordinates.
     */
    void init(const GlobalPosition& lowerLeft,
              const GlobalPosition& upperRight,
              const std::array<int, dim>& cells,
              const std::string& modelParamGroup = "",
              const int overlap = 1,
              const std::bitset<dim> periodic = std::bitset<dim>{})
    {
        static_assert(std::is_same_v<Dune::EquidistantOffsetCoordinates<ct, dim>, Coordinates>,
                      "LowerLeft can only be specified with EquidistantOffsetCoordinates");

        if (!hasParamInGroup(modelParamGroup, "Grid.Partitioning"))
        {
            // construct using default load balancing
            ParentType::gridPtr() = std::make_unique<Grid>(lowerLeft, upperRight, cells, periodic, overlap);
        }
        else
        {
            // construct using user defined partitioning
            const auto partitioning = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Partitioning");
            Dune::YaspFixedSizePartitioner<dim> lb(partitioning);
            ParentType::gridPtr() = std::make_unique<Grid>(lowerLeft, upperRight, cells, periodic, overlap, typename Grid::CollectiveCommunicationType(), &lb);
        }

        postProcessing_(modelParamGroup);
    }

private:

    /*!
     * \brief Postprocessing for YaspGrid
     */
    void postProcessing_(const std::string& modelParamGroup)
    {
        // Check if should refine the grid
        const bool keepPhysicalOverlap = getParamFromGroup<bool>(modelParamGroup, "Grid.KeepPhysicalOverlap", true);
        ParentType::grid().refineOptions(keepPhysicalOverlap);
        ParentType::maybeRefineGrid(modelParamGroup);
        ParentType::loadBalance();
    }
};

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for YaspGrids with different zones and grading
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - Positions0 : position array for x-coordinate
 * - Positions1 : position array for y-coordinate
 * - Positions2 : position array for z-coordinate
 * - Cells0 : number of cells array for x-coordinate
 * - Cells1 : number of cells array for y-coordinate
 * - Cells2 : number of cells array for z-coordinate
 * - Grading0 : grading factor array for x-coordinate
 * - Grading1 : grading factor array for y-coordinate
 * - Grading2 : grading factor array for z-coordinate
 * - Verbosity : whether the grid construction should output to standard out
 * - Periodic : true or false for each direction
 * - Overlap : overlap size in cells
 * - Partitioning : a non-standard load-balancing, number of processors per direction
 * - KeepPyhsicalOverlap : whether to keep the physical overlap
 *     in physical size or in number of cells upon refinement
 * - Refinement : the number of global refines to apply initially.
 *
 * The grading factor \f$ g \f$ specifies the ratio between the next and the current cell size:
 * \f$ g = \frac{h_{i+1}}{h_i} \f$.
 * Negative grading factors are converted to
 * \f$ g = -\frac{1}{g_\textrm{negative}} \f$
 * to avoid issues with imprecise fraction numbers.
 */
template<class ctype, int dim>
class GridManager<Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim> >>
: public GridManagerBase<Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim> > >
{
public:
    using Grid = typename Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim> >;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {
        // Only construction from the input file is possible
        // Look for the necessary keys to construct from the input file
        // The positions
        std::array<std::vector<ctype>, dim> positions;
        for (int i = 0; i < dim; ++i)
            positions[i] = getParamFromGroup<std::vector<ctype>>(modelParamGroup, "Grid.Positions" + std::to_string(i));

        // the number of cells (has a default)
        std::array<std::vector<int>, dim> cells;
        for (int i = 0; i < dim; ++i)
        {
            cells[i].resize(positions[i].size()-1, 1.0);
            cells[i] = getParamFromGroup<std::vector<int>>(modelParamGroup, "Grid.Cells" + std::to_string(i), cells[i]);
        }

        // grading factor (has a default)
        std::array<std::vector<ctype>, dim> grading;
        for (int i = 0; i < dim; ++i)
        {
            grading[i].resize(positions[i].size()-1, 1.0);
            grading[i] = getParamFromGroup<std::vector<ctype>>(modelParamGroup, "Grid.Grading" + std::to_string(i), grading[i]);
        }

        // call the generic function
        init(positions, cells, grading, modelParamGroup);
    }

    /*!
     * \brief Make the grid using input data not read from the input file.
     */
    void init(const std::array<std::vector<ctype>, dim>& positions,
              const std::array<std::vector<int>, dim>& cells,
              const std::array<std::vector<ctype>, dim>& grading,
              const std::string& modelParamGroup = "")
    {
        // Additional parameters (they have a default)
        const int overlap = getParamFromGroup<int>(modelParamGroup, "Grid.Overlap", 1);
        const bool verbose = getParamFromGroup<bool>(modelParamGroup, "Grid.Verbosity", false);
        // \todo TODO periodic boundaries with yasp (the periodicity concept of yasp grid is currently not supported, use dune-spgrid)
        // const auto periodic = getParamFromGroup<std::bitset<dim>>(modelParamGroup, "Grid.Periodic", std::bitset<dim>());
        const std::bitset<dim> periodic;

        // Some sanity checks
        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            if (cells[dimIdx].size() + 1 != positions[dimIdx].size())
            {
                DUNE_THROW(Dune::RangeError, "Make sure to specify correct \"Cells\" and \"Positions\" arrays");
            }
            if (grading[dimIdx].size() + 1 != positions[dimIdx].size())
            {
                DUNE_THROW(Dune::RangeError, "Make sure to specify correct \"Grading\" and \"Positions\" arrays");
            }
            ctype temp = std::numeric_limits<ctype>::lowest();
            for (unsigned int posIdx = 0; posIdx < positions[dimIdx].size(); ++posIdx)
            {
                if (temp > positions[dimIdx][posIdx])
                {
                    DUNE_THROW(Dune::RangeError, "Make sure to specify a monotone increasing \"Positions\" array");
                }
                temp = positions[dimIdx][posIdx];
            }
        }

        const auto globalPositions = computeGlobalPositions_(positions, cells, grading, verbose);

        // make the grid
        if (!hasParamInGroup(modelParamGroup, "Grid.Partitioning"))
        {
            // construct using default load balancing
            ParentType::gridPtr() = std::make_shared<Grid>(globalPositions, periodic, overlap);
        }
        else
        {
            // construct using user defined partitioning
            const auto partitioning = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Partitioning");
            Dune::YaspFixedSizePartitioner<dim> lb(partitioning);
            ParentType::gridPtr() = std::make_shared<Grid>(globalPositions, periodic, overlap, typename Grid::CollectiveCommunicationType(), &lb);
        }

        postProcessing_(modelParamGroup);
    }

private:
    /*!
     * \brief Postprocessing for YaspGrid
     */
    void postProcessing_(const std::string& modelParamGroup)
    {
        // Check if should refine the grid
        const bool keepPhysicalOverlap = getParamFromGroup<bool>(modelParamGroup, "Grid.KeepPhysicalOverlap", true);
        ParentType::grid().refineOptions(keepPhysicalOverlap);
        ParentType::maybeRefineGrid(modelParamGroup);
        ParentType::loadBalance();
    }

    //! Compute the global position tensor grid from the given positions, cells, and grading factors
    std::array<std::vector<ctype>, dim>
    computeGlobalPositions_(const std::array<std::vector<ctype>, dim>& positions,
                            const std::array<std::vector<int>, dim>& cells,
                            const std::array<std::vector<ctype>, dim>& grading,
                            bool verbose = false)
    {
        std::array<std::vector<ctype>, dim> globalPositions;
        using std::pow;
        using Dune::power;
        for (int dimIdx = 0; dimIdx < dim; dimIdx++)
        {
            for (int zoneIdx = 0; zoneIdx < cells[dimIdx].size(); ++zoneIdx)
            {
                ctype lower = positions[dimIdx][zoneIdx];
                ctype upper = positions[dimIdx][zoneIdx+1];
                int numCells = cells[dimIdx][zoneIdx];
                ctype gradingFactor = grading[dimIdx][zoneIdx];
                ctype length = upper - lower;
                ctype height = 1.0;
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
                {
                    increasingCellSize = true;
                }

                // take absolute values and reverse cell size increment to achieve
                // reverse behavior for negative values
                if (gradingFactor < 0.0)
                {
                    using std::abs;
                    gradingFactor = abs(gradingFactor);
                    if (gradingFactor < 1.0)
                    {
                        increasingCellSize = true;
                    }
                }

                // if the grading factor is exactly 1.0 do equal spacing
                if (gradingFactor > 1.0 - 1e-7 && gradingFactor < 1.0 + 1e-7)
                {
                    height = 1.0 / numCells;
                    if (verbose)
                    {
                        std::cout << " -> h "  << height * length << std::endl;
                    }
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

                std::vector<ctype> localPositions;
                localPositions.push_back(0);
                for (int i = 0; i < numCells-1; i++)
                {
                    ctype hI = height;
                    if (!(gradingFactor < 1.0 + 1e-7 && gradingFactor > 1.0 - 1e-7))
                    {
                        if (increasingCellSize)
                        {
                            hI *= power(gradingFactor, i);
                        }
                        else
                        {
                            hI *= power(gradingFactor, numCells-i-1);
                        }
                    }
                    localPositions.push_back(localPositions[i] + hI);
                }

                for (int i = 0; i < localPositions.size(); i++)
                {
                    localPositions[i] *= length;
                    localPositions[i] += lower;
                }

                for (unsigned int i = 0; i < localPositions.size(); ++i)
                {
                    globalPositions[dimIdx].push_back(localPositions[i]);
                }
            }
            globalPositions[dimIdx].push_back(positions[dimIdx].back());
        }

        return globalPositions;
    }
};

} // end namespace Dumux

#endif
