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
 * \brief Grid manager specialization for OneDGrid
 */
#ifndef DUMUX_IO_GRID_MANAGER_ONED_HH
#define DUMUX_IO_GRID_MANAGER_ONED_HH

#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfoned.hh>

#ifndef DUMUX_IO_GRID_MANAGER_BASE_HH
#include <dumux/io/grid/gridmanager_base.hh>
#endif

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for OneDGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - LeftBoundary : start coordinate
 * - RightBoundary : end coordinate
 * - Cells : the number of cell
 * - RefinementType : local or copy
 * - Refinement : the number of global refines to apply initially.
 *
 */
template<>
class GridManager<Dune::OneDGrid>
: public GridManagerBase<Dune::OneDGrid>
{
public:
    using Grid = Dune::OneDGrid;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {

        // try to create it from a DGF or msh file in GridParameterGroup.File
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromDgfFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"));
            postProcessing_(modelParamGroup);
            return;
        }

        // Look for the necessary keys to construct from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.RightBoundary"))
        {
            // The required parameters
            using CoordinateType = typename Grid::ctype;
            const auto leftBoundary = getParamFromGroup<CoordinateType>(modelParamGroup, "Grid.LeftBoundary", 0.0);
            const auto rightBoundary = getParamFromGroup<CoordinateType>(modelParamGroup, "Grid.RightBoundary");
            const int cells = getParamFromGroup<int>(modelParamGroup, "Grid.Cells", 1);

            ParentType::gridPtr() = std::make_shared<Grid>(cells, leftBoundary, rightBoundary);
            postProcessing_(modelParamGroup);
            return;
        }

        // Look for the necessary keys to construct from the input file with just a coordinates vector
        else if (hasParamInGroup(modelParamGroup, "Grid.Coordinates"))
        {
            const auto coordinates = getParamFromGroup<std::vector<typename Grid::ctype>>(modelParamGroup, "Grid.Coordinates");
            ParentType::gridPtr() = std::make_shared<Grid>(coordinates);
            postProcessing_(modelParamGroup);
        }

        // Didn't find a way to construct the grid
        else
        {
            const auto prefix = modelParamGroup.empty() ? modelParamGroup : modelParamGroup + ".";
            DUNE_THROW(ParameterException, "Please supply one of the parameters "
                                           << prefix + "Grid.RightBoundary"
                                           << ", or " << prefix + "Grid.Coordinates"
                                           << ", or a grid file in " << prefix + "Grid.File");
        }
    }

    /*!
     * \brief Call loadBalance() function of the grid. OneDGrid is not parallel an thus cannot communicate.
     */
    void loadBalance() {}

private:
    /*!
     * \brief Do some operatrion after making the grid, like global refinement
     */
    void postProcessing_(const std::string& modelParamGroup)
    {
        // Set refinement type
        const auto refType = getParamFromGroup<std::string>(modelParamGroup, "Grid.RefinementType", "Local");
        if (refType == "Local")
            ParentType::grid().setRefinementType(Dune::OneDGrid::RefinementType::LOCAL);
        else if (refType == "Copy")
            ParentType::grid().setRefinementType(Dune::OneDGrid::RefinementType::COPY);
        else
            DUNE_THROW(Dune::IOError, "OneGrid only supports 'Local' or 'Copy' as refinment type. Not '"<< refType<<"'!");

        // Check if should refine the grid
        ParentType::maybeRefineGrid(modelParamGroup);
        loadBalance();
    }
};

} // end namespace Dumux

#endif
