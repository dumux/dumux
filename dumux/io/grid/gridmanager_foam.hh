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
 * \brief Grid manager specialization for FoamGrid
 */
#ifndef DUMUX_IO_GRID_MANAGER_FOAM_HH
#define DUMUX_IO_GRID_MANAGER_FOAM_HH

// FoamGrid specific includes
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>
#endif

#ifndef DUMUX_IO_GRID_MANAGER_BASE_HH
#include <dumux/io/grid/gridmanager_base.hh>
#endif

namespace Dumux {

#if HAVE_DUNE_FOAMGRID

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for FoamGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - Verbosity : whether the grid construction should output to standard out
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 *
 */
template<int dim, int dimworld>
class GridManager<Dune::FoamGrid<dim, dimworld>>
: public GridManagerBase<Dune::FoamGrid<dim, dimworld>>
{
public:
    using Grid = Dune::FoamGrid<dim, dimworld>;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {
        // try to create it from file
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            ParentType::loadBalance();
            return;
        }

        // Then look for the necessary keys to construct a structured grid from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.UpperRight"))
        {
            ParentType::template makeStructuredGrid<dim, dimworld>(ParentType::CellType::Simplex, modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            ParentType::loadBalance();
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
};

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for FoamGrids of dim 1
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - Verbosity : whether the grid construction should output to standard out
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 *
 */
template<int dimworld>
class GridManager<Dune::FoamGrid<1, dimworld>>
: public GridManagerBase<Dune::FoamGrid<1, dimworld>>
{
public:
    using Grid = Dune::FoamGrid<1, dimworld>;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {
        // try to create it from file
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            ParentType::loadBalance();
            return;
        }

        // The required parameters
        using GlobalPosition = Dune::FieldVector<typename Grid::ctype, dimworld>;
        const auto upperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight");
        const auto lowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.LowerLeft", GlobalPosition(0.0));
        using CellArray = std::array<unsigned int, 1>;
        const auto cells = getParamFromGroup<CellArray>(modelParamGroup, "Grid.Cells", std::array<unsigned int, 1>{{1}});

        // make the grid (structured interval grid in dimworld space)
        Dune::GridFactory<Grid> factory;

        constexpr auto geomType = Dune::GeometryTypes::line;

        // create a step vector
        GlobalPosition step = upperRight;
        step -= lowerLeft, step /= cells[0];

        // create the vertices
        GlobalPosition globalPos = lowerLeft;
        for (unsigned int vIdx = 0; vIdx <= cells[0]; vIdx++, globalPos += step)
            factory.insertVertex(globalPos);

        // create the cells
        for(unsigned int eIdx = 0; eIdx < cells[0]; eIdx++)
            factory.insertElement(geomType, {eIdx, eIdx+1});

        ParentType::gridPtr() = std::shared_ptr<Grid>(factory.createGrid());
        ParentType::maybeRefineGrid(modelParamGroup);
        ParentType::loadBalance();
    }
};

#endif // HAVE_DUNE_FOAMGRID

} // end namespace Dumux

#endif
