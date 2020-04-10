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
 * \brief Grid manager specialization for MMesh
 */
#ifndef DUMUX_IO_GRID_MANAGER_MMESH_HH
#define DUMUX_IO_GRID_MANAGER_MMESH_HH

#if HAVE_DUNE_MMESH
#include <dune/mmesh/mmesh.hh>
#endif

#ifndef DUMUX_IO_GRID_MANAGER_BASE_HH
#include <dumux/io/grid/gridmanager_base.hh>
#endif

namespace Dumux {

#if HAVE_DUNE_MMESH

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for Dune MMesh
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 * - BoundarySegments : whether to insert boundary segments into the grid
 *
 */
template<int dim>
class GridManager<Dune::MovingMesh<dim>>
: public GridManagerBase<Dune::MovingMesh<dim>>
{
public:
    using Grid = Dune::MovingMesh<dim>;
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
            using GlobalPosition = Dune::FieldVector<typename Grid::ctype, dim>;
            const auto upperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight");
            const auto lowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.LowerLeft", GlobalPosition(0.0));

            using CellArray = std::array<unsigned int, dim>;
            CellArray numCells; numCells.fill(1);
            numCells = getParamFromGroup<CellArray>(modelParamGroup, "Grid.Cells", numCells);

            // Insert uniformly spaced vertices
            std::array<unsigned int, dim> numVertices = numCells;
            for (int i = 0; i < dim; ++i)
              numVertices[i]++;

            Dune::MMeshImplicitGridFactory<Grid> factory;

            // Insert equally spaced vertices
            Dune::FactoryUtilities::MultiIndex<dim> index(numVertices);
            for (int i = 0; i < index.cycle(); ++i, ++index)
            {
              GlobalPosition pos(0);
              for (int j=0; j<dim; j++)
                pos[j] = lowerLeft[j] + index[j] * (upperRight[j]-lowerLeft[j])/(numVertices[j]-1);

              factory.insertVertex(pos);
            }

            this->gridPtr() = std::unique_ptr<Grid>(factory.createGrid());
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

#endif // HAVE_DUNE_MMESH

} // end namespace Dumux

#endif
