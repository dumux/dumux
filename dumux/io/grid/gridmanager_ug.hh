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
 * \brief Grid manager specialization for UGGrid
 */
#ifndef DUMUX_IO_GRID_MANAGER_UG_HH
#define DUMUX_IO_GRID_MANAGER_UG_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif

#ifndef DUMUX_IO_GRID_MANAGER_BASE_HH
#include <dumux/io/grid/gridmanager_base.hh>
#endif

namespace Dumux {

#if HAVE_UG

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for UGGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 * - CellType : "Cube" or "Simplex" to be used for structured grids
 * - Refinement : the number of global refines to perform
 * - Verbosity : whether the grid construction should output to standard out
 * - BoundarySegments : whether to insert boundary segments into the grid
 *
 */
template<int dim>
class GridManager<Dune::UGGrid<dim>>
: public GridManagerBase<Dune::UGGrid<dim>>
{
public:
    using Grid = typename Dune::UGGrid<dim>;
    using ParentType = GridManagerBase<Grid>;
    using Element = typename Grid::template Codim<0>::Entity;

    /*!
     * \brief Make the UGGrid.
     */
    void init(const std::string& modelParamGroup = "")
    {

        // try to create it from a DGF or msh file in GridParameterGroup.File
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            preProcessing_(modelParamGroup);
            ParentType::makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            postProcessing_(modelParamGroup);
            return;
        }

        // Then look for the necessary keys to construct from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.UpperRight"))
        {
            preProcessing_(modelParamGroup);
            // make the grid
            const auto cellType = getParamFromGroup<std::string>(modelParamGroup, "Grid.CellType", "Cube");
            if (cellType == "Cube")
                ParentType::template makeStructuredGrid<dim, dim>(ParentType::CellType::Cube, modelParamGroup);
            else if (cellType == "Simplex")
                ParentType::template makeStructuredGrid<dim, dim>(ParentType::CellType::Simplex, modelParamGroup);
            else
                DUNE_THROW(Dune::IOError, "UGGrid only supports 'Cube' or 'Simplex' as cell type. Not '"<< cellType<<"'!");
            postProcessing_(modelParamGroup);
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
     * \brief Call loadBalance() function of the grid.
     * \note This overload is necessary because UGGrid cannot load balance element data yet.
     *       For dgf grids using element data in parallel will (hopefully) through an error
     *       For gmsh grids the parameters are read on every process so we use too much memory but
     *       the parameters are available via the insertionIndex of the level 0 element.
     */
    void loadBalance()
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
        {
            // if we may have dgf parameters use load balancing of the dgf pointer
            if(ParentType::enableDgfGridPointer_)
            {
                ParentType::dgfGridPtr().loadBalance();
                // update the grid data
                ParentType::gridData_ = std::make_shared<typename ParentType::GridData>(ParentType::dgfGridPtr());
            }

            // if we have gmsh parameters we have to manually load balance the data
            else if (ParentType::enableGmshDomainMarkers_)
            {
                // element and face markers are communicated during load balance
                auto dh = ParentType::gridData_->createGmshDataHandle();
                ParentType::gridPtr()->loadBalance(dh.interface());
                // Right now, UGGrid cannot communicate element data. If this gets implemented, communicate the data here:
                // ParentType::gridPtr()->communicate(dh.interface(), Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
            }
            else
                ParentType::gridPtr()->loadBalance();
        }
    }

private:
    /*!
     * \brief Do some operations before making the grid
     */
    void preProcessing_(const std::string& modelParamGroup)
    {}

    /*!
     * \brief Do some operatrion after making the grid, like global refinement
     */
    void postProcessing_(const std::string& modelParamGroup)
    {
        // Set refinement type
        const auto refType = getParamFromGroup<std::string>(modelParamGroup, "Grid.RefinementType", "Local");
        if (refType == "Local")
            ParentType::grid().setRefinementType(Dune::UGGrid<dim>::RefinementType::LOCAL);
        else if (refType == "Copy")
            ParentType::grid().setRefinementType(Dune::UGGrid<dim>::RefinementType::COPY);
        else
            DUNE_THROW(Dune::IOError, "UGGrid only supports 'Local' or 'Copy' as refinment type. Not '"<< refType<<"'!");

        // Set closure type
        const auto closureType = getParamFromGroup<std::string>(modelParamGroup, "Grid.ClosureType", "Green");
        if (closureType == "Green")
            ParentType::grid().setClosureType(Dune::UGGrid<dim>::ClosureType::GREEN);
        else if (closureType == "None")
            ParentType::grid().setClosureType(Dune::UGGrid<dim>::ClosureType::NONE);
        else
            DUNE_THROW(Dune::IOError, "UGGrid only supports 'Green' or 'None' as closure type. Not '"<< closureType<<"'!");

        // Check if should refine the grid
        ParentType::maybeRefineGrid(modelParamGroup);
        // do load balancing
        loadBalance();
    }
};

#endif // HAVE_UG

} // end namespace Dumux

#endif
