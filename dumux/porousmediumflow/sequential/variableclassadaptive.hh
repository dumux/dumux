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
#ifndef DUMUX_VARIABLECLASS_ADAPTIVE_HH
#define DUMUX_VARIABLECLASS_ADAPTIVE_HH

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>
#include "variableclass.hh"

/**
 * @file
 * @brief  Base class holding the variables for sequential models.
 */

namespace Dumux
{
/*!
 * \ingroup IMPET
 */
//! Base class holding the variables and discretized data for sequential models.
/*!
 * Stores global information and variables that are common for all sequential models and also functions needed to access these variables.
 * Can be directly used for a single phase model.
 *
 * @tparam TypeTag The Type Tag
 *
 */
template<class TypeTag>
class VariableClassAdaptive: public VariableClass<TypeTag>
{
private:
    using ParentType = VariableClass<TypeTag>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using CellData = GetPropType<TypeTag, Properties::CellData>;
    using AdaptedValues = typename CellData::AdaptedValues;

    using Grid = typename GridView::Grid;
    using LevelGridView = typename Grid::LevelGridView;
    using PersistentContainer = Dune::PersistentContainer<Grid, AdaptedValues>;

private:
    const Grid& grid_;
    PersistentContainer adaptationMap_;

public:
    //! Constructs an adaptive VariableClass object
    /**
     * In addition to providing a storage object for cell-centered Methods, this class provides
     * mapping functionality to adapt the grid.
     *
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */
    VariableClassAdaptive(const GridView& gridView) :
        ParentType(gridView), grid_(gridView.grid()), adaptationMap_(grid_, 0)
    {}


    /*!
     * Store primary variables
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     * From upper level on downwards, the old solution is stored into an container
     * object, before the grid is adapted. Father elements hold averaged information
     * from the son cells for the case of the sons being coarsened.
     *
     * @param problem The current problem
     */
    void storePrimVars(const Problem& problem)
    {
        adaptationMap_.resize();

        // loop over all levels of the grid
        for (int level = grid_.maxLevel(); level >= 0; level--)
        {
            //get grid view on level grid
            LevelGridView levelView = grid_.levelGridView(level);

            for (const auto& element : elements(levelView))
            {
                //get your map entry
                AdaptedValues &adaptedValues = adaptationMap_[element];

                // put your value in the map
                if (element.isLeaf())
                {
                    // get index
                    int indexI = this->index(element);

                    CellData& cellData = this->cellData(indexI);

                    cellData.storeAdaptionValues(adaptedValues, element);

                    adaptedValues.count = 1;
                }
                //Average in father
                if (element.level() > 0)
                {
                    auto father = element.father();
                    AdaptedValues& adaptedValuesFather = adaptationMap_[father];
                    adaptedValuesFather.count += 1;
                    CellData::storeAdaptionValues(adaptedValues, adaptedValuesFather, father);
                }
            }
        }
    }

    /*!
     * Reconstruct missing primary variables (where elements are created/deleted)
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     * Starting from the lowest level, the old solution is mapped on the new grid:
     * Where coarsened, new cells get information from old father element.
     * Where refined, a new solution is reconstructed from the old father cell,
     * and then a new son is created. That is then stored into the general data
     * structure (CellData).
     *
     * @param problem The current problem
     */
    void reconstructPrimVars(const Problem& problem)
    {
        adaptationMap_.resize();

        for (int level = 0; level <= grid_.maxLevel(); level++)
        {
            LevelGridView levelView = grid_.levelGridView(level);

            for (const auto& element : elements(levelView))
            {
                // only treat non-ghosts, ghost data is communicated afterwards
                if (element.partitionType() == Dune::GhostEntity)
                    continue;

                if (!element.isNew())
                {
                    //entry is in map, write in leaf
                    if (element.isLeaf())
                    {
                        AdaptedValues &adaptedValues = adaptationMap_[element];
                        int newIdxI = this->index(element);

                        CellData& cellData = this->cellData(newIdxI);

                        cellData.setAdaptionValues(adaptedValues, element);
                    }
                }
                else
                {
                    // value is not in map, interpolate from father element
                    if (element.level() > 0)
                    {
                        // create new entry: reconstruct from adaptationMap_[father] to a new
                        // adaptationMap_[son]
                        CellData::reconstructAdaptionValues(adaptationMap_, element.father(), element, problem);

                        // access new son
                        AdaptedValues& adaptedValues = adaptationMap_[element];
                        adaptedValues.count = 1;

                        // if we are on leaf, store reconstructed values of son in CellData object
                        if (element.isLeaf())
                        {
                            // acess new CellData object
                            int newIdxI = this->index(element);
                            CellData& cellData = this->cellData(newIdxI);

                            cellData.setAdaptionValues(adaptedValues, element);
                        }
                    }
                }
            }

        }
        // reset entries in restrictionmap
        adaptationMap_.resize( typename PersistentContainer::Value() );
        adaptationMap_.shrinkToFit();
        adaptationMap_.fill( typename PersistentContainer::Value() );

#if HAVE_MPI
        // communicate ghost data
        using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
        using ElementMapper = typename SolutionTypes::ElementMapper;
        using DataHandle = VectorCommDataHandleEqual<ElementMapper, std::vector<CellData>, 0/*elementCodim*/>;
        DataHandle dataHandle(problem.elementMapper(), this->cellDataGlobal());
        problem.gridView().template communicate<DataHandle>(dataHandle,
                                                            Dune::InteriorBorder_All_Interface,
                                                            Dune::ForwardCommunication);
#endif
    }

};
}
#endif
