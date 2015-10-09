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
#ifndef DUMUX_VARIABLECLASS_ADAPTIVE_HH
#define DUMUX_VARIABLECLASS_ADAPTIVE_HH

#include <dune/grid/utility/persistentcontainer.hh>
#include <dumux/linear/vectorexchange.hh>
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
    typedef VariableClass<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename CellData::AdaptedValues AdaptedValues;

    typedef typename GridView::Grid Grid;
    typedef typename Grid::LevelGridView LevelGridView;
    typedef typename LevelGridView::template Codim<0>::Iterator LevelIterator;
    typedef Dune::PersistentContainer<Grid, AdaptedValues> PersistentContainer;

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

            for (LevelIterator eIt = levelView.template begin<0>(); eIt != levelView.template end<0>(); ++eIt)
            {
                //get your map entry
                AdaptedValues &adaptedValues = adaptationMap_[*eIt];

                // put your value in the map
                if (eIt->isLeaf())
                {
                    // get index
                    int indexI = this->index(*eIt);

                    CellData& cellData = this->cellData(indexI);

                    cellData.storeAdaptionValues(adaptedValues, *eIt);

                    adaptedValues.count = 1;
                }
                //Average in father
                if (eIt->level() > 0)
                {
                    auto father = eIt->father();
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

            for (LevelIterator eIt = levelView.template begin<0>(); eIt != levelView.template end<0>(); ++eIt)
            {
                // only treat non-ghosts, ghost data is communicated afterwards
                if (eIt->partitionType() == Dune::GhostEntity)
                    continue;

                if (!eIt->isNew())
                {
                    //entry is in map, write in leaf
                    if (eIt->isLeaf())
                    {
                        AdaptedValues &adaptedValues = adaptationMap_[*eIt];
                        int newIdxI = this->index(*eIt);

                        CellData& cellData = this->cellData(newIdxI);

                        cellData.setAdaptionValues(adaptedValues, *eIt);
                    }
                }
                else
                {
                    // value is not in map, interpolate from father element
                    if (eIt->level() > 0)
                    {
                        // create new entry: reconstruct from adaptationMap_[father] to a new
                        // adaptationMap_[son]
                        CellData::reconstructAdaptionValues(adaptationMap_, eIt->father(), *eIt, problem);

                        // access new son
                        AdaptedValues& adaptedValues = adaptationMap_[*eIt];
                        adaptedValues.count = 1;

                        // if we are on leaf, store reconstructed values of son in CellData object
                        if (eIt->isLeaf())
                        {
                            // acess new CellData object
                            int newIdxI = this->index(*eIt);
                            CellData& cellData = this->cellData(newIdxI);

                            cellData.setAdaptionValues(adaptedValues, *eIt);
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
        typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
        typedef typename SolutionTypes::ElementMapper ElementMapper;
        typedef VectorExchange<ElementMapper, std::vector<CellData> > DataHandle;
        DataHandle dataHandle(problem.elementMapper(), this->cellDataGlobal());
        problem.gridView().template communicate<DataHandle>(dataHandle,
                                                            Dune::InteriorBorder_All_Interface,
                                                            Dune::ForwardCommunication);
#endif
    }

};
}
#endif
