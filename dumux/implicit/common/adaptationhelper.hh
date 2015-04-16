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
#ifndef DUMUX_ADAPTATIONHELPER_HH
#define DUMUX_ADAPTATIONHELPER_HH

#include <dune/grid/utility/persistentcontainer.hh>
//#include <dumux/linear/vectorexchange.hh>

/**
 * @file
 * @brief  Base class holding the variables for implicit models.
 */

namespace Dumux
{

template<class TypeTag>
class AdaptationHelper
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    struct AdaptedValues
    {
    	PrimaryVariables u;
        int count;
        AdaptedValues()
        {
            count = 0;
        }
    };

    typedef typename GridView::Grid Grid;
    typedef typename Grid::LevelGridView LevelGridView;
    typedef typename LevelGridView::template Codim<0>::Iterator LevelIterator;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::PersistentContainer<Grid, AdaptedValues> PersistentContainer;

private:
    const GridView gridView_;
    const Grid& grid_;
    PersistentContainer adaptationMap_;

public:
    //! Constructs an adaptive helper object
    /**
     * In addition to providing a storage object for cell-centered Methods, this class provides
     * mapping functionality to adapt the grid.
     *
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */
    AdaptationHelper(const GridView& gridView) :
    	gridView_(gridView), grid_(gridView.grid()), adaptationMap_(grid_, 0)
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
    void storePrimVars(Problem& problem)
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
                    int indexI = this->elementIndex(*eIt, problem);

                    storeAdaptionValues(adaptedValues, problem.model().curSol()[indexI]);

                    adaptedValues.count = 1;
                }
                //Average in father
                if (eIt->level() > 0)
                {
                    ElementPointer epFather = eIt->father();
                    int indexI = this->elementIndex(*epFather, problem);
                    AdaptedValues& adaptedValuesFather = adaptationMap_[*epFather];
                    adaptedValuesFather.count += 1;
                    storeAdaptionValues(adaptedValues, adaptedValuesFather,
                    					problem.model().curSol()[indexI]);
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
    void reconstructPrimVars(Problem& problem)
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
                        int newIdxI = this->elementIndex(*eIt, problem);

                        setAdaptionValues(adaptedValues, problem.model().curSol()[newIdxI]);
                    }
                }
                else
                {
                    // value is not in map, interpolate from father element
                    if (eIt->level() > 0)
                    {
                        ElementPointer epFather = eIt->father();

                        // create new entry: reconstruct from adaptationMap_[*father] to a new
                        // adaptationMap_[*son]
                        reconstructAdaptionValues(adaptationMap_, *epFather, *eIt, problem);

                        // access new son
                        AdaptedValues& adaptedValues = adaptationMap_[*eIt];
                        adaptedValues.count = 1;

                        // if we are on leaf, store reconstructed values of son in CellData object
                        if (eIt->isLeaf())
                        {
                            // acess new CellData object
                            int newIdxI = this->elementIndex(*eIt, problem);

                            setAdaptionValues(adaptedValues, problem.model().curSol()[newIdxI]);
                        }
                    }
                }
            }

        }
        // reset entries in restrictionmap
        adaptationMap_.resize( typename PersistentContainer::Value() );
        adaptationMap_.shrinkToFit();
        adaptationMap_.fill( typename PersistentContainer::Value() );

//#if HAVE_MPI
//        // communicate ghost data
//        typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
//        typedef typename SolutionTypes::ElementMapper ElementMapper;
//        typedef VectorExchange<ElementMapper, std::vector<CellData> > DataHandle;
//        DataHandle dataHandle(problem.elementMapper(), this->cellDataGlobal());
//        problem.gridView().template communicate<DataHandle>(dataHandle,
//                                                            Dune::InteriorBorder_All_Interface,
//                                                            Dune::ForwardCommunication);
//#endif
    }

    //! Stores values to be adapted in an adaptedValues container
    /**
     * Stores values to be adapted from the current CellData objects into
     * the adaptation container in order to be mapped on a new grid.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param element The element to be stored
     */
    static void storeAdaptionValues(AdaptedValues& adaptedValues, const PrimaryVariables& u)
    {
        adaptedValues.u = u;
    }
    //! Stores sons entries into father element for averaging
    /**
     * Sum up the adaptedValues (sons values) into father element. We store from leaf
     * upwards, so sons are stored first, then cells on the next leaf (=fathers)
     * can be averaged.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param adaptedValuesFather Values to be adapted of father cell
     * \param fatherElement The element of the father
     */
    static void storeAdaptionValues(AdaptedValues& adaptedValues,
                                    AdaptedValues& adaptedValuesFather,
									const PrimaryVariables& u)
    {
    	adaptedValuesFather.u += adaptedValues.u;
    	adaptedValuesFather.u /= adaptedValues.count;
    }
    //! Set adapted values in CellData
    /**
     * This methods stores reconstructed values into the cellData object, by
     * this setting a newly mapped solution to the storage container of the
     * decoupled models.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param element The element where things are stored.
     */
    static void setAdaptionValues(AdaptedValues& adaptedValues, PrimaryVariables& u)
    {
    	PrimaryVariables uNew = adaptedValues.u;
    	uNew /= adaptedValues.count;

    	u = uNew;
    }

    //! Reconstructs sons entries from data of father cell
    /**
     * Reconstructs a new solution from a father cell into a newly
     * generated son cell. New cell is stored into the global
     * adaptationMap.
     *
     * \param adaptionMap Global map storing all values to be adapted
     * \param father Entity Pointer to the father cell
     * \param son Entity Pointer to the newly created son cell
     * \param problem The problem
     */
    static void reconstructAdaptionValues(Dune::PersistentContainer<Grid, AdaptedValues>& adaptionMap,
            const Element& father, const Element& son, const Problem& problem)
    {
        AdaptedValues& adaptedValues = adaptionMap[son];
        AdaptedValues& adaptedValuesFather = adaptionMap[father];

        adaptedValues.u = adaptedValuesFather.u;
        adaptedValues.u /= adaptedValuesFather.count;

    }

    int elementIndex(const Element& element, const Problem& problem) const
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        return problem.elementMapper().index(element);
#else
        return problem.elementMapper().map(element);
#endif
    }

};
}
#endif
