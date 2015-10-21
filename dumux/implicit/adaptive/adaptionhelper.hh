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
#ifndef DUMUX_ADAPTIONHELPER_HH
#define DUMUX_ADAPTIONHELPER_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/propertysystem.hh>

/**
 * @file
 * @brief  Base class holding the variables for implicit models.
 */

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ImplicitIsBox);
NEW_PROP_TAG(PrimaryVariables);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(Scalar);
}

template<class TypeTag>
class AdaptionHelper
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    struct AdaptedValues
    {
    	PrimaryVariables u;
        int count;
        AdaptedValues()
        {
            count = 0;
        }
    };

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    typedef typename GridView::Grid Grid;
    typedef typename Grid::LevelGridView LevelGridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<dofCodim>::Entity DofEntity;
    typedef Dune::PersistentContainer<Grid, AdaptedValues> PersistentContainer;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;

    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;

private:
    const GridView gridView_;
    const Grid& grid_;
    PersistentContainer adaptionMap_;

public:
    //! Constructs an adaptive helper object
    /**
     * In addition to providing a storage object for cell-centered Methods, this class provides
     * mapping functionality to adapt the grid.
     *
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */
    AdaptionHelper(const GridView& gridView) :
    	gridView_(gridView), grid_(gridView.grid()), adaptionMap_(grid_, dofCodim)
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
        adaptionMap_.resize();

        // loop over all levels of the grid
        for (int level = grid_.maxLevel(); level >= 0; level--)
        {
            //get grid view on level grid
            LevelGridView levelView = grid_.levelGridView(level);

            if(!isBox)
            {
				for (const auto& element : Dune::elements(levelView))
				{
					//get your map entry
					AdaptedValues &adaptedValues = adaptionMap_[element];

					// put your value in the map
					if (element.isLeaf())
					{
						// get index
						int indexI = this->elementIndex(problem, element);

						storeAdaptionValues(adaptedValues, problem.model().curSol()[indexI]);

						adaptedValues.count = 1;
					}
					//Average in father
					if (element.level() > 0)
					{
						AdaptedValues& adaptedValuesFather = adaptionMap_[element.father()];
						adaptedValuesFather.count += 1;
						storeAdaptionValues(adaptedValues, adaptedValuesFather);
					}

				}
        	}
            else
            {
				for (const auto& entity : Dune::entities(levelView, Dune::Codim<dofCodim>()))
				{
					//get your map entry
					AdaptedValues &adaptedValues = adaptionMap_[entity];

					// put your value in the map
					int indexI = this->dofIndex(problem, entity);

					storeAdaptionValues(adaptedValues, problem.model().curSol()[indexI]);

					adaptedValues.count = 1;

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
        adaptionMap_.resize();

        for (int level = 0; level <= grid_.maxLevel(); level++)
        {
            LevelGridView levelView = grid_.levelGridView(level);

            for (const auto& element : Dune::elements(levelView))
            {
                // only treat non-ghosts, ghost data is communicated afterwards
                if (element.partitionType() == Dune::GhostEntity)
                    continue;

                if (!element.isNew())
                {
                    //entry is in map, write in leaf
                    if (element.isLeaf())
                    {
                    	if(!isBox)
                    	{
							AdaptedValues &adaptedValues = adaptionMap_[element];
							int newIdxI = this->elementIndex(problem, element);

							setAdaptionValues(adaptedValues, problem.model().curSol()[newIdxI]);
						}
                    	else
                    	{
                            unsigned int numSubEntities = element.subEntities(dofCodim);

                        	for(unsigned int i = 0; i < numSubEntities; i++)
                        	{
                                auto subEntity = element.template subEntity <dofCodim>(i);
                                AdaptedValues &adaptedValues = adaptionMap_[subEntity];
    							int newIdxI = this->dofIndex(problem, subEntity);

    							setAdaptionValues(adaptedValues, problem.model().curSol()[newIdxI]);

                        	}
                    	}
                    }
                }
                else
                {
                    // value is not in map, interpolate from father element
                    if (element.level() > 0 && element.hasFather())
                    {
                        auto epFather = element.father();

                        if(!isBox)
                        {
							// create new entry: reconstruct from adaptionMap_[*father] to a new
							// adaptionMap_[*son]
							reconstructAdaptionValues(adaptionMap_, epFather, element, problem);

							// access new son
							AdaptedValues& adaptedValues = adaptionMap_[element];
							adaptedValues.count = 1;

							// if we are on leaf, store reconstructed values of son in CellData object
							if (element.isLeaf())
							{
								// acess new CellData object
								int newIdxI = this->elementIndex(problem, element);

								setAdaptionValues(adaptedValues, problem.model().curSol()[newIdxI]);
							}
                        }
                        else
                        {
                            unsigned int numSubEntities= element.subEntities(dofCodim);
                    		const auto geometryI = element.geometry();

                        	for(unsigned int i = 0; i < numSubEntities; i++)
                        	{
                        		auto subEntity = element.template subEntity <dofCodim>(i);
    							AdaptedValues &adaptedValues = adaptionMap_[subEntity];

    							if(adaptedValues.count == 0){
									LocalPosition dofCenterPos = geometryI.local(subEntity.geometry().center());
									const LocalFiniteElementCache feCache;
									Dune::GeometryType geomType = epFather.geometry().type();

									const LocalFiniteElement &localFiniteElement = feCache.get(geomType);
									std::vector<Dune::FieldVector<Scalar, 1> > shapeVal;
									localFiniteElement.localBasis().evaluateFunction(dofCenterPos, shapeVal);
									PrimaryVariables u(0);
									for (int j = 0; j < shapeVal.size(); ++j)
									{
										AdaptedValues & adaptedValuesFather = adaptionMap_[epFather.template subEntity <dofCodim>(j)];
										u.axpy(shapeVal[j], adaptedValuesFather.u);
									}

									adaptedValues.u = u;
									adaptedValues.count = 1;
    							}

    							if (element.isLeaf())
    							{
        							int newIdxI = this->dofIndex(problem, subEntity);
    								setAdaptionValues(adaptedValues, problem.model().curSol()[newIdxI]);
    							}

                        	}
                        }
                    }

                }
            }

        }
        // reset entries in restrictionmap
        adaptionMap_.resize( typename PersistentContainer::Value() );
        adaptionMap_.shrinkToFit();
        adaptionMap_.fill( typename PersistentContainer::Value() );

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
     * the adaption container in order to be mapped on a new grid.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param u The variables to be stored
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
     */
    static void storeAdaptionValues(AdaptedValues& adaptedValues,
                                    AdaptedValues& adaptedValuesFather)
    {
    	if(!isBox)
    	{
			adaptedValuesFather.u += adaptedValues.u;
			adaptedValuesFather.u /= adaptedValues.count;
    	}
    	else
    	{
    		adaptedValuesFather.u = adaptedValues.u;
    	}
    }
    //! Set adapted values in CellData
    /**
     * This methods stores reconstructed values into the cellData object, by
     * this setting a newly mapped solution to the storage container of the
     * decoupled models.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param u The variables to be stored
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
     * adaptionMap.
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

    int dofIndex(const Problem& problem, const DofEntity& entity) const
    {
                return problem.model().dofMapper().index(entity);
    }

    int elementIndex(const Problem& problem, const Element& element) const
    {
                return problem.elementMapper().index(element);
    }

};
}
#endif
