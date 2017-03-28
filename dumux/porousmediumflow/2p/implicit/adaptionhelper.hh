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
#ifndef DUMUX_TWOP_ADAPTIONHELPER_HH
#define DUMUX_TWOP_ADAPTIONHELPER_HH

#include <dumux/implicit/adaptive/adaptionhelper.hh>
#include "properties.hh"

namespace Dumux {

namespace Properties
{
NEW_PROP_TAG(PrimaryVariables);
NEW_PROP_TAG(Scalar);
}

/*!
 * \brief Base class holding the variables for implicit models.
 */
template<class TypeTag>
class TwoPAdaptionHelper : public ImplicitAdaptionHelper<TypeTag>
{
private:
    typedef ImplicitAdaptionHelper<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        //indices
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pwsn = Indices::pwsn,
        pnsw = Indices::pnsw,
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    typedef typename GridView::Grid Grid;
    typedef typename Grid::LevelGridView LevelGridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;

    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;

    struct AdaptedValues
    {
        std::vector<PrimaryVariables> u;
        int count;
        PrimaryVariables associatedMass;
        AdaptedValues(): associatedMass(0)
        {
            count = 0;
        }
    };

    typedef Dune::PersistentContainer<Grid, AdaptedValues> PersistentContainer;
    PersistentContainer adaptionMap_;

public:
    //! Constructs an adaption helper object
    /**
     *  @param gridView a DUNE gridview object
     */
    TwoPAdaptionHelper(Problem& problem) : ParentType(problem), adaptionMap_(problem.grid(), 0)
    {
        if(FluidSystem::isCompressible(wPhaseIdx) || FluidSystem::isCompressible(nPhaseIdx))
            DUNE_THROW(Dune::InvalidStateException, "Adaptionhelper is only for incompressible fluids mass-conservative!");
    }

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
        for (int level = problem.grid().maxLevel(); level >= 0; level--)
        {
            //get grid view on level grid
            LevelGridView levelView = problem.grid().levelGridView(level);

            for (const auto& element : elements(levelView))
            {
                //get your map entry
                AdaptedValues &adaptedValues = adaptionMap_[element];

                // put values in the map
                if (element.isLeaf())
                {
                    FVElementGeometry fvGeometry;
                    fvGeometry.update(problem.gridView(), element);

                    for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                    {
                        // get index
                        int dofIdx = this->dofIndex(problem, element, scvIdx);

                        adaptedValues.u.push_back(problem.model().curSol()[dofIdx]);

                        VolumeVariables volVars;
                        volVars.update(adaptedValues.u[scvIdx],
                                      problem,
                                      element,
                                      fvGeometry,
                                      scvIdx,
                                      false);

                        Scalar volume = fvGeometry.subContVol[scvIdx].volume;

                        adaptedValues.associatedMass[nPhaseIdx] += volume*volVars.density(nPhaseIdx)
                                            * volVars.porosity() * volVars.saturation(nPhaseIdx);
                        adaptedValues.associatedMass[wPhaseIdx] += volume*volVars.density(wPhaseIdx)
                                            * volVars.porosity() * volVars.saturation(wPhaseIdx);
                    }
                    adaptedValues.count = 1;
                }
                //Average in father
                if (element.level() > 0)
                {
                    if(!element.hasFather())
                        DUNE_THROW(Dune::InvalidStateException, "Element on level > 0 has no father element!");

                    AdaptedValues& adaptedValuesFather = adaptionMap_[element.father()];
                    //For some grids the father element is identical to the son element.
                    //For that case averaging is not necessary.
                    if(&adaptedValues != &adaptedValuesFather)
                    {
                        adaptedValuesFather.count += 1;
                        storeAdaptionValues(adaptedValues, adaptedValuesFather);
                    }
                }

                if(isBox && !element.isLeaf())
                {
                    FVElementGeometry fvGeometry;
                    fvGeometry.update(problem.gridView(), element);

                    for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                    {
                        int dofIdx = this->dofIndex(problem, element, scvIdx);

                        adaptedValues.u.push_back(problem.model().curSol()[dofIdx]);
                    }
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
        //vectors storing the mass associated with each vertex, when using the box method
        std::vector<Scalar> massCoeff;
        std::vector<Scalar> associatedMass;

        if(isBox)
        {
            massCoeff.resize(problem.model().numDofs(),0.0);
            associatedMass.resize(problem.model().numDofs(),0.0);
        }

        for (int level = 0; level <= problem.grid().maxLevel(); level++)
        {
            LevelGridView levelView = problem.grid().levelGridView(level);

            for (const auto& element : elements(levelView))
            {
                // only treat non-ghosts, ghost data is communicated afterwards
                if (element.partitionType() == Dune::GhostEntity)
                    continue;

                if (!element.isNew() || element.level() == 0)
                {
                    //entry is in map, write in leaf
                    if (element.isLeaf())
                    {
                        AdaptedValues &adaptedValues = adaptionMap_[element];

                        FVElementGeometry fvGeometry;
                        fvGeometry.update(problem.gridView(), element);

                        for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                        {
                            // get index
                            int dofIdx = this->dofIndex(problem, element, scvIdx);

                            VolumeVariables volVars;
                            volVars.update(adaptedValues.u[scvIdx],
                                          problem,
                                          element,
                                          fvGeometry,
                                          scvIdx,
                                          false);

                            Scalar volumeElement = fvGeometry.elementVolume;

                            this->setAdaptionValues(adaptedValues, problem.model().curSol()[dofIdx],scvIdx);

                            if (int(formulation) == pwsn)
                            {
                                problem.model().curSol()[dofIdx][saturationIdx] = adaptedValues.associatedMass[nPhaseIdx];
                                problem.model().curSol()[dofIdx][saturationIdx] /= volumeElement * volVars.density(nPhaseIdx) * volVars.porosity();
                            }
                            else if (int(formulation) == pnsw)
                            {
                                problem.model().curSol()[dofIdx][saturationIdx] = adaptedValues.associatedMass[wPhaseIdx];
                                problem.model().curSol()[dofIdx][saturationIdx] /= volumeElement * volVars.density(wPhaseIdx) * volVars.porosity();
                            }

                            if(isBox)
                            {
                                Scalar volume = fvGeometry.subContVol[scvIdx].volume;
                                if (int(formulation) == pwsn)
                                {
                                    massCoeff[dofIdx] += volume * volVars.density(nPhaseIdx) * volVars.porosity();
                                    associatedMass[dofIdx] += volume/volumeElement*adaptedValues.associatedMass[nPhaseIdx];
                                }
                                else if (int(formulation) == pnsw)
                                {
                                    massCoeff[dofIdx] += volume * volVars.density(wPhaseIdx) * volVars.porosity();
                                    associatedMass[dofIdx] += volume/volumeElement*adaptedValues.associatedMass[wPhaseIdx];
                                }
                            }

                        }

                    }
                }
                else
                {
                    // value is not in map, interpolate from father element
                    if (element.hasFather())
                    {
                        auto eFather = element.father();
                        while(eFather.isNew() && eFather.level() > 0)
                            eFather = eFather.father();

                        Scalar massFather = 0.0;

                        if(!isBox)
                        {
                            AdaptedValues& adaptedValuesFather = adaptionMap_[eFather];

                            if (int(formulation) == pwsn)
                            {
                                massFather = adaptedValuesFather.associatedMass[nPhaseIdx];
                            }
                            else if (int(formulation) == pnsw)
                            {
                                massFather = adaptedValuesFather.associatedMass[wPhaseIdx];
                            }

                            // access new son
                            AdaptedValues& adaptedValues = adaptionMap_[element];
                            adaptedValues.count = 1;

                            FVElementGeometry fvGeometry;
                            fvGeometry.update(problem.gridView(), element);

                            adaptedValues.u.push_back(adaptedValuesFather.u[0]);

                            VolumeVariables volVars;
                            volVars.update(adaptedValues.u[0],
                                          problem,
                                          element,
                                          fvGeometry,
                                          /*scvIdx=*/0,
                                          false);

                            Scalar volume = fvGeometry.subContVol[0].volume;
                            Scalar massCoeffSon = 0.0;
                            if (int(formulation) == pwsn)
                            {
                                massCoeffSon = volume * volVars.density(nPhaseIdx) * volVars.porosity();
                            }
                            else if (int(formulation) == pnsw)
                            {
                                massCoeffSon = volume * volVars.density(wPhaseIdx) * volVars.porosity();
                            }
                            Scalar volumeFather = eFather.geometry().volume();
                            adaptedValues.u[0][saturationIdx] = (volume/volumeFather*massFather)/massCoeffSon;

                            // if we are on leaf, store reconstructed values of son in CellData object
                            if (element.isLeaf())
                            {
                                // access new CellData object
                                int newIdxI = this->elementIndex(problem, element);

                                this->setAdaptionValues(adaptedValues, problem.model().curSol()[newIdxI],0);
                            }
                        }
                        else
                        {
                            AdaptedValues& adaptedValuesFather = adaptionMap_[eFather];
                            // access new son
                            AdaptedValues& adaptedValues = adaptionMap_[element];
                            adaptedValues.u.clear();
                            adaptedValues.count = 1;

                            const auto geometryI = element.geometry();

                            FVElementGeometry fvGeometry;
                            fvGeometry.update(problem.gridView(), element);

                            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                            {
                                auto subEntity = element.template subEntity <dofCodim>(scvIdx);

                                LocalPosition dofCenterPos = geometryI.local(subEntity.geometry().center());
                                const LocalFiniteElementCache feCache;
                                Dune::GeometryType geomType = eFather.geometry().type();

                                // Interpolate values from father element by using ansatz functions
                                const LocalFiniteElement &localFiniteElement = feCache.get(geomType);
                                std::vector<Dune::FieldVector<Scalar, 1> > shapeVal;
                                localFiniteElement.localBasis().evaluateFunction(dofCenterPos, shapeVal);
                                PrimaryVariables u(0);
                                for (int j = 0; j < shapeVal.size(); ++j)
                                {
                                    u.axpy(shapeVal[j], adaptedValuesFather.u[j]);
                                }

                                adaptedValues.u.push_back(u);
                                adaptedValues.count = 1;

                                if (element.isLeaf())
                                {
                                    VolumeVariables volVars;
                                    volVars.update(adaptedValues.u[scvIdx],
                                                  problem,
                                                  element,
                                                  fvGeometry,
                                                  scvIdx,
                                                  false);

                                    Scalar volume = fvGeometry.subContVol[scvIdx].volume;
                                    Scalar volumeFather = eFather.geometry().volume();

                                    int dofIdx = this->dofIndex(problem, element, scvIdx);
                                    if (int(formulation) == pwsn)
                                    {
                                        massCoeff[dofIdx] += volume * volVars.density(nPhaseIdx) * volVars.porosity();
                                        associatedMass[dofIdx] += volume/volumeFather*adaptedValuesFather.associatedMass[nPhaseIdx];
                                    }
                                    else if (int(formulation) == pnsw)
                                    {
                                        massCoeff[dofIdx] += volume * volVars.density(wPhaseIdx) * volVars.porosity();
                                        associatedMass[dofIdx] += volume/volumeFather*adaptedValuesFather.associatedMass[wPhaseIdx];
                                    }

                                    this->setAdaptionValues(adaptedValues, problem.model().curSol()[dofIdx],scvIdx);
                                }

                            }
                        }
                    }
                    else
                    {
                        DUNE_THROW(Dune::InvalidStateException, "Element is new but has no father element!");
                    }

                }
            }

        }

        if(isBox)
        {
            for(int dofIdx = 0; dofIdx < problem.model().numDofs(); dofIdx++)
            {
                problem.model().curSol()[dofIdx][saturationIdx] = associatedMass[dofIdx]/massCoeff[dofIdx];
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
            if(adaptedValuesFather.u.size() == 0)
                adaptedValuesFather.u.resize(1);

            adaptedValuesFather.u[0] += adaptedValues.u[0];
            adaptedValuesFather.u[0] /= adaptedValues.count;
            adaptedValuesFather.associatedMass += adaptedValues.associatedMass;
        }
        else
        {
            adaptedValuesFather.associatedMass += adaptedValues.associatedMass;
        }
    }
    //! Set adapted values in CellData
    /**
     * This methods stores reconstructed values into the cellData object, by
     * this setting a newly mapped solution to the storage container of the
     * sequential models.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param u The variables to be stored
     */
    static void setAdaptionValues(AdaptedValues& adaptedValues, PrimaryVariables& u, int scvIdx)
    {
        PrimaryVariables uNew(0);
        if(!isBox)
        {
            uNew = adaptedValues.u[scvIdx];
            uNew /= adaptedValues.count;
        }
        else
        {
            uNew = adaptedValues.u[scvIdx];
        }

        u = uNew;
    }
};
}
#endif
