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

#include <dumux/common/properties.hh>
#include <dumux/implicit/adaptive/adaptionhelper.hh>

namespace Dumux {

/*!
 * \brief Base class holding the variables for implicit models.
 */
template<class TypeTag>
class TwoPAdaptionHelper : public ImplicitAdaptionHelper<TypeTag>
{
private:
    using ParentType = ImplicitAdaptionHelper<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

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

    using Grid = typename GridView::Grid;
    using Element = typename GridView::Traits::template Codim<0>::Entity;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using CoordScalar = typename GridView::ctype;
    using LocalPosition = Dune::FieldVector<CoordScalar, dim>;

    using LocalFiniteElementCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using LocalFiniteElement = typename LocalFiniteElementCache::FiniteElementType;

    struct AdaptedValues
    {
        ElementSolutionVector u;
        int count;
        PrimaryVariables associatedMass;
        AdaptedValues(): count(0), associatedMass(0.0) {}
    };

    using PersistentContainer = Dune::PersistentContainer<Grid, AdaptedValues>;
    PersistentContainer adaptionMap_;

public:
    //! Constructs an adaption helper object
    /**
     *  @param gridView a DUNE gridview object
     */
    TwoPAdaptionHelper(Problem& problem) : ParentType(problem), adaptionMap_(problem.grid(), 0)
    {
        if(FluidSystem::isCompressible(wPhaseIdx) || FluidSystem::isCompressible(nPhaseIdx))
            DUNE_THROW(Dune::InvalidStateException, "This adaption helper is only mass conservative for incompressible fluids!");
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
            // get grid view on level grid
            auto levelView = problem.grid().levelGridView(level);

            for (const auto& element : elements(levelView))
            {
                //get your map entry
                auto& adaptedValues = adaptionMap_[element];

                // put values in the map
                if (element.isLeaf())
                {
                    auto fvGeometry = localView(problem.model().globalFvGeometry());
                    fvGeometry.bindElement(element);

                    auto elemVolVars = localView(problem.model().curGlobalVolVars());
                    elemVolVars.bindElement(element, fvGeometry, problem.model().curSol());

                    for (auto&& scv : scvs(fvGeometry))
                    {
                        adaptedValues.u = problem.model().elementSolution(element, problem.model().curSol());

                        VolumeVariables volVars;
                        volVars.update(adaptedValues.u, problem, element, scv);

                        adaptedValues.associatedMass[nPhaseIdx] += scv.volume() * volVars.density(nPhaseIdx)
                                                                   * volVars.porosity() * volVars.saturation(nPhaseIdx);
                        adaptedValues.associatedMass[wPhaseIdx] += scv.volume() * volVars.density(wPhaseIdx)
                                                                   * volVars.porosity() * volVars.saturation(wPhaseIdx);
                    }
                    adaptedValues.count = 1;
                }
                // Average in father
                if (element.level() > 0)
                {
                    auto& adaptedValuesFather = adaptionMap_[element.father()];
                    // For some grids the father element is identical to the son element.
                    // For that case averaging is not necessary.
                    if(&adaptedValues != &adaptedValuesFather)
                    {
                        adaptedValuesFather.count += 1;
                        storeAdaptionValues(adaptedValues, adaptedValuesFather);
                    }
                }

                if(isBox && !element.isLeaf())
                {
                    adaptedValues.u = problem.model().elementSolution(element, problem.model().curSol());
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
        // vectors storing the mass associated with each vertex, when using the box method
        std::vector<Scalar> massCoeff;
        std::vector<Scalar> associatedMass;

        if(isBox)
        {
            massCoeff.resize(problem.model().numDofs(), 0.0);
            associatedMass.resize(problem.model().numDofs(), 0.0);
        }

        for (int level = 0; level <= problem.grid().maxLevel(); level++)
        {
            auto levelView = problem.grid().levelGridView(level);

            for (const auto& element : elements(levelView))
            {
                // only treat non-ghosts, ghost data is communicated afterwards
                if (element.partitionType() == Dune::GhostEntity)
                    continue;

                if (!element.isNew() || element.level() == 0)
                {
                    // entry is in map, write in leaf
                    if (element.isLeaf())
                    {
                        auto& adaptedValues = adaptionMap_[element];

                        auto fvGeometry = localView(problem.model().globalFvGeometry());
                        fvGeometry.bindElement(element);

                        auto elementVolume = element.geometry().volume();

                        for (auto&& scv : scvs(fvGeometry))
                        {
                            VolumeVariables volVars;
                            volVars.update(adaptedValues.u, problem, element, scv);

                            problem.model().curSol()[scv.dofIndex()] = getPriVars(adaptedValues, scv);

                            if (formulation == pwsn)
                            {
                                problem.model().curSol()[scv.dofIndex()][saturationIdx] = adaptedValues.associatedMass[nPhaseIdx];
                                problem.model().curSol()[scv.dofIndex()][saturationIdx] /= elementVolume * volVars.density(nPhaseIdx) * volVars.porosity();
                            }
                            else if (formulation == pnsw)
                            {
                                problem.model().curSol()[scv.dofIndex()][saturationIdx] = adaptedValues.associatedMass[wPhaseIdx];
                                problem.model().curSol()[scv.dofIndex()][saturationIdx] /= elementVolume * volVars.density(wPhaseIdx) * volVars.porosity();
                            }

                            if(isBox)
                            {
                                if (int(formulation) == pwsn)
                                {
                                    massCoeff[scv.dofIndex()] += scv.volume() * volVars.density(nPhaseIdx) * volVars.porosity();
                                    associatedMass[scv.dofIndex()] += scv.volume() / elementVolume * adaptedValues.associatedMass[nPhaseIdx];
                                }
                                else if (int(formulation) == pnsw)
                                {
                                    massCoeff[scv.dofIndex()] += scv.volume() * volVars.density(wPhaseIdx) * volVars.porosity();
                                    associatedMass[scv.dofIndex()] += scv.volume() / elementVolume * adaptedValues.associatedMass[wPhaseIdx];
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
                        auto fatherElement = element.father();
                        while(fatherElement.isNew() && fatherElement.level() > 0)
                            fatherElement = fatherElement.father();

                        Scalar massFather = 0.0;

                        if(!isBox)
                        {
                            auto& adaptedValuesFather = adaptionMap_[fatherElement];

                            if (formulation == pwsn)
                                massFather = adaptedValuesFather.associatedMass[nPhaseIdx];

                            else if (formulation == pnsw)
                                massFather = adaptedValuesFather.associatedMass[wPhaseIdx];

                            // access new son
                            auto& adaptedValues = adaptionMap_[element];
                            adaptedValues.count = 1;

                            auto fvGeometry = localView(problem.model().globalFvGeometry());
                            fvGeometry.bindElement(element);

                            adaptedValues.u = adaptedValuesFather.u;

                            for (auto&& scv : scvs(fvGeometry))
                            {
                                VolumeVariables volVars;
                                volVars.update(adaptedValues.u, problem, element, scv);

                                Scalar massCoeffSon = 0.0;
                                if (int(formulation) == pwsn)
                                    massCoeffSon = scv.volume() * volVars.density(nPhaseIdx) * volVars.porosity();

                                else if (int(formulation) == pnsw)
                                    massCoeffSon = scv.volume() * volVars.density(wPhaseIdx) * volVars.porosity();

                                auto fatherElementVolume = fatherElement.geometry().volume();
                                adaptedValues.u[0][saturationIdx] = (scv.volume() / fatherElementVolume * massFather)/massCoeffSon;
                            }

                            // if we are on leaf, store reconstructed values of son in CellData object
                            if (element.isLeaf())
                            {
                                // access new CellData object
                                auto newIdxI = problem.elementMapper().index(element);
                                problem.model().curSol()[newIdxI] = adaptedValues.u[0];
                            }
                        }
                        else
                        {
                            // auto& adaptedValuesFather = adaptionMap_[fatherElement];
                            // // access new son
                            // auto& adaptedValues = adaptionMap_[element];
                            // adaptedValues.u.clear();
                            // adaptedValues.count = 1;
                            //
                            // const auto geometryI = element.geometry();
                            //
                            // FVElementGeometry fvGeometry;
                            // fvGeometry.update(problem.gridView(), element);
                            //
                            // for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                            // {
                            //     auto subEntity = element.template subEntity <dofCodim>(scvIdx);
                            //
                            //     LocalPosition dofCenterPos = geometryI.local(subEntity.geometry().center());
                            //     const LocalFiniteElementCache feCache;
                            //     Dune::GeometryType geomType = fatherElement.geometry().type();
                            //
                            //     // Interpolate values from father element by using ansatz functions
                            //     const LocalFiniteElement &localFiniteElement = feCache.get(geomType);
                            //     std::vector<Dune::FieldVector<Scalar, 1> > shapeVal;
                            //     localFiniteElement.localBasis().evaluateFunction(dofCenterPos, shapeVal);
                            //     PrimaryVariables u(0);
                            //     for (int j = 0; j < shapeVal.size(); ++j)
                            //     {
                            //         u.axpy(shapeVal[j], adaptedValuesFather.u[j]);
                            //     }
                            //
                            //     adaptedValues.u.push_back(u);
                            //     adaptedValues.count = 1;
                            //
                            //     if (element.isLeaf())
                            //     {
                            //         VolumeVariables volVars;
                            //         volVars.update(adaptedValues.u[scvIdx],
                            //                       problem,
                            //                       element,
                            //                       fvGeometry,
                            //                       scvIdx,
                            //                       false);
                            //
                            //         Scalar volume = fvGeometry.subContVol[scvIdx].volume;
                            //         Scalar fatherElementVolume = fatherElement.geometry().volume();
                            //
                            //         int dofIdxGlobal = this->dofIndex(problem, element, scvIdx);
                            //         if (int(formulation) == pwsn)
                            //         {
                            //             massCoeff[dofIdxGlobal] += volume * volVars.density(nPhaseIdx) * volVars.porosity();
                            //             associatedMass[dofIdxGlobal] += volume/fatherElementVolume*adaptedValuesFather.associatedMass[nPhaseIdx];
                            //         }
                            //         else if (int(formulation) == pnsw)
                            //         {
                            //             massCoeff[dofIdxGlobal] += volume * volVars.density(wPhaseIdx) * volVars.porosity();
                            //             associatedMass[dofIdxGlobal] += volume/fatherElementVolume*adaptedValuesFather.associatedMass[wPhaseIdx];
                            //         }
                            //
                            //         this->setAdaptionValues(adaptedValues, problem.model().curSol()[dofIdxGlobal],scvIdx);
                            //     }
                            //
                            // }
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
            for(int dofIdxGlobal = 0; dofIdxGlobal < problem.model().numDofs(); dofIdxGlobal++)
                problem.model().curSol()[dofIdxGlobal][saturationIdx] = associatedMass[dofIdxGlobal] / massCoeff[dofIdxGlobal];
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
    static PrimaryVariables getPriVars(const AdaptedValues& adaptedValues, const SubControlVolume& scv)
    {
        if(!isBox)
        {
            auto uNew = adaptedValues.u[0];
            uNew /= adaptedValues.count;
            return uNew;
        }
        else
        {
            return adaptedValues.u[scv.index()];
        }
    }
};

} // end namespace Dumux

#endif
