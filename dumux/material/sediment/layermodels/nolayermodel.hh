// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LayerModels
 * \copydoc Dumux::NoLayerModel
 */
#ifndef DUMUX_MATERIAL_NOLAYERMODEL_HH
#define DUMUX_MATERIAL_NOLAYERMODEL_HH

#include "layermodel.hh"

namespace Dumux {
/*!
 * \ingroup LayerModels
 * \brief Implementation of the default layer model, which doesn't use different layers.
 */

template <class GridGeometry, class SpatialParams, class VolumeVariables, class SolutionVector, class GridVariables>
class NoLayerModel : public LayerModel<GridGeometry, SpatialParams, VolumeVariables, SolutionVector, GridVariables>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
    using ParentType = LayerModel<GridGeometry, SpatialParams, VolumeVariables, SolutionVector, GridVariables>;
public:
    /*!
    * \brief Fill the layer model with data
    *
    * \param gridGeometry The grid geometry
    * \param spatialParams The spatial parameters contain the initial data of the layer model for each element.
      *                    Theses are on the one hand "bedSurface" and "fixedGroundLevel". On the other hand the upper limit
    *                      of each layer ("upperLimitLayer1", "upperLimitLayer2"...) and the mass fractions for each layer
    *                      and grain class ("massFractionLayer1GrainClass1", "massFractionLayer1GrainClass2", ...).
    *                      The length of the corresponding vector is the number of elements (per rank).
    */
    NoLayerModel(std::shared_ptr<const GridGeometry> gridGeometry,
                 std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        if (this->nLayer_ != 1) {
            DUNE_THROW(Dune::InvalidStateException, "'Sediment.NumberLayers' in params.input is set to '" + std::to_string(this->nLayer_)
                                                    + "', which is invalid, since 'NoLayerModel' is used!"
                                                    + "'Sediment.NumberLayers' must be set to 1, when 'NoLayerModel' is used!");
        }

        createActiveLayer_();
    }

    //! Update the layer model
    void update(SolutionVector& curSol,
                const GridVariables& gridVariables)
    {
        Dune::MPIGuard guard(this->gridGeometry_->gridView().comm());
        bool correct = true;
        for (const auto& element : elements(this->gridGeometry_->gridView()))
        {
            auto fvGeometry = localView(*this->gridGeometry_);
            fvGeometry.bindElement(element);
            auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);
            auto elemId = this->gridGeometry_->elementMapper().index(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                Scalar layerThickness = elemVolVars[scv].layerThickness();
                // Check for negative layer thicknesses
                if (layerThickness < 0.0) {
                    std::cout<<"Error: Negative layer thickness at position x="+std::to_string(scv.center()[0])+", y="+std::to_string(scv.center()[1])
                             <<". Possible reasons are a to small choice of the parameters `FluxLimiterLET.LowerSedimentThickness` and `FluxLimiterLET.UpperSedimentThickness`"
                               " or to much bedload transport."<<std::endl;
                    DUNE_THROW(Dune::InvalidStateException, "Negative layer thickness at position x="+std::to_string(scv.center()[0])+", y="+std::to_string(scv.center()[1])+
                                                            ". Possible reasons are a to small choice of the parameters `FluxLimiterLET.LowerSedimentThickness` and "
                                                            "`FluxLimiterLET.UpperSedimentThickness` or to much bedload transport.");
                }
                // Check for negative sediment masses
                for (int i=0; i<this->nGrainClasses_; i++) {
                    if (curSol[elemId][i] < 0) {
                        std::cout<<"Error: Negative sediment mass ("+std::to_string(curSol[elemId][i])+" kg) in the active layer (grain class "+std::to_string(i+1)
                                   + ") at position x="+std::to_string(scv.center()[0])+", y="+std::to_string(scv.center()[1])<<std::endl;
                        DUNE_THROW(Dune::InvalidStateException, "Negative sediment mass ("+std::to_string(curSol[elemId][i])+" kg) int the active layer (grain class "+std::to_string(i+1)
                                                                 + ") at position x="+std::to_string(scv.center()[0])+", y="+std::to_string(scv.center()[1]));
                    }
                }
            }
        }
        guard.finalize(correct);
    }

private:
    /*!
    * \brief Create the active layer
    *
    * Initialize the following data member:
    *   - initialMassActiveLayer_
    *   - initialMassFractionsActiveLayer_
    *   - initialActiveLayerThickness_
    *   - sedimentMasses_
    */
    void createActiveLayer_() const
    {
        int nElems = this->gridGeometry_->gridView().size(0);
        this->sedimentMasses_.resize(nElems);
        this->initialActiveLayerThickness_.resize(nElems);
        this->initialMassFractionsActiveLayer_.resize(nElems, std::vector<Scalar>(this->nGrainClasses_, 0.0));
        this->initialMassActiveLayer_.resize(nElems);
        for (const auto& element : elements(this->gridGeometry_->gridView()))
        {
            auto fvGeometry = localView(*this->gridGeometry_);
            fvGeometry.bindElement(element);
            auto elemId = this->gridGeometry_->elementMapper().index(element);
            this->initialActiveLayerThickness_[elemId] = this->initialBedSurface_[elemId] - this->fixedGroundLevel_[elemId];
            for (auto&& scv : scvs(fvGeometry))
            {
                // We don't have to loop over the layer since there is just one layer
                int layerId = 1;
                this->initialMassActiveLayer_[elemId] = this->initialActiveLayerThickness_[elemId] * element.geometry().volume()
                               * (1 - this->spatialParams_->porosity()) * this->harmonicAverageGrainDensities_[elemId][layerId];
                for (int grainClassId=0; grainClassId<this->nGrainClasses_; grainClassId++) {
                    this->initialMassFractionsActiveLayer_[elemId][grainClassId] = this->initialMassFractions_[elemId][layerId][grainClassId];
                }
                this->upperLayerBoundaries_[elemId].erase(layerId);
            }
        }
    }
};
} // end namespace Dumux

#endif // DUMUX_MATERIAL_NOLAYERMODEL_HH
