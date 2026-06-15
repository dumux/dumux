// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LayerModels
 * \copydoc Dumux::LayerModelHirano
 */
#ifndef DUMUX_MATERIAL_LAYERMODEL_HIRANO_HH
#define DUMUX_MATERIAL_LAYERMODEL_HIRANO_HH

#include "layermodel.hh"

namespace Dumux {
/*!
 * \ingroup LayerModels
 * \brief Implementation of the Hirano layer model.
 */

template <class GridGeometry, class SpatialParams, class VolumeVariables, class SolutionVector, class GridVariables>
class LayerModelHirano : public LayerModel<GridGeometry, SpatialParams, VolumeVariables, SolutionVector, GridVariables>
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
    LayerModelHirano(std::shared_ptr<const GridGeometry> gridGeometry,
                     std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        maxActiveLayerThickness_ = 0.1;
        eps_ = 1e-14;

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
            int indexUppermostLayer = this->indexUppermostLayer(element);
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
                        DUNE_THROW(Dune::InvalidStateException, "Negative sediment mass ("+std::to_string(curSol[elemId][i])+" kg) in the active layer (grain class "+std::to_string(i+1)
                                                                 + ") at position x="+std::to_string(scv.center()[0])+", y="+std::to_string(scv.center()[1]));
                    }
                }
                if (layerThickness<maxActiveLayerThickness_) {
                    // **** erosion ****
                    if (indexUppermostLayer != 0){
                        // there is at least one sublayer left
                        // -> prepare the sediment exchange between the first sublayer and the active layer
                        Scalar diff = maxActiveLayerThickness_ - layerThickness;
                        Scalar massDiff = diff * element.geometry().volume() * (1 - this->spatialParams_->porosity())
                                          * this->harmonicAverageGrainDensities_[elemId][indexUppermostLayer];
                        Scalar massFirstSubLayer = std::accumulate(this->sedimentMasses_[elemId][indexUppermostLayer].begin(),
                                                                   this->sedimentMasses_[elemId][indexUppermostLayer].end(), 0.0);
                        // Check if the first sublayer vanishes. The while loop is necessary, because theoreticaly
                        // its possible, that more than one layer vanishes in one time step
                        while (massDiff > massFirstSubLayer) {
                            // first sublayer vanishes
                            // -> add the rest of the first sublayer to the active layer
                            for (int i=0; i<this->nGrainClasses_; i++) {
                                curSol[elemId][i] += this->sedimentMasses_[elemId][indexUppermostLayer][i];
                            }
                            massDiff -= massFirstSubLayer;
                            diff -= massFirstSubLayer / element.geometry().volume() / (1 - this->spatialParams_->porosity())
                                   / this->harmonicAverageGrainDensities_[elemId][indexUppermostLayer];
                            // remove the (now empty) first sublayer
                            this->sedimentMasses_[elemId].erase(indexUppermostLayer);
                            this->upperLayerBoundaries_[elemId].erase(indexUppermostLayer);
                            // get the index of the new first sublayer
                            indexUppermostLayer = this->indexUppermostLayer(element);
                            if (indexUppermostLayer == 0){
                                // there is no layer left, but the active layer -> no further update necessary
                                break;
                            }
                            // calculate the mass of the new first sublayer
                            massFirstSubLayer = std::accumulate(this->sedimentMasses_[elemId][indexUppermostLayer].begin(),
                                                                this->sedimentMasses_[elemId][indexUppermostLayer].end(), 0.0);
                        }
                        if (indexUppermostLayer != 0){
                        // there is still at least one sublayer left
                        // -> exchange sediment between the first sublayer and the active layer
                            // Modify the uppermost sub layer boundary
                            this->upperLayerBoundaries_[elemId][indexUppermostLayer] -= diff;

                            Scalar massFirstSubLayer = std::accumulate(this->sedimentMasses_[elemId][indexUppermostLayer].begin(),
                                                                       this->sedimentMasses_[elemId][indexUppermostLayer].end(), 0.0);
                            for (int i=0; i<this->nGrainClasses_; i++) {
                                Scalar massChange = this->sedimentMasses_[elemId][indexUppermostLayer][i]/massFirstSubLayer * massDiff;
                                // Modify mass in active layer (primary variable)
                                curSol[elemId][i] += massChange;
                                // Modify the mass in the uppermost sub layer
                                this->sedimentMasses_[elemId][indexUppermostLayer][i] -= massChange;
                                // The average grain density does not change in case of erosion, since the mass fraction in the uppermost sublayer remain the same
                            }
                        }
                    }
                    else {
                        // There is only the active layer left and it does not exceed the maximum layer thickness.
                        // Therefore, the primary variables does not have to be modified.
                    }
                }
                else if (layerThickness>maxActiveLayerThickness_) {
                    // **** deposition ****
                    if (indexUppermostLayer == 0){
                        // there is no layer left, therefore we have to create one to put the surplus sediment into
                        this->upperLayerBoundaries_[elemId][this->nLayer_] = this->fixedGroundLevel_[elemId];
                        this->sedimentMasses_[elemId][this->nLayer_] = std::vector<Scalar>(this->nGrainClasses_, 0.0);
                        indexUppermostLayer = this->nLayer_;
                    }

                    Scalar diff = layerThickness - maxActiveLayerThickness_;
                    Scalar massDiff = diff * element.geometry().volume() * (1 - this->spatialParams_->porosity()) * elemVolVars[scv].harmonicAverageGrainDensity();

                    // Modify the uppermost sub layer boundary
                    this->upperLayerBoundaries_[elemId][indexUppermostLayer] += diff;

                    Scalar massActiveLayer = std::accumulate(curSol[elemId].begin(), curSol[elemId].end(), 0.0);
                    for (int i=0; i<this->nGrainClasses_; i++) {
                        Scalar massChange = curSol[elemId][i] / massActiveLayer * massDiff;
                        // Modify mass in active layer (primary variable)
                        curSol[elemId][i] -= massChange;
                        // Modify the mass in the uppermost sub layer
                        this->sedimentMasses_[elemId][indexUppermostLayer][i] += massChange;
                    }
                    // Recalculate the average grain density
                    Scalar temp = 0.0;
                    Scalar massFirstSubLayer = std::accumulate(this->sedimentMasses_[elemId][indexUppermostLayer].begin(),
                                                               this->sedimentMasses_[elemId][indexUppermostLayer].end(), 0.0);
                    for (int i=0; i<this->nGrainClasses_; i++) {
                        temp += this->sedimentMasses_[elemId][indexUppermostLayer][i] / massFirstSubLayer / this->spatialParams_->grainDensity(i);
                    }
                    this->harmonicAverageGrainDensities_[elemId][indexUppermostLayer] = 1/temp;
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
        this->initialMassActiveLayer_.resize(nElems, 0.0);
        for (const auto& element : elements(this->gridGeometry_->gridView()))
        {
            auto fvGeometry = localView(*this->gridGeometry_);
            fvGeometry.bindElement(element);
            auto elemId = this->gridGeometry_->elementMapper().index(element);
            bool settingActiveLayer = true;
            Scalar missingActiveLayerThickness;
            if (this->initialBedSurface_[elemId]-this->fixedGroundLevel_[elemId] > maxActiveLayerThickness_) {
                this->initialActiveLayerThickness_[elemId] = maxActiveLayerThickness_;
                missingActiveLayerThickness = maxActiveLayerThickness_;
            }
            else {
                this->initialActiveLayerThickness_[elemId] = this->initialBedSurface_[elemId]-this->fixedGroundLevel_[elemId];
                missingActiveLayerThickness = this->initialBedSurface_[elemId]-this->fixedGroundLevel_[elemId];
            }
            for (int layerId=1; layerId<=this->nLayer_; layerId++) {
                Scalar layerThickness;
                if (layerId < this->nLayer_) {
                    layerThickness = this->upperLayerBoundaries_[elemId][layerId] - this->upperLayerBoundaries_[elemId][layerId+1];
                }
                else if(layerId == this->nLayer_) {
                    layerThickness = this->upperLayerBoundaries_[elemId][layerId] - this->fixedGroundLevel_[elemId];
                }
                if (settingActiveLayer) { // The active layer is not yet initialized
                    if (layerThickness > missingActiveLayerThickness+eps_) { // The current layer has enough sediment to build the active layer
                        settingActiveLayer = false;
                        this->upperLayerBoundaries_[elemId][layerId] -= missingActiveLayerThickness;
                        layerThickness -= missingActiveLayerThickness;
                        this->initialMassActiveLayer_[elemId] += missingActiveLayerThickness * element.geometry().volume() * (1 - this->spatialParams_->porosity())
                                                                    * this->harmonicAverageGrainDensities_[elemId][layerId];
                        for (int grainClassId=0; grainClassId<this->nGrainClasses_; grainClassId++) {
                            this->initialMassFractionsActiveLayer_[elemId][grainClassId] += this->initialMassFractions_[elemId][layerId][grainClassId]
                                                            // scale the mass fractions with the layer thickness
                                                            * missingActiveLayerThickness / this->initialActiveLayerThickness_[elemId];
                        }
                        // Set the (remaining) sediment mass of the current layer
                        Scalar sedimentMassLayer = layerThickness * element.geometry().volume()
                                        * (1 - this->spatialParams_->porosity()) * this->harmonicAverageGrainDensities_[elemId][layerId];
                        this->sedimentMasses_[elemId][layerId].resize(this->nGrainClasses_);
                        for (int grainClassId=0; grainClassId<this->nGrainClasses_; grainClassId++) {
                            this->sedimentMasses_[elemId][layerId][grainClassId] = sedimentMassLayer * this->initialMassFractions_[elemId][layerId][grainClassId];
                        }
                    }
                    else { // The current layer has exactly as much sediment as needed or less to build the active layer
                        // move the sediment of the current layer to the active layer
                        missingActiveLayerThickness -= layerThickness;
                        this->initialMassActiveLayer_[elemId] += layerThickness * element.geometry().volume() * (1 - this->spatialParams_->porosity())
                                                                    * this->harmonicAverageGrainDensities_[elemId][layerId];
                        for (int grainClassId=0; grainClassId<this->nGrainClasses_; grainClassId++) {
                            this->initialMassFractionsActiveLayer_[elemId][grainClassId] += this->initialMassFractions_[elemId][layerId][grainClassId]
                                                            // scale the mass fractions with the layer thickness
                                                            * layerThickness / this->initialActiveLayerThickness_[elemId];
                        }
                        // remove the layer
                        this->upperLayerBoundaries_[elemId].erase(layerId);
                        // Since the layer vanished it is not necessary to calculate the sediment mass
                    }
                }
                else { // The active layer is already initialized
                    // Set the sediment mass of the current layer
                    Scalar sedimentMassLayer = layerThickness * element.geometry().volume()
                                    * (1 - this->spatialParams_->porosity()) * this->harmonicAverageGrainDensities_[elemId][layerId];
                    this->sedimentMasses_[elemId][layerId].resize(this->nGrainClasses_);
                    for (int grainClassId=0; grainClassId<this->nGrainClasses_; grainClassId++) {
                        this->sedimentMasses_[elemId][layerId][grainClassId] = sedimentMassLayer * this->initialMassFractions_[elemId][layerId][grainClassId];
                    }
                }
            }
        }
    }

    Scalar eps_;
    Scalar maxActiveLayerThickness_;
};
} // end namespace Dumux

#endif // DUMUX_MATERIAL_LAYERMODEL_HIRANO_HH
