// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LayerModels
 * \copydoc Dumux::LayerModel
 */
#ifndef DUMUX_MATERIAL_LAYERMODEL_HH
#define DUMUX_MATERIAL_LAYERMODEL_HH

namespace Dumux {
/*!
 * \ingroup LayerModels
 * \brief Abstract base class for layer models.
 */

template <class GridGeometry, class SpatialParams, class VolumeVariables, class SolutionVector, class GridVariables>
class LayerModel
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
    using VertexMapper = typename GridGeometry::VertexMapper;
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
    LayerModel(std::shared_ptr<const GridGeometry> gridGeometry,
               std::shared_ptr<SpatialParams> spatialParams)
    {
        // read the params
        nLayer_ = Dumux::getParam<int>("Sediment.NumberLayers");
        nGrainClasses_ = Dumux::getParam<int>("Sediment.NumberGrainClasses");
        gridGeometry_ = gridGeometry;
        spatialParams_ = spatialParams;

        Dune::MPIGuard guard(gridGeometry_->gridView().comm());
        bool correct = true;

        // prepare data input
        int nElems = gridGeometry_->gridView().size(0);
        initialBedSurface_.resize(nElems);
        fixedGroundLevel_.resize(nElems);

        for (int layerId=1; layerId<=nLayer_; layerId++) {
            if (layerId>1) {
                auto varName = "upperLimitLayer" + std::to_string(layerId);
                elementdata_[varName] = std::vector<Scalar>(nElems, 0.0);
            }
            for (int grainClassId=0; grainClassId<nGrainClasses_; grainClassId++) {
                auto varName = "massFractionLayer" + std::to_string(layerId) + "GrainClass" + std::to_string(grainClassId+1);
                elementdata_[varName] = std::vector<Scalar>(nElems, 0.0);
            }
        }

        // get the data input
        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            auto fvGeometry = localView(*gridGeometry_);
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto eIdx = scv.elementIndex();
                initialBedSurface_[eIdx] = spatialParams_->initialBedSurface(element, scv);
                fixedGroundLevel_[eIdx] = spatialParams_->fixedGroundLevel(element, scv);
                for (int layerId=1; layerId<=nLayer_; layerId++) {
                    if (layerId>1) {
                        auto varName = "upperLimitLayer" + std::to_string(layerId);
                        elementdata_[varName][eIdx] = spatialParams_->initialUpperLayerLimit(element, scv, layerId);
                    }
                    for (int grainClassId=0; grainClassId<nGrainClasses_; grainClassId++) {
                        auto varName = "massFractionLayer" + std::to_string(layerId) + "GrainClass" + std::to_string(grainClassId+1);
                        elementdata_[varName][eIdx] = spatialParams_->initialMassFraction(element, scv, layerId, grainClassId);
                    }
                }
            }
        }

        // fill the layer model with the input data
        upperLayerBoundaries_.resize(nElems);
        harmonicAverageGrainDensities_.resize(nElems);
        initialMassFractions_.resize(nElems);
        for (int layerId=1; layerId<=nLayer_; layerId++) {
            // initialize upperLayerBoundaries_
            if (layerId==1) {
                for (int elemId = 0; elemId<nElems; elemId++) {
                    upperLayerBoundaries_[elemId][layerId] = initialBedSurface_[elemId];
                }
            }
            else {
                auto varName = "upperLimitLayer" + std::to_string(layerId);
                if(elementdata_.find(varName) != elementdata_.end()){
                    for (int elemId = 0; elemId<nElems; elemId++) {
                        upperLayerBoundaries_[elemId][layerId] = elementdata_.at(varName)[elemId];
                    }
                }
                else {
                    DUNE_THROW(Dune::InvalidStateException, "Can not find " + varName + " in input data");
                }
            }
            // initialize initialMassFractions_
            for (int elemId = 0; elemId<nElems; elemId++) {
                initialMassFractions_[elemId][layerId].resize(nGrainClasses_);
            }
            // check if massFraction data exist
            for (int grainClassId=0; grainClassId<nGrainClasses_; grainClassId++) {
                auto varName = "massFractionLayer" + std::to_string(layerId) + "GrainClass" + std::to_string(grainClassId+1);
                if (elementdata_.find(varName) == elementdata_.end()) { DUNE_THROW(Dune::InvalidStateException, "Can not find " + varName + " in input data"); }
            }
            // fill initialMassFractions_ and harmonicAverageGrainDensities_ with values
            for (int elemId=0; elemId<nElems; elemId++) {
                Scalar temp = 0.0;
                for (int grainClassId=0; grainClassId<nGrainClasses_; grainClassId++) {
                    auto varName = "massFractionLayer" + std::to_string(layerId) + "GrainClass" + std::to_string(grainClassId+1);
                    initialMassFractions_[elemId][layerId][grainClassId] = elementdata_.at(varName)[elemId];
                    temp += initialMassFractions_[elemId][layerId][grainClassId] / spatialParams_->grainDensity(grainClassId);
                }
                harmonicAverageGrainDensities_[elemId][layerId] = 1/temp;
            }
        }

        // check the input
        for (int elemId=0; elemId<nElems; elemId++) {
            for (int layerId=1; layerId<=nLayer_; layerId++) {
                // check massfractions sum
                Scalar massFractionSum = 0.0;
                for (int grainClassId=0; grainClassId<nGrainClasses_; grainClassId++) {
                    massFractionSum += initialMassFractions_[elemId][layerId][grainClassId];
                }
                if (massFractionSum > 1.0+eps_ || massFractionSum < 1.0-eps_) {
                    DUNE_THROW(Dune::InvalidStateException, "Invalid input data: Sum of mass fractions in layer "
                                                             + std::to_string(layerId) + " is different to 1.0");
                }
                // check layer order
                if (layerId == 2) {
                    if (upperLayerBoundaries_[elemId][layerId] > initialBedSurface_[elemId]) {
                        DUNE_THROW(Dune::InvalidStateException, "Invalid input data: 'upperLimitLayer" + std::to_string(layerId) +
                                                                "' is higher than 'bedSurface'");
                    }
                }
                else if (layerId > 2 && upperLayerBoundaries_[elemId][layerId] > upperLayerBoundaries_[elemId][layerId-1]) {
                    DUNE_THROW(Dune::InvalidStateException, "Invalid input data: 'upperLimitLayer" + std::to_string(layerId) +
                                                        "' is higher than 'upperLimitLayer" + std::to_string(layerId-1) + "'");
                }
            }
            if (fixedGroundLevel_[elemId] > initialBedSurface_[elemId]) {
                    DUNE_THROW(Dune::InvalidStateException, "Invalid input data: 'fixedGroundLevel' is higher than 'bedSurface'");
            }
            if (fixedGroundLevel_[elemId] > upperLayerBoundaries_[elemId][nLayer_]) {
                    DUNE_THROW(Dune::InvalidStateException, "Invalid input data: 'fixedGroundLevel'"
                                                            " is higher than 'upperLimitLayer" + std::to_string(nLayer_) + "'");
            }
            // check distance between bedSurface and fixedGroundLevel
            if (initialBedSurface_[elemId] - fixedGroundLevel_[elemId] < 0.001-eps_) {
                    DUNE_THROW(Dune::InvalidStateException, "Invalid input data: The minimal layer thickness of 1e-3 m of the"
                                                            " erodible layer is undercut.");
            }
        }
        guard.finalize(correct);

    }

    //! Update the layer model
    virtual void update(SolutionVector& curSol, const GridVariables& gridVariables) = 0;

    /*!
    * \brief Provide the initial sediment mass of a specific grain class for the active layer
    *
    * This value is used as primary variable within the simulation.
    *
    * \param element
    * \param grainClassId Index of the grain class. Starts with zero.
    */
    Scalar initialSedimentMassActiveLayer(Element element, const int grainClassId) const
    {
        auto elemId = gridGeometry_->elementMapper().index(element);
        return initialMassActiveLayer_[elemId] * initialMassFractionsActiveLayer_[elemId][grainClassId];
    }

    /*!
    * \brief Return the bottom of the active layer for a specific element.
    *
    * \param element
    */
    Scalar bottomActiveLayer(const Element& element) const
    {
        auto elemId = gridGeometry_->elementMapper().index(element);

        if (indexUppermostLayer(element)==0) {
            return fixedGroundLevel_[elemId];
        }
        else {
            return upperLayerBoundaries_[elemId][indexUppermostLayer(element)];
        }
    }

    /*!
    * \brief Return the fixed ground level for a specific element.
    *
    * \param element
    */
    Scalar fixedGroundLevel(const Element& element) const
    {
        auto elemId = gridGeometry_->elementMapper().index(element);
        return fixedGroundLevel_[elemId];
    }

    /*!
    * \brief Return the sediment masses for a given element
    *
    * \param element
    */
    std::map<int, std::vector<Scalar>>& sedimentMasses(const Element& element) const
    {
        auto elemId = gridGeometry_->elementMapper().index(element);
        return sedimentMasses_[elemId];
    }

    /*!
    * \brief Return the upper layer boundaries of a given element
    *
    * \param element
    */
    std::map<int, Scalar>& upperLayerBoundaries(const Element& element) const
    {
        auto elemId = gridGeometry_->elementMapper().index(element);
        return upperLayerBoundaries_[elemId];
    }

    /*!
    * \brief Return the initial mass fractions of the active layer
    */
    Scalar &initialMassFractionsActiveLayer(Element element, const int grainClass) const
    {
        auto elemId = gridGeometry_->elementMapper().index(element);
        return initialMassFractionsActiveLayer_[elemId][grainClass];
    }

    /*!
    * \brief Return the initial bed surface
    */
    Scalar& initialBedSurface(Element element) const
    {
        auto elemId = gridGeometry_->elementMapper().index(element);
        return initialBedSurface_[elemId];
    }

    /*!
    * \brief Return the key of the uppermost (still existing) layer
    *
    * Returns 0 if there is no uppermost layer because all layers are erroded.
    *
    * \param element
    */
    int indexUppermostLayer(const Element& element) const
    {
        auto elemId = gridGeometry_->elementMapper().index(element);
        std::vector<int> key;
        for(auto it = upperLayerBoundaries_[elemId].begin(); it != upperLayerBoundaries_[elemId].end(); ++it) {
            key.push_back(it->first);
        }
        if (key.size()==0) {
            return 0;
        }
        else {
            sort(key.begin(), key.end());
            return key[0];
        }
    }

    virtual ~LayerModel() {}

protected:
    /*!
    * \brief Create the active layer
    *
    * Initialize the following data member:
    *   - initialMassActiveLayer_
    *   - initialMassFractionsActiveLayer_
    *   - initialActiveLayerThickness_
    *   - sedimentMasses_
    */
    virtual void createActiveLayer_() const = 0;

    Scalar eps_ = 1e-9;
    int nLayer_;
    int nGrainClasses_;
    std::shared_ptr<const GridGeometry> gridGeometry_;
    std::shared_ptr<SpatialParams> spatialParams_;

    /*Initial data */
    inline static std::vector<Scalar> fixedGroundLevel_;
    inline static std::vector<Scalar> initialBedSurface_;
    inline static std::vector<Scalar> initialMassActiveLayer_;
    inline static std::map<std::string, std::vector<Scalar>> elementdata_;
    // initialMassFractions_.size() -> nElems
    // initialMassFractions_[i]::iterator->first -> layer index
    // initialMassFractions_[i]::iterator->second.size() -> nGrainClasses_
    inline static std::vector<std::map<int, std::vector<Scalar>>> initialMassFractions_;
    // initialMassFractionsActiveLayer_.size() -> nElems
    // initialMassFractionsActiveLayer_[i].size() -> nGrainClasses_
    inline static std::vector<std::vector<Scalar>> initialMassFractionsActiveLayer_;
    // initialActiveLayerThickness_.size() -> nElems
    inline static std::vector<Scalar> initialActiveLayerThickness_;
    // sedimentMasses_.size() -> nElems
    // sedimentMasses_[i]::iterator->first -> layer index
    // sedimentMasses_[i]::iterator->second.size() -> nGrainClasses_
    inline static std::vector<std::map<int, std::vector<Scalar>>> sedimentMasses_;
    // upperLayerBoundaries_.size() -> nElems
    // upperLayerBoundaries_[i]::iterator->first -> layer index
    inline static std::vector<std::map<int, Scalar>> upperLayerBoundaries_;
    // upperLayerBoundaries_.size() -> nElems
    // upperLayerBoundaries_[i]::iterator->first -> layer index
    inline static std::vector<std::map<int, Scalar>> harmonicAverageGrainDensities_;
};

} // end namespace Dumux

#endif // DUMUX_MATERIAL_LAYERMODEL_HH
