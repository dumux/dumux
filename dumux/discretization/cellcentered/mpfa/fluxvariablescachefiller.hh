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
/*!
 * \file
 * \brief The flux variables cache filler class
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_FLUXVARSCACHE_FILLER_HH
#define DUMUX_DISCRETIZATION_CCMPFA_FLUXVARSCACHE_FILLER_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/tensorlambdafactory.hh>

namespace Dumux
{

//! forward declaration of properties
namespace Properties
{
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(NumComponents);
};

/*!
 * \ingroup ImplicitModel
 * \brief Helper class to fill the flux var caches
 */
template<class TypeTag>
class CCMpfaFluxVariablesCacheFiller
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using DataHandle = typename BoundaryInteractionVolume::Traits::DataHandle;

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr bool doAdvection = GET_PROP_VALUE(TypeTag, EnableAdvection);
    static constexpr bool doDiffusion = GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion);
    static constexpr bool doHeatConduction = GET_PROP_VALUE(TypeTag, EnableEnergyBalance);

    static constexpr bool soldependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static constexpr bool soldependentDiffusion = GET_PROP_VALUE(TypeTag, SolutionDependentMolecularDiffusion);
    static constexpr bool soldependentHeatConduction = GET_PROP_VALUE(TypeTag, SolutionDependentHeatConduction);

public:
    //! The constructor. Sets the problem pointer
    CCMpfaFluxVariablesCacheFiller(const Problem& problem) : problemPtr_(&problem) {}

    /*!
     * \brief function to fill the flux variables caches
     *
     * \param fluxVarsCacheContainer Either the element or global flux variables cache
     * \param scvfFluxVarsCache The flux var cache to be updated corresponding to the given scvf
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables
     * \param scvf The corresponding sub-control volume face
     * \param doSubCaches Array of bools indicating which sub caches have to be updated
     */
    template<class FluxVariablesCacheContainer>
    void fill(FluxVariablesCacheContainer& fluxVarsCacheContainer,
              FluxVariablesCache& scvfFluxVarsCache,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              bool isUpdate = false)
    {
        // Set pointers
        elementPtr_ = &element;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;
        scvfPtr_ = &scvf;

        // prepare interaction volume and fill caches of all the scvfs connected to it
        const auto& globalFvGeometry = problem().model().globalFvGeometry();
        if (globalFvGeometry.isInBoundaryInteractionVolume(scvf))
        {
            if (!isUpdate)
            {
                // prepare the locally cached boundary interaction volume
                fluxVarsCacheContainer.boundaryInteractionVolumes_.emplace_back();
                bIv_ = &fluxVarsCacheContainer.boundaryInteractionVolumes_.back();

                // prepare the corresponding data handle
                fluxVarsCacheContainer.boundaryIvDataHandles_.emplace_back();
                ivDataHandle_ = &fluxVarsCacheContainer.boundaryIvDataHandles_.back();

                // the local index of the actual interaction volume in its container
                const auto ivIndexInContainer = fluxVarsCacheContainer.boundaryInteractionVolumes_.size()-1;

                bIv_->bind(globalFvGeometry.boundaryInteractionVolumeSeed(scvf),
                           problem(),
                           fvGeometry,
                           elemVolVars,
                           *ivDataHandle_);

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, boundaryInteractionVolume(), *ivDataHandle_, ivIndexInContainer);
            }
            else
            {
                const auto ivIndexInContainer = scvfFluxVarsCache.ivIndexInContainer();
                bIv_ = &fluxVarsCacheContainer.boundaryInteractionVolumes_[ivIndexInContainer];
                ivDataHandle_ = &fluxVarsCacheContainer.boundaryIvDataHandles_[ivIndexInContainer];

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, boundaryInteractionVolume(), *ivDataHandle_, ivIndexInContainer, true);
            }
        }
        else
        {
            if (!isUpdate)
            {
                // prepare the locally cached boundary interaction volume
                fluxVarsCacheContainer.interactionVolumes_.emplace_back();
                iv_ = &fluxVarsCacheContainer.interactionVolumes_.back();

                // prepare the corresponding data handle
                fluxVarsCacheContainer.ivDataHandles_.emplace_back();
                ivDataHandle_ = &fluxVarsCacheContainer.ivDataHandles_.back();

                // the local index of the actual interaction volume in its container
                const auto ivIndexInContainer = fluxVarsCacheContainer.interactionVolumes_.size()-1;

                iv_->bind(globalFvGeometry.interactionVolumeSeed(scvf),
                          problem(),
                          fvGeometry,
                          elemVolVars,
                          *ivDataHandle_);

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, interactionVolume(), *ivDataHandle_, ivIndexInContainer);
            }
            else
            {
                const auto ivIndexInContainer = scvfFluxVarsCache.ivIndexInContainer();
                iv_ = &fluxVarsCacheContainer.interactionVolumes_[ivIndexInContainer];
                ivDataHandle_ = &fluxVarsCacheContainer.ivDataHandles_[ivIndexInContainer];

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, interactionVolume(), *ivDataHandle_, ivIndexInContainer, true);
            }
        }
    }

    /*!
     * \brief function to update the flux variables caches during derivative calculation
     *
     * \copydoc fill
     */
    template<class FluxVariablesCacheContainer>
    void update(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                FluxVariablesCache& scvfFluxVarsCache,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf)
    {
        // forward to fill routine
        fill(fluxVarsCacheContainer, scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf, true);
    }

    static bool isSolutionDependent()
    {
        static const bool isSolDependent = (doAdvection && soldependentAdvection) ||
                                           (doDiffusion && soldependentDiffusion) ||
                                           (doHeatConduction && soldependentHeatConduction);
        return isSolDependent;
    }

    const InteractionVolume& interactionVolume() const
    { return *iv_; }

    const BoundaryInteractionVolume& boundaryInteractionVolume() const
    { return *bIv_; }

    const DataHandle& dataHandle() const
    { return *ivDataHandle_; }

private:

    const Problem& problem() const
    { return *problemPtr_; }

    const Element& element() const
    { return *elementPtr_; }

    const FVElementGeometry& fvGeometry() const
    { return *fvGeometryPtr_; }

    const ElementVolumeVariables& elemVolVars() const
    { return *elemVolVarsPtr_; }

    const SubControlVolumeFace& scvFace() const
    { return *scvfPtr_; }

    InteractionVolume& interactionVolume()
    { return *iv_; }

    BoundaryInteractionVolume& boundaryInteractionVolume()
    { return *bIv_; }

    //! Method to fill the flux var caches within an interaction volume
    template<class FluxVariablesCacheContainer, class InteractionVolumeType>
    void fillCachesInInteractionVolume_(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                                        InteractionVolumeType& iv,
                                        DataHandle& handle,
                                        unsigned int ivIndexInContainer,
                                        bool isUpdate = false)
    {
        // First we upate data which are not dependent on the physical processes.
        // We store pointers to the other flux var caches, so that we have to obtain
        // this data only once and can use it again in the sub-cache fillers.
        if (!isUpdate)
        {
            std::vector<FluxVariablesCache*> ivFluxVarCaches(iv.globalLocalScvfPairedData().size());
            unsigned int i = 0;
            for (const auto& dataPair : iv.globalLocalScvfPairedData())
            {
                // obtain the scvf
                const auto& scvfJ = *dataPair.first;
                ivFluxVarCaches[i] = &fluxVarsCacheContainer[scvfJ];
                ivFluxVarCaches[i]->updateInteriorBoundaryData(iv, scvfJ);
                ivFluxVarCaches[i]->setIvIndexInContainer(ivIndexInContainer);
                ivFluxVarCaches[i]->setUpdateStatus(true);
                i++;
            }

            fillAdvection(fluxVarsCacheContainer, iv, handle, ivFluxVarCaches);
            fillDiffusion(fluxVarsCacheContainer, iv, handle, ivFluxVarCaches);
            fillHeatConduction(fluxVarsCacheContainer, iv, handle, ivFluxVarCaches);
        }
        else
        {
            std::vector<FluxVariablesCache*> ivFluxVarCaches(iv.globalLocalScvfPairedData().size());
            unsigned int i = 0;
            for (const auto& dataPair : iv.globalLocalScvfPairedData())
            {
                // the interior boundary data and the iv index have been set already
                ivFluxVarCaches[i] = &fluxVarsCacheContainer[*dataPair.first];
                ivFluxVarCaches[i]->setUpdateStatus(true);
                i++;
            }

            if (doAdvection && soldependentAdvection)
                fillAdvection(fluxVarsCacheContainer, iv, handle, ivFluxVarCaches);
            if (doDiffusion && soldependentDiffusion)
                fillDiffusion(fluxVarsCacheContainer, iv, handle, ivFluxVarCaches);
            if (doHeatConduction && soldependentHeatConduction)
                fillHeatConduction(fluxVarsCacheContainer, iv, handle, ivFluxVarCaches);
        }
    }

    //! method to fill the advective quantities
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool advectionEnabled = doAdvection>
    typename std::enable_if<advectionEnabled>::type
    fillAdvection(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  InteractionVolumeType& iv,
                  DataHandle& handle,
                  const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
        using AdvectionFiller = typename AdvectionType::CacheFiller;

        static constexpr auto AdvectionMethod = AdvectionType::myDiscretizationMethod;
        using LambdaFactory = TensorLambdaFactory<TypeTag, AdvectionMethod>;

        // set the advection context in the data handle
        handle.setAdvectionContext();

        // maybe solve the local system subject to K (if AdvectionType uses mpfa)
        if (AdvectionMethod == DiscretizationMethods::CCMpfa)
            iv.solveLocalSystem(LambdaFactory::getAdvectionLambda(), handle);

        // fill advection caches
        unsigned int i = 0;
        for (const auto& dataPair : iv.globalLocalScvfPairedData())
            AdvectionFiller::fill(*ivFluxVarCaches[i++],
                                  problem(),
                                  iv.localElement(dataPair.second.localScvIndex),
                                  fvGeometry(),
                                  elemVolVars(),
                                  *dataPair.first,
                                  *this);
    }

    //! do nothing if advection is not enabled
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool advectionEnabled = doAdvection>
    typename std::enable_if<!advectionEnabled>::type
    fillAdvection(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  InteractionVolumeType& iv,
                  DataHandle& handle,
                  const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {}

    //! method to fill the diffusive quantities
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool diffusionEnabled = doDiffusion>
    typename std::enable_if<diffusionEnabled>::type
    fillDiffusion(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  InteractionVolumeType& iv,
                  DataHandle& handle,
                  const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {
        using DiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
        using DiffusionFiller = typename DiffusionType::CacheFiller;

        static constexpr auto DiffusionMethod = DiffusionType::myDiscretizationMethod;
        using LambdaFactory = TensorLambdaFactory<TypeTag, DiffusionMethod>;

        static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
        static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                if (phaseIdx == compIdx)
                    continue;

                // set the diffusion context in the data handle
                handle.setDiffusionContext(phaseIdx, compIdx);

                // solve the local system subject to the diffusion tensor (if uses mpfa)
                if (DiffusionMethod == DiscretizationMethods::CCMpfa)
                    iv.solveLocalSystem(LambdaFactory::getDiffusionLambda(phaseIdx, compIdx), handle);

                // fill diffusion caches
                unsigned int i = 0;
                for (const auto& dataPair : iv.globalLocalScvfPairedData())
                    DiffusionFiller::fill(*ivFluxVarCaches[i++],
                                          phaseIdx,
                                          compIdx,
                                          problem(),
                                          iv.localElement(dataPair.second.localScvIndex),
                                          fvGeometry(),
                                          elemVolVars(),
                                          *dataPair.first,
                                          *this);
            }
        }
    }

    //! do nothing if diffusion is not enabled
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool diffusionEnabled = doDiffusion>
    typename std::enable_if<!diffusionEnabled>::type
    fillDiffusion(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  InteractionVolumeType& iv,
                  DataHandle& handle,
                  const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {}

    //! method to fill the quantities related to heat conduction
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool heatConductionEnabled = doHeatConduction>
    typename std::enable_if<heatConductionEnabled>::type
    fillHeatConduction(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                       InteractionVolumeType& iv,
                       DataHandle& handle,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {
        using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);
        using HeatConductionFiller = typename HeatConductionType::CacheFiller;

        static constexpr auto HeatConductionMethod = HeatConductionType::myDiscretizationMethod;
        using LambdaFactory = TensorLambdaFactory<TypeTag, HeatConductionMethod>;

        // set the advection context in the data handle
        handle.setHeatConductionContext();

        // maybe solve the local system subject to fourier coefficient
        if (HeatConductionMethod == DiscretizationMethods::CCMpfa)
            iv.solveLocalSystem(LambdaFactory::getHeatConductionLambda(), handle);

        // fill heat conduction caches
        unsigned int i = 0;
        for (const auto& dataPair : iv.globalLocalScvfPairedData())
            HeatConductionFiller::fill(*ivFluxVarCaches[i++],
                                       problem(),
                                       iv.localElement(dataPair.second.localScvIndex),
                                       fvGeometry(),
                                       elemVolVars(),
                                       *dataPair.first,
                                       *this);
    }

    //! do nothing if heat conduction is disabled
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool heatConductionEnabled = doHeatConduction>
    typename std::enable_if<!heatConductionEnabled>::type
    fillHeatConduction(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                       InteractionVolumeType& iv,
                       DataHandle& handle,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {}

    const Problem* problemPtr_;
    const Element* elementPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;
    const SubControlVolumeFace* scvfPtr_;

    // We store pointers to an inner and boundary interaction volume
    // these are updated during the filling of the caches and the
    // physics-related caches have access to them
    InteractionVolume* iv_;
    BoundaryInteractionVolume* bIv_;

    // pointer to the current interaction volume data handle
    DataHandle* ivDataHandle_;
};

} // end namespace

#endif
