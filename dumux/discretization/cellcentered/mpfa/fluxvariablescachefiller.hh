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
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
    using DataHandle = typename PrimaryInteractionVolume::Traits::DataHandle;

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr bool doAdvection = GET_PROP_VALUE(TypeTag, EnableAdvection);
    static constexpr bool doDiffusion = GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion);
    static constexpr bool doHeatConduction = GET_PROP_VALUE(TypeTag, EnableEnergyBalance);

    static constexpr bool soldependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static constexpr bool soldependentDiffusion = GET_PROP_VALUE(TypeTag, SolutionDependentMolecularDiffusion);
    static constexpr bool soldependentHeatConduction = GET_PROP_VALUE(TypeTag, SolutionDependentHeatConduction);

public:
    static constexpr bool isSolDependent = (doAdvection && soldependentAdvection) ||
                                           (doDiffusion && soldependentDiffusion) ||
                                           (doHeatConduction && soldependentHeatConduction);

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
     * \param forceUpdateAll if true, forces all caches to be updated (even the solution-independent ones)
     */
    template<class FluxVariablesCacheContainer>
    void fill(FluxVariablesCacheContainer& fluxVarsCacheContainer,
              FluxVariablesCache& scvfFluxVarsCache,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              bool forceUpdateAll = false)
    {
        // Set pointers
        elementPtr_ = &element;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;

        // prepare interaction volume and fill caches of all the scvfs connected to it
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        if (fvGridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
        {
            if (forceUpdateAll)
            {
                // the local index of the interaction volume to be created in its container
                const auto ivIndexInContainer = fluxVarsCacheContainer.secondaryInteractionVolumes_.size();

                // prepare the locally cached boundary interaction volume
                const auto& indexSet = fvGridGeometry.gridInteractionVolumeIndexSets().secondaryIndexSet(scvf);
                fluxVarsCacheContainer.secondaryInteractionVolumes_.emplace_back();
                secondaryIv_ = &fluxVarsCacheContainer.secondaryInteractionVolumes_.back();
                secondaryIv_->setUpLocalScope(indexSet, problem(), fvGeometry);

                // prepare the corresponding data handle
                fluxVarsCacheContainer.secondaryIvDataHandles_.emplace_back();
                ivDataHandle_ = &fluxVarsCacheContainer.secondaryIvDataHandles_.back();
                secondaryIv_->prepareDataHandle(*ivDataHandle_);

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, *secondaryIv_, *ivDataHandle_, ivIndexInContainer, true);
            }
            else
            {
                const auto ivIndexInContainer = scvfFluxVarsCache.ivIndexInContainer();
                secondaryIv_ = &fluxVarsCacheContainer.secondaryInteractionVolumes_[ivIndexInContainer];
                ivDataHandle_ = &fluxVarsCacheContainer.secondaryIvDataHandles_[ivIndexInContainer];

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, *secondaryIv_, *ivDataHandle_, ivIndexInContainer);
            }
        }
        else
        {
            if (forceUpdateAll)
            {
                // the local index of the interaction volume to be created in its container
                const auto ivIndexInContainer = fluxVarsCacheContainer.primaryInteractionVolumes_.size();

                // prepare the locally cached boundary interaction volume
                const auto& indexSet = fvGridGeometry.gridInteractionVolumeIndexSets().primaryIndexSet(scvf);
                fluxVarsCacheContainer.primaryInteractionVolumes_.emplace_back();
                primaryIv_ = &fluxVarsCacheContainer.primaryInteractionVolumes_.back();
                primaryIv_->setUpLocalScope(indexSet, problem(), fvGeometry);

                // prepare the corresponding data handle
                fluxVarsCacheContainer.primaryIvDataHandles_.emplace_back();
                ivDataHandle_ = &fluxVarsCacheContainer.primaryIvDataHandles_.back();
                primaryIv_->prepareDataHandle(*ivDataHandle_);

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, *primaryIv_, *ivDataHandle_, ivIndexInContainer, true);
            }
            else
            {
                const auto ivIndexInContainer = scvfFluxVarsCache.ivIndexInContainer();
                primaryIv_ = &fluxVarsCacheContainer.primaryInteractionVolumes_[ivIndexInContainer];
                ivDataHandle_ = &fluxVarsCacheContainer.primaryIvDataHandles_[ivIndexInContainer];

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, *primaryIv_, *ivDataHandle_, ivIndexInContainer);
            }
        }
    }

    const PrimaryInteractionVolume& primaryInteractionVolume() const
    { return *primaryIv_; }

    const SecondaryInteractionVolume& secondaryInteractionVolume() const
    { return *secondaryIv_; }

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

    //! Method to fill the flux var caches within an interaction volume
    template<class FluxVariablesCacheContainer, class InteractionVolumeType>
    void fillCachesInInteractionVolume_(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                                        InteractionVolumeType& iv,
                                        DataHandle& handle,
                                        unsigned int ivIndexInContainer,
                                        bool forceUpdateAll = false)
    {
        // First we upate data which are not dependent on the physical processes.
        // We store pointers to the other flux var caches, so that we have to obtain
        // this data only once and can use it again in the sub-cache fillers.
        if (forceUpdateAll)
        {
            const auto numGlobalScvfs = iv.localFaceData().size();
            std::vector<const SubControlVolumeFace*> ivScvfs(numGlobalScvfs);
            std::vector<FluxVariablesCache*> ivFluxVarCaches(numGlobalScvfs);

            unsigned int i = 0;
            for (const auto& d : iv.localFaceData())
            {
                // obtain the scvf
                const auto& scvfJ = fvGeometry().scvf(d.globalScvfIndex());
                ivScvfs[i] = &scvfJ;
                ivFluxVarCaches[i] = &fluxVarsCacheContainer[scvfJ];
                ivFluxVarCaches[i]->setIvIndexInContainer(ivIndexInContainer);
                ivFluxVarCaches[i]->setUpdateStatus(true);
                i++;
            }

            fillAdvection(fluxVarsCacheContainer, iv, handle, ivScvfs, ivFluxVarCaches);
            fillDiffusion(fluxVarsCacheContainer, iv, handle, ivScvfs, ivFluxVarCaches);
            fillHeatConduction(fluxVarsCacheContainer, iv, handle, ivScvfs, ivFluxVarCaches);
        }
        else
        {
            const auto numGlobalScvfs = iv.localFaceData().size();
            std::vector<const SubControlVolumeFace*> ivScvfs(numGlobalScvfs);
            std::vector<FluxVariablesCache*> ivFluxVarCaches(numGlobalScvfs);

            unsigned int i = 0;
            for (const auto& d : iv.localFaceData())
            {
                // the iv index has been set already
                const auto& scvfJ = fvGeometry().scvf(d.globalScvfIndex());
                ivScvfs[i] = &scvfJ;
                ivFluxVarCaches[i] = &fluxVarsCacheContainer[scvfJ];
                ivFluxVarCaches[i]->setUpdateStatus(true);
                i++;
            }

            if (doAdvection && soldependentAdvection)
                fillAdvection(fluxVarsCacheContainer, iv, handle, ivScvfs, ivFluxVarCaches);
            if (doDiffusion && soldependentDiffusion)
                fillDiffusion(fluxVarsCacheContainer, iv, handle, ivScvfs, ivFluxVarCaches);
            if (doHeatConduction && soldependentHeatConduction)
                fillHeatConduction(fluxVarsCacheContainer, iv, handle, ivScvfs, ivFluxVarCaches);
        }
    }

    //! method to fill the advective quantities
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool advectionEnabled = doAdvection>
    typename std::enable_if<advectionEnabled>::type
    fillAdvection(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  InteractionVolumeType& iv,
                  DataHandle& handle,
                  const std::vector<const SubControlVolumeFace*>& ivScvfs,
                  const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
        using AdvectionFiller = typename AdvectionType::Cache::Filler;

        static constexpr auto AdvectionMethod = AdvectionType::myDiscretizationMethod;
        using LambdaFactory = TensorLambdaFactory<TypeTag, AdvectionMethod>;

        // set the advection context in the data handle
        handle.setAdvectionContext();

        // maybe solve the local system subject to K (if AdvectionType uses mpfa)
        if (AdvectionMethod == DiscretizationMethods::CCMpfa)
            iv.solveLocalSystem(LambdaFactory::getAdvectionLambda(), problem(), fvGeometry(), elemVolVars(), handle);

        // fill advection caches
        for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
            AdvectionFiller::fill(*ivFluxVarCaches[i],
                                  problem(),
                                  iv.element(iv.localFaceData()[i].ivLocalInsideScvIndex()),
                                  fvGeometry(),
                                  elemVolVars(),
                                  *ivScvfs[i],
                                  *this);
    }

    //! do nothing if advection is not enabled
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool advectionEnabled = doAdvection>
    typename std::enable_if<!advectionEnabled>::type
    fillAdvection(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  InteractionVolumeType& iv,
                  DataHandle& handle,
                  const std::vector<const SubControlVolumeFace*>& ivScvfs,
                  const std::vector<FluxVariablesCache*>& ivFluxVarCaches) {}

    //! method to fill the diffusive quantities
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool diffusionEnabled = doDiffusion>
    typename std::enable_if<diffusionEnabled>::type
    fillDiffusion(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  InteractionVolumeType& iv,
                  DataHandle& handle,
                  const std::vector<const SubControlVolumeFace*>& ivScvfs,
                  const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {
        using DiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
        using DiffusionFiller = typename DiffusionType::Cache::Filler;

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
                    iv.solveLocalSystem(LambdaFactory::getDiffusionLambda(phaseIdx, compIdx), problem(), fvGeometry(), elemVolVars(), handle);

                // fill diffusion caches
                for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
                    DiffusionFiller::fill(*ivFluxVarCaches[i],
                                          phaseIdx,
                                          compIdx,
                                          problem(),
                                          iv.element(iv.localFaceData()[i].ivLocalInsideScvIndex()),
                                          fvGeometry(),
                                          elemVolVars(),
                                          *ivScvfs[i],
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
                  const std::vector<const SubControlVolumeFace*>& ivScvfs,
                  const std::vector<FluxVariablesCache*>& ivFluxVarCaches) {}

    //! method to fill the quantities related to heat conduction
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool heatConductionEnabled = doHeatConduction>
    typename std::enable_if<heatConductionEnabled>::type
    fillHeatConduction(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                       InteractionVolumeType& iv,
                       DataHandle& handle,
                       const std::vector<const SubControlVolumeFace*>& ivScvfs,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {
        using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);
        using HeatConductionFiller = typename HeatConductionType::Cache::Filler;

        static constexpr auto HeatConductionMethod = HeatConductionType::myDiscretizationMethod;
        using LambdaFactory = TensorLambdaFactory<TypeTag, HeatConductionMethod>;

        // set the advection context in the data handle
        handle.setHeatConductionContext();

        // maybe solve the local system subject to fourier coefficient
        if (HeatConductionMethod == DiscretizationMethods::CCMpfa)
            iv.solveLocalSystem(LambdaFactory::getHeatConductionLambda(), problem(), fvGeometry(), elemVolVars(), handle);

        // fill heat conduction caches
        for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
            HeatConductionFiller::fill(*ivFluxVarCaches[i],
                                       problem(),
                                       iv.element(iv.localFaceData()[i].ivLocalInsideScvIndex()),
                                       fvGeometry(),
                                       elemVolVars(),
                                       *ivScvfs[i],
                                       *this);
    }

    //! do nothing if heat conduction is disabled
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool heatConductionEnabled = doHeatConduction>
    typename std::enable_if<!heatConductionEnabled>::type
    fillHeatConduction(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                       InteractionVolumeType& iv,
                       DataHandle& handle,
                       const std::vector<const SubControlVolumeFace*>& ivScvfs,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches) {}

    const Problem* problemPtr_;
    const Element* elementPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;

    // We store pointers to an inner and a boundary interaction volume.
    // These are updated during the filling of the caches and the
    // physics-related caches have access to them
    PrimaryInteractionVolume* primaryIv_;
    SecondaryInteractionVolume* secondaryIv_;

    // pointer to the current interaction volume data handle
    DataHandle* ivDataHandle_;
};

} // end namespace

#endif
